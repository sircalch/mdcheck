"""
Validation of MDCheck on real molecular dynamics generated with OpenMM.

Systems
-------
  lj     : 864-atom Lennard-Jones argon (sigma = 3.405 A, epsilon/kB = 119.8 K), rho* = 0.8442,
           T = 94.4 K (T* = 0.788), started from an FCC lattice so every replica has a genuine
           melting / equilibration transient.
  water  : TIP3P water box (~3 nm, PME, NPT at 300 K / 1 bar), started from a Modeller box
           heated from 200 K velocities.

Protocol
--------
  * `n_long` long reference runs give the reference ensemble mean of each observable
    (first `ref_discard_ps` discarded).
  * `n_short` independent short replicas are analysed exactly as a user would: MDCheck
    detects t_eq, estimates g / N_eff and reports a 95% CI of the production mean.
    pymbar.timeseries.detect_equilibration is run on the same data.
  * Reported: t_eq (MDCheck vs pymbar), g, and whether each CI covers the reference mean
    (combining the reference standard error in quadrature is unnecessary when it is <10% of the
    replica CI half-width; this ratio is reported).

Usage:
    python validation/validate_md_openmm.py --system lj --out validation/results
    python validation/validate_md_openmm.py --system water --out validation/results
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

import openmm as mm
from openmm import app, unit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import mdcheck  # noqa: E402
from mdcheck.core.autocorrelation import integrated_autocorrelation_time  # noqa: E402
from mdcheck.core.bootstrap import inefficiency_corrected_ci  # noqa: E402
from mdcheck.core.equilibration import detect_equilibration  # noqa: E402
from mdcheck.core.replicas import assess_replica_consistency  # noqa: E402
from pymbar import timeseries  # noqa: E402


def pick_platform():
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    for pref in ("CUDA", "OpenCL", "CPU"):
        if pref in names:
            return mm.Platform.getPlatformByName(pref)
    return mm.Platform.getPlatformByName("Reference")


# ----------------------------------------------------------------------------------------------
# Lennard-Jones argon
# ----------------------------------------------------------------------------------------------
def build_lj():
    sigma = 0.3405  # nm
    epsilon = (119.8 * unit.kelvin * unit.MOLAR_GAS_CONSTANT_R).value_in_unit(unit.kilojoule_per_mole)
    mass = 39.948
    ncell = 6
    rho_star = 0.8442
    n = 4 * ncell ** 3
    box = (n / rho_star) ** (1.0 / 3.0) * sigma
    a = box / ncell
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    pos = np.array([(np.array([i, j, k]) + b) * a
                    for i in range(ncell) for j in range(ncell) for k in range(ncell) for b in basis])

    system = mm.System()
    system.setDefaultPeriodicBoxVectors(mm.Vec3(box, 0, 0), mm.Vec3(0, box, 0), mm.Vec3(0, 0, box))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(3.0 * sigma)
    nb.setUseDispersionCorrection(True)
    for _ in range(n):
        system.addParticle(mass)
        nb.addParticle(0.0, sigma, epsilon)
    system.addForce(nb)
    return system, pos, {"temperature_K": 94.4, "dt_fs": 10.0, "report_steps": 25,
                         "observables": ["potential_energy_kJmol"], "barostat": False}


# ----------------------------------------------------------------------------------------------
# TIP3P water
# ----------------------------------------------------------------------------------------------
def build_water():
    ff = app.ForceField("amber14/tip3p.xml")
    modeller = app.Modeller(app.Topology(), [])
    modeller.addSolvent(ff, model="tip3p", boxSize=mm.Vec3(3.0, 3.0, 3.0) * unit.nanometer)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=0.9 * unit.nanometer, constraints=app.HBonds,
                             rigidWater=True)
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 25))
    pos = np.array(modeller.positions.value_in_unit(unit.nanometer))
    return system, pos, {"temperature_K": 300.0, "dt_fs": 2.0, "report_steps": 100,
                         "observables": ["potential_energy_kJmol", "density_gcm3"], "barostat": True,
                         "start_velocity_K": 200.0}


def run_replica(system, pos, cfg, n_samples, seed, platform):
    integ = mm.LangevinMiddleIntegrator(cfg["temperature_K"] * unit.kelvin, 1.0 / unit.picosecond,
                                        cfg["dt_fs"] * unit.femtosecond)
    integ.setRandomNumberSeed(seed)
    if cfg["barostat"]:
        for f in system.getForces():
            if isinstance(f, mm.MonteCarloBarostat):
                f.setRandomNumberSeed(seed)
    ctx = mm.Context(system, integ, platform)
    ctx.setPositions(pos * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, maxIterations=200)
    ctx.setVelocitiesToTemperature(cfg.get("start_velocity_K", cfg["temperature_K"]) * unit.kelvin, seed)
    total_mass = sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))

    out = {k: np.empty(n_samples) for k in cfg["observables"]}
    for i in range(n_samples):
        integ.step(cfg["report_steps"])
        st = ctx.getState(getEnergy=True)
        out["potential_energy_kJmol"][i] = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if "density_gcm3" in out:
            vol_nm3 = st.getPeriodicBoxVolume().value_in_unit(unit.nanometer ** 3)
            out["density_gcm3"][i] = total_mass * 1.66053906660e-3 / vol_nm3  # 1 amu/nm^3 = 1.66054e-3 g/cm^3
    del ctx, integ
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", choices=["lj", "water"], default="lj")
    ap.add_argument("--n-short", type=int, default=24)
    ap.add_argument("--n-long", type=int, default=4)
    ap.add_argument("--short-ps", type=float, default=None)
    ap.add_argument("--long-ps", type=float, default=None)
    ap.add_argument("--ref-discard-ps", type=float, default=None)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    defaults = {"lj": (150.0, 1500.0, 100.0), "water": (100.0, 600.0, 50.0)}[args.system]
    short_ps = args.short_ps or defaults[0]
    long_ps = args.long_ps or defaults[1]
    discard_ps = args.ref_discard_ps or defaults[2]

    system, pos, cfg = (build_lj if args.system == "lj" else build_water)()
    platform = pick_platform()
    dt_ps = cfg["dt_fs"] * 1e-3 * cfg["report_steps"]
    n_short = int(round(short_ps / dt_ps))
    n_long = int(round(long_ps / dt_ps))
    n_disc = int(round(discard_ps / dt_ps))
    t_start = time.time()

    # Reference runs
    ref = {k: [] for k in cfg["observables"]}
    raw_long, raw_short = [], []
    for r in range(args.n_long):
        o = run_replica(system, pos, cfg, n_long, args.seed + 1000 + r, platform)
        raw_long.append(o)
        for k in ref:
            ref[k].append(o[k][n_disc:])
    ref_stats = {}
    for k, runs in ref.items():
        means = np.array([x.mean() for x in runs])
        # SE from independent long runs (each also checked with g below)
        gs = [integrated_autocorrelation_time(x)[1] for x in runs]
        pooled = np.concatenate(runs)
        se_g = pooled.std(ddof=1) * np.sqrt(np.mean(gs) / len(pooled))
        ref_stats[k] = {"mean": float(means.mean()), "se_between_runs": float(means.std(ddof=1) / np.sqrt(len(means))),
                        "se_g": float(se_g), "g_mean": float(np.mean(gs))}

    # Short replicas analysed as a user would
    rows = []
    prod_segments = {k: [] for k in cfg["observables"]}
    for r in range(args.n_short):
        o = run_replica(system, pos, cfg, n_short, args.seed + r, platform)
        raw_short.append(o)
        for k in cfg["observables"]:
            x = o[k]
            md = detect_equilibration(x)
            t0_pm, g_pm, neff_pm = timeseries.detect_equilibration(x)
            prod = x[md["t_eq_index"]:]
            prod_segments[k].append(prod)
            m, lo, hi = inefficiency_corrected_ci(prod, g=md["g"])
            ref_m = ref_stats[k]["mean"]
            rows.append({
                "system": args.system, "observable": k, "replica": r,
                "t_eq_ps_mdcheck": md["t_eq_index"] * dt_ps, "t_eq_ps_pymbar": int(t0_pm) * dt_ps,
                "g_mdcheck": md["g"], "g_pymbar": float(g_pm),
                "neff_mdcheck": md["n_eff"], "neff_pymbar": float(neff_pm),
                "mean_prod": m, "ci_lo": lo, "ci_hi": hi,
                "ref_mean": ref_m, "covers_ref": lo <= ref_m <= hi,
                "naive_covers_ref": abs(prod.mean() - ref_m) <= 1.96 * prod.std(ddof=1) / np.sqrt(len(prod)),
                "no_discard_covers_ref": abs(x.mean() - ref_m) <= 1.96 * x.std(ddof=1) * np.sqrt(md["g"] / len(x)),
                "ref_se_over_halfwidth": ref_stats[k]["se_between_runs"] / max(1e-300, (hi - lo) / 2.0),
            })
    # Raw series, so alternative estimators can be evaluated on exactly the same trajectories
    np.savez_compressed(
        os.path.join(args.out, f"md_{args.system}_series.npz"),
        **{f"long_{k}": np.array([o[k] for o in raw_long]) for k in cfg["observables"]},
        **{f"short_{k}": np.array([o[k] for o in raw_short]) for k in cfg["observables"]},
        sample_interval_ps=dt_ps, ref_discard_samples=n_disc)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out, f"md_{args.system}_raw.csv"), index=False)
    summ = df.groupby("observable").agg(
        n_replicas=("replica", "count"),
        t_eq_ps_mdcheck_median=("t_eq_ps_mdcheck", "median"),
        t_eq_ps_pymbar_median=("t_eq_ps_pymbar", "median"),
        g_mdcheck_median=("g_mdcheck", "median"), g_pymbar_median=("g_pymbar", "median"),
        coverage_mdcheck=("covers_ref", "mean"),
        coverage_naive_se=("naive_covers_ref", "mean"),
        coverage_without_discard=("no_discard_covers_ref", "mean"),
        ref_se_over_halfwidth_median=("ref_se_over_halfwidth", "median"),
    ).reset_index()
    # Replica-consistency false-alarm rate: groups of 3 replicas of the same ensemble should pass.
    rep_rows = []
    for k, segs in prod_segments.items():
        for gi in range(len(segs) // 3):
            res = assess_replica_consistency(segs[3 * gi:3 * gi + 3])
            rep_rows.append({"observable": k, "group": gi, "status": res["status"],
                             "max_jsd": res["max_jsd"], "mean_jsd": res["mean_jsd"]})
    rep = pd.DataFrame(rep_rows)
    rep.to_csv(os.path.join(args.out, f"md_{args.system}_replica_consistency.csv"), index=False)
    rep_s = rep.groupby("observable").agg(
        n_groups=("group", "count"),
        replica_false_alarm_rate=("status", lambda s: float((s != "PASS").mean())),
        max_jsd_median=("max_jsd", "median"),
    ).reset_index()
    summ = summ.merge(rep_s, on="observable", how="left")
    summ.to_csv(os.path.join(args.out, f"md_{args.system}_summary.csv"), index=False)
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    meta = {"system": args.system, "platform": platform.getName(), "openmm": mm.__version__,
            "mdcheck": mdcheck.__version__, "short_ps": short_ps, "long_ps": long_ps,
            "ref_discard_ps": discard_ps, "sample_interval_ps": dt_ps, "n_short": args.n_short,
            "n_long": args.n_long, "seed": args.seed, "reference": ref_stats,
            "runtime_s": round(time.time() - t_start, 1), "config": cfg}
    with open(os.path.join(args.out, f"md_{args.system}_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(json.dumps(meta["reference"], indent=2), "\nruntime_s", meta["runtime_s"], "platform", meta["platform"])


if __name__ == "__main__":
    main()
