"""
Builds every figure and table of the MDCheck validation manuscript from validation/results/*.csv.

    python validation/make_figures.py [--results validation/results] [--out validation/figures]

Figure 3 re-runs one short Lennard-Jones replica with OpenMM (seed fixed) to show a raw trace.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

COLORS = {
    "naive": "#9e9e9e",
    "bootstrap": "#d95f02",
    "mdcheck": "#1b6ca8",
    "pymbar": "#2ca25f",
    "pyblock": "#7b3294",
}


def fig_inefficiency(summ, out):
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    s = summ[summ.phi > 0]
    lim = [0.8, s.g_true.max() * 1.4]
    ax.plot(lim, lim, color="black", lw=0.8, ls="--", label="exact")
    ax.errorbar(s.g_true * 0.97, s.g_mdcheck_mean, yerr=s.g_mdcheck_sd, fmt="o", ms=4,
                color=COLORS["mdcheck"], label="MDCheck")
    ax.errorbar(s.g_true * 1.03, s.g_pymbar_mean, yerr=s.g_pymbar_sd, fmt="s", ms=4,
                color=COLORS["pymbar"], label="pymbar")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"exact $g = (1+\phi)/(1-\phi)$")
    ax.set_ylabel(r"estimated $g$ (mean $\pm$ SD)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"fig1_inefficiency.{ext}"), dpi=300)
    plt.close(fig)


def fig_coverage(summ, reps, out):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    series = [
        ("cov_naive", "naive $s/\\sqrt{N}$", "naive", "x"),
        ("cov_mdcheck_bootstrap_v100", "MDCheck 1.0.0 block bootstrap", "bootstrap", "v"),
        ("cov_pyblock", "pyblock (optimal block)", "pyblock", "D"),
        ("cov_g_pymbar", "$g$ from pymbar", "pymbar", "s"),
        ("cov_mdcheck_ci", "MDCheck $t$ interval", "mdcheck", "o"),
    ]
    x = np.arange(len(summ))
    band = 1.96 * np.sqrt(0.95 * 0.05 / reps)
    ax.axhspan(0.95 - band, 0.95 + band, color="black", alpha=0.08, lw=0)
    ax.axhline(0.95, color="black", lw=0.8, ls="--")
    for col, lab, c, mk in series:
        if col in summ:
            ax.plot(x, summ[col], marker=mk, ms=4, lw=1, color=COLORS[c], label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p:g}" for p in summ.phi])
    ax.set_xlabel(r"AR(1) coefficient $\phi$")
    ax.set_ylabel("coverage of nominal 95% CI")
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"fig2_coverage.{ext}"), dpi=300)
    plt.close(fig)


def fig_lj_trace(out):
    try:
        from validation.validate_md_openmm import build_lj, run_replica, pick_platform
    except Exception:
        sys.path.insert(0, os.path.dirname(__file__))
        from validate_md_openmm import build_lj, run_replica, pick_platform
    from mdcheck.core.equilibration import detect_equilibration
    from pymbar import timeseries

    system, pos, cfg = build_lj()
    n = int(round(150.0 / (cfg["dt_fs"] * 1e-3 * cfg["report_steps"])))
    x = run_replica(system, pos, cfg, n, 1234, pick_platform())["potential_energy_kJmol"]
    dt = cfg["dt_fs"] * 1e-3 * cfg["report_steps"]
    t = np.arange(n) * dt
    md = detect_equilibration(x)
    t0_pm, _, _ = timeseries.detect_equilibration(x)
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    ax.plot(t, x / 864.0, lw=0.6, color="#444444")
    ax.axvline(md["t_eq_index"] * dt, color=COLORS["mdcheck"], lw=1.2, label=f"MDCheck $t_{{eq}}$ = {md['t_eq_index'] * dt:.1f} ps")
    ax.axvline(int(t0_pm) * dt, color=COLORS["pymbar"], lw=1.2, ls=":", label=f"pymbar $t_{{eq}}$ = {int(t0_pm) * dt:.1f} ps")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("potential energy per atom (kJ/mol)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"fig3_lj_trace.{ext}"), dpi=300)
    plt.close(fig)


def fig_md_coverage(res_dir, out):
    frames = []
    for sub, tag in (("results", ""), ("results_water300", " 300 ps")):
        p = os.path.join(os.path.dirname(res_dir.rstrip("/\\")), sub, "md_final_summary.csv")
        if os.path.exists(p):
            part = pd.read_csv(p)
            part["tag"] = tag
            frames.append(part)
    if not frames:
        return None
    d = pd.concat(frames, ignore_index=True)
    names = {"potential_energy_kJmol": "U", "density_gcm3": r"$\rho$"}
    labels = [f"{'LJ' if r.system == 'lj' else 'TIP3P'} {names.get(r.observable, r.observable)}{r.tag}\n(n = {r.n_replicas})"
              for r in d.itertuples()]
    x = np.arange(len(d))
    w = 0.26
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.bar(x - w, d.coverage_naive, w, color=COLORS["naive"], label="naive SE")
    ax.bar(x, d.coverage_mdcheck, w, color=COLORS["mdcheck"], label="MDCheck")
    ax.bar(x + w, d.coverage_extra_discard, w, color=COLORS["pymbar"], label="MDCheck + 20 ps discard")
    for xi, n in zip(x, d.n_replicas):
        band = 1.96 * np.sqrt(0.95 * 0.05 / n)
        ax.plot([xi - 1.5 * w, xi + 1.5 * w], [0.95 - band] * 2, color="black", lw=0.6, ls=":")
    ax.axhline(0.95, color="black", lw=0.8, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("coverage of reference mean")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"fig4_md_coverage.{ext}"), dpi=300)
    plt.close(fig)
    return d


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(__file__)
    ap.add_argument("--results", default=os.path.join(here, "results"))
    ap.add_argument("--out", default=os.path.join(here, "figures"))
    ap.add_argument("--skip-trace", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    summ = pd.read_csv(os.path.join(args.results, "ar1_inefficiency_coverage.csv"))
    reps = json.load(open(os.path.join(args.results, "versions.json")))["reps"]
    fig_inefficiency(summ, args.out)
    fig_coverage(summ, reps, args.out)
    fig_md_coverage(args.results, args.out)
    if not args.skip_trace:
        fig_lj_trace(args.out)
    print("figures written to", args.out)


if __name__ == "__main__":
    main()
