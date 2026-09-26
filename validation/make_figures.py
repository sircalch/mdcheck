"""
Publication figures for the MDCheck validation manuscript, built only from validation/results*/.

    python validation/make_figures.py [--out validation/figures]

Style: double-column width 174 mm (Taylor & Francis), 8 pt sans-serif text, one fixed colour per
method across all figures (validated categorical palette; every series also has its own marker
and a legend entry, so identity never relies on colour alone).
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

HERE = os.path.dirname(__file__)
MM = 1 / 25.4
DOUBLE = 174 * MM

# Fixed method -> style mapping (colour, marker, label); the same in every figure.
STYLE = {
    "naive": ("#8a8984", "x", "naive $s/\\sqrt{N}$"),
    "bootstrap": ("#eda100", "v", "MDCheck 1.0.0 block bootstrap"),
    "pyblock": ("#1baf7a", "D", "pyblock (optimal block)"),
    "pymbar": ("#eb6834", "s", "pymbar $g$"),
    "mdcheck": ("#2a78d6", "o", "MDCheck 1.1.0"),
    "rhat101": ("#e87ba4", "^", r"$\hat{R} > 1.01$"),
    "rhat11": ("#4a3aa7", "P", r"$\hat{R} > 1.1$"),
}
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"


def setup():
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "legend.frameon": False,
        "lines.linewidth": 1.2, "lines.markersize": 4,
        "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel(ax, letter):
    ax.text(-0.14, 1.04, f"({letter})", transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=INK)


def mc_band(ax, n):
    band = 1.96 * np.sqrt(0.95 * 0.05 / n)
    ax.axhspan(0.95 - band, 0.95 + band, color=INK2, alpha=0.08, lw=0)
    ax.axhline(0.95, color=INK2, lw=0.7, ls="--")


def save(fig, name, out):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def fig_ar1(res, out):
    s = pd.read_csv(os.path.join(res, "ar1_inefficiency_coverage.csv"))
    reps = json.load(open(os.path.join(res, "versions.json")))["reps"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 62 * MM), gridspec_kw={"width_ratios": [1, 1.35]})
    t = s[s.phi > 0]
    lim = [0.9, 380]
    a.plot(lim, lim, color=INK2, lw=0.7, ls="--", zorder=1)
    for key, col, sd, dx in (("mdcheck", "g_mdcheck_mean", "g_mdcheck_sd", 0.96),
                             ("pymbar", "g_pymbar_mean", "g_pymbar_sd", 1.04)):
        c, m, lab = STYLE[key]
        a.errorbar(t.g_true * dx, t[col], yerr=t[sd], fmt=m, color=c, ms=4, elinewidth=0.8,
                   capsize=1.5, label=lab.replace(" $g$", ""), zorder=3)
    a.set(xscale="log", yscale="log", xlim=lim, ylim=lim,
          xlabel="exact $g=(1+\\phi)/(1-\\phi)$", ylabel="estimated $g$ (mean $\\pm$ SD)")
    a.legend(loc="upper left")
    panel(a, "a")
    x = np.arange(len(s))
    mc_band(b, reps)
    for key, col in (("naive", "cov_naive"), ("bootstrap", "cov_mdcheck_bootstrap_v100"),
                     ("pyblock", "cov_pyblock"), ("pymbar", "cov_g_pymbar"), ("mdcheck", "cov_mdcheck_ci")):
        c, m, lab = STYLE[key]
        b.plot(x, s[col], marker=m, color=c, label=lab, zorder=3 if key == "mdcheck" else 2)
    b.set_xticks(x, [f"{p:g}\n({g:.0f})" for p, g in zip(s.phi, s.g_true)])
    b.set(xlabel="AR(1) coefficient $\\phi$ (exact $g$)", ylabel="coverage of nominal 95% CI", ylim=(0, 1.02))
    b.legend(loc="lower left", ncol=1)
    panel(b, "b")
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig2_ar1", out)


def fig_length_and_trace(res, out):
    s = pd.read_csv(os.path.join(res, "ar1_length_coverage.csv"))
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 62 * MM), gridspec_kw={"width_ratios": [1, 1.2]})
    mc_band(a, 500)
    markers = {0.5: "o", 0.9: "s", 0.95: "D"}
    for phi, grp in s.groupby("phi"):
        grp = grp.sort_values("neff_true")
        a.plot(grp.neff_true, grp.t_c6, marker=markers[phi], color=STYLE["mdcheck"][0],
               label=f"MDCheck, $\\phi$ = {phi:g}", lw=1.0)
        a.plot(grp.neff_true, grp.naive, marker=markers[phi], color=STYLE["naive"][0],
               label=f"naive, $\\phi$ = {phi:g}", lw=0.8, mfc="none")
    a.axvline(50, color=INK2, lw=0.6, ls=":")
    a.text(55, 0.05, "$N_{\\mathrm{eff}}=50$", fontsize=6.5, color=INK2)
    a.set(xscale="log", xlabel="true effective sample size $N/g$", ylabel="coverage of nominal 95% CI",
          ylim=(0, 1.02))
    a.legend(loc="center left", bbox_to_anchor=(0.07, 0.55), ncol=3, fontsize=5.6, columnspacing=0.8,
             handlelength=1.6)
    panel(a, "a")
    from mdcheck.core.equilibration import detect_equilibration
    from pymbar import timeseries
    d = np.load(os.path.join(res, "md_lj_series.npz"))
    x = d["short_potential_energy_kJmol"][0]
    dt = float(d["sample_interval_ps"])
    t = np.arange(len(x)) * dt
    md = detect_equilibration(x)
    t0_pm, _, _ = timeseries.detect_equilibration(x)
    b.plot(t, x / 864.0, color=INK2, lw=0.6)
    b.axvline(md["t_eq_index"] * dt, color=STYLE["mdcheck"][0], lw=1.2,
              label=f"MDCheck $t_{{\\mathrm{{eq}}}}$ = {md['t_eq_index'] * dt:.1f} ps")
    b.axvline(int(t0_pm) * dt, color=STYLE["pymbar"][0], lw=1.2, ls=(0, (3, 2)),
              label=f"pymbar $t_{{\\mathrm{{eq}}}}$ = {int(t0_pm) * dt:.1f} ps")
    b.set(xlabel="time (ps)", ylabel="$U$ per atom (kJ mol$^{-1}$)")
    b.legend(loc="lower right")
    panel(b, "b")
    fig.tight_layout(w_pad=2.0)
    save(fig, "fig3_length_trace", out)


def fig_md(out):
    frames = []
    for sub in ("results", "results_water300"):
        p = os.path.join(HERE, sub, "md_final_summary.csv")
        if os.path.exists(p):
            part = pd.read_csv(p)
            part["is300"] = sub.endswith("300")
            frames.append(part)
    d = pd.concat(frames, ignore_index=True)
    names = {"potential_energy_kJmol": "$U$", "density_gcm3": r"$\rho$"}
    labels = []
    for r in d.itertuples():
        sysname = "LJ" if r.system == "lj" else "TIP3P"
        length = "150 ps" if r.system == "lj" else ("300 ps" if r.is300 else "100 ps")
        labels.append(f"{sysname} {names[r.observable]}\n{length}")
    rt = pd.read_csv(os.path.join(HERE, "results", "replica_tests_summary.csv"))

    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 66 * MM), gridspec_kw={"width_ratios": [1.25, 1]})
    x = np.arange(len(d))
    w = 0.27
    for off, key, col, lab in ((-w, "naive", "coverage_naive", "naive SE"),
                               (0, "mdcheck", "coverage_mdcheck", "MDCheck"),
                               (w, "pyblock", "coverage_extra_discard", "MDCheck + 20 ps discard")):
        a.bar(x + off, d[col], w * 0.92, color=STYLE[key][0], label=lab, zorder=2)
    for xi, n in zip(x, d.n_replicas):
        band = 1.96 * np.sqrt(0.95 * 0.05 / n)
        a.plot([xi - 1.5 * w, xi + 1.5 * w], [0.95 - band] * 2, color=INK2, lw=0.6, ls=":")
    a.axhline(0.95, color=INK2, lw=0.7, ls="--")
    a.set_xticks(x, labels)
    a.set(ylabel="coverage of reference mean", ylim=(0, 1.05))
    a.grid(axis="x", visible=False)
    a.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.13))
    panel(a, "a")
    from matplotlib.lines import Line2D
    tests = (("mdcheck", "q_test", "Cochran $Q$ (MDCheck)"), ("rhat101", "rhat_gt_1_01", STYLE["rhat101"][2]),
             ("rhat11", "rhat_gt_1_1", STYLE["rhat11"][2]))
    for ds, ls in (("LJ U (150 ps)", "-"), ("TIP3P U (300 ps)", "--")):
        g = rt[rt.dataset == ds].sort_values("shift_se")
        for key, col, _ in tests:
            c, m, _ = STYLE[key]
            b.plot(g.shift_se, g[col], color=c, marker=m, ls=ls, lw=1.0, ms=3.5)
    b.axhline(0.05, color=INK2, lw=0.6, ls=":")
    b.text(6.1, 0.05, "5%", fontsize=6, color=INK2, va="center")
    handles = [Line2D([], [], color=STYLE[k][0], marker=STYLE[k][1], lw=1.0, ms=3.5, label=lab) for k, _, lab in tests]
    handles += [Line2D([], [], color=INK2, ls="-", lw=1.0, label="LJ $U$, 150 ps"),
                Line2D([], [], color=INK2, ls="--", lw=1.0, label="TIP3P $U$, 300 ps")]
    b.legend(handles=handles, loc="upper left", ncol=1, fontsize=6, handlelength=2.2,
             bbox_to_anchor=(0.0, 1.0))
    b.set(xlabel="shift of one replica (standard errors)", ylabel="alarm rate", ylim=(-0.02, 1.02),
          xlim=(-0.3, 6.5))
    panel(b, "b")
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig4_md", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "figures"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for f in os.listdir(args.out):
        os.remove(os.path.join(args.out, f))
    setup()
    res = os.path.join(HERE, "results")
    fig_ar1(res, args.out)
    fig_length_and_trace(res, args.out)
    fig_md(args.out)
    print("figures written to", args.out)


if __name__ == "__main__":
    main()
