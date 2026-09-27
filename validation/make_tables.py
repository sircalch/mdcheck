"""
Writes the LaTeX tables of the MDCheck validation manuscript from validation/results*/.

    python validation/make_tables.py [--out validation/tables]
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(__file__)


def f(v, nd=2):
    return f"{v:.{nd}f}"


def table_ar1(res, out):
    s = pd.read_csv(os.path.join(res, "ar1_inefficiency_coverage.csv"))
    meta = json.load(open(os.path.join(res, "versions.json")))
    lines = [
        r"\begin{tabular}{rrrrrrrr}", r"\toprule",
        r"$\phi$ & $g$ exact & $g$ MDCheck & $g$ pymbar & naive & MDCheck 1.0.0 & MDCheck 1.1.0 & pyblock \\",
        r" & & mean (SD) & mean (SD) & \multicolumn{4}{c}{coverage of nominal 95\% CI} \\",
        r"\midrule",
    ]
    for r in s.itertuples():
        lines.append(
            f"{r.phi:g} & {r.g_true:.0f} & {f(r.g_mdcheck_mean)} ({f(r.g_mdcheck_sd)}) & "
            f"{f(r.g_pymbar_mean)} ({f(r.g_pymbar_sd)}) & {f(r.cov_naive)} & {f(r.cov_mdcheck_bootstrap_v100)} & "
            f"{f(r.cov_mdcheck_ci)} & {f(r.cov_pyblock)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(out, "table_ar1.tex"), "w").write("\n".join(lines) + "\n")
    return meta


def table_eq(res, out):
    s = pd.read_csv(os.path.join(res, "equilibration_summary.csv"))
    lines = [
        r"\begin{tabular}{rrrrrrr}", r"\toprule",
        r"$\phi$ & $\tau_{\mathrm{relax}}$ & $t_{\mathrm{eq}}$ MDCheck & $t_{\mathrm{eq}}$ pymbar & "
        r"$|\Delta t_{\mathrm{eq}}|$ & \multicolumn{2}{c}{bias of mean} \\",
        r" & (steps) & (median) & (median) & (median) & no discard & MDCheck \\",
        r"\midrule",
    ]
    for r in s.itertuples():
        lines.append(f"{r.phi:g} & {r.tau_relax:.0f} & {r.t_eq_mdcheck:.0f} & {r.t_eq_pymbar:.0f} & "
                     f"{r.median_abs_diff_t_eq:.0f} & {f(r.bias_no_discard, 3)} & {f(r.bias_mdcheck, 3)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(out, "table_equilibration.tex"), "w").write("\n".join(lines) + "\n")


def table_md(out, sources):
    frames = []
    for label, res in sources:
        p = os.path.join(res, "md_final_summary.csv")
        if os.path.exists(p):
            d = pd.read_csv(p)
            d["label"] = d.apply(lambda r: label(r), axis=1)
            frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    lines = [
        r"\begin{tabular}{lrrrrrrrrrr}", r"\toprule",
        r"System & $n$ & $t_{\mathrm{eq}}$ (ps) & $g$ & $N_{\mathrm{eff}}$ & "
        r"\multicolumn{3}{c}{coverage} & no-discard & \multicolumn{2}{c}{replica alarm rate} \\",
        r"\cmidrule(lr){6-8}\cmidrule(lr){10-11}",
        r" & & MDCheck / pymbar & & & naive & MDCheck & extra discard & bias (SE) & $Q$ test & JSD $>0.15$ \\",
        r"\midrule",
    ]
    for r in d.itertuples():
        lines.append(
            f"{r.label} & {r.n_replicas} & {r.t_eq_ps_mdcheck:.1f} / {r.t_eq_ps_pymbar:.1f} & {r.g_median:.1f} & "
            f"{r.n_eff_median:.0f} & {f(r.coverage_naive)} & {f(r.coverage_mdcheck)} & {f(r.coverage_extra_discard)} & "
            f"{r.bias_no_discard_over_se:+.1f} & {f(r.replica_alarm_rate)} & {f(r.legacy_jsd_alarm_rate)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(out, "table_md.tex"), "w").write("\n".join(lines) + "\n")
    d.to_csv(os.path.join(out, "table_md.csv"), index=False)


def table_estimators(res, out):
    s = pd.read_csv(os.path.join(res, "estimator_comparison.csv"))
    cov = s.pivot_table(index=["data", "observable"], columns="estimator", values="coverage")
    g = s.pivot_table(index=["data", "observable"], columns="estimator", values="g_median")
    names = {"ar1_phi0.5": r"AR(1) $\phi=0.5$", "ar1_phi0.9": r"AR(1) $\phi=0.9$",
             "ar1_phi0.99": r"AR(1) $\phi=0.99$", "lj": "LJ", "water": "TIP3P"}
    obs = {"x": "", "potential_energy_kJmol": " $U$", "density_gcm3": r" $\rho$"}
    cols = ["sokal_c6", "first_negative", "geyer_imse", "pymbar"]
    lines = [r"\begin{tabular}{lrrrr}", r"\toprule",
             r"Data & Sokal ($c=6$) & first negative & Geyer IMSE & pymbar \\", r"\midrule"]
    for idx in cov.index:
        lines.append(f"{names.get(idx[0], idx[0])}{obs.get(idx[1], idx[1])} & " +
                     " & ".join(f"{cov.loc[idx, c]:.3f} ({g.loc[idx, c]:.1f})" for c in cols) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(out, "table_estimators.tex"), "w").write("\n".join(lines) + "\n")


def table_replica_tests(res, out):
    s = pd.read_csv(os.path.join(res, "replica_tests_summary.csv"))
    order = [d for d in ["LJ U (150 ps)", "TIP3P rho (100 ps)", "TIP3P U (100 ps)", "TIP3P rho (300 ps)",
                         "TIP3P U (300 ps)", "Villin RMSD (2 ns)", "Villin Rg (2 ns)"] if d in set(s.dataset)]
    names = {"LJ U (150 ps)": "LJ $U$ (150 ps)", "TIP3P rho (100 ps)": r"TIP3P $\rho$ (100 ps)",
             "TIP3P U (100 ps)": "TIP3P $U$ (100 ps)", "TIP3P rho (300 ps)": r"TIP3P $\rho$ (300 ps)",
             "TIP3P U (300 ps)": "TIP3P $U$ (300 ps)", "Villin RMSD (2 ns)": r"Villin C$\alpha$ RMSD (2 ns)",
             "Villin Rg (2 ns)": r"Villin $R_g$ (2 ns)"}
    lines = [r"\begin{tabular}{lrrrrrrr}", r"\toprule",
             r" & \multicolumn{3}{c}{no injected shift} & \multicolumn{2}{c}{3 SE shift} & \multicolumn{2}{c}{6 SE shift} \\",
             r"\cmidrule(lr){2-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}",
             r"Data & $Q$ test & $\hat{R}>1.01$ & JSD$>0.15$ & $Q$ test & $\hat{R}>1.01$ & $Q$ test & $\hat{R}>1.1$ \\",
             r"\midrule"]
    for ds in order:
        g = s[s.dataset == ds].set_index("shift_se")
        lines.append(f"{names[ds]} & {g.loc[0, 'q_test']:.2f} & {g.loc[0, 'rhat_gt_1_01']:.2f} & "
                     f"{g.loc[0, 'jsd_gt_0_15']:.2f} & {g.loc[3, 'q_test']:.2f} & {g.loc[3, 'rhat_gt_1_01']:.2f} & "
                     f"{g.loc[6, 'q_test']:.2f} & {g.loc[6, 'rhat_gt_1_1']:.2f} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(out, "table_replica_tests.tex"), "w").write("\n".join(lines) + "\n")


def table_length(res, out):
    s = pd.read_csv(os.path.join(res, "ar1_length_coverage.csv"))
    lines = [r"\begin{tabular}{rrrrrrrr}", r"\toprule",
             r"$\phi$ & $N$ & $N/g$ & median $\hat g$ & naive & $z$, $c=6$ & $t$, $c=6$ & $t$, $c=4$ / $c=8$ \\",
             r"\midrule"]
    for r in s.itertuples():
        lines.append(f"{r.phi:g} & {r.n} & {r.neff_true:.0f} & {r.g_c6_median:.1f} & {r.naive:.2f} & {r.z_c6:.2f} & "
                     f"{r.t_c6:.2f} & {r.t_c4:.2f} / {r.t_c8:.2f} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(out, "table_length.tex"), "w").write("\n".join(lines) + "\n")


def table_runtime(res, out):
    p = os.path.join(res, "runtime_benchmark.csv")
    if not os.path.exists(p):
        return
    s = pd.read_csv(p)

    def cell(t, t0):
        return "--" if pd.isna(t) else f"{t:.2f} ({int(t0)})"

    lines = [r"\begin{tabular}{rrrrr}", r"\toprule",
             r"$N$ & MDCheck 1.1 & MDCheck 1.0 grid & pymbar & pymbar \texttt{fast} \\",
             r" & s ($t_{\mathrm{eq}}$) & s ($t_{\mathrm{eq}}$) & s ($t_{\mathrm{eq}}$) & s ($t_{\mathrm{eq}}$) \\",
             r"\midrule"]
    for _, r in s.iterrows():
        lines.append(f"{int(r['n']):,} & {cell(r['mdcheck_s'], r['mdcheck_t0'])} & "
                     f"{cell(r.get('mdcheck10_grid_s', np.nan), r.get('mdcheck10_grid_t0', np.nan))} & "
                     f"{cell(r.get('pymbar_s', np.nan), r.get('pymbar_t0', np.nan))} & "
                     f"{cell(r['pymbar_fast_s'], r['pymbar_fast_t0'])} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(out, "table_runtime.tex"), "w").write("\n".join(lines) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "tables"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    res = os.path.join(HERE, "results")
    table_ar1(res, args.out)
    table_eq(res, args.out)
    table_estimators(res, args.out)
    table_replica_tests(res, args.out)
    table_length(res, args.out)
    table_runtime(res, args.out)
    obs = {"potential_energy_kJmol": "$U$", "density_gcm3": r"$\rho$", "rmsd_ca_nm": r"C$\alpha$ RMSD",
           "rg_protein_nm": "$R_g$"}
    table_md(args.out, [
        (lambda r: f"{'LJ' if r.system == 'lj' else 'TIP3P'} {obs[r.observable]} "
                   f"({'150' if r.system == 'lj' else '100'} ps)", res),
        (lambda r: f"TIP3P {obs[r.observable]} (300 ps)", os.path.join(HERE, "results_water300")),
        (lambda r: f"Villin {obs[r.observable]} (2 ns)", os.path.join(HERE, "results_villin")),
    ])
    print("tables written to", args.out)


if __name__ == "__main__":
    main()
