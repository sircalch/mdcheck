# Worked example: three replicas of the radius of gyration

This example uses `sample_data/multi_replica_gyrate.csv`, a file with a time column (ns) and the radius of gyration of three replicas of 1,500 frames each. Run it from the repository root:

```bash
python examples/replicas_to_report.py
```

## 1. Reading the file

```text
Replicas: ['Replica_1_Rg', 'Replica_2_Rg', 'Replica_3_Rg']  (1500 frames each)
```

The first column is taken as time and every other column as one observable. Here each column is one replica of the same observable.

## 2. One trajectory

The first replica is assessed as the primary trajectory:

```text
Radius of gyration, replica 1: t_eq = 6.35 ns, N_eff = 178, mean = 1.451 [1.448, 1.453], status WARNING
```

- The first 6.35 ns are discarded as equilibration.
- Of the 1,500 frames, the production region carries the information of about 178 independent samples. That is above 100, so the observable is not a `FAIL`, but below 500, so it is flagged `WARNING` for **marginal sampling**.
- The mean radius of gyration is 1.451 with a 95 % confidence interval of [1.448, 1.453] that already accounts for the autocorrelation.

## 3. Comparing the three replicas

```text
Replica comparison: status PASS, mean JSD = 0.114
Overall status: WARNING
```

The three replica means (1.4536, 1.4521, 1.4533) are statistically consistent: Cochran's \(Q = 0.13\), \(p = 0.94\), \(I^2 = 0\). The mean Jensen–Shannon distance of 0.114 is reported as `LOW`.

The overall status is `WARNING` because of the marginal sampling of the primary trajectory, not because of the replicas. The practical reading is: the replicas agree, and the average is probably right, but a longer run would give a tighter and better-supported estimate.

## 4. The report

The script writes the files described in [Outputs](../outputs.md). In `report.html` each observable has its diagnostic messages, the convergence overview, the block-averaging curve and the replica distributions. `methods_snippet.txt` contains a draft methods paragraph; edit it so that it describes your study.
