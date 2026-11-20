# Outputs

`mdcheck assess` and `mdcheck demo` write these files to the output directory:

| File | Content |
|---|---|
| `report.html` | Self-contained report with every check, its status and diagnostics |
| `mdcheck_convergence_overview.svg`, `.pdf`, `.png` | Each observable with the equilibration cut-off and the production region |
| `mdcheck_block_averaging.svg`, `.pdf`, `.png` | Standard error as a function of block size |
| `mdcheck_replica_distributions.svg`, `.pdf`, `.png` | Distributions of each replica (when replicas are given) |
| `mdcheck_summary_table.csv`, `.tex` | Summary table, ready for a spreadsheet or a LaTeX manuscript |
| `methods_snippet.txt` | Draft methods paragraph describing the checks and their results |
| `citation.bib` | BibTeX entry for MDCheck |

The methods paragraph is a starting point: edit it so that it describes your study accurately.
