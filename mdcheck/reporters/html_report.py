"""
Interactive, self-contained HTML report generator for MDCheck.
"""

from typing import Dict, Any, Optional
import os
import jinja2
from mdcheck.core.scoring import SimulationQualityReport

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MDCheck Simulation Quality Report</title>
    <style>
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --pass-color: #10b981;
            --pass-bg: rgba(16, 185, 129, 0.15);
            --warn-color: #f59e0b;
            --warn-bg: rgba(245, 158, 11, 0.15);
            --fail-color: #ef4444;
            --fail-bg: rgba(239, 68, 68, 0.15);
            --accent-blue: #38bdf8;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem 1rem;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 1.5rem;
            margin-bottom: 2rem;
            flex-wrap: wrap;
            gap: 1rem;
        }
        .title-group h1 {
            font-size: 2rem;
            font-weight: 700;
            letter-spacing: -0.025em;
            color: var(--text-primary);
        }
        .title-group p {
            color: var(--text-secondary);
            font-size: 0.95rem;
        }
        .status-badge {
            display: inline-flex;
            align-items: center;
            padding: 0.5rem 1.25rem;
            border-radius: 9999px;
            font-weight: 700;
            font-size: 1.1rem;
            letter-spacing: 0.05em;
            text-transform: uppercase;
        }
        .badge-pass { background-color: var(--pass-bg); color: var(--pass-color); border: 1px solid var(--pass-color); }
        .badge-warning { background-color: var(--warn-bg); color: var(--warn-color); border: 1px solid var(--warn-color); }
        .badge-fail { background-color: var(--fail-bg); color: var(--fail-color); border: 1px solid var(--fail-color); }
        
        .grid-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: 1.25rem;
            margin-bottom: 2rem;
        }
        .card {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.5rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
        }
        .card-label {
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-secondary);
            margin-bottom: 0.5rem;
        }
        .card-value {
            font-size: 1.8rem;
            font-weight: 700;
            color: var(--text-primary);
        }
        .card-subtext {
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-top: 0.25rem;
        }

        .section-title {
            font-size: 1.4rem;
            font-weight: 600;
            margin-bottom: 1rem;
            color: var(--accent-blue);
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 2rem;
            font-size: 0.95rem;
            background-color: var(--card-bg);
            border-radius: 0.75rem;
            overflow: hidden;
            border: 1px solid var(--card-border);
        }
        th, td {
            padding: 0.85rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--card-border);
        }
        th {
            background-color: rgba(255, 255, 255, 0.03);
            color: var(--text-secondary);
            font-weight: 600;
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }
        tr:last-child td { border-bottom: none; }
        tr:hover td { background-color: rgba(255, 255, 255, 0.02); }

        .tag {
            display: inline-block;
            padding: 0.2rem 0.6rem;
            border-radius: 0.375rem;
            font-size: 0.75rem;
            font-weight: 700;
        }
        .tag-pass { background-color: var(--pass-bg); color: var(--pass-color); }
        .tag-warning { background-color: var(--warn-bg); color: var(--warn-color); }
        .tag-fail { background-color: var(--fail-bg); color: var(--fail-color); }

        .box {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.5rem;
            margin-bottom: 2rem;
        }
        pre {
            background-color: rgba(0, 0, 0, 0.4);
            padding: 1rem;
            border-radius: 0.5rem;
            overflow-x: auto;
            color: #38bdf8;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
            font-size: 0.9rem;
            border: 1px solid rgba(255, 255, 255, 0.05);
            white-space: pre-wrap;
        }
        .btn-copy {
            background-color: #2563eb;
            color: white;
            border: none;
            padding: 0.4rem 0.8rem;
            border-radius: 0.375rem;
            cursor: pointer;
            font-size: 0.8rem;
            margin-top: 0.5rem;
            transition: background 0.2s;
        }
        .btn-copy:hover { background-color: #1d4ed8; }

        .figure-container {
            text-align: center;
            margin-bottom: 2rem;
        }
        .figure-container img {
            max-width: 100%;
            height: auto;
            border-radius: 0.5rem;
            border: 1px solid var(--card-border);
        }
        footer {
            text-align: center;
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-top: 3rem;
            padding-top: 1.5rem;
            border-top: 1px solid var(--card-border);
        }
    </style>
</head>
<body>
    <div class="container">
        <header class="header">
            <div class="title-group">
                <h1>MDCheck Trajectory Quality Report</h1>
                <p>Automated Convergence, Equilibration & Statistical Inefficiency Analysis</p>
            </div>
            <div>
                <span class="status-badge badge-{{ report.overall_status.lower() }}">
                    {{ report.overall_status }}
                </span>
            </div>
        </header>

        <div class="grid-cards">
            <div class="card">
                <div class="card-label">Overall Assessment</div>
                <div class="card-value" style="font-size: 1.3rem;">{{ report.score_summary }}</div>
                <div class="card-subtext">{{ report.key_metrics.n_observables_evaluated }} observables evaluated</div>
            </div>
            <div class="card">
                <div class="card-label">Max Equilibration Time (t_eq)</div>
                <div class="card-value">{{ "%.2f"|format(report.key_metrics.max_equilibration_time) }} <span style="font-size: 1rem; color: var(--text-secondary)">ns</span></div>
                <div class="card-subtext">Production region automatically sliced</div>
            </div>
            <div class="card">
                <div class="card-label">Min Effective Samples (N_eff)</div>
                <div class="card-value">{{ "%.0f"|format(report.key_metrics.min_effective_sample_size) }}</div>
                <div class="card-subtext">Independent statistical observations</div>
            </div>
            <div class="card">
                <div class="card-label">Replicas Evaluated</div>
                <div class="card-value">{{ report.key_metrics.n_replicas_evaluated }}</div>
                <div class="card-subtext">
                    {% if report.replica_assessment %}
                    JSD: {{ "%.3f"|format(report.replica_assessment.mean_jsd) }} ({{ report.replica_assessment.status }})
                    {% else %}
                    Single trajectory mode
                    {% endif %}
                </div>
            </div>
        </div>

        <h2 class="section-title">Observable Quality Breakdown</h2>
        <table>
            <thead>
                <tr>
                    <th>Observable</th>
                    <th>Total Frames</th>
                    <th>t_eq (ns)</th>
                    <th>Production Frames</th>
                    <th>&tau;<sub>int</sub> (frames)</th>
                    <th>N<sub>eff</sub></th>
                    <th>Production Mean [95% CI]</th>
                    <th>Drift</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
                {% for name, obs in report.observables.items() %}
                <tr>
                    <td><strong>{{ name }}</strong></td>
                    <td>{{ obs.n_total }}</td>
                    <td>{{ "%.2f"|format(obs.t_eq_time) }} ({{ "%.1f"|format(obs.fraction_discarded * 100) }}%)</td>
                    <td>{{ obs.n_prod }}</td>
                    <td>{{ "%.1f"|format(obs.tau_int) }} (g={{ "%.2f"|format(obs.g_inefficiency) }})</td>
                    <td><strong>{{ "%.0f"|format(obs.n_eff) }}</strong></td>
                    <td>{{ "%.3f"|format(obs.mean_prod) }} [{{ "%.3f"|format(obs.ci_lower_95) }}, {{ "%.3f"|format(obs.ci_upper_95) }}]</td>
                    <td><span class="tag tag-{{ obs.drift_status.lower() }}">{{ obs.drift_status }}</span></td>
                    <td><span class="tag tag-{{ obs.status.lower() }}">{{ obs.status }}</span></td>
                </tr>
                {% endfor %}
            </tbody>
        </table>

        {% if report.recommendations %}
        <div class="box" style="border-left: 4px solid var(--warn-color);">
            <h3 style="color: var(--warn-color); margin-bottom: 0.5rem;">Diagnostic Notes & Recommendations</h3>
            <ul style="padding-left: 1.25rem;">
                {% for rec in report.recommendations %}
                <li style="margin-bottom: 0.25rem; color: var(--text-secondary);">{{ rec }}</li>
                {% endfor %}
            </ul>
        </div>
        {% endif %}

        <h2 class="section-title">Methods Section Snippet (Publication Ready)</h2>
        <div class="box">
            <p style="color: var(--text-secondary); margin-bottom: 0.5rem;">Copy and paste directly into your manuscript's Methods or Computational Details section:</p>
            <pre id="methodsSnippet">{{ methods_text }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('methodsSnippet')">Copy Methods Paragraph</button>
        </div>

        <h2 class="section-title">BibTeX Citation</h2>
        <div class="box">
            <pre id="bibSnippet">{{ citation_bib }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('bibSnippet')">Copy BibTeX</button>
        </div>

        <footer>
            Generated automatically by <strong>MDCheck v1.0.0</strong> &bull; Scientific Simulation Quality & Convergence Toolkit &bull; Monreal-Hernández, 2026.
        </footer>
    </div>

    <script>
        function copyToClipboard(elementId) {
            const text = document.getElementById(elementId).innerText;
            navigator.clipboard.writeText(text).then(() => {
                alert('Copied to clipboard!');
            }).catch(err => {
                console.error('Error copying: ', err);
            });
        }
    </script>
</body>
</html>
"""


def generate_html_report(
    report: SimulationQualityReport,
    output_path: str,
    methods_text: str = "",
    citation_bib: str = ""
) -> str:
    """
    Renders the HTML report template and writes it to disk.

    Parameters
    ----------
    report : SimulationQualityReport
        Quality assessment report.
    output_path : str
        Path to output HTML file.
    methods_text : str
        Generated Methods paragraph.
    citation_bib : str
        BibTeX string.

    Returns
    -------
    output_path : str
        Path to generated HTML file.
    """
    template = jinja2.Template(HTML_TEMPLATE)
    rendered = template.render(
        report=report,
        methods_text=methods_text,
        citation_bib=citation_bib
    )
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered)
        
    return output_path
