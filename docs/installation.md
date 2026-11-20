# Installation

MDCheck requires Python 3.9 or newer. Its dependencies are NumPy, SciPy, pandas, Matplotlib and Jinja2.

## From PyPI

```bash
pip install mdcheck
mdcheck --version
```

## From source

```bash
git clone https://github.com/sircalch/mdcheck.git
cd mdcheck
pip install -e ".[dev]"
pytest
```

## Building this documentation locally

```bash
pip install -e ".[docs]"
mkdocs serve
```

Then open <http://127.0.0.1:8000>.
