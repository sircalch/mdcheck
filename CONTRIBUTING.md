# Contributing to MDCheck

Thank you for your interest in MDCheck (convergence and reproducibility checks for molecular-dynamics simulations). Contributions are welcome: bug reports, documentation fixes, new tests, and code.

## Before you start
- **Open an issue first** for anything beyond a small fix, so we can agree on scope. Describe the problem, what you expected, and a minimal example (ideally with open data).
- Please be clear about what a change claims. MDCheck is a research software tool; it does not replace independent validation or expert judgement, and contributions should not describe it as if it did.

## Development setup
```bash
git clone https://github.com/sircalch/mdcheck.git
cd mdcheck
python -m pip install --upgrade pip
pip install -e .[dev]
pytest tests/ --cov=mdcheck --cov-report=term-missing
```
The test suite runs in CI on Linux, Windows and macOS with several Python versions (see `.github/workflows/test.yml`).

## What a good pull request includes
1. A focused change with a clear description of the problem it solves.
2. Tests for new behaviour and for boundary cases (the existing tests in `tests/` show the style).
3. Updated documentation or README text when behaviour changes.
4. A line in `CHANGELOG.md` under an "Unreleased" heading for user-visible changes.
5. Code formatted with `black` (included in the `dev` extras).

## Reporting bugs and asking questions
Use the issue templates. For security problems, do **not** open a public issue; see `SECURITY.md`.

## Support and governance

- **Getting help.** Ask questions and report bugs in the [issue tracker](https://github.com/sircalch/mdcheck/issues). Include the MDCheck version (`mdcheck --version`), the command or code you ran, the full error message and, when possible, a small input file that reproduces the problem.
- **What to expect.** MDCheck is maintained by one person alongside research and teaching. Issues and pull requests are answered on a best-effort basis; the aim is to acknowledge new reports within two weeks. Security problems follow `SECURITY.md`.
- **How decisions are made.** The maintainer, Andrés Monreal Hernández ([@sircalch](https://github.com/sircalch)), reviews and merges changes. Changes that alter a scientific result, a default threshold or the meaning of a status are discussed in an issue first and are documented in `CHANGELOG.md` and in the [documentation](https://sircalch.github.io/mdcheck/).
- **Becoming a maintainer.** Regular contributors who have had several pull requests merged may be invited to help maintain the project.

## Citing
If you use MDCheck in your work, please cite it using `CITATION.cff` (GitHub shows a "Cite this repository" button).

## Code of conduct
By participating you agree to follow `CODE_OF_CONDUCT.md`.
