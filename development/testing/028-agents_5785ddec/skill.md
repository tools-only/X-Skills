# Repository Guidelines

## Project Structure & Module Organization
Source lives under `feishu_docx/`, with subpackages for `auth/` (OAuth), `core/` (API + parsing), `schema/` (Pydantic models), `cli/` (Typer commands), `tui/` (Textual UI), and `utils/` (config/temp helpers). Tests are in `tests/` with `test_*.py` files and sample outputs under `tests/output/`. User-facing docs and assets live in `docs/` (e.g., screenshots). The CLI entry point is exposed via the `feishu-docx` script, with `main.py` as a local runner.

## Build, Test, and Development Commands
Use editable installs for development and keep Python 3.11+.
```bash
pip install -e ".[dev]"   # install dev deps (pytest/ruff/mypy)
pytest tests/ -v          # run the test suite
feishu-docx --help        # list CLI commands
```
Packaging uses Hatchling; build artifacts with `python -m build` (see `pyproject.toml`). The Makefile contains a release upload target (`make upload`) for maintainers.

## Coding Style & Naming Conventions
Follow standard Python conventions: 4-space indentation, `snake_case` for functions/variables, `PascalCase` for classes, and type hints where practical. Static checks use Ruff and MyPy. Ruff is configured for `line-length = 120` and ignores `E501`; prefer Ruff fixes over manual formatting. Keep modules focused on a single responsibility (e.g., CLI wiring stays in `feishu_docx/cli`).

## Testing Guidelines
Tests run with `pytest`. Name new tests as `tests/test_<area>.py` and keep fixtures close to the modules they validate. When adding new CLI behavior, include at least one test in `tests/test_cli.py`.

## Commit & Pull Request Guidelines
Recent history uses short, imperative subjects with conventional prefixes like `fix:` or emoji tags such as `:sparkles:` and `:recycle:`. Keep titles concise and scoped to one change. PRs should include a clear description, test results (commands run), and links to related issues or tickets. Add screenshots for TUI/CLI output changes when relevant.

## Configuration & Auth Notes
Authentication is handled via `feishu-docx config set` and `feishu-docx auth` (OAuth flow). Avoid committing tokens or generated output; keep samples in `docs/` or `tests/output/`.
