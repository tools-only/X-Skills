# OST Entry Flow

Goal: Start a session, pick a workspace, and select or create an outcome.

## Steps

1) Ensure DB exists
- If `.agr/ost.db` is missing, initialize it:
  - `uv run python scripts/ost.py init --path .agr/ost.db`

2) Select workspace
- List workspaces:
  - `uv run python scripts/ost.py workspace list`
- If needed, create one:
  - `uv run python scripts/ost.py workspace create "<name>"`

3) Select outcome
- List outcomes:
  - `uv run python scripts/ost.py outcome list --workspace <name>`
- If needed, create one:
  - `uv run python scripts/ost.py outcome add --workspace <name> "<title>"`

4) Route to phase
- If user wants to define outcomes: `references/outcomes.md`
- If user wants opportunities: `references/opportunities.md`
- If user wants solutions: `references/solutions.md`
- If user wants assumptions/experiments: `references/assumptions.md`
