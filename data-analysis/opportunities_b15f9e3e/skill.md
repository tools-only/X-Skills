# Opportunities Phase

Use when discovering and ranking opportunities tied to outcomes.

## Inputs
- Selected outcome
- Target users / JTBD
- Evidence of frictions and alternatives
- Constraints

## Process
1) Write 5-10 opportunity statements.
2) Score and rank (impact, urgency, underservedness, feasibility, WTP).
3) Save top opportunities to the OST graph, linked to the outcome.

## CLI Actions
- Add opportunity: `uv run python scripts/ost.py opportunity add --outcome <id> "<title>" --data <json>`
- List opportunities: `uv run python scripts/ost.py opportunity list --outcome <id>`

## Output
- Ranked opportunity nodes linked to the outcome.
- Next step: solutions.
