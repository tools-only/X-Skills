# Outcomes Phase

Use when defining or refining outcomes for an OST.

## Inputs
- Business goal / north star
- Target segment or market
- Baseline metrics
- Time horizon
- Constraints

## Process
1) Define outcome statements (actor + behavior change + context + metric).
2) Ladder outcomes (top outcome → supporting outcomes).
3) Add metrics (baseline, target, time window).
4) Save outcomes to the OST graph.

## CLI Actions
- List outcomes: `uv run python scripts/ost.py outcome list --workspace <name>`
- Add outcome: `uv run python scripts/ost.py outcome add --workspace <name> "<title>" --data <json>`
- Update outcome: `uv run python scripts/ost.py outcome update --id <id> --data <json>`

## Output
- Updated outcome nodes with metrics and ladder links.
- Next step: opportunities.
