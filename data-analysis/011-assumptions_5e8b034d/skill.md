# Assumptions and Experiments Phase

Use when validating solution concepts with the riskiest assumptions first.

## Inputs
- Selected solution
- Evidence so far
- Constraints (time, budget, ethics)

## Process
1) List assumptions (desirability, usability, feasibility, viability, risk).
2) Prioritize by impact x uncertainty.
3) Design smallest experiments with decision thresholds.
4) Save assumptions and experiment plans to the OST graph linked to the solution.

## CLI Actions
- Add assumption: `uv run python scripts/ost.py assumption add --solution <id> "<title>" --data <json>`
- List assumptions: `uv run python scripts/ost.py assumption list --solution <id>`

## Output
- Assumption nodes linked to the solution.
- Next step: decide whether to proceed, pivot, or stop.
