# Solutions Phase

Use when ideating solution concepts tied to top opportunities.

## Inputs
- Selected opportunity
- Constraints and priorities
- Differentiation requirements

## Process
1) Generate 3-5 distinct solution concepts.
2) Evaluate trade-offs (pros/cons, feasibility, differentiation).
3) Select a leading concept.
4) Save solutions to the OST graph linked to the opportunity.

## CLI Actions
- Add solution: `uv run python scripts/ost.py solution add --opportunity <id> "<title>" --data <json>`
- List solutions: `uv run python scripts/ost.py solution list --opportunity <id>`

## Output
- Solution nodes linked to the opportunity.
- Next step: assumptions/experiments.
