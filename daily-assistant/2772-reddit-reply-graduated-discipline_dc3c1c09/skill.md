# Graduated Discipline

Yes — corrections are logged, but per-mission, not in a persistent cross-mission store.

The escalation ladder has three levels:

1. **Signal** — first occurrence, admiral references the relevant standing order in a coordination message.
2. **Standing Order Remedy** — repeated or moderate impact. Apply the formal remedy and log it in the quarterdeck report.
3. **Damage Control** — mission-threatening. Invoke a full procedure (replace agent, partial rollback, abort, or escalate to human).

Levels aren't skipped unless the issue is immediately critical.

At mission end the Captain's Log captures decisions, failure modes, and reusable patterns. The partial rollback procedure also feeds the failure mode back as a constraint when re-executing the task. So the tuning signal exists — it just lives in structured logs rather than an automated feedback loop across missions. A natural next step would be tracking correction frequency per anti-pattern and using that to adjust defaults, but that's not built yet.
