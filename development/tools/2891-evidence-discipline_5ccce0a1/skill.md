# Evidence Discipline

Primary-source verification. Training data recall is explicitly rejected as evidence.

---

## fact-check

- **Evidence rules**: WebFetch, WebSearch, CLI output, repo source, MCP tools = valid. Training data recall = invalid.
- **Verdicts**: VERIFIED / REFUTED / INCONCLUSIVE
- **Integration**: REFUTED → RT-ICA treats as MISSING

---

## find-cause

Evidence-chain protocol for investigations:

- **Valid**: Command output, file content, direct observation
- **Invalid**: Docs intent, training recall, inference from absence

Steps: Disambiguate → Reproduce → Read source → Build evidence chain → Present findings.

---

## scientific-thinking

Hypothesis-driven reasoning:

- Observation (factual only) → Hypothesis (H₀, Hₐ) → Prediction → Experiment Plan → Execute → Conclusion

Use for: bugs, architecture, refactor, strange behavior. Not for trivial tasks.

---

## Sources

- [fact-check SKILL.md](../../.claude/skills/fact-check/SKILL.md)
- [find-cause SKILL.md](../../.claude/skills/find-cause/SKILL.md)
- [scientific-thinking SKILL.md](../../.claude/skills/scientific-thinking/SKILL.md)
