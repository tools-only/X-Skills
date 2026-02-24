# Agent Quality Standards & Validation Checklist

Load after Phase 4 improvements are applied. Use quality standards to verify the overall
repair meets the bar, then run the checklist item by item.

---

## Quality Standards

A fully improved agent satisfies all of the following:

**Voice:**
- System prompt is second-person throughout — no first-person, no third-person
- Persona established in the first sentence
- Process steps are numbered, sequential, and imperative

**Description:**
- Starts with "Use this agent when..." (exact routing pattern)
- 2–4 `<example>` blocks, each fully formed
- Commentary explains routing reasoning, not just restates user intent
- Proactive trigger examples present for proactive agents

**Tools and modifiers:**
- `tools` restricted to minimum needed (least-privilege for autonomous execution)
- No unscoped `Bash` unless full shell access is genuinely required
- `color` set and semantically meaningful
- `skills:` preloads listed match what is actually used in the process

**Efficiency:**
- System prompt under 400 lines; domain reference data deferred to `skills:` preloads
- No routing guidance embedded in the system prompt body
- No hedging language; direct imperatives throughout

---

## Validation Checklist

**Description:**
- [ ] Starts with "Use this agent when..."
- [ ] Uses `|` literal scalar (not `>` — XML `<example>` blocks require literal newlines)
- [ ] 2–4 `<example>` blocks present
- [ ] Each example has: `Context:`, `user:`, `assistant:`, `<commentary>`
- [ ] `<commentary>` explains routing reasoning (not just restates the user message)
- [ ] Proactive two-turn assistant pattern included if agent fires after events

**System Prompt:**
- [ ] Written entirely in second person
- [ ] No first-person language ("I will", "I'll", "I am")
- [ ] No bare imperatives without "you" (risk of reading as skill instruction)
- [ ] First sentence is a persona statement ("You are a [role] specializing in [domain]")
- [ ] Process steps are numbered and each has a clear completion criterion
- [ ] Output format section is present and explicit
- [ ] Edge cases section is present

**Frontmatter:**
- [ ] `name` is lowercase-hyphens, 3–50 chars, not a generic term
- [ ] `color` is set
- [ ] `model` is set or intentionally omitted (inherit)
- [ ] `tools` is restricted for analysis-only agents
- [ ] No unscoped `Bash` unless required
- [ ] `skills:` preloads listed are actually used in the process steps

**Efficiency:**
- [ ] System prompt is under 400 lines
- [ ] No embedded reference tables > 100 lines (use `skills:` instead)
- [ ] No "When to trigger this agent" language in the body
- [ ] No hedging ("you might", "generally", "consider possibly")

**Script Opportunities:**
- [ ] No code blocks re-generated identically across invocations
- [ ] Any referenced scripts have trigger condition, exact invocation, and output handling
- [ ] Steps that must produce consistent output are scripted, not left to LLM generation
