# Quality Standards & Validation Checklist

Load this after Phase 4 improvements are applied. Use the quality standards to verify
the overall repair meets the bar, then run the checklist item by item.

---

## Quality Standards

A fully improved skill satisfies all of the following:

**Anatomy:**
- Directory structure matches complexity tier (see `skill-anatomy.md`)
- Every resource file in `scripts/`, `references/`, `examples/` is referenced in SKILL.md
- No extraneous documentation files in the skill directory

**Format Economy:**
- Simple instruction → single imperative, no surrounding prose
- Repeated operation → script in `scripts/`, not inlined code block
- Invocation-selective detail → deferred to `references/`

**Intensional Instruction:**
- Every rule states the *why* alongside the *what*
- Examples confirm stated principles; they do not carry the instructional weight alone
- Degrees of freedom match task fragility (see `skill-anatomy.md`)

**Balance Flexibility with Precision:**
- Agentic steps are loose enough for judgment with explicit outcome criteria
- Deterministic steps are scripted, not reproduced by the model each time

**Remove ruthlessly:** Filler phrases, headers that restate their content, hedging
language, routing guidance in the body, extraneous documentation files.

---

## Validation Checklist

Run after all improvements are applied:

**Structure:**
- [ ] SKILL.md exists with valid YAML frontmatter
- [ ] Frontmatter has `name` and `description` fields
- [ ] Markdown body is present and substantial
- [ ] Directory structure matches skill complexity tier
- [ ] Every resource file is referenced from SKILL.md

**Frontmatter Quality:**
- [ ] Description uses third-person ("This skill should be used when...")
- [ ] 3–5+ varied trigger phrases; includes naive user phrasing
- [ ] Uses `>` scalar, not `|`
- [ ] `argument-hint` present if skill reads `$ARGUMENTS`/`$1`

**Content Quality:**
- [ ] Body uses imperative voice; no first-person, no second-person
- [ ] No "When to Use This Skill" section in the body
- [ ] No headers deeper than H3
- [ ] No extraneous files (`README.md`, `CHANGELOG.md`, etc.)
- [ ] Instructions are intensional (rule + reasoning), not purely extensional
- [ ] Agentic steps have explicit outcome criteria
- [ ] Deterministic operations are scripted, not inlined

**Progressive Disclosure:**
- [ ] SKILL.md under 500 lines; invocation-selective detail in `references/`
- [ ] `scripts/` contains deterministic operations for low-freedom steps
- [ ] `examples/` exists if skill produces user-adaptable output
- [ ] `references/` defers topic-specific detail not needed every invocation

**Tool Selection:**
- [ ] `AskUserQuestion` present if skill needs user decisions mid-workflow
- [ ] `Skill` present if skill invokes other skills
- [ ] `Bash` scoped or absent; unrestricted `Bash` is flagged
- [ ] No dead tool entries (tools listed but never used)

**Script Opportunities:**
- [ ] No code blocks that would be re-generated identically across invocations
- [ ] Scripts in `scripts/` are referenced with trigger condition and invocation
- [ ] No vague script references ("run if needed" without specifying when/how)
- [ ] Deterministic steps with consistency requirements are scripted, not LLM-generated
