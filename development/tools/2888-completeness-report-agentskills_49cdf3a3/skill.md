# Skill Completeness Report: agentskills

**Evaluated:** 2026-02-17
**Skill Path:** plugins/plugin-creator/skills/agentskills/

## Overall Score: 50% (12/24)

| Category | Score | Label | Findings |
|----------|-------|-------|----------|
| 1. Preparation | 0 | None | No environment verification, no input inspection, no metadata extraction scripts |
| 2. Progression | 2 | Adequate | Clear structure with field reference table, validation rules, progressive disclosure patterns. No decision tree or deterministic scripts. |
| 3. Verification | 1 | Minimal | Validation CLI documented but not bundled. No automated checks or error-correction loops. |
| 4. Scripts | 0 | None | No scripts provided |
| 5. Examples | 3 | Exemplary | Extensive concrete examples: valid/invalid names, good/bad descriptions, frontmatter with all fields, 3 progressive disclosure patterns, directory structures, validation CLI + Python API |
| 6. Anti-Patterns | 2 | Adequate | Anti-patterns shown with corrections (bad descriptions, invalid names, what NOT to include). Missing explicit side-by-side WRONG/CORRECT formatting and reasoning for why each fails. |
| 7. References | 3 | Exemplary | 3 reference files organized by topic (specification, best-practices, integration), all linked from SKILL.md workflow sections, each with table of contents |
| 8. Assets | 1 | Minimal | No assets/ directory. However, this is a reference/knowledge skill — assets are not strongly expected. External links to example skills repo partially compensate. |

## Category Details

### 1. Preparation (0/3 — None)

**What was evaluated:**
- Environment verification before starting
- Input inspection before acting
- Metadata extraction scripts

**Evidence found:**
- ❌ No environment checks (no verification of `skills-ref` installation, Python version, or filesystem state)
- ❌ No input inspection (no step to analyze the user's target skill directory before advising)
- ❌ No metadata extraction scripts (no script to parse/validate existing SKILL.md frontmatter)

**Recommendation:**
This is a reference/knowledge skill, not an operational skill, so preparation is less critical than for document-processing skills. However, adding a `scripts/validate_frontmatter.py` that validates a SKILL.md against the open standard schema would significantly increase utility. It would let the agent verify a skill's compliance before advising on fixes.

---

### 2. Progression (2/3 — Adequate)

**What was evaluated:**
- Clear sequence of steps
- Fragile operations handled by deterministic scripts
- Working code examples with imports
- WRONG/CORRECT contrast pairs
- Decision tree for multiple paths

**Evidence found:**
- ✅ Clear structure: SKILL.md:16-86 covers format → frontmatter → fields → name rules → description guidelines in logical sequence
- ✅ Working code examples: SKILL.md:199-206 shows Python API with imports (`from skills_ref import validate, read_properties, to_prompt`)
- ✅ Multiple valid/invalid examples for name field at SKILL.md:71-72
- ✅ Good/bad description contrast at SKILL.md:78-84
- ❌ No deterministic scripts for fragile operations (validation is documented but not bundled)
- ❌ No explicit decision tree (e.g., "Creating portable skill? → use only standard fields. Claude Code only? → see claude-skills-overview-2026")

**Recommendation:**
Add a decision tree at the top of SKILL.md that routes the agent based on the task: creating a new portable skill vs. auditing an existing skill vs. understanding the spec for integration. Also consider bundling a lightweight validation script.

---

### 3. Verification (1/3 — Minimal)

**What was evaluated:**
- Explicit verification steps
- Automated verification scripts
- Error-correction loops
- Concrete acceptance criteria

**Evidence found:**
- ✅ Validation CLI documented at SKILL.md:184-208 (`skills-ref validate`, Python API)
- ❌ No bundled verification scripts — relies on external `skills-ref` library
- ❌ No error-correction loop (no "validate → fix → re-validate" workflow)
- ❌ No concrete acceptance criteria for "this skill is spec-compliant"

**Recommendation:**
Add a verification section with a concrete checklist: "A spec-compliant skill MUST have: (1) name field matching directory, (2) description 1-1024 chars, (3) no unknown frontmatter fields, (4) SKILL.md under 500 lines." Bundle a lightweight validation script or at minimum define a manual verification workflow.

---

### 4. Scripts (0/3 — None)

**What was evaluated:**
- Repetitive operations captured in scripts
- Scripts with --help and self-documenting
- Edge case handling
- Tested output

**Evidence found:**
- ❌ No `scripts/` directory exists
- ❌ No validation script bundled
- ❌ No scaffolding script (the skill-creator has `init_skill.py` but this skill doesn't reference or bundle it)

**Recommendation:**
Consider adding:
1. `scripts/validate_agentskills.py` — validates a skill directory against the open standard (name rules, description length, no unknown fields, directory name match)
2. Alternatively, reference the existing `skill-creator` scripts with a clear pointer: "For scaffolding, use `init_skill.py` from the skill-creator skill"

For a pure reference skill, scripts are less critical, but a validation script would be the single highest-value addition.

---

### 5. Examples (3/3 — Exemplary)

**What was evaluated:**
- Working code with imports and realistic data
- Exact input→output pairs
- Common cases covered
- Edge cases demonstrated

**Evidence found:**
- ✅ Complete frontmatter examples: minimal (SKILL.md:30-35), full (SKILL.md:39-50)
- ✅ Field reference table with constraints (SKILL.md:54-61)
- ✅ Valid/invalid name examples with specific failure reasons (SKILL.md:71-72)
- ✅ Good/bad description examples (SKILL.md:78-84)
- ✅ Three progressive disclosure patterns with directory trees (SKILL.md:102-131)
- ✅ Python API with imports (SKILL.md:199-206)
- ✅ CLI examples (SKILL.md:186-195)
- ✅ Portable vs Claude Code comparison table (SKILL.md:216-231)
- ✅ Reference files include additional examples: specification.md has 4 valid + 4 invalid name examples, best-practices.md has template/examples/conditional patterns with code

**Recommendation:**
No improvement needed. Examples are comprehensive and concrete.

---

### 6. Anti-Patterns (2/3 — Adequate)

**What was evaluated:**
- Known failure modes documented
- Bad output shown (not just described)
- Corrections shown side-by-side

**Evidence found:**
- ✅ Invalid name examples with reasons at SKILL.md:72 and specification.md:103-110
- ✅ Bad description example at SKILL.md:83-84
- ✅ "What NOT to Include" section at SKILL.md:173-178 (README.md, CHANGELOG.md, etc.)
- ✅ best-practices.md:453-461 lists anti-patterns (Windows paths, too many options, assuming tools installed, vague descriptions, nested references, time-sensitive logic)
- ❌ Missing explicit ❌/✅ formatting for contrast pairs (the WRONG/CORRECT pattern from Anthropic's DOCX/XLSX skills)
- ❌ Missing reasoning for WHY each anti-pattern fails (e.g., "vague descriptions cause the agent to activate on wrong tasks")

**Recommendation:**
Add explicit ❌ WRONG / ✅ CORRECT contrast pairs for the most common mistakes:
- ❌ `name: My-Skill` → ✅ `name: my-skill` (uppercase breaks cross-agent discovery)
- ❌ Using YAML multiline `>-` → ✅ Single-line quoted string (broken in Claude Code indexer)
- ❌ SKILL.md >800 lines with everything inline → ✅ Split to references/ (context window bloat)

---

### 7. References (3/3 — Exemplary)

**What was evaluated:**
- Reference documentation for APIs, schemas, formats
- Organized by topic
- Linked from workflow steps

**Evidence found:**
- ✅ 3 reference files: specification.md (322 lines), best-practices.md (461 lines), integration.md (131 lines)
- ✅ Each organized by topic with table of contents (specification.md:9-24, best-practices.md:9-18, integration.md:9-17)
- ✅ All linked from SKILL.md:240-242 ("Detailed References" section)
- ✅ Best practices also linked from SKILL.md:162 ("Authoring Best Practices" section)
- ✅ Source attribution on each reference file (specification.md:3, best-practices.md:3, integration.md:3)
- ✅ References are one level deep from SKILL.md (no nesting)

**Recommendation:**
No improvement needed. Reference structure follows the progressive disclosure pattern the skill itself documents.

---

### 8. Assets (1/3 — Minimal)

**What was evaluated:**
- Templates, fonts, images, boilerplate bundled
- Assets the AI uses (not reads into context)

**Evidence found:**
- ❌ No `assets/` directory
- ✅ External links to example skills repo (SKILL.md:248) and reference library (SKILL.md:249) partially compensate
- Mitigating factor: This is a reference/knowledge skill. Unlike document-processing skills (PPTX with fonts, XLSX with templates), a specification reference skill has no natural output assets.

**Recommendation:**
Consider adding an `assets/skill-template/` directory with a minimal spec-compliant SKILL.md template that agents can copy when creating new portable skills. This would be the single asset that makes sense for this skill type. Alternatively, a frontmatter template YAML snippet file.

---

## Recommendations for Improvement

1. **High Priority:** Add a lightweight `scripts/validate_agentskills.py` that validates a skill directory against the open standard schema (name rules, description length, unknown fields, directory name match). This is the single highest-value addition — it turns a passive reference into an active tool. (Preparation +1, Scripts +2, Verification +1 = **+4 points**)

2. **Medium Priority:** Add a decision tree at the top of SKILL.md routing the agent: "Creating portable skill?" / "Auditing existing skill?" / "Understanding spec for integration?" Each path points to the relevant section. (Progression +1 = **+1 point**)

3. **Medium Priority:** Add explicit ❌/✅ contrast pairs for the 3-4 most common spec violations, with reasoning for why each fails. (Anti-Patterns +1 = **+1 point**)

4. **Low Priority:** Add `assets/skill-template/SKILL.md` with a minimal spec-compliant template agents can copy. (Assets +1 = **+1 point**)

**Projected score after improvements:** 19/24 (79%)

## Reference

This audit follows patterns from Anthropic's official skills repository:
- <https://github.com/anthropics/skills>

Checklist: [Skill Completeness Checklist](../../plugins/plugin-creator/skills/audit-skill-completeness/references/skill-completeness-checklist.md)
