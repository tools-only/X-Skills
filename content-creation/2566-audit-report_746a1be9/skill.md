# Quiz Question Audit Report

**Generated**: 2026-02-04

---

## Executive Summary

**Total Questions Reviewed**: 256
**Pass**: 231 (90.2%)
**Issues Found**: 25 (9.8%)

### Issue Breakdown

- **Critical**: 6 (wrong answer, major factual error)
- **Warning**: 16 (ambiguous, outdated, misleading)
- **Info**: 3 (minor wording, trivial)

---

## Critical Issues (Immediate Fix Required)

### Q01-001

**Type**: CORRECT_ANSWER
**Issue**: Guide shows npm as Universal Method, not curl

### Q03-011

**Type**: CORRECT_ANSWER
**Issue**: Guide states CLAUDE.md in .claude/ should be committed (project memory), not gitignored

### Q08-019

**Type**: CORRECT_ANSWER
**Issue**: The explanation states "There is no configurable 'auto:N' parameter" yet the question claims auto:N controls lazy loading. The guide (architecture.md:996) shows ENABLE_TOOL_SEARCH=auto:N sets thresholds (5%/10%/20% context), NOT max tools. The 10,000 token threshold is automatic. The question confuses threshold configuration with tool count.

### Q09-003

**Type**: CORRECT_ANSWER
**Issue**: Guide shows no CLI flag for headless mode, only mentions "Headless Mode" as section title with no content

### Q09-029

**Type**: CORRECT_ANSWER
**Issue**: Guide quote says "I treat Claude.md as compounding memory: every mistake becomes a durable rule for the team" (line 12856) - the 4-step cycle explanation is NOT from Boris Cherny, it's an interpretation. The correct answer explains a process not explicitly stated in the guide.

### Q12-012

**Type**: CORRECT_ANSWER
**Issue**: Guide shows 3 sub-agent types (Explore, Plan, general-purpose) but option b incorrectly lists 4 types including Bash as a sub-agent type. Bash is a TOOL not a sub-agent type.

---

## Warnings (Review & Consider Fixing)

### AMBIGUITY (5 questions)

- **Q01-014**: Guide context doesn't clearly list what's preserved vs not preserved
- **Q02-007**: Guide context shows generic section header, not specific content about context poisoning
- **Q02-015**: Guide context points to Fresh Context Pattern section, not XML prompts usage
- **Q04-011**: The guide context (line 17354) is irrelevant to multi-agent orchestration pattern; the correct context appears at lines 5280-5310
- **Q09-005**: "Rev the Engine" describes multiple rounds of PLANNING (think → plan → think harder → refine), not write-critique-improve cycles as stated in explanation

### CORRECT_ANSWER (2 questions)

- **Q10-001**: Shift+Tab does NOT toggle plan/execute; it cycles through permission modes (default→auto→plan). Use /plan to explicitly enter Plan Mode.
- **Q14-011**: Both interpretations are technically valid in guide context (ai-traceability.md lines 114-138), but the nuance is slightly different. Guide shows Assisted-by is for when you're the primary author with AI help (LLVM standard), Co-Authored-By is Claude's default (shared authorship). The answer is correct in principle but could be more precise.

### FACTUAL_ACCURACY (7 questions)

- **Q02-018**: Explanation says "76% fewer tokens with better results" but this specific metric is not in guide context provided
- **Q03-018**: Explanation mixes guide's 8 domains with Boris Cherny's 4 methods without clear distinction; potential confusion
- **Q04-018**: Stats cited as "53-79%" but guide (line 6259) shows "~56%" auto-invocation rate from Gao 2026; also "100% reliable" for CLAUDE.md is overstated (no source confirms 100%)
- **Q06-003**: Guide shows both $ARGUMENTS[0] and $0 as valid syntax, explanation incomplete
- **Q09-006**: Guide context excerpt shows generic "## Output Format" header (line 4849) unrelated to CLI flags; actual flag exists but wrong context provided
- **Q09-026**: Guide says ">10 occurrences = established" (line 5227) not ">10 occurrences", threshold should be "10+" or "≥10"
- **Q10-014**: Guide context shows "nano ~/.claude.json" (line 17354) which is NOT about .gitignore patterns. Correct info is in the explanation but context snippet is wrong file location.

### OUTDATED (2 questions)

- **Q10-004**: Guide shows 75-90% for /compact (line 1449: "🔴 Red | 75-90% | Use /compact or /clear"). Explanation says 70-90% which conflicts. Threshold updated from 70% to 75% in recent versions.
- **Q15-011**: Bridge script exists in examples/scripts/bridge.py, not unresolved. Guide context was incorrectly marked as unresolved.

---

## Info (Minor Issues)

- **Q09-028** (FACTUAL_ACCURACY): Guide references Osmani's article which mentions "comprehension debt" but doesn't explicitly define it as "Code you shipped but don't fully understand" - explanation is correct but attribution could be clearer
- **Q10-002** (FACTUAL_ACCURACY): Explanation correct but context snippet does not show the Esc×2 shortcut (line 15862 only shows "Esc: Dismiss current suggestion", not Esc×2). Guide context incomplete.
- **Q10-006** (TRIVIAL): Question shows answer in option text: "Bash(git *)" is the only option with wildcard syntax matching the question "allow ALL git commands".

---

## Health by Category

| Category | Pass | Issues | Pass Rate |
|----------|------|--------|-----------|
| Category Q01 | 16 | 2 | 88.9% |
| Category Q02 | 15 | 3 | 83.3% |
| Category Q03 | 17 | 2 | 89.5% |
| Category Q04 | 16 | 2 | 88.9% |
| Category Q05 | 18 | 0 | 100.0% |
| Category Q06 | 11 | 1 | 91.7% |
| Category Q07 | 16 | 0 | 100.0% |
| Category Q08 | 19 | 1 | 95.0% |
| Category Q09 | 23 | 6 | 79.3% |
| Category Q10 | 15 | 5 | 75.0% |
| Category Q11 | 17 | 0 | 100.0% |
| Category Q12 | 14 | 1 | 93.3% |
| Category Q13 | 12 | 0 | 100.0% |
| Category Q14 | 10 | 1 | 90.9% |
| Category Q15 | 12 | 1 | 92.3% |

---

## Recommended Actions

1. **Fix Critical Issues** (Priority 1)
   - Review each critical issue
   - Fix question/answer or update explanation
   - Rebuild: `python3 scripts/build-questions.py`

2. **Review Warnings** (Priority 2)
   - Evaluate ambiguities and outdated info
   - Decide: fix, clarify, or accept

3. **Consider Info Issues** (Priority 3)
   - Minor improvements for quality
