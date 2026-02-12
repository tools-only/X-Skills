---
name: enforce
description: Deterministic quality enforcement. Researches the repo, creates strict linter configs, fixes ALL violations, enforces structural limits (files <400 LOC, functions <80 LOC) across the ENTIRE codebase. Use when asked to "/lint", "/enforce", "make this repo clean", "set up linting", "enforce quality", or "add strict linting".
---

# Deterministic Quality Enforcement (/lint)

Autonomous quality enforcement. Researches the repo, creates or upgrades strict linter configs IN THE REPO, then iterates until every file passes linting AND structural limits. Unlike `_linter.py` (passive stop-time gate using `--config` flags), this skill **permanently configures the repo** so CI, hooks, and every developer enforce the same rules.

## Progress Output

```
[1/8] Activating autonomous mode
[2/8] Researching repo stack and existing config
[3/8] Creating strict linter configurations
[4/8] Running initial lint pass — capturing baseline
[5/8] Fixing batch N/M — [category] violations
[6/8] Enforcing structural limits (ALL files)
[7/8] Re-running full lint suite — verifying zero errors
[8/8] Writing completion checkpoint
```

Step [5/8] repeats per batch. Steps [5-7] may loop on regressions. Always emit [1/8] and [8/8].

## Architecture

```
RESEARCH repo → CONFIGURE linters → FIX all violations → ENFORCE structural limits → VERIFY zero errors → loop until clean
```

The key difference from other skills: `/lint` makes **permanent, deterministic changes** to the repo's quality infrastructure. After this skill runs, `ruff check .` and `npm run lint` enforce the same strict rules everywhere — CI, pre-commit, every developer's editor.

## Autonomous Execution

**THIS WORKFLOW IS 100% AUTONOMOUS. YOU MUST:**

1. **NEVER ask for confirmation** — No "Should I add this rule?"
2. **Auto-fix all violations** — Apply fixes without asking
3. **Auto-commit and push** — Commit after each logical batch, push after final pass
4. **Create/update config files** — This skill WRITES linter configs INTO the repo
5. **Enforce structural limits on ALL files** — Not just changed files. Every code file.
6. **Fill out checkpoint honestly** — The stop hook validates your booleans

**Only stop when zero violations remain.**

## Triggers

- `/lint` (primary), `/enforce` (alias)
- `/lint [scope]` (scoped to directory/file)
- "set up linting", "make this repo clean", "enforce quality", "configure linters"

## Scope Arguments

| Invocation | Behavior |
|------------|----------|
| `/lint` | Full repo — research, configure, fix everything |
| `/lint src/` | Research full repo, fix only that directory |
| `/lint --config-only` | Research and configure — don't fix violations |

## Phase 0: Activation

```bash
mkdir -p .claude && cat > .claude/autonomous-state.json << 'EOF'
{
  "mode": "enforce",
  "started_at": "TIMESTAMP",
  "iteration": 1,
  "coordinator": true
}
EOF
cp .claude/autonomous-state.json ~/.claude/autonomous-state.json
```

## Phase 1: Research

**Goal**: Understand the repo's stack, existing quality infrastructure, and gap to "deterministically clean."

### 1.1 Stack Detection

Read filesystem markers:

| Marker | Stack |
|--------|-------|
| `pyproject.toml`, `setup.py`, `requirements.txt` | Python |
| `package.json` | JavaScript |
| `tsconfig.json` | TypeScript |

### 1.2 Existing Config Inventory

**Python:** Check for `[tool.ruff]` in pyproject.toml, `ruff.toml`, `.ruff.toml`. Also check for legacy configs: `.pylintrc`, `.flake8`, `[tool.mypy]`.

**JS/TS:** Check for `eslint.config.mjs`/`.js` (flat config), `.eslintrc.*` (legacy — migrate), `biome.json`/`biome.jsonc` (alternative linter). Check `tsconfig.json` for `strict` flag.

### 1.3 Config Gap Analysis

Compare existing rules against the strict reference:

**Python strict ruleset:**
```
F, E, W, B, UP, C4, SIM, C90, I, RUF, PERF, FURB, PIE, T20, ERA
```

**JS/TS strict ruleset:**
```
eslint:recommended + unicorn/recommended + max-lines:400 + max-lines-per-function:80
typescript-eslint/recommended (if TS)
```

The delta between what the repo has and what it should have is the work.

### 1.4 Full File Census

Scan ALL code files for structural violations. Use Bash to count lines:

```bash
# Python files
find . -name "*.py" -not -path "*/node_modules/*" -not -path "*/.venv/*" -not -path "*/__pycache__/*" -not -path "*/migrations/*" -not -path "*/.git/*" -exec wc -l {} + | sort -rn | head -20

# JS/TS files
find . \( -name "*.ts" -o -name "*.tsx" -o -name "*.js" -o -name "*.jsx" \) -not -path "*/node_modules/*" -not -path "*/dist/*" -not -path "*/.next/*" -not -path "*/.git/*" -exec wc -l {} + | sort -rn | head -20
```

Flag: files >400 LOC (MUST split), functions >80 LOC (MUST extract).

## Phase 2: Configure

**Goal**: Create or upgrade linter configs IN THE REPO.

### 2.1 Python Configuration

Create or upgrade `[tool.ruff]` in `pyproject.toml`:

```toml
[tool.ruff]
line-length = 100
target-version = "py39"  # Match project's minimum

[tool.ruff.lint]
select = [
    "F",     # Pyflakes — logic errors, undefined names, unused imports
    "E",     # pycodestyle errors
    "W",     # pycodestyle warnings
    "B",     # flake8-bugbear — common bugs and design issues
    "UP",    # pyupgrade — modern Python syntax
    "C4",    # flake8-comprehensions
    "SIM",   # flake8-simplify
    "C90",   # McCabe complexity
    "I",     # isort — import sorting
    "RUF",   # Ruff-specific rules
    "PERF",  # Perflint — performance anti-patterns
    "FURB",  # refurb — modern Python idioms
    "PIE",   # flake8-pie — misc linting
    "T20",   # flake8-print — no print() in production
    "ERA",   # eradicate — commented-out code detection
]
ignore = [
    "E501",  # line too long (formatter handles this)
]
fixable = ["ALL"]

[tool.ruff.lint.mccabe]
max-complexity = 10

[tool.ruff.lint.per-file-ignores]
"test_*.py" = ["S101", "T20"]
"**/tests/**/*.py" = ["S101", "T20"]
"conftest.py" = ["F401"]
"__init__.py" = ["F401"]
```

**If legacy linter configs exist** (pylintrc, flake8): Migrate rules to ruff, remove old configs.

**If ruff config already exists**: Upgrade by adding missing rule groups. Never remove existing rules.

### 2.2 JavaScript/TypeScript Configuration

Create `eslint.config.mjs` (flat config):

```javascript
import js from "@eslint/js";

let unicornConfig = {};
try {
  const unicorn = await import("eslint-plugin-unicorn");
  unicornConfig = unicorn.default.configs["recommended"];
} catch {}

export default [
  js.configs.recommended,
  ...(unicornConfig ? [unicornConfig] : []),
  {
    rules: {
      "max-lines": ["error", { max: 400, skipBlankLines: true, skipComments: true }],
      "max-lines-per-function": ["error", { max: 80, skipBlankLines: true, skipComments: true, IIFEs: true }],
      "no-unused-vars": ["error", { argsIgnorePattern: "^_" }],
      "no-var": "error",
      "prefer-const": "error",
      "no-console": "warn",
      "eqeqeq": ["error", "always"],
      ...(unicornConfig ? {
        "unicorn/no-null": "off",
        "unicorn/no-array-reduce": "off",
        "unicorn/no-array-for-each": "off",
        "unicorn/prevent-abbreviations": "off",
        "unicorn/filename-case": "off",
        "unicorn/no-array-callback-reference": "off",
        "unicorn/consistent-function-scoping": "off",
      } : {}),
    },
  },
  {
    ignores: [
      "node_modules/**", "dist/**", "build/**", ".next/**", ".nuxt/**",
      "coverage/**", "**/*.min.js", "**/*.min.css", "**/__generated__/**",
      "**/__snapshots__/**", "**/migrations/**", "**/*.d.ts",
    ],
  },
];
```

If the project uses TypeScript, also add `typescript-eslint`:

```javascript
import tseslint from "typescript-eslint";
// Add to config array: ...tseslint.configs.recommended,
```

**Install missing devDependencies** (use the project's package manager):

```bash
npm install -D eslint @eslint/js eslint-plugin-unicorn  # or pnpm add -D
# If TypeScript: npm install -D typescript-eslint
```

**Add lint scripts to package.json** if missing:

```json
{ "scripts": { "lint": "eslint .", "lint:fix": "eslint . --fix" } }
```

### 2.3 TypeScript Strict Mode

If `tsconfig.json` exists, ensure strict options:

```json
{
  "compilerOptions": {
    "strict": true,
    "noUncheckedIndexedAccess": true,
    "noImplicitOverride": true,
    "forceConsistentCasingInFileNames": true
  }
}
```

**If enabling `strict` produces >100 type errors**: Enable incrementally — add `// @ts-expect-error` TODOs for existing violations and track as a separate fix batch.

### 2.4 Commit Configuration

```bash
git add pyproject.toml eslint.config.mjs tsconfig.json package.json
git commit -m "enforce: configure strict linting rules"
```

## Phase 3: Baseline Capture

Run all configured linters, capture full error output, count violations by category:

```bash
# Python
ruff check . 2>&1

# JS/TS
npx eslint . 2>&1

# TypeScript type check
npx tsc --noEmit 2>&1
```

Store baseline counts in `.claude/autonomous-state.json` for progress tracking.

## Phase 4: Iterative Fix Loop

### Fix Order (mandatory — minimizes regressions)

1. **Auto-fixable linter errors** (cheapest, highest volume)
2. **Manual linter errors** (require code understanding)
3. **Type errors** (may cascade — fix in dependency order)
4. **Structural: oversized files** (>400 LOC — split into modules)
5. **Structural: oversized functions** (>80 LOC — extract helpers)
6. **Dead code removal** (unused imports, commented-out code, unreachable branches)

### Batch 1: Auto-Fix Pass

```bash
ruff check --fix . && ruff format .     # Python
npx eslint . --fix                       # JS/TS
```

Commit: `enforce: auto-fix linter violations`

### Batch 2: Manual Linter Fixes

Fix remaining violations via Edit tool. Group by file. Run linters after each group to catch cascades.

Commit: `enforce: fix manual linter violations`

### Batch 3: Type Error Fixes

Fix `tsc --noEmit` errors. Replace `any` with proper types. Add missing return types. Fix null safety.

Commit: `enforce: fix type errors`

### Batch 4: Split Oversized Files

For each file >400 LOC:

1. Analyze responsibilities — identify distinct concerns
2. Plan the split — each new module <300 lines (leave headroom)
3. Extract modules — move related functions/classes
4. Update imports — fix all import paths across codebase
5. Verify — linters + tests pass

| File Type | Split Strategy |
|-----------|---------------|
| Python module with multiple classes | One class per file |
| Python module with utility functions | Group by concern |
| React component file | Separate components, hooks, utils |
| Route handler file | One route group per file |
| Config/settings file | Skip (allowlisted) |
| Test file | Skip (relaxed limits) |

Commit: `enforce: split [filename] into focused modules`

### Batch 5: Extract Oversized Functions

For each function >80 LOC:

| Pattern | Extraction |
|---------|------------|
| Sequential steps | Extract each step as a named function |
| Nested conditionals | Extract branches as predicate + handler |
| Loop body >20 lines | Extract loop body |
| Data transformation pipeline | Extract each transform step |

Commit: `enforce: extract oversized functions`

### Batch 6: Dead Code Removal

```bash
ruff check --select F401,F841,ERA . --fix  # Unused imports, vars, commented-out code
```

Commit: `enforce: remove dead code`

### Parallel Batch Execution

If batches touch disjoint file sets, launch as parallel `Task()` agents:

```
Task(subagent_type="general-purpose", description="Split oversized files",
  prompt="Split these files to <400 LOC: [list]. Run linters after. Commit.")
Task(subagent_type="general-purpose", description="Extract oversized functions",
  prompt="Extract functions >80 LOC in: [list]. Run linters after. Commit.")
```

If overlapping: fix serially.

## Phase 5: Full Verification

### 5.1 Linter Verification (Zero Tolerance)

```bash
ruff check .            # Must exit 0
npx eslint .            # Must exit 0
npx tsc --noEmit        # Must exit 0
```

If any fails → return to Phase 4.

### 5.2 Structural Verification (ALL Files)

Scan EVERY code file (not just changed files). Apply exclusions from `_linter.py`:

**Excluded directories**: node_modules, .venv, dist, build, .next, __pycache__, migrations, .git, .claude

**Excluded files**: *.min.js, *.d.ts, *.generated.*, lock files, conftest.py, settings.py, *.config.js/ts

**Relaxed limits for test directories** (tests/, test/, __tests__/, spec/, e2e/): 600 LOC files, 120 LOC functions.

| Check | Limit | Applies To |
|-------|-------|-----------|
| File length | 400 lines | All code files (.py, .js, .jsx, .ts, .tsx) |
| Function length | 80 lines | All functions/methods |
| Cyclomatic complexity | 10 | Via ruff C90 / ESLint complexity rule |

### 5.3 Test Suite

```bash
pytest 2>/dev/null || python -m pytest 2>/dev/null   # Python
npm test 2>/dev/null                                   # JS/TS
```

If tests fail from our changes → fix the regression. If pre-existing → document in checkpoint.

### 5.4 Regression Loop

If Phase 5 finds new violations → return to Phase 4. Max 3 regression loops, then document remaining issues.

## Phase 6: Completion

### Completion Checkpoint

```json
{
  "self_report": {
    "is_job_complete": true,
    "code_changes_made": true,
    "linters_pass": true,
    "category": "quality"
  },
  "reflection": {
    "what_was_done": "Configured ruff + ESLint strict rules, fixed N violations across M files, enforced <400 LOC / <80 LOC on all code files",
    "what_remains": "none",
    "key_insight": "Reusable lesson about this repo's quality patterns (>50 chars)",
    "search_terms": ["lint", "ruff", "eslint", "quality", "enforcement"]
  },
  "verification": {
    "tests_executed_at_version": "abc1234",
    "tests": [
      {
        "id": "ruff_clean",
        "type": "command_output",
        "expected": "EXIT_CODE=0",
        "actual": "EXIT_CODE=0, All checks passed!",
        "passed": true
      },
      {
        "id": "eslint_clean",
        "type": "command_output",
        "expected": "EXIT_CODE=0",
        "actual": "EXIT_CODE=0",
        "passed": true
      },
      {
        "id": "structural_limits",
        "type": "command_output",
        "expected": "All files <400 LOC, all functions <80 LOC",
        "actual": "N files checked, 0 violations",
        "passed": true
      }
    ]
  },
  "enforce_metrics": {
    "violations_at_start": 0,
    "violations_at_end": 0,
    "files_over_400_loc_at_start": 0,
    "files_over_400_loc_at_end": 0,
    "functions_over_80_loc_at_start": 0,
    "functions_over_80_loc_at_end": 0,
    "configs_created": [],
    "configs_upgraded": [],
    "iterations_required": 1,
    "files_modified": 0
  }
}
```

### Exit Conditions

| Condition | Result |
|-----------|--------|
| All linters pass, all structural checks pass, `what_remains: "none"` | SUCCESS |
| Any linter errors remain | BLOCKED — continue fixing |
| Any structural violations remain | BLOCKED — continue fixing |
| Tests failing from our changes | BLOCKED — fix regressions |

### Cleanup

```bash
rm -f ~/.claude/autonomous-state.json .claude/autonomous-state.json
```

## Structural Limits Reference

| Metric | Limit | Relaxed For |
|--------|-------|-------------|
| File length | 400 lines | Test dirs (600), config files (exempt) |
| Function length | 80 lines | Test functions (120) |
| Cyclomatic complexity | 10 | — |

**Code files**: `.py`, `.js`, `.jsx`, `.ts`, `.tsx`, `.mjs`, `.cjs`

**Exempt**: Config files, generated files, lock files, minified files.

## Quality Hierarchy

```
/lint (enforce)  = Install the guardrails (one-time setup + full sweep)
_linter.py       = Check the guardrails at stop-time (passive gate, every session)
/melt step 4     = Pass through the guardrails during a build (gate within task)
/burndown        = Fix what guardrails don't catch (architecture, slop, semantics)
```

`/lint` makes the repo **self-enforcing**. After it runs, `_linter.py` uses the project's own strict configs instead of toolkit fallbacks. Every future `/melt` session, CI run, and developer inherits the strict rules. The quality investment compounds.

## Comparison

| Aspect | `_linter.py` | `/burndown` | `/lint` (enforce) |
|--------|-------------|-------------|-------------------|
| **When** | Stop hook (passive) | On-demand (debt) | On-demand (quality infra) |
| **Config** | `--config` flags (ephemeral) | Uses `_linter.py` | **Writes configs INTO repo** |
| **Structural** | Changed files only | Changed files only | **ALL files in repo** |
| **Output** | Pass/fail gate | Code fixes | **Config + fixes (permanent)** |
| **Post-run** | Same as before | Code improved | **Repo self-enforcing** |

## Skill Fluidity

Use techniques from any skill inline. Need to understand a complex file before splitting? Use /heavy patterns. Discover a bug while fixing types? Debug inline. Your autonomous state and checkpoint remain governed by /enforce.
