# Input Resolution and Docs Discovery

Covers Phase 0 of `user-docs-to-ai-skill` — resolving the `source` input to a local directory, deriving `output_skill` when not provided, and locating documentation within the resolved directory.

## Table of Contents

1. [Source Type Resolution](#source-type-resolution)
2. [Project Name Derivation](#project-name-derivation)
3. [output_skill Derivation](#output_skill-derivation)
4. [Docs Discovery](#docs-discovery)
5. [Anti-Patterns](#anti-patterns)

---

## Source Type Resolution

```mermaid
flowchart TD
    Input([source input received]) --> Q{Starts with https://github.com/?}
    Q -->|Yes — GitHub URL| Clone["Run: git clone source .claude/worktrees/project-name/\nproject-name = last URL path segment"]
    Q -->|No — local path| UseLocal[Use source value as docs_root directly]
    Clone --> Cloned{Clone succeeded?}
    Cloned -->|Yes| SetRootClone[docs_root = .claude/worktrees/project-name/]
    Cloned -->|No — network or auth error| Block([BLOCKED — report clone failure with command output])
    UseLocal --> Exists{Directory exists?}
    Exists -->|Yes| SetRootLocal[docs_root = source]
    Exists -->|No| Block2([BLOCKED — report path not found])
    SetRootClone --> Done([docs_root resolved])
    SetRootLocal --> Done
```

### Git Clone Command

```bash
git clone <source> .clone/worktrees/<project-name>/
```

- Path is relative to the project root — do not use absolute paths
- `project-name` is derived from the URL before cloning (see next section)
- If `.clone/worktrees/<project-name>/` already exists, skip the clone and use the existing directory

---

## Project Name Derivation

Extract `project-name` from the last non-empty path segment of the input.

| Input | project-name |
|-------|-------------|
| `https://github.com/astral-sh/ty` | `ty` |
| `https://github.com/anthropics/anthropic-sdk-python` | `anthropic-sdk-python` |
| `/home/user/repos/my-tool` | `my-tool` |
| `.claude/worktrees/fastmcp` | `fastmcp` |

Strip trailing slashes before extracting the segment.

---

## output_skill Derivation

```mermaid
flowchart TD
    Q{output_skill provided by caller?}
    Q -->|Yes| Use[Use provided value as-is]
    Q -->|No| Derive[Derive from project-name]
    Derive --> Strip[Strip common suffixes: -tool, -sdk, -py, -python, -lib]
    Strip --> Lower[Lowercase, replace non-alphanumeric with hyphens]
    Lower --> Done([output_skill set])
    Use --> Done
```

**Examples:**

| project-name | output_skill |
|-------------|-------------|
| `ty` | `ty` |
| `anthropic-sdk-python` | `anthropic-sdk` |
| `fastmcp` | `fastmcp` |
| `httpx` | `httpx` |

---

## Docs Discovery

After `docs_root` is set, locate where the documentation lives.

```mermaid
flowchart TD
    Start([docs_root resolved]) --> Q1{docs_root/docs/ exists?}
    Q1 -->|Yes| UseDocs[docs_path = docs_root/docs/\nProceed to inventory]
    Q1 -->|No| Q2{docs_root/doc/ exists?}
    Q2 -->|Yes| UseDoc[docs_path = docs_root/doc/\nProceed to inventory]
    Q2 -->|No| Scan["Delegate to Explore subagent:\nGlob all *.md files under docs_root\nReturn list of file paths"]
    Scan --> Q3{Any .md files found?}
    Q3 -->|Yes — include README.md| UseList[docs_path = list of discovered .md files\nProceed to inventory treating each as a source file]
    Q3 -->|No markdown found| CheckRST["Glob *.rst files under docs_root"]
    CheckRST --> Q4{Any .rst files found?}
    Q4 -->|Yes| UseRST[docs_path = list of .rst files]
    Q4 -->|No| Block([BLOCKED — no documentation found in docs_root\nReport what was searched])
    UseDocs --> Inventory([Phase 0d — Inventory])
    UseDoc --> Inventory
    UseList --> Inventory
    UseRST --> Inventory
```

### Explore Subagent Delegation

When no `docs/` or `doc/` directory exists, delegate discovery:

```text
Task: subagent_type="general-purpose"
Prompt: Glob all *.md and *.rst files under <docs_root>.
        Also check for a README.md at docs_root root level.
        Return the full list of file paths found.
        Do not read file contents — paths only.
Output: flat list of file paths
```

Use `general-purpose` (not `Explore`) because the Explore agent has a ~50% hallucination rate on pattern-matching tasks.

---

## Anti-Patterns

**Using absolute paths in clone destination:**

```bash
# WRONG
git clone https://github.com/astral-sh/ty /home/user/repos/.clone/worktrees/ty/

# CORRECT
git clone https://github.com/astral-sh/ty .clone/worktrees/ty/
```

**Hardcoding docs_path before checking:**

```text
# WRONG — assumes docs/ always exists
docs_path = docs_root + "/docs/"

# CORRECT — check first, fall back to discovery
```

**Fabricating output_skill from partial URL parsing:**

```text
# WRONG — astral-sh/ty parsed as "astral-sh"
project-name = first URL segment after github.com

# CORRECT — use the last segment
project-name = last URL path segment = "ty"
```
