# Field Conversion Mappings

Complete field-by-field conversion tables for all resource types.

## Subagent → Skill Conversion

| Subagent Field | Skill Field | Conversion Notes |
|----------------|-------------|------------------|
| `name` | `name` | Normalize to lowercase-hyphenated |
| `description` | `description` | Enhance with trigger keywords |
| `model: claude-opus-4-*` | `model: opus` | Simplify model names |
| `model: claude-sonnet-4-*` | `model: sonnet` | |
| `model: claude-3-5-haiku-*` | `model: haiku` | |
| `model: inherit` | *(omit field)* | Inherit is default |
| `allowed-tools` / `tools` | `allowed-tools` | Keep unchanged |
| `disallowedTools` | *(document in body)* | Not supported in spec |
| `readonly: true` | `allowed-tools: Read Grep Glob` | Restrict to read-only |
| `permissionMode` | *(document in body)* | Not supported in spec |
| `skills` | *(document in body)* | Import referenced skills |
| `hooks` | `hooks` | Keep if Claude Code 2.1+ |
| *(body)* | *(body)* | Keep as skill instructions |
| *(new)* | `context: fork` | **Always add** for subagent behavior |

### Subagent Template

```yaml
---
name: {normalized-name}
description: |
  {enhanced-description}
  Use when [trigger conditions]. Identifies [what it finds/does].
context: fork
model: {simplified-model}
allowed-tools: {original-tools}
---

{original-body-content}
```

## Slash Command → Skill Conversion

| Command Field | Skill Field | Conversion Notes |
|---------------|-------------|------------------|
| `name` | `name` | Normalize to lowercase-hyphenated |
| `description` | `description` | Keep unchanged |
| `argument-hint` / `arguments` | `argument-hint` | Keep unchanged |
| `agent` | `agent` | Map to Explore/Plan/general-purpose |
| `model` | `model` | Simplify model names |
| *(body)* | *(body)* | Keep as skill instructions |
| *(new)* | `disable-model-invocation: true` | **Add** - user-invoked only |

### Command Template

```yaml
---
name: {normalized-name}
description: {original-description}
disable-model-invocation: true
argument-hint: {original-arguments}
---

{original-body-content}

## Arguments

- `$ARGUMENTS` - All arguments as single string
- `$1`, `$2`, etc. - Positional arguments
```

## Rule Conversion (Three Paths)

### Path A: Always-Apply Rules → CLAUDE.md @import

For rules with `alwaysApply: true` or high-priority project rules:

1. Keep rule file or move to `.claude/rules/` or `references/`
2. Add `@import` to CLAUDE.md:

```markdown
@import .claude/rules/security-rules.md
@import ./shared-rules/code-style.md
```

### Path B: Glob-Scoped Rules → Keep as Rules

Rules with `globs` patterns should remain as rules (skills don't support glob-based auto-activation):

```yaml
# Keep in .cursor/rules/ or .claude/rules/
---
description: React component guidelines
globs: "src/components/**/*.tsx"
alwaysApply: false
---
```

### Path C: Background Knowledge Rules → Skills

General knowledge rules without globs:

```yaml
---
name: security-guidelines
description: |
  Security best practices for this codebase. 
  Auto-activates when writing forms, controllers, queries, or handling user input.
# No disable-model-invocation - Claude can use automatically
---

{original-rule-content}
```

## Custom Prompt (Codex) → Skill Conversion

| Custom Prompt Field | Skill Field | Notes |
|--------------------|-------------|-------|
| `name` / filename | `name` | Normalize format |
| `description` | `description` | Keep unchanged |
| `model` | `model` | Map to sonnet/opus/haiku |
| `tools` | `allowed-tools` | Keep unchanged |
| *(content)* | *(body)* | Keep as instructions |

## AGENTS.md → CLAUDE.md or Skills

AGENTS.md files are simple markdown alternatives to rules. Convert by:
1. Moving content to CLAUDE.md as project memory
2. Or splitting into topic-specific skills with appropriate descriptions

## Model Name Standardization

| Legacy Format | Skills Format |
|---------------|---------------|
| `claude-opus-4-20250514` | `opus` |
| `claude-opus-4-5-20251101` | `opus` |
| `claude-sonnet-4-20250514` | `sonnet` |
| `claude-sonnet-4-5-20250514` | `sonnet` |
| `claude-3-5-haiku-20241022` | `haiku` |
| `claude-3-5-haiku-latest` | `haiku` |
| `inherit` | *(omit - inherit is default)* |
| `fast` | `haiku` |
| `default` | *(omit)* |
| `gpt-4o` | `metadata.original-model: gpt-4o` |
| `o1` | `metadata.original-model: o1` |

## Tool Name Standardization

| Variations | Standard Name |
|------------|---------------|
| `read`, `Read`, `ReadFile`, `file_read` | `Read` |
| `write`, `Write`, `WriteFile`, `file_write` | `Write` |
| `edit`, `Edit`, `EditFile`, `file_edit` | `Edit` |
| `bash`, `Bash`, `Terminal`, `Shell`, `shell` | `Bash` |
| `grep`, `Grep`, `Search`, `search` | `Grep` |
| `glob`, `Glob`, `Find`, `find_files` | `Glob` |
| `web_search`, `WebSearch`, `search_web` | `WebSearch` |
| `web_fetch`, `WebFetch`, `fetch_url` | `WebFetch` |
| `list_files`, `ListFiles`, `ls` | `Glob` |
| `ask_human`, `AskHuman`, `human_input` | `AskHuman` |

## Name Normalization

Convert names to valid skill names:

| Original | Normalized |
|----------|------------|
| `Code Reviewer` | `code-reviewer` |
| `runTests` | `run-tests` |
| `API_Builder` | `api-builder` |
| `my--skill` | `my-skill` |
| `-leadinghyphen` | `leadinghyphen` |

Rules:
- Lowercase only
- Replace spaces/underscores with hyphens
- Remove consecutive hyphens
- Remove leading/trailing hyphens
- Max 64 characters
