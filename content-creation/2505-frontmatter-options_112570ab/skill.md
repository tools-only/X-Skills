# OpenClaw Frontmatter Options & Patterns Reference

Load before writing frontmatter in Phase 2, Step 2. Contains the full OpenClaw field catalog,
description patterns, the `metadata` JSON constraint, and path conventions. OpenClaw uses the
AgentSkills spec (pi-coding-agent) — its frontmatter fields are substantially different from
Claude Code's.

---

## Valid Frontmatter Fields

```yaml
---
name: identifier                      # Required. Hyphen-case, lowercase, digits, hyphens only.
description: >                        # Required. See description patterns below.
  [See patterns below]

# User-facing behavior
user-invocable: true                  # Show in / command menu (default true). Set false for
                                      # background-knowledge skills that auto-trigger only.
disable-model-invocation: true        # Prevent auto-triggering; user must type slash command.
argument-hint: "[arg1] [arg2]"        # Document expected arguments. Quote if value contains [...].

# Command dispatch (bypass model entirely)
command-dispatch: tool                # Route slash command directly to a tool without model.
command-tool: web_search              # Tool name for dispatch (used with command-dispatch).
command-arg-mode: raw                 # Argument passing mode. "raw" is default and most common.

# Gating and ecosystem metadata (MUST be single-line JSON on one line)
metadata: '{"openclaw": {"emoji": "🔧", "requires": {"bins": ["jq"]}}}'
---
```

---

## Fields NOT Valid in OpenClaw

These are Claude Code-only fields. Do not include them in OpenClaw skills — they will be
ignored or cause parse errors:

| Field | Claude Code purpose | OpenClaw equivalent |
|-------|---------------------|---------------------|
| `model` | Select haiku/sonnet/opus per skill | Set at agent level in `openclaw.json` |
| `context: fork` | Run in isolated sub-agent | Use `sessions_spawn` in body |
| `agent` | Route to specialized agent type | All sub-agents are general purpose |
| `allowed-tools` | Restrict tools per skill | Tool policy is gateway config, not frontmatter |
| `hooks` | Pre/PostToolUse scripts | System-level hooks, not skill-level |
| `license` | SPDX license identifier | Not a skill field |

If any of these appear in a skill being validated, `validate_claw_skill.py` will report them
as invalid with the explanation above.

---

## The `metadata` Field: Single-Line JSON Constraint

**Critical:** `metadata` must be a single-line JSON string on one line. Multi-line YAML
mappings under `metadata` are not parsed correctly by the pi-coding-agent skill loader.

```yaml
# Correct — single-line JSON
metadata: '{"openclaw": {"emoji": "🔍", "requires": {"bins": ["rg"]}}}'

# Wrong — multi-line YAML mapping
metadata:
  openclaw:
    emoji: "🔍"
    requires:
      bins:
        - rg
```

### `metadata.openclaw` Subfield Catalog

All subfields are nested inside `{"openclaw": { ... }}`:

| Field | Type | Purpose |
|-------|------|---------|
| `always` | bool | Load skill body into every session (use sparingly — costs tokens always) |
| `emoji` | string | Emoji shown in macOS OpenClaw UI |
| `homepage` | string | URL for skill documentation or repo |
| `os` | string | Platform filter: `"macos"`, `"linux"`, `"windows"`, or omit for all |
| `requires.bins` | array | Binary gating — skill hidden if any listed binary not in PATH |
| `requires.anyBins` | array | Binary gating — skill hidden if NONE of the listed binaries are in PATH |
| `requires.env` | array | Env var gating — skill hidden if any listed var not set |
| `requires.config` | array | Config key gating — skill hidden if any listed key absent |
| `primaryEnv` | string | Primary env var driving the skill (informational) |
| `skillKey` | string | Stable identifier for programmatic cross-referencing |
| `install` | object | Installer spec: `{"brew": "pkg"}`, `{"node": "pkg"}`, `{"go": "pkg"}`, `{"uv": "pkg"}` |

Examples:

```yaml
# Skill requires jq and only runs on macOS
metadata: '{"openclaw": {"emoji": "📋", "os": "macos", "requires": {"bins": ["jq"]}}}'

# Skill needs either rg or grep (anyBins — at least one must exist)
metadata: '{"openclaw": {"requires": {"anyBins": ["rg", "grep"]}}}'

# Skill needs OPENAI_API_KEY env var to be set
metadata: '{"openclaw": {"requires": {"env": ["OPENAI_API_KEY"]}}}'

# Skill with installer and homepage
metadata: '{"openclaw": {"emoji": "🌐", "homepage": "https://example.com", "install": {"brew": "my-tool"}}}'
```

---

## Description Patterns

### For Skills (auto-triggered)

- **Third-person framing is a routing signal, not a stylistic choice.** The routing model evaluates the description as a triggering condition. Third-person ("This skill should be used when...") reads as a condition to test. First-person reads as an instruction to execute.
- **Quoted phrases must be verbatim user speech.** Routing matches on literal token patterns. Write the exact words a user would type, not paraphrases.
- **The description is always in context, even when the skill isn't active.** Every session pays the token cost of every skill's name, description, and location. Density matters. Avoid restating the skill name or explaining what skills are.
- **Cover the naive phrasing.** A user who doesn't know this skill exists won't search for it by name — they'll describe their problem in plain language.
- **3–5 trigger phrases minimum.** Single-phrase descriptions have high miss rates.
- **Use `>` scalar, not `|`.** Folded scalar (`>`) collapses newlines to spaces, producing a single continuous string. Literal scalar (`|`) preserves newlines, which creates unexpected whitespace.

```yaml
# Correct — third-person, verbatim phrases, folded scalar
description: >
  This skill should be used when the user asks to "create an OpenClaw skill",
  "build a claw skill", "write a SKILL.md for openclaw", or wants to author
  a skill for the pi-coding-agent ecosystem.

# Wrong — vague, no trigger phrases, not third-person
description: Helps with skill creation.
```

### For Commands (user-invoked)

- **Verb-first, under 60 chars.** The description appears as a single scannable line in the `/` menu — treat it as a menu label.
- **Describe the action, not the tool.** "Search web for topic" orients by outcome.

```yaml
description: Search web for topic
description: Fetch and summarize a URL
description: Publish skill to ClawHub
```

---

## Token Cost Formula

Each skill costs: **195 base + 97/skill + field lengths** tokens in the system prompt
(name + description + location only — body is lazy-loaded on demand via `read`).

Keep descriptions dense and informative but not verbose. The body is free until triggered.

---

## `{baseDir}` Path Substitution

`{baseDir}` is substituted before the model sees the skill body — it becomes the absolute path
to the skill directory at load time. Use it for all skill-relative file references:

```markdown
Read `{baseDir}/references/my-reference.md` for the full option catalog.
```

```bash
python3 {baseDir}/scripts/my-script.py --output json
```

Do not use `$CLAUDE_PLUGIN_ROOT` (Claude Code) or hardcoded absolute paths in OpenClaw skills.
