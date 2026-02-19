# Memory and Rules — Detailed Reference

Complete reference for Claude Code memory system configuration, behavior, and edge cases.

---

## Table of Contents

1. [Memory Type Details](#memory-type-details)
2. [CLAUDE.md Import System](#claudemd-import-system)
3. [Path-Specific Rules Reference](#path-specific-rules-reference)
4. [Auto Memory Internals](#auto-memory-internals)
5. [Organization Deployment](#organization-deployment)
6. [Environment Variables](#environment-variables)
7. [Troubleshooting](#troubleshooting)

---

## Memory Type Details

### Managed Policy

System-level instructions deployed by IT/DevOps. Cannot be overridden by users.

**Locations by platform**:

- macOS: `/Library/Application Support/ClaudeCode/CLAUDE.md`
- Linux: `/etc/claude-code/CLAUDE.md`
- Windows: `C:\Program Files\ClaudeCode\CLAUDE.md`

**Use cases**: Company coding standards, security policies, compliance requirements.

**Deployment**: Via MDM (macOS), Group Policy (Windows), Ansible/Chef/Puppet (Linux), or any configuration management system.

### User Memory

Personal preferences that apply to all projects.

**Location**: `~/.claude/CLAUDE.md`

**Shared with**: No one — personal to the user.

**Examples**: Code styling preferences, personal tooling shortcuts, preferred frameworks.

### Project Memory

Team-shared instructions for a specific project.

**Location**: `./CLAUDE.md` or `./.claude/CLAUDE.md` (either works, `.claude/` keeps root clean).

**Shared with**: Team members via source control.

**Examples**: Project architecture, coding standards, common workflows, build commands.

### Project Rules

Modular, topic-specific project instructions. Automatically discovered and loaded.

**Location**: `./.claude/rules/*.md` (recursive, supports subdirectories).

**Shared with**: Team members via source control.

**Priority**: Same level as `.claude/CLAUDE.md` — project-level.

### Project Local Memory

Personal project-specific preferences that should NOT be shared.

**Location**: `./CLAUDE.local.md`

**Shared with**: No one — automatically added to `.gitignore`.

**Use cases**: Personal sandbox URLs, preferred test data, local environment specifics.

**Worktree note**: `CLAUDE.local.md` only exists in one worktree. For cross-worktree personal instructions, use a home-directory import: `@~/.claude/my-project-instructions.md`

### Auto Memory

Claude's automatic notes and learnings.

**Location**: `~/.claude/projects/<project>/memory/`

**Project path derivation**: Based on git repository root. All subdirectories within the same repo share one auto memory directory. Git worktrees get separate directories. Outside git repos, the working directory path is used.

**Loading**: Only first 200 lines of `MEMORY.md` are loaded into system prompt. Topic files load on demand.

**What Claude saves**:
- Project patterns (build commands, test conventions, code style)
- Debugging insights (solutions to tricky problems, common error causes)
- Architecture notes (key files, module relationships, abstractions)
- User preferences (communication style, workflow habits, tool choices)

---

## CLAUDE.md Import System

### Syntax

Use `@path/to/file` anywhere in CLAUDE.md content:

```markdown
See @README for project overview.
Follow @docs/git-workflow.md for branching strategy.
Personal prefs: @~/.claude/my-prefs.md
```

### Path Resolution

- Relative paths resolve relative to the **file containing the import** (not the working directory)
- Absolute paths supported
- Home directory (`~`) paths supported

### Recursion

Imported files can recursively import additional files. Maximum depth: 5 hops.

### Code Block Protection

Imports are NOT evaluated inside code spans or code blocks:

````markdown
This is NOT imported: `@anthropic-ai/claude-code`

```text
This is NOT imported: @some/file.md
```
````

### Approval Dialog

The first time Claude Code encounters external imports in a project, it shows an approval dialog listing the specific files. This is a one-time decision per project. Declined imports remain disabled.

---

## Path-Specific Rules Reference

### Frontmatter Format

```yaml
---
paths:
  - "pattern1"
  - "pattern2"
---
```

Rules with `paths` only apply when Claude works with files matching those patterns. Rules without `paths` apply unconditionally.

### Glob Pattern Reference

**Basic patterns**:

| Pattern | Matches |
|---------|---------|
| `**/*.ts` | All TypeScript files in any directory |
| `src/**/*` | All files under `src/` directory |
| `*.md` | Markdown files in project root only |
| `src/components/*.tsx` | React components in specific directory |

**Brace expansion**:

| Pattern | Expands to |
|---------|------------|
| `src/**/*.{ts,tsx}` | `.ts` and `.tsx` files under src/ |
| `{src,lib}/**/*.ts` | TypeScript files in either src/ or lib/ |

**Multiple patterns** — specify as array:

```yaml
---
paths:
  - "src/**/*.ts"
  - "lib/**/*.ts"
  - "tests/**/*.test.ts"
---
```

### Priority Ordering

1. User-level rules (`~/.claude/rules/`) load first (lowest priority)
2. Project rules (`.claude/rules/`) load second (higher priority)
3. Within same level, more specific paths take precedence

---

## Auto Memory Internals

### Directory Structure

```text
~/.claude/projects/<project>/memory/
├── MEMORY.md          # Index — first 200 lines in system prompt
├── debugging.md       # On-demand topic file
├── patterns.md        # On-demand topic file
└── ...                # Any other topic files
```

### MEMORY.md Constraints

- First 200 lines loaded into system prompt every session
- Lines beyond 200 are NOT loaded automatically
- Keep index concise — link to topic files for detailed content
- Claude reads and writes this file during sessions

### Topic Files

- NOT loaded at startup
- Claude reads on demand using standard file tools
- No size constraint
- Use descriptive names (e.g., `debugging.md`, `api-conventions.md`)
- `MEMORY.md` should reference topic files so Claude knows they exist

### User Interaction Commands

- `/memory` — Open file selector for all memory files
- "Remember that we use pnpm" — Direct instruction to Claude to save
- "Save to memory that API tests require Redis" — Explicit save request
- "Stop remembering X" or "Forget X" — Claude finds and removes entries

---

## Organization Deployment

### Setup

1. Create managed memory file at the platform-specific location
2. Deploy via configuration management (MDM, Group Policy, Ansible, etc.)
3. File is read-only for users — enforced by the system

### Content Guidelines

Organization-level memory is for:
- Security policies that must be enforced organization-wide
- Compliance requirements that cannot be overridden
- Standardized configurations (approved tools, deployment targets)

---

## Environment Variables

| Variable | Values | Effect |
|----------|--------|--------|
| `CLAUDE_CODE_DISABLE_AUTO_MEMORY` | `1` | Force auto memory OFF |
| `CLAUDE_CODE_DISABLE_AUTO_MEMORY` | `0` | Force auto memory ON |
| (unset) | — | Follow gradual rollout |
| `CLAUDE_CODE_ADDITIONAL_DIRECTORIES_CLAUDE_MD` | `1` | Load CLAUDE.md from `--add-dir` directories |

### Additional Directory Memory

```bash
# Add directory without loading its memory
claude --add-dir ../shared-config

# Add directory AND load its CLAUDE.md, .claude/CLAUDE.md, .claude/rules/*.md
CLAUDE_CODE_ADDITIONAL_DIRECTORIES_CLAUDE_MD=1 claude --add-dir ../shared-config
```

---

## Troubleshooting

### Memory Not Loading

1. **Run `/memory`** to see which files are detected and loaded
2. **Check file location** against the hierarchy — files must be in expected paths
3. **Verify file exists** at the exact path (case-sensitive on Linux)
4. **Check `.gitignore`** — `CLAUDE.local.md` is auto-gitignored, other files should not be

### Import Not Resolving

1. **Check path is relative to the importing file**, not the working directory
2. **Verify the file exists** at the resolved path
3. **Check for code blocks** — imports inside `` ` `` or ` ``` ` are not evaluated
4. **Check approval status** — first-time imports require user approval; declined imports stay disabled
5. **Check recursion depth** — max 5 hops for recursive imports

### Path-Specific Rules Not Applying

1. **Check glob pattern** matches actual file paths
2. **Test patterns** against specific files you expect to match
3. **Verify YAML frontmatter** is properly formatted with `---` delimiters
4. **Check file is in `.claude/rules/`** directory (not somewhere else)

### Auto Memory Not Working

1. **Check environment variable**: `echo $CLAUDE_CODE_DISABLE_AUTO_MEMORY`
2. **Set explicitly**: `export CLAUDE_CODE_DISABLE_AUTO_MEMORY=0` to force on
3. **Check MEMORY.md size**: Content beyond line 200 is not auto-loaded
4. **Verify project directory**: Auto memory is per-project based on git root

---

## Sources

- [Claude Code Memory Documentation](https://code.claude.com/docs/en/memory.md) (accessed 2026-02-17)
- [Claude Code Settings Documentation](https://code.claude.com/docs/en/settings.md) (accessed 2026-02-17)
