# Auto-Sync Manifests - Automatic Version Bumping and Manifest Maintenance

## Overview

`auto-sync-manifests.py` automatically maintains `plugin.json` and `marketplace.json` files based on git changes detected during pre-commit. It detects CRUD operations on plugins and components, updates manifests accordingly, and bumps semantic versions.

## How It Works

### CRUD Detection

The script analyzes staged git changes to detect:

**Plugin-level changes:**

- **Created**: New `plugins/{name}/.claude-plugin/plugin.json` file added
- **Deleted**: `plugins/{name}/` directory removed

**Component-level changes:**

- **Skill created/deleted**: New/removed `plugins/{name}/skills/{skill-name}/SKILL.md`
- **Agent created/deleted**: New/removed `plugins/{name}/agents/{agent-name}.md`
- **Command created/deleted**: New/removed `plugins/{name}/commands/{command-name}.md`
- **Hook created/deleted**: New/removed `plugins/{name}/hooks/*.json`
- **MCP server created/deleted**: New/removed MCP configuration files

### Automatic Updates

**Plugin.json updates:**

1. Detects added/removed components
2. Updates `skills`, `agents`, `commands` arrays with proper `./` paths
3. Bumps plugin version:
   - **Major** (X.0.0): Component deleted (breaking change)
   - **Minor** (0.X.0): Component added (new feature)
   - **Patch** (0.0.X): Component modified
4. Stages updated `plugin.json` automatically

**Marketplace.json updates:**

1. Adds new plugins to `plugins` array
2. Removes deleted plugins from `plugins` array
3. Bumps marketplace version:
   - **Major** (X.0.0): Plugin deleted
   - **Minor** (0.X.0): Plugin added
   - **Patch** (0.0.X): Plugin component changed
4. Stages updated `marketplace.json` automatically

### Version Bumping Rules

```
Component Change → Plugin Version Bump → Marketplace Version Bump
───────────────────────────────────────────────────────────────────
+ New skill      → Minor (0.1.0 → 0.2.0) → Patch (1.0.0 → 1.0.1)
- Delete agent   → Major (0.1.0 → 1.0.0) → Patch (1.0.0 → 1.0.1)
~ Modify command → Patch (0.1.0 → 0.1.1) → Patch (1.0.0 → 1.0.1)

+ New plugin     → N/A                   → Minor (1.0.0 → 1.1.0)
- Delete plugin  → N/A                   → Major (1.0.0 → 2.0.0)
```

## Usage

### Pre-Commit Hook (Automatic)

The script runs automatically during `git commit` via `.pre-commit-config.yaml`:

```yaml
- repo: local
  hooks:
    - id: auto-sync-manifests
      name: Auto-sync plugin and marketplace manifests
      entry: uv run -q --no-sync plugins/plugin-creator/scripts/auto-sync-manifests.py
      language: system
      stages: [pre-commit]
      pass_filenames: false
      files: ^plugins/
```

### Manual Execution

Run manually to preview changes without committing:

```bash
# From repository root
uv run -q --no-sync plugins/plugin-creator/scripts/auto-sync-manifests.py
```

The script only operates on **staged changes** (via `git diff --cached`), so manually running it requires staged changes.

## Example Workflows

### Adding a New Skill to Existing Plugin

```bash
# 1. Create skill file
mkdir -p plugins/my-plugin/skills/new-feature
touch plugins/my-plugin/skills/new-feature/SKILL.md

# 2. Stage changes
git add plugins/my-plugin/skills/new-feature/SKILL.md

# 3. Commit (hook runs automatically)
git commit -m "feat(my-plugin): add new-feature skill"

# Output:
# ✅ Updated my-plugin → 1.1.0 (+1 components)
# ✅ Updated marketplace → 2.0.1
```

**What happened:**

- `plugins/my-plugin/.claude-plugin/plugin.json` updated:
  - `skills` array gained `./skills/new-feature`
  - `version` bumped from `1.0.0` → `1.1.0` (minor)
- `marketplace.json` updated:
  - `metadata.version` bumped from `2.0.0` → `2.0.1` (patch)

### Deleting an Agent (Breaking Change)

```bash
# 1. Remove agent file
git rm plugins/my-plugin/agents/deprecated-agent.md

# 2. Commit
git commit -m "feat(my-plugin)!: remove deprecated agent"

# Output:
# ✅ Updated my-plugin → 2.0.0 (-1 components)
# ✅ Updated marketplace → 2.0.1
```

**What happened:**

- `plugin.json` updated:
  - `agents` array no longer contains `./agents/deprecated-agent.md`
  - `version` bumped from `1.1.0` → `2.0.0` (major - breaking)
- `marketplace.json` bumped `2.0.0` → `2.0.1` (patch)

### Creating a New Plugin

```bash
# 1. Scaffold plugin
uv run -q --no-sync plugins/plugin-creator/scripts/create_plugin.py

# 2. Stage plugin files
git add plugins/new-plugin/

# 3. Commit
git commit -m "feat: add new-plugin"

# Output:
# ✅ Updated marketplace → 2.1.0
#    • Added plugins: new-plugin
```

**What happened:**

- `marketplace.json` updated:
  - `plugins` array gained `{"name": "new-plugin", "source": "./plugins/new-plugin"}`
  - `metadata.version` bumped from `2.0.1` → `2.1.0` (minor)

### Deleting a Plugin (Breaking Change)

```bash
# 1. Remove plugin directory
git rm -r plugins/obsolete-plugin/

# 2. Commit
git commit -m "feat!: remove obsolete-plugin"

# Output:
# ✅ Updated marketplace → 3.0.0
#    • Removed plugins: obsolete-plugin
```

**What happened:**

- `marketplace.json` updated:
  - `plugins` array no longer contains `obsolete-plugin`
  - `metadata.version` bumped from `2.1.0` → `3.0.0` (major - breaking)

## Edge Cases Handled

### Renamed Files

Git renames are treated as delete + add:

```bash
git mv plugins/my-plugin/agents/old-name.md plugins/my-plugin/agents/new-name.md
git commit -m "refactor(my-plugin): rename agent"

# Treated as:
# - Delete: agents/old-name.md
# + Add: agents/new-name.md
# Result: Major version bump (breaking change)
```

### Multiple Plugins Changed in One Commit

Each plugin version bumps independently:

```bash
# Changed files:
# - plugins/plugin-a/skills/new-skill/SKILL.md (added)
# - plugins/plugin-b/agents/agent.md (modified)

git commit -m "feat: update multiple plugins"

# Output:
# ✅ Updated plugin-a → 0.2.0 (+1 components)
# ✅ Updated plugin-b → 0.1.1 (~1 components)
# ✅ Updated marketplace → 2.0.1
```

### String to Array Conversion

If `plugin.json` uses string format for component paths, the script converts to array automatically:

```json
{
  "skills": "./skills/main-skill"
}
```

After adding a second skill, becomes:

```json
{
  "skills": ["./skills/main-skill", "./skills/new-skill"]
}
```

### Non-Plugin Changes

Changes outside `plugins/` directory are ignored:

```bash
git add README.md CONTRIBUTING.md
git commit -m "docs: update documentation"

# Output:
# ℹ️  No manifest updates needed
```

## Hook Behavior Details

### Silent Success (Exit Code 0)

The hook always returns exit code 0, even when it modifies files. This means:

- **Commit proceeds without interruption** - You won't see "files modified, please review" messages
- **Modified manifests auto-stage** - Updated `plugin.json` and `marketplace.json` are added to your commit automatically
- **Convenient workflow** - No need to re-run `git commit` after version bumps

This differs from typical pre-commit hooks that return exit code 1 when modifying files.

### Failure Recovery and Double-Bump Protection

**Scenario: Later hook fails (e.g., mypy, linting)**

1. You stage `SKILL.md` and run `git commit`
2. `auto-sync-manifests` runs successfully:
   - Detects modified skill
   - Bumps plugin version `2.14.0` → `2.14.1`
   - Stages `plugin.json` and `marketplace.json`
   - Returns exit code 0
3. `mypy` runs and fails
4. Commit is aborted
5. **Current state:** All three files remain staged with version `2.14.1`

**Re-running the commit:**

1. You run `git commit` again (after fixing mypy errors)
2. `auto-sync-manifests` runs again
3. **Protection activates** - Script detects `plugin.json` already staged (lines 276-280 of script)
4. **Skips version bumping** - Returns early without modification
5. Commit proceeds with version `2.14.1` (no double-bump to `2.14.2`)

**Code reference:**

```python
# Check if plugin.json is already staged - if so, skip modifying it
staged_status = run_git_command(["diff", "--cached", "--name-only"])
if str(plugin_json_path) in staged_status:
    # File is already staged, don't modify it
    return False, "0.0.0"
```

Same protection exists for `marketplace.json`.

**Key insight:** The script only bumps versions when manifests are NOT already staged. This prevents version inflation from repeated commit attempts.

## Bypassing the Hook

If you need to commit without automatic manifest updates:

```bash
git commit --no-verify -m "message"
```

**Warning:** This can create inconsistency between actual files and manifests. Only use when you have a specific reason to skip validation.

## Integration with Conventional Commits

The script is **independent** of commit message content. Version bumping is based solely on file changes detected in git, not on commit message parsing.

This means:

- ✅ Works with any commit message format
- ✅ Version bumps reflect actual code changes
- ✅ No dependency on conventional commit enforcement

However, we **recommend** using conventional commits for consistency:

```bash
# Good commit messages
feat(plugin-name): add new feature
fix(plugin-name): resolve bug
feat(plugin-name)!: breaking change
chore(plugin-name): update dependencies
```

## Troubleshooting

### Hook Not Running

**Symptom:** Commit succeeds but no version bumps occur

**Diagnosis:**

```bash
# Check if hook is registered
prek run --hook-stage pre-commit --verbose

# Check if files match pattern
git diff --cached --name-only | grep '^plugins/'
```

**Solution:** Ensure `.pre-commit-config.yaml` contains the hook and changes are staged.

### Version Not Bumping

**Symptom:** Hook runs but version unchanged

**Diagnosis:**

```bash
# Manually run script to see output
uv run -q --no-sync plugins/plugin-creator/scripts/auto-sync-manifests.py
```

**Common causes:**

- Changes to non-component files (README, CLAUDE.md, etc.)
- Changes already reflected in manifest
- `plugin.json` missing or malformed
- **Manifest already staged** - Protection prevents double-bumping (see Failure Recovery section)

### Merge Conflicts in manifests

**Symptom:** Git reports conflicts in `plugin.json` or `marketplace.json`

**Solution:**

1. Resolve conflicts manually in the JSON files
2. Validate JSON syntax: `python3 -m json.tool <file>`
3. Re-stage: `git add <file>`
4. Continue: `git commit` or `git rebase --continue`

## Performance

**Typical execution time:** <100ms

The script only processes staged changes, not the entire repository, ensuring fast pre-commit execution.

## Security Considerations

**Sandbox:** The script only modifies files within `plugins/` and `.claude-plugin/` directories.

**Git operations:** Only reads staged changes via `git diff --cached`, never modifies git history.

**File writes:** Only writes to `plugin.json` and `marketplace.json`, never creates new files.

## Source Code

Location: `plugins/plugin-creator/scripts/auto-sync-manifests.py`

Dependencies:

- Python 3.11+ standard library only
- `git` command-line tool

No external Python packages required.
