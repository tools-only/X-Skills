---
name: lint
description: Run linting and formatting on files or discover project linters. Usage with /lint path or /lint init with optional --force flag
color: orange
---

# Linting Command

This command provides manual control over linting workflows. It can run linters on specific files/directories or discover and document project linters.

## Usage Modes

### 1. Lint Files

Lint one or more files or directories:

```bash
/lint path/to/file.py
/lint path/to/directory/
/lint file1.py file2.py file3.py
```

**Behavior**:

1. Read `## LINTERS` section from project's `CLAUDE.md` to identify configured linters
2. Run formatters first (auto-fix trivial issues)
3. Run linters second (report substantive issues)
4. If errors found, use linting-root-cause-resolver agent to fix systematically
5. Re-run linters to verify resolution

### 2. Discover Project Linters

Scan the project and generate the `## LINTERS` section for `CLAUDE.md`:

```bash
/lint init
/lint init --force  # Overwrite existing LINTERS section
```

**Behavior**:

1. Scan for linting configuration files:
   - `.pre-commit-config.yaml`
   - `pyproject.toml` (ruff, mypy, pyright, bandit)
   - `package.json` (eslint, prettier)
   - `.husky/` directory
   - Root `.eslintrc*`, `.prettierrc*`, `.markdownlint*` files
2. Identify configured formatters and linters
3. Generate `## LINTERS` section in standard format
4. Append to `CLAUDE.md` (or update if `--force` provided)

## Implementation

When this command is invoked, perform the following steps based on the arguments:

### For `/lint <path>` (Lint Mode)

1. **Read project linting configuration**:

   ```claude
   Grep(pattern="^## LINTERS", path="CLAUDE.md", output_mode="content", -A=50)
   ```

2. **If LINTERS section not found**:

   - Inform user: "No ## LINTERS section found in CLAUDE.md. Run `/lint init` first to discover project linters."
   - Exit

3. **Parse LINTERS section** to identify:

   - Formatters list with file patterns
   - Linters list with file patterns

4. **Match file paths to tools**:

   - For each provided path, determine which formatters/linters apply based on file extension
   - Example: `*.py` files → ruff format, ruff check, mypy, pyright

5. **Run formatters first** (auto-fix phase):

   ```bash
   # Python
   uv run ruff format <file.py>

   # JavaScript/TypeScript
   npx prettier --write <file.ts>

   # Markdown
   npx markdownlint-cli2 --fix <file.md>

   # Shell
   shfmt -w <script.sh>
   ```

6. **Run linters second** (validation phase):

   ```bash
   # Python
   uv run ruff check <file.py>
   uv run mypy <file.py>
   uv run pyright <file.py>

   # JavaScript/TypeScript
   npx eslint <file.ts>

   # Shell
   shellcheck <script.sh>
   ```

7. **If linting errors found**:

   - For each file with errors, launch a linting-root-cause-resolver agent:
     ```claude
     Task(subagent_type="linting-root-cause-resolver",
          description="Fix linting errors in <filename>",
          prompt="...")
     ```

8. **Verify resolution**:
   - Re-run linters on all files
   - Confirm all issues resolved

### For `/lint init` (Discovery Mode)

1. **Check for existing LINTERS section**:

   ```claude
   Grep(pattern="^## LINTERS", path="CLAUDE.md", output_mode="content")
   ```

2. **If section exists and `--force` not provided**:

   - Inform user: "## LINTERS section already exists in CLAUDE.md. Use `/lint init --force` to overwrite."
   - Show existing configuration
   - Exit

3. **Scan for git pre-commit hooks**:

   ```bash
   test -d .git && echo "Git repository: yes" || echo "Git repository: no"
   test -f .pre-commit-config.yaml && echo "pre-commit config: found" || echo "pre-commit config: not found"
   test -d .husky && echo "husky: found" || echo "husky: not found"
   ```

4. **Scan for Python linting config** (pyproject.toml):

   ```claude
   Read(file_path="pyproject.toml")
   ```

   - Look for `[tool.ruff]`, `[tool.mypy]`, `[tool.pyright]`, `[tool.bandit]` sections
   - Identify which tools are configured

5. **Scan for JavaScript/TypeScript config** (package.json):

   ```claude
   Read(file_path="package.json")
   ```

   - Look for `eslint`, `prettier`, `@typescript-eslint/*` in devDependencies
   - Check for `.eslintrc*`, `.prettierrc*` config files

6. **Scan for Markdown linting**:

   ```bash
   test -f .markdownlint.json && echo "markdownlint config: found"
   test -f .markdownlint.yaml && echo "markdownlint config: found"
   ```

7. **Scan for Shell linting**:

   - Look for shellcheck in `.pre-commit-config.yaml` or installed globally
   - Look for shfmt in `.pre-commit-config.yaml` or installed globally

8. **Generate LINTERS section**:

   ```markdown
   ## LINTERS

   git pre-commit hooks: [enabled|disabled] pre-commit tool: [husky|pre-commit|manual]

   ### Formatters

   - [tool] [file patterns] ...

   ### Static Checking and Linting

   - [tool] [file patterns] ...
   ```

9. **Append or update CLAUDE.md**:
   - If `--force` provided, remove existing section first
   - Append generated section to `CLAUDE.md`
   - Confirm success: "✓ LINTERS section written to CLAUDE.md"

## File Pattern Matching

When determining which linters apply to files, use these standard patterns:

- **Python**: `*.py` → ruff format, ruff check, mypy, pyright, bandit
- **JavaScript/TypeScript**: `*.{js,ts,jsx,tsx}` → prettier, eslint
- **Markdown**: `*.{md,markdown}` → markdownlint-cli2
- **Shell**: `*.{sh,bash,zsh,fish}` → shfmt, shellcheck
- **JSON**: `*.json` → prettier
- **YAML**: `*.{yml,yaml}` → prettier (if configured)

## Error Handling

**If CLAUDE.md doesn't exist**:

- In lint mode: Warn and suggest running `/lint init`
- In init mode: Create CLAUDE.md with LINTERS section

**If tools aren't installed**:

- Show which tools are missing
- Suggest installation commands (e.g., `uv add --dev ruff mypy pyright`)

**If linting errors persist after agent resolution**:

- Show remaining errors
- Ask user if they want to continue investigation or accept current state

## Examples

### Example 1: Lint a single file

```text
User: /lint src/auth.py

Claude:
1. Reading ## LINTERS section from CLAUDE.md...
2. Found formatters: ruff format [*.py]
3. Found linters: ruff check [*.py], mypy [*.py], pyright [*.py]
4. Running formatter: uv run ruff format src/auth.py
   ✓ Formatted 1 file
5. Running linters:
   - uv run ruff check src/auth.py
     ✗ Found 2 errors (E501, F401)
   - uv run mypy src/auth.py
     ✗ Found 1 error (arg-type)
6. Launching linting-root-cause-resolver agent for src/auth.py...
7. Agent resolved all 3 issues
8. Re-running linters:
   - uv run ruff check src/auth.py ✓
   - uv run mypy src/auth.py ✓
   - uv run pyright src/auth.py ✓
9. All linting errors resolved ✓
```

### Example 2: Lint multiple files

```text
User: /lint src/auth.py src/models.py tests/test_auth.py

Claude:
1. Reading ## LINTERS section from CLAUDE.md...
2. Running formatters on 3 files...
   ✓ All files formatted
3. Running linters on 3 files...
   ✗ src/auth.py: 2 errors
   ✗ src/models.py: 5 errors
   ✓ tests/test_auth.py: no errors
4. Launching concurrent linting-root-cause-resolver agents:
   - Agent 1: Resolving src/auth.py
   - Agent 2: Resolving src/models.py
5. Both agents completed successfully
6. Re-running linters: All files clean ✓
```

### Example 3: Discover project linters

```text
User: /lint init

Claude:
1. Checking for existing ## LINTERS section in CLAUDE.md...
   Not found - proceeding with discovery
2. Scanning project configuration:
   ✓ Git repository detected
   ✓ Found .pre-commit-config.yaml (6 hooks)
   ✓ Found pyproject.toml with [tool.ruff], [tool.mypy], [tool.pyright]
   ✓ Found package.json with prettier, eslint
   ✓ Found .markdownlint.json
3. Generating LINTERS section...

## LINTERS

git pre-commit hooks: enabled
pre-commit tool: pre-commit

### Formatters

- ruff format [*.py]
- prettier [*.{ts,tsx,json,md}]
- markdownlint-cli2 [*.md]

### Static Checking and Linting

- ruff check [*.py]
- mypy [*.py]
- pyright [*.py]
- eslint [*.{ts,tsx}]
- markdownlint-cli2 [*.md]

4. Appended to CLAUDE.md ✓
```

### Example 4: Re-discover linters (force mode)

```text
User: /lint init --force

Claude:
1. Checking for existing ## LINTERS section in CLAUDE.md...
   Found existing section
2. --force flag provided, removing existing section
3. Re-scanning project configuration...
   [Same discovery process as Example 3]
4. Updated CLAUDE.md with new configuration ✓
```

## Notes

- The `/lint` command respects the holistic-linting skill philosophy: format first (auto-fix trivial issues), lint second (report substantive issues), resolve systematically (understand root causes)
- Orchestrators should launch concurrent linting-root-cause-resolver agents when multiple files have errors
- Sub-agents should use this command on their touched files before completing tasks
- The init mode creates a "cache" in CLAUDE.md to avoid repeated configuration discovery overhead
