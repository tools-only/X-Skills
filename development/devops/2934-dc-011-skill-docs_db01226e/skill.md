# DC-011: Update Skill Files with Container Mode Docs

**Level**: 5 | **Critical Path**: No | **Estimate**: 20 min
**Dependencies**: DC-010

## Objective

Update `zerg:init.md` and `zerg:rush.md` skill files to document:
- Multi-language project detection
- Container mode options
- Examples of container-based execution

## Files Owned

- `.claude/commands/zerg:init.md` (modify)
- `.claude/commands/zerg:rush.md` (modify)

## Implementation

### 1. Update `.claude/commands/zerg:init.md`

Add section about multi-language detection:

```markdown
## Multi-Language Detection

ZERG automatically detects all languages in your project:

- **Python**: pyproject.toml, requirements.txt, setup.py, *.py
- **JavaScript/TypeScript**: package.json, tsconfig.json, *.js, *.ts
- **Go**: go.mod, *.go
- **Rust**: Cargo.toml, *.rs
- **Java**: pom.xml, build.gradle, *.java
- **Ruby**: Gemfile, *.rb
- **C++**: CMakeLists.txt, *.cpp, *.hpp
- **R**: *.R, *.r
- **Julia**: *.jl

### Multi-Language Devcontainer

For projects with multiple languages, ZERG generates a devcontainer with multiple runtime features:

```json
{
  "image": "mcr.microsoft.com/devcontainers/base:ubuntu",
  "features": {
    "ghcr.io/devcontainers/features/python:1": {"version": "3.12"},
    "ghcr.io/devcontainers/features/node:1": {"version": "20"},
    "ghcr.io/devcontainers/features/go:1": {"version": "1.22"}
  }
}
```

### Building Container Image

To build the devcontainer image during init:

```bash
zerg init --with-containers
```

This enables container mode for `zerg rush`.
```

### 2. Update `.claude/commands/zerg:rush.md`

Add section about execution modes:

```markdown
## Execution Modes

ZERG supports three execution modes:

### Subprocess Mode (Default)

Workers run as local Python subprocesses:

```bash
zerg rush --mode subprocess --workers 5
```

- No Docker required
- Uses local environment
- Good for development

### Container Mode

Workers run in isolated Docker containers:

```bash
zerg rush --mode container --workers 5
```

Requirements:
- Docker installed and running
- Devcontainer image built (`zerg init --with-containers`)
- ANTHROPIC_API_KEY in environment

Benefits:
- Isolated environments
- Consistent dependencies
- Resource limits per worker

### Auto Mode (Default)

Automatically selects the best mode:

```bash
zerg rush --mode auto --workers 5
# or simply:
zerg rush --workers 5
```

Auto-detection logic:
1. If `.devcontainer/devcontainer.json` exists AND
2. Docker is available AND
3. `zerg-worker` image is built
→ Uses container mode

Otherwise → Uses subprocess mode

## Examples

```bash
# Quick local execution
zerg rush --workers 3

# Force container mode
zerg rush --mode container --workers 5

# Dry run to check mode selection
zerg rush --dry-run
# Shows: "Auto-detected mode: container" or "subprocess"

# Build containers first, then run
zerg init --with-containers
zerg rush --mode container --workers 5
```
```

## Verification

```bash
# Check zerg:init.md has multi-language docs
grep -l 'Multi-Language' .claude/commands/zerg:init.md

# Check zerg:rush.md has container mode docs
grep -l 'Container Mode' .claude/commands/zerg:rush.md
grep -l '\-\-mode' .claude/commands/zerg:rush.md

# Count occurrences
echo "Container mentions:"
grep -c 'container' .claude/commands/zerg:*.md
```

## Acceptance Criteria

- [ ] zerg:init.md documents multi-language detection
- [ ] zerg:init.md shows example multi-feature devcontainer.json
- [ ] zerg:init.md documents --with-containers flag
- [ ] zerg:rush.md documents subprocess, container, auto modes
- [ ] zerg:rush.md shows --mode flag usage
- [ ] zerg:rush.md explains auto-detection logic
- [ ] Both files have usage examples
- [ ] No broken markdown formatting
