# Repository Guidelines

## Project Structure & Module Organization
- The runtime package lives in `hodor/`: `cli.py` exposes the Click entrypoint, `agent.py` manages the review loop using OpenHands SDK
- `workspace.py` handles repo cloning, CI detection, and PR branch checkout for GitHub and GitLab
- `llm/openhands_client.py` wraps the OpenHands SDK with proper model configuration and LLM provider setup
- `hodor/prompts/pr_review_prompt.py` contains the PR review prompt (get diff → review only changed files → analyze)
- `skills.py` implements skill discovery from `.cursorrules`, `agents.md`, and `.hodor/skills/` for repository-specific guidelines
- `tools/` hosts integration helpers (GitHub/GitLab CLI wrappers)
- Keep credentials and local overrides out of git by editing `.env` (mirrors `.env.example`)
- New tests and fixtures should land in `tests/`, which pytest auto-discovers via `test_*.py`

## OpenHands SDK Architecture
Hodor uses the [OpenHands Agent SDK](https://github.com/OpenHands/agent-sdk) to run AI-powered code reviews:
- **Agent Creation**: `create_hodor_agent()` configures the LLM provider, model, and terminal settings
- **Workspace**: Creates temporary directories or reuses CI workspace (`$CI_PROJECT_DIR`, `$GITHUB_WORKSPACE`)
- **Conversation**: Sends the review prompt and streams events back for monitoring
- **Terminal Type**: Uses `subprocess` PTY (not tmux) to avoid environment variable length limits
- **Event Callbacks**: Monitors agent progress with `on_event()` for real-time logging in verbose mode

### Workspace Setup
The workspace module (`workspace.py`) handles three scenarios:
1. **CI Environment**: Auto-detects GitLab CI (`$GITLAB_CI=true`) or GitHub Actions (`$GITHUB_ACTIONS=true`) and skips cloning
2. **Local Review**: Clones repo via `gh` or `glab` CLI tools and checks out PR branch
3. **Reusable Workspace**: If `--workspace-dir` specified, reuses existing clone and fetches latest changes

**Self-hosted GitLab Support**: Extracts host from URL (e.g., `gitlab.example.com`) and uses it for cloning. Falls back to `$GITLAB_HOST` env var.

### PR Review Process (3-Step Prompt)
The agent follows a strict 3-step process to avoid reviewing the entire codebase:
1. **Get Diff**: Run `gh pr diff` or `glab mr diff` to see ONLY changed files
2. **Review Changed Files**: ONLY analyze files that appear in the diff (no exploring other files)
3. **Analyze Each Change**: Look for bugs in new/modified code, ignore pre-existing code unless PR breaks it

**Critical Rules**:
- Use three-dot diff (`origin/main...HEAD`) to show only PR changes, not entire branch state
- Disable git pager with `git --no-pager` to avoid interactive prompts
- Never flag "files will be deleted" or "dependency downgrade" from stale branches

## Build, Test, and Development Commands
Use `uv` via Just to guarantee the locked toolchain:
```bash
just sync        # create/update the uv environment
just check       # run formatting, lint, and type checks
just test-cov    # execute pytest with HTML + terminal coverage
just review URL  # shortcut for `uv run hodor URL [flags]`
```
Docker workflows: `docker-build` for local testing, `docker-push REGISTRY` to publish (amd64 only).

## Coding Style & Naming Conventions
Target Python 3.13, Black formatting (120-char lines, 4-space indents), and Ruff linting with the same width. Prefer descriptive module-level names (`github_tools.py`, `workspace.py`) and snake_case for functions, UPPER_SNAKE_CASE for constants, and CapWords for classes. Run `just fmt` before committing to avoid noisy diffs.

## Testing Guidelines
Pytest is configured via `pyproject.toml` with `tests` as the root and files/functions named `test_*`. Aim to keep coverage near the default `pytest --cov=hodor` output (term-missing must stay clean). Reach for `just test-watch` when iterating locally, and regenerate the HTML report with `just test-cov` when validating larger refactors.

## Commit & Pull Request Guidelines
Follow the existing history: short, imperative subjects (`Add Hodor - …`) under 72 chars plus optional body. Squash unrelated changes, reference issues (`Fixes #123`) when applicable, and attach before/after screenshots for CLI output tweaks. Every PR description should list the motivation, testing proof (commands run), and any follow-up TODOs so reviewers can reproduce your setup.

## Security & Configuration Tips
Store API keys (`ANTHROPIC_API_KEY`, `OPENAI_API_KEY`, `GITHUB_TOKEN`, `GITLAB_TOKEN`) in `.env` and pass them through uv or Docker (`just docker-run` respects exported variables with `${VAR:-}` syntax for optional vars). Never commit credentials, lockfiles containing secrets, or sanitized logs.

### CI/CD Environment Variables
When running in GitLab CI or GitHub Actions, Hodor automatically detects:
- **GitLab CI**: `$GITLAB_CI`, `$CI_PROJECT_DIR`, `$CI_PROJECT_PATH`, `$CI_MERGE_REQUEST_IID`
- **GitHub Actions**: `$GITHUB_ACTIONS`, `$GITHUB_WORKSPACE`, `$GITHUB_REPOSITORY`

If detected, workspace setup skips cloning and uses the existing checkout. For self-hosted GitLab, the host is extracted from the MR URL automatically.

### Event Monitoring
Use `--verbose` flag to see real-time agent progress:
- Bash commands being executed
- File operations (reads/edits)
- Command exit codes (✓/✗)
- Token usage and cost estimates

Always-on metrics show: input/output tokens, cache hits, cost estimate, and review time.
