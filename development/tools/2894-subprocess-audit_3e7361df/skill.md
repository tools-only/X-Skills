# Subprocess Audit: Python Scripts Using CLI Tools via subprocess

Scripts that call `gh`, `glab`, or `git` via `subprocess` instead of using Python libraries.
All identified instances require migration to library equivalents.

## Files Requiring Migration

### 1. `.claude/skills/create-merge-request-changelog/scripts/fetch_gitlab_mr.py`

**Current pattern:** `subprocess.run(["glab", ...])` for MR metadata, commits, and diffs

**Replace with:** `python-gitlab` library

**Reference example:** `plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py` (line 4+)

---

### 2. `.claude/skills/create-merge-request-changelog/scripts/analyze_git_changes.py`

**Current pattern:** `subprocess.run(["git", ...])` for `log`, `diff`, `rev-parse`, `merge-base`, `rev-list`

**Replace with:** `GitPython` library

**Reference example:** `plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py` (lines 31, 48–50, 64)

---

### 3. `.claude/skills/daily-releases/scripts/publish_daily_release.py`

**Current pattern:** `subprocess` calls for both `git` and `gh` CLI operations

**Replace git calls with:** `GitPython`

**Replace gh calls with:** `PyGithub`

**Reference examples:**
- PyGithub: `.claude/skills/gh/scripts/github_project_setup.py` (line 6+)
- PyGithub: `.claude/skills/backlog/scripts/backlog.py` (line 55+)
- GitPython: `plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py` (lines 31, 48–50, 64)

---

### 4. `.claude/skills/daily-releases/scripts/list_daily_ranges.py`

**Current pattern:** `subprocess` calls for `git` operations

**Replace with:** `GitPython`

**Reference example:** `plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py` (lines 31, 48–50, 64)

---

### 5. `plugins/holistic-linting/skills/holistic-linting/scripts/install_agents.py`

**Current pattern:** `subprocess.run(["git", "rev-parse", "--show-toplevel"])` — single call to find git root

**Replace with:**

```python
import git
git.Repo(search_parent_directories=True).git.rev_parse("--show-toplevel")
```

**Reference example:** `plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py`

---

## Borderline Cases

### 6. `plugins/holistic-linting/skills/holistic-linting/scripts/discover_linters.py`

**Current pattern:** `subprocess` to `git` — used only to detect if running inside a git repo and to find `hooksPath` config

**Replace with:** `GitPython`

**Replacement pattern:**

```python
import git
try:
    repo = git.Repo(search_parent_directories=True)
    hooks_path = repo.config_reader().get_value("core", "hooksPath", None)
except git.InvalidGitRepositoryError:
    # not in a git repo
    hooks_path = None
```

---

## Available Libraries

Libraries already present in `pyproject.toml`:

| Library | Version constraint | Purpose |
|---|---|---|
| `gitpython` | `>=3.1.45` | Git operations (log, diff, rev-parse, merge-base, rev-list, repo root) |
| `pygithub` | `>=2.1.1` | GitHub API (gh CLI replacement) |

`python-gitlab` is available via PEP 723 inline metadata only — not in `pyproject.toml`. Add it as a dev dependency if scripts using it need IDE tooling support.

---

## Reference File Locations

- **PyGithub pattern:** `.claude/skills/backlog/scripts/backlog.py` line 55+
- **PyGithub pattern:** `.claude/skills/gh/scripts/github_project_setup.py` line 6+
- **python-gitlab pattern:** `plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py` line 4+
- **GitPython pattern:** `plugins/gitlab-skill/skills/gitlab-skill/scripts/gitlab_context.py` lines 31, 48–50, 64
