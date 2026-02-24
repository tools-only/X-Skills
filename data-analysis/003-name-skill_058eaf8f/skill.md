---
name: daily-releases
description: "Create GitHub Releases with AI-analyzed changelogs for every calendar day with commits on origin/main. Uses the same analyze → AI-categorize → format pipeline as /create-merge-request-changelog for rich, structured output. Idempotent: skips days that are already up to date, updates releases where new commits have been added. Automatically invoked when command is run; accepts optional --start-date, --end-date, --branch, --dry-run arguments."
argument-hint: '[--start-date YYYY-MM-DD] [--end-date YYYY-MM-DD] [--branch BRANCH] [--dry-run]'
---
# Daily Releases

Create GitHub Releases with AI-categorized changelogs for every day that had commits. Uses the same pipeline as `/create-merge-request-changelog` — real AI analysis, not template substitution.

## Automatic Invocation

When this skill is activated, immediately begin processing without asking the user. Parse any arguments from `$ARGUMENTS`:

```text
--start-date YYYY-MM-DD   Only process days on or after this date
--end-date YYYY-MM-DD     Only process days on or before this date (default: today)
--branch BRANCH           Git branch (default: origin/main)
--dry-run                 Preview without creating releases
```

## Process

Requires `GITHUB_TOKEN` for release status checks (list) and publishing.

**Working directory:** Run all commands from the repository root. Paths below assume cwd is the repo root.

### Step 1: List days to process

```bash
uv run .claude/skills/daily-releases/scripts/list_daily_ranges.py [--branch BRANCH] [--start-date ...] [--end-date ...] [-R OWNER/REPO]
```

This outputs a JSON array. Each entry has:

```json
{
  "date": "2026-02-21",
  "tag": "v2026.02.21",
  "base_ref": "<parent-commit-hash>",
  "head_ref": "<last-commit-hash-of-day>",
  "commit_count": 12,
  "release_exists": true,
  "needs_update": false
}
```

Skip entries where `release_exists: true` and `needs_update: false` — those are up to date.

For `--dry-run`, print the list and stop.

### Step 2: For each day that needs a release

Work through days chronologically. For each day:

#### 2a. Extract git data

```bash
uv run .claude/skills/create-merge-request-changelog/scripts/analyze_git_changes.py \
  <base_ref> <head_ref> ./daily-releases/<date>/
```

This writes to `./daily-releases/<date>/`:

- `commits_oneline.txt` — one-line commit list
- `commits_detailed.txt` — full commit messages with metadata
- `changes.diff` — unified diff
- `changes_stat.txt` — diffstat summary
- `changed_files.txt` — file list with A/M/D status
- `changes_numstat.txt` — per-file line counts
- `summary.json` — machine-readable stats

#### 2b. AI analysis

Read the extracted data and apply the primary analysis prompt from:

```text
.claude/skills/create-merge-request-changelog/references/analysis_prompts.md
```

Use the **Primary Analysis Prompt** section with the day's data substituted in:

- `{commit_details}` → contents of `commits_detailed.txt`
- `{changes_diff}` → contents of `changes.diff` (truncate to ~4000 chars if very large; use `changes_stat.txt` summary instead)
- `{changed_files}` → contents of `changed_files.txt`
- `{commit_count}`, `{files_changed}`, `{lines_added}`, `{lines_deleted}` → from `summary.json`

Produce the structured JSON analysis output defined in the prompt's `<output_format>` section.

Also generate a title using the **Title Generation Prompt**:

```text
Daily Release - <date>
```

Save the complete analysis JSON to `./daily-releases/<date>/analysis.json`.

#### 2c. Format into release notes

```bash
uv run .claude/skills/create-merge-request-changelog/scripts/format_mr_description.py \
  ./daily-releases/<date>/analysis.json \
  --no-preview \
  --output ./daily-releases/<date>/description.md
```

#### 2d. Publish the release

```bash
uv run .claude/skills/daily-releases/scripts/publish_daily_release.py \
  --date <date> \
  --tag <tag> \
  --head-ref <head_ref> \
  --notes-file ./daily-releases/<date>/description.md
```

Add `--keep-existing-tag=false` if updating a release that already has the correct tag commit.

### Step 3: Report

After processing all days, print a summary:

```text
Processed N days:
  - Created: X new releases
  - Updated: Y existing releases
  - Skipped: Z already up to date
```

## Reference files

- [./scripts/list_daily_ranges.py](./scripts/list_daily_ranges.py) — list days + commit ranges
- [./scripts/publish_daily_release.py](./scripts/publish_daily_release.py) — create/update git tag + GitHub release
- [../create-merge-request-changelog/scripts/analyze_git_changes.py](../create-merge-request-changelog/scripts/analyze_git_changes.py) — extract git data per day
- [../create-merge-request-changelog/references/analysis_prompts.md](../create-merge-request-changelog/references/analysis_prompts.md) — AI analysis prompts
- [../create-merge-request-changelog/scripts/format_mr_description.py](../create-merge-request-changelog/scripts/format_mr_description.py) — render AI analysis to markdown

Reference paths above are relative to this skill directory; CLI commands use repo-root paths.
