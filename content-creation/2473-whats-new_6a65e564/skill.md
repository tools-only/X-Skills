---
description: Review Claude Code updates and surface relevant new features
---

Analyze what changed in Claude Code since you last checked. Follow these steps exactly:

## Step 1: Get Versions

Run `claude --version` to get the current version.

Read `~/.claude/.claude-code-last-seen-version` to get the last-seen version (may not exist).

## Step 2: Fetch Changelog

Use WebFetch to fetch `https://raw.githubusercontent.com/anthropics/claude-code/main/CHANGELOG.md` and extract entries between the last-seen version and the current version.

If no last-seen version exists, extract the last 10 version entries.

## Step 3: Fetch Blog Announcements

Use WebFetch to fetch `https://claude.com/blog/category/announcements`.

Extract the **5 most recent** blog post titles and dates from the page. For each post whose title contains any of these keywords: "Claude Code", "CLI", "agent", "hooks", "MCP", "tools", "coding" — fetch the full post URL and extract a 1-2 sentence summary.

If the last-seen version exists, only include posts published after the date of that version's changelog entry (found in Step 2). If no last-seen version exists, include all matching posts from the 5 most recent.

If the fetch fails or no posts match the keywords, output: "No Claude Code blog announcements found."

## Step 4: Analyze Relevance

Compare changes against the current project's setup:
- Read CLAUDE.md for project conventions and tools
- Check `.claude/settings.json` for current configuration
- Check for installed plugins, skills, hooks, and MCP servers
- Identify which new features could improve the current workflow

## Step 5: Present Findings

Use this format:

```
## Claude Code Updates: vX.Y.Z → vA.B.C

### Relevant Changes
- [Feature]: [Why it matters for this project] → [Recommended action]

### Other Notable Changes
- [Feature]: [Brief description]

### Blog Announcements
- [Title] ([Date]): [1-2 sentence summary]

(or "No Claude Code blog announcements found." if none)

✅ Version tracking updated. Next nudge on next Claude Code update.
```

If no last-seen version existed, adjust the header:
```
## Claude Code Features: vA.B.C (initial review)
```

## Step 6: Update Tracking

After presenting findings, write the current version to `~/.claude/.claude-code-last-seen-version`:

```bash
mkdir -p ~/.claude && claude --version > ~/.claude/.claude-code-last-seen-version
```

This ensures the SessionStart hook won't nudge again until the next actual update.
