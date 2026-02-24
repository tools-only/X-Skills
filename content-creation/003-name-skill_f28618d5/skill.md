---
name: backlog
description: "Single interface for BACKLOG.md and GitHub Issues. All backlog CRUD goes through the backlog script — no direct edits. Use when create-backlog-item or work-backlog-item invoke it, or when syncing P0/P1 to GitHub."
---
# Backlog Script

The backlog script is the **sole interface** for `.claude/BACKLOG.md` and GitHub Issues. Skills and agents invoke it; no direct Write/Grep edits to BACKLOG.md.

## Invocation

```bash
uv run .claude/skills/backlog/scripts/backlog.py <subcommand> [options]
```

## Subcommands

| Subcommand | Purpose |
|------------|---------|
| `add` | Create item in BACKLOG; create GitHub issue if P0/P1 |
| `list` | List items (for interactive browser) |
| `sync` | Create issues for P0/P1 items missing them |
| `close` | Mark DONE, close issue (requires `--checklist-pass`) |
| `resolve` | Mark RESOLVED, close issue |
| `update` | Add Plan, set status:in-progress, or create issue |

## Usage

### add

```bash
backlog add --title "Item title" --priority P1 --description "Description" \
  [--source "Source"] [--type Feature] [--research-first "questions"] [--create-issue]
```

### list

```bash
backlog list [--format text|json]
```

### sync

```bash
backlog sync [--dry-run]
```

### close

Skill must verify checklist first, then:

```bash
backlog close "<title or #N>" --plan <path> --checklist-pass
```

### resolve

```bash
backlog resolve "<title or #N>" --reason "reason text"
```

### update

```bash
backlog update "<title or #N>" [--plan <path>] [--status in-progress] [--create-issue]
```

## Environment

- `GITHUB_TOKEN` — Required for issue operations (add, sync, close, resolve, update --create-issue)

## Integration

- **create-backlog-item** — invokes `backlog add` instead of writing BACKLOG.md
- **work-backlog-item** — invokes `backlog list`, `backlog close`, `backlog resolve`, `backlog update`
- **GitHub Action** — invokes `backlog sync` on BACKLOG.md change
