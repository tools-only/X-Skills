---
description: Manage development sessions — start, checkpoint, wrap, resume
argument-hint: "[start | checkpoint | wrap | resume]"
allowed-tools: "*"
---

Load the `dev-session` skill and execute based on $ARGUMENTS:

| Argument | Action |
|----------|--------|
| `start` | Create SESSION.md, begin tracking progress |
| `checkpoint` | WIP commit + update SESSION.md with current progress |
| `wrap` | Final SESSION.md update, capture learnings, commit |
| `resume` | Read SESSION.md and continue where we left off |
| (none) | Check if SESSION.md exists — resume if yes, start if no |
