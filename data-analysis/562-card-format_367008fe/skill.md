# Card Format

Bugs, docs, and module cards all use this unified structure.

## Frontmatter

```yaml
---
title: <Short human-readable title>
link: <stable-slug-for-links>
type: <delta | doc | module>
path: <relative path, for modules/docs>
depth: <nesting level, 0 = root>
seams: <[A] architecture | [E] entry | [M] module | [S] state | [D] data>
ontological_relations:
  - relates_to: [[<primary-system-or-area>]]
  - affects: [[<component-or-module>]]
  - fixes: [[<bug-or-failure-mode>]] # for deltas
tags:
  - <area>
  - <symptom>
  - <tool-or-feature>
created_at: <ISO-8601 timestamp>
updated_at: <ISO-8601 timestamp>
uuid: <generated-uuid>
---
```

## Body Sections

- **Summary** - One paragraph, what and why, readable in a year
- **Context** - Optional, where it lives, how it surfaced
- **Root Cause** - For bugs: mechanism, not blame
- **Changes** - Bullet list of concrete changes
- **Behavioral Impact** - What users notice, what didn't change
- **Related Cards** - `[[links]]` to other cards

## Bug Rules

- Fix it, we don't care who caused it
- Document: what? where? when? why?
- No shims, fix at root
- Document how we missed it and how to prevent it
