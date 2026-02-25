---
description: "Documentation and authoring workflow router — audit docs vs code drift, sync docs after code changes, optimize prompts and SKILL.md files for AI consumption, validate GLFM and Markdown, summarize files/URLs/images with fidelity enforcement. Routes each task to the right specialist agent automatically."
allowed-tools: Read, Grep, Glob, Bash, Task, Write, Edit
---
# The Rewrite Room

Routes documentation, authoring, and optimization tasks to the correct specialist agents. Does not rewrite source agents or skills — orchestrates them. Governs authoring, docs, prompts, and summaries — not product code.

## Quick Start

```text
/rwr:audit "check if kaizen plugin docs match the code"
/rwr:optimize "plugins/plugin-creator/skills/add-doc-updater/SKILL.md"
/rwr:author "summarize plugins/summarizer/skills/summarizer/SKILL.md"
```

## Command Reference

| Command | Entry Agent | Use When |
|---------|-------------|----------|
| `/rwr:audit <task>` | rewrite-room-auditor | Docs vs code drift, doc sync after changes, freshness tracking |
| `/rwr:optimize <file>` | rewrite-room-optimizer | CLAUDE.md, SKILL.md, agent .md improvement |
| `/rwr:author <task>` | rewrite-room-author | User-facing docs, GLFM validation, summarization |

Each command loads the corresponding workflow file and follows its numbered steps.

## Workflow Index

```mermaid
flowchart TD
    User([User invokes /rwr:*]) --> Q{Which command?}
    Q -->|/rwr:audit| Audit[rewrite-room-auditor\nLoads: plugins/the-rewrite-room/the-rewrite-room/workflows/audit.md]
    Q -->|/rwr:optimize| Opt[rewrite-room-optimizer\nLoads: plugins/the-rewrite-room/the-rewrite-room/workflows/optimize.md]
    Q -->|/rwr:author| Auth[rewrite-room-author\nLoads: plugins/the-rewrite-room/the-rewrite-room/workflows/author.md]
    Audit --> A1[development-harness:doc-drift-auditor]
    Audit --> A2[development-harness:service-docs-maintainer]
    Audit --> A3[doc-freshness-guardian]
    Opt --> O1[plugin-creator:contextual-ai-documentation-optimizer]
    Opt --> O2[plugin-creator:subagent-refactorer]
    Auth --> B1[gitlab-docs-expert]
    Auth --> B2[documentation-expert]
    Auth --> B3[summarizer:file-summarizer / summarizer:url-summarizer / summarizer:image-summarizer]
```

## Workflow Files

Each command agent loads the corresponding workflow file at runtime:

- Audit workflow: `plugins/the-rewrite-room/the-rewrite-room/workflows/audit.md`
- Optimize workflow: `plugins/the-rewrite-room/the-rewrite-room/workflows/optimize.md`
- Author workflow: `plugins/the-rewrite-room/the-rewrite-room/workflows/author.md`

Workflow files contain numbered steps, conditional branching, explicit agent spawn instructions, structured return handling, and output contracts.

## Adding New Workflows

1. Create workflow file in `plugins/the-rewrite-room/the-rewrite-room/workflows/<name>.md`
2. Create command file in the plugin commands directory referencing the workflow
3. Add agent in `plugins/the-rewrite-room/agents/<name>.md` if a new routing agent is needed
4. Register new agent path in `.claude-plugin/plugin.json` agents array

## Source Components

This plugin routes to these specialist agents and scripts (not copied — referenced by path):

**Audit agents:**

- `plugins/development-harness/agents/doc-drift-auditor.md` — evidence-based drift audit with file:line citations
- `plugins/development-harness/agents/service-docs-maintainer.md` — post-implementation doc sync via git diff
- `/home/ubuntulinuxqa2/.claude/agents/doc-freshness-guardian.md` — freshness headers and staleness alerts

**Optimize agents:**

- `plugins/plugin-creator/agents/contextual-ai-documentation-optimizer.md` — RT-ICA + CoVe prompt optimization with token impact reporting
- `plugins/plugin-creator/agents/subagent-refactorer.md` — Anthropic official best practices refactoring with mandatory research phase

**Author agents:**

- `gitlab-docs-expert` — GitLab Wiki, MR descriptions, GitLab README authoring
- `documentation-expert` — general README, tutorials, API docs, user-facing docs

**Summarizer agents:**

- `plugins/summarizer/agents/file-summarizer.md` — file content summarization with fidelity enforcement
- `plugins/summarizer/agents/url-summarizer.md` — URL content summarization
- `plugins/summarizer/agents/image-summarizer.md` — image/screenshot description

**Validation scripts:**

- `plugins/gitlab-skill/skills/gitlab-skill/scripts/validate_glfm.py` — GitLab Flavored Markdown validation via GitLab API
- `plugins/plugin-creator/scripts/validate_frontmatter.py` — YAML frontmatter schema validation

**Reference files consulted by workflows:**

- `plugins/summarizer/skills/summarizer/references/fidelity-rules.md` — summarizer fidelity rules
- `plugins/gitlab-skill/skills/gitlab-skill/references/glfm-syntax.md` — GLFM syntax reference
- `plugins/prompt-optimization-claude-45/skills/prompt-optimization-claude-45/SKILL.md` — prompt optimization principles
