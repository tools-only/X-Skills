---
phase: hooks-guide-skill
task: 0
total_tasks: 14
status: in_progress
last_updated: 2026-02-27T20:03:18.401Z
---

<current_state>
Pre-execution state. Plan written and approved. Pre-flight passed. Ready to dispatch Task 1
(scaffold hooks-guide skill directory). No code written yet. Context ran out before first
subagent was dispatched.
</current_state>

<completed_work>

- Design doc written: `docs/plans/2026-02-27-hooks-guide-design.md`
- Implementation plan written: `docs/plans/2026-02-27-hooks-guide-implementation.md`
- Raw source docs fetched to `.claude/plan/plugin-creator-hooks/`:
  - `hooks-doc.md` — Claude Code hooks reference (1744 lines, from code.claude.com/docs/en/hooks.md)
  - `sub-agents-doc.md` — Claude Code sub-agents reference (814 lines, inline hooks/MCP section)
  - `github-copilot-hooks-doc.md` — GitHub Copilot coding agent hooks
  - `claude-hooks-doc.md` — (also present, likely duplicate of hooks-doc.md)
  - `vscode-copilot-hooks-doc.md` — VS Code Copilot hooks (bonus — not in original plan)
- Pre-flight verified: all 4+ raw docs confirmed present
- Skill not yet scaffolded

</completed_work>

<remaining_work>

All 14 tasks from `docs/plans/2026-02-27-hooks-guide-implementation.md` are pending:

- Task 1: Scaffold hooks-guide skill directory (init_skill.py)
- Task 2: Transform hooks-doc.md → claude-code.md (general-purpose agent, rwr:doc-to-skill)
- Task 3: Transform sub-agents-doc.md → inline-agent-hooks.md (general-purpose agent)
- Task 4: Transform github-copilot-hooks-doc.md → github-copilot.md (general-purpose agent)
- Task 5: Author hooks-cjs.md (contextual-ai-documentation-optimizer agent)
- Task 6: Author hooks-python.md (contextual-ai-documentation-optimizer agent)
- Task 7: Author common-schema.md and best-practices.md (two separate agent spawns)
- Task 8: Write platform-coverage.md (direct write — static registry)
- Task 9: Write SKILL.md (contextual-ai-documentation-optimizer agent)
- Task 10: Write fetch-and-transform-hooks-docs.sh (bash script, verify with bash -n)
- Task 11: Update agents/hook-creator.md (4 surgical changes)
- Task 12: Update skills/claude-hooks-reference-2026/SKILL.md (add hooks-guide to umbrella)
- Task 13: Full plugin validation (uv run plugin_validator.py plugins/plugin-creator)
- Task 14: Lint all modified files (uv run prek run --files ...)

</remaining_work>

<decisions_made>

- Option A chosen: single hooks-guide skill with references/ per platform/topic
- Scope: Claude Code + GitHub Copilot verified; Cursor/Windsurf/Amp attempted by fetch script
  (graceful skip if no hooks content found)
- No new agent — existing hook-creator agent updated in place
- claude-hooks-reference-2026 umbrella gains hooks-guide in its load list
- Fetch script uses `CLAUDECODE= claude -p` (not bare `claude -p`) — linter enforced this
- vscode-copilot-hooks-doc.md exists in .claude/plan/ — not in original plan but could be
  added as bonus Task 4b (transform to references/vscode-copilot.md)
- Execution method: subagent-driven-development (same session, fresh subagent per task)

</decisions_made>

<blockers>

- None. All raw docs present. Plan approved. Ready to execute.

</blockers>

<context>
Building a cross-platform hooks reference skill for the plugin-creator plugin. The goal is
an AI-facing skill that knows about hooks for Claude Code (plugin-level, settings-level,
AND inline-agent frontmatter), GitHub Copilot (.github/hooks/), and other platforms.

The key insight: existing skills only cover Claude Code + Node.js CJS. This adds Python
authoring, inline-agent hooks (new Claude Code sub-agents convention), GitHub Copilot,
and a self-updating fetch pipeline via rwr:doc-to-skill transformation.

The transform pipeline converts human-facing official docs into AI-optimised reference files.
</context>

<next_action>
Resume with: /plugin-creator:plugin-creator or load superpowers:subagent-driven-development

Then immediately dispatch Task 1 implementer subagent with full task text from the plan.

Task 1 text (exact):
  Run: uv run plugins/plugin-creator/skills/skill-creator/scripts/init_skill.py hooks-guide --path plugins/plugin-creator/skills
  Remove example files.
  Commit: "feat(hooks-guide): scaffold hooks-guide skill directory"

Note: vscode-copilot-hooks-doc.md also exists in .claude/plan/plugin-creator-hooks/ —
consider adding as Task 4b after Task 4 (same transform pattern, output to references/vscode-copilot.md).
</next_action>
