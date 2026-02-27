# P2 Validate Agent Browser For Web Automation

| Property | Value |
|----------|-------|
| **Name** | P2 Validate Agent Browser For Web Automation |
| **Repository** | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-validate-agent-browser-for-web-automation.md) (⭐ 20) |
| **Original Path** | `.claude/backlog/p2-validate-agent-browser-for-web-automation.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-23 |
| **Updated** | 2026-02-27 |
| **File Hash** | `300bab583efd6e72...` |

## Description

Test agent-browser (Playwright-based) on host with unrestricted network and Playwright browsers installed.\n**Validation steps**:\n- Install browsers: `npx playwright install`\n- Test: `npx agent-browser open https://code.claude.com/docs/en/skills`\n- Test: `npx agent-browser snapshot -i` (get element refs)\n- Test: `npx agent-browser get text body` (extract page text)\n- Verify snapshot/interact/re-snapshot workflow works\n- Document prerequisites for skill to function\n**Blocked on 2026-02-05**: Could not download Playwright browsers (DNS resolution failed, missing system libs)\n**Skill location**: `.claude/skills/agent-browser/SKILL.md`

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-validate-agent-browser-for-web-automation.md)*
