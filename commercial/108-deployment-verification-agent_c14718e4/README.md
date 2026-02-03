# Deployment Verification Agent

| Property | Value |
|----------|-------|
| **Name** | Deployment Verification Agent |
| **Repository** | [davekilleen/Dex](https://raw.githubusercontent.com/davekilleen/Dex/main/.claude/plugins/compound-engineering/agents/review/deployment-verification-agent.md) (⭐ 58) |
| **Original Path** | `.claude/plugins/compound-engineering/agents/review/deployment-verification-agent.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-01-30 |
| **Updated** | 2026-01-30 |
| **File Hash** | `c14718e4eb4be95c...` |

## Description

Use this agent when a PR touches production data, migrations, or any behavior that could silently discard or duplicate records. Produces a concrete pre/post-deploy checklist with SQL verification queries, rollback procedures, and monitoring plans. Essential for risky data changes where you need a Go/No-Go decision. <example>Context: The user has a PR that modifies how emails are classified. user: \

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [davekilleen/Dex](https://raw.githubusercontent.com/davekilleen/Dex/main/.claude/plugins/compound-engineering/agents/review/deployment-verification-agent.md)*
