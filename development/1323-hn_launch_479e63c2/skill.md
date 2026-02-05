# Hacker News Launch: skill-gen

## Title (79 chars)
Show HN: skill-gen – A skill that helps you create agent skills

## Body

I built a skill that teaches Claude (and other AI agents) how to create skills.

If you've tried writing custom skills for your team's APIs, databases, or workflows, you know the challenge: skills aren't just docs. They're context engineering artifacts that need to balance token efficiency, progressive disclosure, and bundled resources.

skill-gen provides a complete workflow:

1. **Initialize** - Scaffold with proper SKILL.md template + example files
2. **Validate** - Check frontmatter, naming, description requirements
3. **Package** - Create distributable .skill files (zip format)

Includes Python scripts (init_skill.py, quick_validate.py, package_skill.py) and reference docs for workflow patterns and output templates.

Install:
```
npx skills add crafter-station/skills --skill skill-gen -g
```

Works with Claude Code, Cursor, Copilot, and 10+ agents.

I've used this to create 7+ production skills for Clerk (webhooks, orgs, setup) and validated the patterns across multiple projects.

This is the second skill in our marketplace after intent-layer (AGENTS.md infrastructure). Both are open source, battle-tested from real codebases.

Built in Peru 🇵🇪 as part of Crafter Station's mission to build tools for AI-first development.

GitHub: https://github.com/crafter-station/skills
Blog: https://railly.dev/blog/skill-gen
Demo: Run `/skill-gen` after installing

---

## Alternative Shorter Version (if HN prefers brevity)

I built a skill that teaches AI agents how to create skills.

If you've tried writing custom skills for Claude/Cursor/Copilot, you know it's more than docs—it's context engineering. skill-gen provides init, validate, and package scripts plus guided workflows.

```
npx skills add crafter-station/skills --skill skill-gen -g
```

Used to create 7+ production Clerk skills. Open source, works across 10+ agents.

GitHub: https://github.com/crafter-station/skills
Blog: https://railly.dev/blog/skill-gen

---

## Reddit r/ClaudeAI Version

Just released skill-gen - a skill that helps you create skills

Been building custom Claude skills for my team's APIs and workflows. Realized the hardest part isn't knowing what to build, it's structuring it properly.

skill-gen teaches Claude how to:
- Initialize skill scaffolds
- Validate frontmatter and naming
- Package for distribution
- Follow progressive disclosure patterns

Includes Python scripts and reference docs. Works with Claude Code, Cursor, Copilot.

Install: `npx skills add crafter-station/skills --skill skill-gen -g`

Used it to create 7 production skills for Clerk auth. Open source.

https://github.com/crafter-station/skills
https://railly.dev/blog/skill-gen

---

## Twitter/X Thread

🧵 Releasing skill-gen: a skill that helps you create skills

1/ If you've tried writing custom skills for Claude, you know it's not just docs—it's context engineering

Token efficiency, progressive disclosure, bundled resources... easy to mess up

2/ skill-gen provides the complete workflow:
- init_skill.py: scaffold proper structure
- quick_validate.py: catch errors early
- package_skill.py: ship .skill files

Plus guided workflows for patterns

3/ Works with Claude Code, Cursor, Copilot, 10+ agents

Install: npx skills add crafter-station/skills --skill skill-gen -g

Then run /skill-gen when creating a skill

4/ Used to create 7+ production Clerk skills (webhooks, orgs, setup)

Part of Crafter Station marketplace alongside intent-layer (AGENTS.md infrastructure)

5/ Both skills open source, battle-tested from real codebases

Built in Peru 🇵🇪

📦 github.com/crafter-station/skills
📝 railly.dev/blog/skill-gen

---

## LinkedIn Version

Excited to release skill-gen, a tool for creating high-quality agent skills.

As AI agents become standard in development workflows, the need for custom skills grows. But writing effective skills requires understanding context engineering principles that most developers haven't encountered.

skill-gen provides:
✅ Scaffolding with proper structure
✅ Validation for common errors
✅ Packaging for distribution
✅ Battle-tested patterns from production

Used internally to create 7+ Clerk authentication skills, now open source for the community.

Works with Claude Code, Cursor, GitHub Copilot, and 10+ other agents following the Agent Skills specification.

Part of Crafter Station's mission to build tools for AI-first development.

🔗 https://github.com/crafter-station/skills
📖 https://railly.dev/blog/skill-gen

#AI #DeveloperTools #OpenSource #ContextEngineering
