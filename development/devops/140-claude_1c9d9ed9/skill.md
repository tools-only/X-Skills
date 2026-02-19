# Claude Skills

**Repository**: https://github.com/jezweb/claude-skills
**Owner**: Jeremy Dawes (Jez) | Jezweb

Production workflow skills for Claude Code CLI. Each skill guides Claude through a recipe to produce tangible output — not knowledge dumps, but working deliverables.

## Philosophy

- Every skill must produce visible output (files, configurations, deployable projects)
- "The context window is a public good" — only include what Claude doesn't already know
- Follow the Agent Skills spec: https://agentskills.io/specification

## Directory Structure

```
claude-skills/
├── skills/                         # All skills (24)
│   ├── skill-creator/              # Foundation: create new skills (Anthropic official)
│   ├── cloudflare-worker-builder/  # Scaffold Cloudflare Worker projects
│   ├── vite-flare-starter/         # Full-stack Cloudflare app from starter template
│   ├── hono-api-scaffolder/        # Scaffold Hono API routes, middleware, endpoint docs
│   ├── d1-drizzle-schema/          # Drizzle ORM schemas for Cloudflare D1
│   ├── tailwind-theme-builder/     # Set up Tailwind v4 + shadcn/ui themes
│   ├── shadcn-ui/                  # shadcn/ui component installation, recipes, customisation
│   ├── color-palette/              # Generate colour palettes from brand hex
│   ├── favicon-gen/                # Generate favicon packages
│   ├── icon-set-generator/         # Generate custom SVG icon sets
│   ├── web-design-methodology/     # Universal web design implementation (BEM, responsive, a11y)
│   ├── web-design-patterns/        # Heroes, cards, CTAs, trust signals, testimonials
│   ├── seo-local-business/         # SEO setup for local businesses (JSON-LD, meta, sitemap)
│   ├── google-chat-messages/       # Google Chat webhooks (text, cards, threads)
│   ├── google-apps-script/         # Google Sheets Apps Script automation
│   ├── elevenlabs-agents/          # Build ElevenLabs voice agents
│   ├── mcp-builder/                # Build MCP servers with FastMCP
│   ├── context-manager/             # Audit project context: memory, docs, footprint, overlap
│   ├── dev-session/                # Multi-session progress tracking and handoff
│   ├── ux-audit/                   # UX walkthroughs and QA sweeps via browser automation
│   ├── github-release/             # Sanitize and publish GitHub releases
│   ├── gemini-peer-review/         # Second opinion from Gemini on code/architecture
│   ├── claude-capabilities/        # Current Claude AI & Code capabilities reference
│   └── aussie-business-english/    # Australian business English writing style
├── CLAUDE.md                       # This file
├── README.md                       # Public-facing overview
└── LICENSE                         # MIT
```

## Skill Anatomy (Anthropic Spec)

```
skill-name/
├── SKILL.md (required)     # Frontmatter + instructions, under 500 lines
├── scripts/                # Executable code (run directly)
├── references/             # Docs loaded on demand by Claude
└── assets/                 # Files used in output (templates, images)
```

No README.md, no CHANGELOG.md, no rules/ — just what the AI agent needs.

## Creating a Skill

Use the skill-creator skill:

```bash
python3 skills/skill-creator/scripts/init_skill.py my-skill --path skills/
```

Or ask Claude: "Create a new skill for [use case]"

Key principle: **every skill must produce something.** If it's just reference material Claude already knows, it doesn't earn a place here.

## Installing Skills

```bash
# Add marketplace (one-time)
/plugin marketplace add jezweb/claude-skills

# Install all skills
/plugin install all@jezweb-skills

# Or by category: design, cloudflare, frontend, web, integrations, ai, mcp, development, testing, writing
/plugin install design@jezweb-skills

# Local dev (loads without install)
claude --plugin-dir ./skills/cloudflare-worker-builder
```

After installing, restart Claude Code to load new skills.

## Quality Bar

Before committing a skill:
- [ ] SKILL.md has valid YAML frontmatter (name + description)
- [ ] Under 500 lines
- [ ] Produces tangible output (not just reference material)
- [ ] Passes validation: `python3 skills/skill-creator/scripts/quick_validate.py skills/my-skill`
- [ ] Tested by actually using it on a real task

## Git History

All 105 skills from the v1 era are preserved:
- Tag `v1-final` — the complete 105-skill collection
- Branch `archive/low-priority-skills` — 13 previously archived skills
- Full git history available via `git log v1-final`
