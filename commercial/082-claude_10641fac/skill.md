# Claude Skills

**Repository**: https://github.com/jezweb/claude-skills
**Owner**: Jeremy Dawes (Jez) | Jezweb

Production workflow skills for Claude Code CLI. Each skill guides Claude through a recipe to produce tangible output — not knowledge dumps, but working deliverables.

## Philosophy

- Every skill must produce visible output (files, configurations, deployable projects)
- "The context window is a public good" — only include what Claude doesn't already know
- Follows the official Claude Code plugin spec

## Directory Structure

```
claude-skills/
├── plugins/                                # 9 installable plugins (34 skills)
│   ├── cloudflare/                         # Cloudflare Workers, Hono, D1/Drizzle, Vite, TanStack Start
│   │   └── skills/
│   │       ├── cloudflare-worker-builder/
│   │       ├── vite-flare-starter/
│   │       ├── tanstack-start/
│   │       ├── hono-api-scaffolder/
│   │       └── d1-drizzle-schema/
│   ├── web-design/                         # Web design methodology, patterns, SEO
│   │   └── skills/
│   │       ├── web-design-methodology/
│   │       ├── web-design-patterns/
│   │       └── seo-local-business/
│   ├── frontend/                           # Tailwind v4 + shadcn/ui
│   │   └── skills/
│   │       ├── tailwind-theme-builder/
│   │       └── shadcn-ui/
│   ├── design-assets/                      # Colour palettes, favicons, icons, image gen/processing
│   │   └── skills/
│   │       ├── color-palette/
│   │       ├── favicon-gen/
│   │       ├── icon-set-generator/
│   │       ├── gemini-image-gen/
│   │       └── image-processing/
│   ├── integrations/                       # Google, ElevenLabs, MCP
│   │   └── skills/
│   │       ├── google-chat-messages/
│   │       ├── google-apps-script/
│   │       ├── elevenlabs-agents/
│   │       └── mcp-builder/
│   ├── dev-tools/                          # Skill creation, context, sessions, releases
│   │   └── skills/
│   │       ├── skill-creator/
│   │       ├── project-kickoff/
│   │       ├── context-manager/
│   │       ├── dev-session/
│   │       ├── github-release/
│   │       ├── gemini-peer-review/
│   │       ├── gemini-guide/
│   │       ├── claude-capabilities/
│   │       ├── ux-audit/
│   │       └── responsiveness-check/
│   ├── shopify/                            # Shopify store management
│   │   └── skills/
│   │       ├── shopify-setup/
│   │       ├── shopify-products/
│   │       └── shopify-content/
│   ├── wordpress/                          # WordPress content & Elementor
│   │   └── skills/
│   │       ├── wordpress-setup/
│   │       ├── wordpress-content/
│   │       └── wordpress-elementor/
│   └── writing/                            # Australian business English
│       └── skills/
│           └── aussie-business-english/
├── .claude-plugin/                         # Marketplace + plugin config
│   ├── marketplace.json
│   └── plugin.json
├── CLAUDE.md                               # This file
├── README.md                               # Public-facing overview
└── LICENSE                                 # MIT
```

## Plugin Anatomy (Anthropic Spec)

Each plugin contains one or more skills, auto-discovered from `skills/`:

```
plugin-name/
├── .claude-plugin/
│   └── plugin.json        # name, description, author
└── skills/
    └── skill-name/
        ├── SKILL.md       # Frontmatter + instructions, under 500 lines
        ├── scripts/       # Executable code (run directly)
        ├── references/    # Docs loaded on demand by Claude
        └── assets/        # Files used in output (templates, images)
```

## Adding a New Plugin

1. Create the plugin directory:
   ```bash
   mkdir -p plugins/my-plugin/{.claude-plugin,skills}
   ```

2. Create `.claude-plugin/plugin.json`:
   ```json
   {
     "name": "my-plugin",
     "description": "What this plugin does.",
     "author": { "name": "Jeremy Dawes / Jezweb", "email": "jeremy@jezweb.net" }
   }
   ```

3. Add skills inside `plugins/my-plugin/skills/` (each with SKILL.md)

4. Add an entry to `.claude-plugin/marketplace.json`:
   ```json
   { "name": "my-plugin", "description": "...", "source": "./plugins/my-plugin", "category": "development" }
   ```

5. Update the directory tree in this file and the table in README.md

**Categories**: `development`, `design`, `productivity`, `testing`, `security`, `database`, `monitoring`, `deployment`

## Creating a Skill

Use the skill-creator skill:

```bash
python3 plugins/dev-tools/skills/skill-creator/scripts/init_skill.py my-skill --path plugins/my-plugin/skills/
```

Or ask Claude: "Create a new skill for [use case]"

Key principle: **every skill must produce something.** If it's just reference material Claude already knows, it doesn't earn a place here.

## Installing Plugins

```bash
# Add marketplace (one-time)
/plugin marketplace add jezweb/claude-skills

# Install individual plugins
/plugin install cloudflare@jezweb-skills
/plugin install dev-tools@jezweb-skills
/plugin install frontend@jezweb-skills

# Local dev (loads a single plugin without install)
claude --plugin-dir ./plugins/cloudflare
```

After installing, restart Claude Code to load new plugins.

## Quality Bar

Before committing a skill:
- [ ] SKILL.md has valid YAML frontmatter (name + description)
- [ ] Under 500 lines
- [ ] Produces tangible output (not just reference material)
- [ ] Passes validation: `python3 plugins/dev-tools/skills/skill-creator/scripts/quick_validate.py plugins/category/skills/my-skill`
- [ ] Tested by actually using it on a real task

## Git History

All 105 skills from the v1 era are preserved:
- Tag `v1-final` — the complete 105-skill collection
- Branch `archive/low-priority-skills` — 13 previously archived skills
- Full git history available via `git log v1-final`
