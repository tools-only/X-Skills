# X-Skills

A curated collection of **442 AI-powered skills** organized into 14 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (3 skills)
- **Automation/Workflow** (13 skills)
- **Commercial** (18 skills)
- **Communication** (12 skills)
- **Content Creation** (44 skills)
- **Daily Assistant** (5 skills)
- **Data Analysis** (54 skills)
- **Development** (104 skills)
- **Development/Devops** (102 skills)
- **Development/Testing** (5 skills)
- **Development/Tools** (55 skills)
- **Investment** (1 skill)
- **Productivity** (3 skills)
- **Research** (23 skills)

## Patches - Curated Skill Bundles

Pre-configured skill bundles for common AI agent use cases.

### Quick Start with Claude Code

```bash
# Install X-Skills CLI
pip install -e skillflow_repos/X-Skills

# List available patches
xskills patches list

# Install a patch (e.g., research-agent)
xskills patch install research-agent

# Skills are now available in Claude Code!
```

### Available Patches

| Patch | Skills | Categories | Install |
|-------|--------|------------|---------|
| [Automation Agent](automation-agent/) | 0 | automation | `xskills patch install automation-agent` |
| [Communication Agent](communication-agent/) | 40 | communication | `xskills patch install communication-agent` |
| [Content Creator](content-creator/) | 0 | content-creation | `xskills patch install content-creator` |
| [Data Analyst](data-analyst/) | 0 | data-analysis | `xskills patch install data-analyst` |
| [DevOps Engineer](devops-engineer/) | 60 | development | `xskills patch install devops-engineer` |
| [Productivity Assistant](productivity-assistant/) | 50 | productivity, daily-assistant | `xskills patch install productivity-assistant` |
| [Python Developer](python-dev/) | 0 | development | `xskills patch install python-dev` |
| [Research Agent](research-agent/) | 0 | research | `xskills patch install research-agent` |
| [Web Development Agent](web-dev-agent/) | 0 | development | `xskills patch install web-dev-agent` |

### Installation Methods

**Method 1: Using X-Skills CLI (Recommended)**

```bash
# Install to Claude Code skills directory
xskills patch install <patch-id>

# Skills install to: ~/.claude/skills/patch-<patch-id>/
# Immediately available in Claude Code
```

**Method 2: Direct Installation**

```bash
# Clone this repository
git clone https://github.com/tools-only/X-Skills.git
cd X-Skills

# Install a patch
python -m src.patch_installer install <patch-id>
```

### Usage Examples

```bash
# Research Agent - Academic papers and literature review
xskills patch install research-agent

# Web Development Agent - Full-stack web development
xskills patch install web-dev-agent

# Content Creator - Writing and content generation
xskills patch install content-creator

# Data Analyst - Data analysis and visualization
xskills patch install data-analyst

# View installed patches
xskills patches list

# Uninstall a patch
xskills patch uninstall <patch-id>
```

### Browse Skills

```bash
# Browse all skills
xskills browse

# Browse by category
xskills browse --category research

# Search skills
xskills search "web development"
```

**[View all patches and documentation](INDEX.md)**


## Claude Code Integration

### Quick Start

X-Skills integrates seamlessly with Claude Code through patches - curated skill bundles for common use cases.

```bash
# 1. Clone this repository
git clone https://github.com/tools-only/X-Skills.git
cd X-Skills

# 2. Install the X-Skills CLI
pip install -e .

# 3. List available patches
xskills patches list

# 4. Install a patch
xskills patch install research-agent

# 5. Skills are now available in Claude Code!
```

### Available Commands

```bash
# Patch Management
xskills patches list              # List all patches
xskills patch install <patch>     # Install a patch
xskills patch uninstall <patch>   # Uninstall a patch

# Skill Browsing
xskills browse                    # Browse all skills
xskills browse --category research # Browse by category
xskills search "web dev"          # Search skills

# System Status
xskills status                    # Show system status
```

### How It Works

When you install a patch:

1. Skills are symlinked to `~/.claude/skills/patch-<patch-id>/`
2. Skills become immediately available in Claude Code
3. No manual configuration needed

### Example Workflow

```bash
# Install research capabilities
xskills patch install research-agent

# Browse what was installed
ls ~/.claude/skills/patch-research-agent/

# Use in Claude Code - skills are now available!
```

## Skills Directory


### Automation/Scripting (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agent Calibration](automation/scripting/085-agent-calibration_a87c080f/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/agent-calibration.md) | ⭐ 61 | `automation` |
| [Skill](automation/scripting/003-name-skill_570982e2/) | [ALBEDO-TABAI/lets-go-rss](https://raw.githubusercontent.com/ALBEDO-TABAI/lets-go-rss/main/SKILL.md) | ⭐ 31 | `automation` |
| [Skill](automation/scripting/003-name-skill_28ff4f04/) | [artwist-polyakov/polyakov-claude-skills](https://raw.githubusercontent.com/artwist-polyakov/polyakov-claude-skills/main/plugins/yandex-search-api/skills/yandex-search-api/SKILL.md) | ⭐ 38 | `automation` |

### Automation/Workflow (13 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Vfunc Sig](automation/workflow/136-vfunc_sig_9dfecfa2/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.serena/memories/vfunc_sig.md) | ⭐ 19 | `automation` |
| [Skill](automation/workflow/002-name-skill_2854fe0b/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-deploy/SKILL.md) | ⭐ 80 | `automation` |
| [Cli Guide](automation/workflow/137-cli-guide_69513c0f/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-ai-agentscript/references/cli-guide.md) | ⭐ 80 | `automation` |
| [Export Import Architecture](automation/workflow/148-export-import-architecture_1fc68b51/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/export-import-architecture.md) | ⭐ 3.3k | `automation` |
| [Export Import Tutorial](automation/workflow/149-export-import-tutorial_372778e0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import-tutorial.md) | ⭐ 3.3k | `automation` |
| [Export Import](automation/workflow/150-export-import_af809418/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import.md) | ⭐ 3.3k | `automation` |
| [Mcp Cli](automation/workflow/151-mcp-cli_dbea27bc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/mcp-cli.md) | ⭐ 3.3k | `automation` |
| [Libreoffice Server](automation/workflow/152-libreoffice-server_2c0ac43b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/libreoffice-server.md) | ⭐ 3.3k | `automation` |
| [Skill](automation/workflow/002-name-skill_15b7286f/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/gas-fee-optimizer/skills/optimizing-gas-fees/SKILL.md) | ⭐ 1.4k | `automation` |
| [Readme.Ja](automation/workflow/040-readmeja_c8bfa185/) | [japan1988/multi-agent-mediation](https://raw.githubusercontent.com/japan1988/multi-agent-mediation/main/README.ja.md) | ⭐ 28 | `automation` |
| [Custom Workflows](automation/workflow/133-custom-workflows_4c3bac3c/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/advanced/custom-workflows.md) | 🔥 9.7k | `automation` |
| [04 Next Steps](automation/workflow/134-04-next-steps_12451d6e/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/04-next-steps.md) | 🔥 9.7k | `automation` |
| [04 Packaging](automation/workflow/135-04-packaging_49de9086/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/04-packaging.md) | 🔥 9.7k | `automation` |

### Commercial (18 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Express Adapter](commercial/387-express-adapter_6db5888c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/express-adapter.md) | ⭐ 102 | `commercial` |
| [Passkey](commercial/388-passkey_b16def64/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/passkey.md) | ⭐ 102 | `commercial` |
| [Routing Patterns](commercial/389-routing-patterns_959a3406/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/routing-patterns.md) | ⭐ 102 | `commercial` |
| [Caching Strategies](commercial/390-caching-strategies_3decca32/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/caching-strategies.md) | ⭐ 102 | `products` |
| [Server Components](commercial/179-server-components_249b33ee/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/server-components.md) | ⭐ 102 | `commercial` |
| [Advanced](commercial/002-advanced_87d29bd8/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/advanced.md) | ⭐ 102 | `commercial` |
| [Migration 0.7.0](commercial/391-migration-070_0d56b710/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/MIGRATION-0.7.0.md) | ⭐ 3.3k | `commercial` |
| [Kubernetes](commercial/344-kubernetes_0493cab9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/kubernetes.md) | ⭐ 3.3k | `commercial` |
| [Developer Workstation](commercial/060-developer-workstation_2f62c709/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/developer-workstation.md) | ⭐ 3.3k | `commercial` |
| [Oauth Troubleshooting](commercial/392-oauth-troubleshooting_2786d167/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/oauth-troubleshooting.md) | ⭐ 3.3k | `commercial` |
| [Supported Databases](commercial/393-supported-databases_d8a56a4d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/supported-databases.md) | ⭐ 3.3k | `commercial` |
| [Index](commercial/102-index_833de47b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/index.md) | ⭐ 3.3k | `commercial` |
| [Grpc Services](commercial/394-grpc-services_a969a94f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/grpc-services.md) | ⭐ 3.3k | `commercial` |
| [Query Param Auth](commercial/395-query-param-auth_e000b191/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/query-param-auth.md) | ⭐ 3.3k | `commercial` |
| [Cards](commercial/378-cards_6c04640b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/cards.md) | ⭐ 117 | `commercial` |
| [Skill](commercial/210-name-skill_1b74b4fe/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/stored-procedure-generator/skills/generating-stored-procedures/SKILL.md) | ⭐ 1.4k | `commercial` |
| [Permissions](commercial/permissions_1f64e2d9/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/permissions.md) | ⭐ 36 | `commercial` |
| [Multi User](commercial/multi-user_505b1d3c/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/examples/multi-user.md) | ⭐ 36 | `commercial` |

### Communication (12 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Schema](communication/267-schema_737de31c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/SCHEMA.md) | ⭐ 102 | `communication` |
| [Nestjs Setup](communication/268-nestjs-setup_476bb950/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/nestjs-setup.md) | ⭐ 102 | `communication` |
| [Nextjs Setup](communication/269-nextjs-setup_c674b3e1/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/nextjs-setup.md) | ⭐ 102 | `communication` |
| [Server Actions](communication/270-server-actions_4b85acf1/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/server-actions.md) | ⭐ 102 | `communication` |
| [Database Adapter](communication/271-database-adapter_0b083935/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-authentication/references/database-adapter.md) | ⭐ 102 | `communication` |
| [Rfc9728 Compliance](communication/272-rfc9728-compliance_0ee721f6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/rfc9728-compliance.md) | ⭐ 3.3k | `communication` |
| [Password Management](communication/273-password-management_bc5d8451/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/password-management.md) | ⭐ 3.3k | `communication` |
| [Inputs](communication/269-inputs_0e8fceb4/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/inputs.md) | ⭐ 117 | `communication` |
| [Toasts](communication/270-toasts_88d5e31e/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/toasts.md) | ⭐ 117 | `communication` |
| [Skill](communication/127-name-skill_5efa5b4b/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/SKILL.md) | ⭐ 61 | `team:support` |
| [Import Export](communication/251-import-export_9220633d/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/import-export.md) | ⭐ 61 | `communication` |
| [Sdk Usage](communication/252-sdk-usage_484ce649/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/sdk-usage.md) | ⭐ 61 | `tenant:acme` |

### Content Creation (44 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](content-creation/049-name-skill_df94dda4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/drizzle-orm-patterns/SKILL.md) | ⭐ 102 | `drizzle` `orm` `database` |
| [Skill](content-creation/049-name-skill_feaebcdb/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/SKILL.md) | ⭐ 102 | `nextjs` `next.js` `app-router` |
| [Skill](content-creation/049-name-skill_2cfa63f4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-data-fetching/SKILL.md) | ⭐ 102 | `posts` |
| [App Router Fundamentals](content-creation/405-app-router-fundamentals_5db2d98d/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/app-router-fundamentals.md) | ⭐ 102 | `content creation` |
| [Metadata Api](content-creation/406-metadata-api_464de7c9/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/metadata-api.md) | ⭐ 102 | `content creation` |
| [Nextjs16 Migration](content-creation/407-nextjs16-migration_807e73c6/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/nextjs16-migration.md) | ⭐ 102 | `content creation` |
| [React Query](content-creation/408-react-query_66bf81c6/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-data-fetching/references/REACT-QUERY.md) | ⭐ 102 | `content creation` |
| [Core Web Vitals](content-creation/409-core-web-vitals_c726a13b/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/core-web-vitals.md) | ⭐ 102 | `content creation` |
| [Image Optimization](content-creation/410-image-optimization_8a4da78f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/image-optimization.md) | ⭐ 102 | `content creation` |
| [Metadata Seo](content-creation/411-metadata-seo_eea42aa7/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/metadata-seo.md) | ⭐ 102 | `content creation` |
| [Nextjs 16 Patterns](content-creation/412-nextjs-16-patterns_1faa29b3/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/nextjs-16-patterns.md) | ⭐ 102 | `content creation` |
| [Nextjs Config](content-creation/413-nextjs-config_d4dbc259/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/nextjs-config.md) | ⭐ 102 | `content creation` |
| [Troubleshooting](content-creation/110-troubleshooting_7c93e51b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/troubleshooting.md) | ⭐ 3.3k | `content creation` |
| [Passthrough](content-creation/414-passthrough_6dcf2b83/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/passthrough.md) | ⭐ 3.3k | `content creation` |
| [Tool Annotations](content-creation/415-tool-annotations_9f3db504/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/tool-annotations.md) | ⭐ 3.3k | `content creation` |
| [038 Experimental Rust Transport Backend](content-creation/416-038-experimental-rust-transport-backend_0eab4f20/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/038-experimental-rust-transport-backend.md) | ⭐ 3.3k | `content creation` |
| [Index](content-creation/019-index_c4edf73d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/kit/index.md) | ⭐ 3.3k | `content creation` |
| [Index](content-creation/019-index_1e66c9e7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/press/index.md) | ⭐ 3.3k | `content creation` |
| [Chunker Server](content-creation/417-chunker-server_334bf908/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/chunker-server.md) | ⭐ 3.3k | `content creation` |
| [Eval Server](content-creation/418-eval-server_35d409dc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/eval-server.md) | ⭐ 3.3k | `content creation` |
| [Url To Markdown Server](content-creation/419-url-to-markdown-server_55cb28c7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/url-to-markdown-server.md) | ⭐ 3.3k | `content creation` |
| [Skill](content-creation/049-name-skill_6b88bec6/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/SKILL.md) | ⭐ 117 | `package` |
| [Skill](content-creation/049-name-skill_e7bf1e39/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/SKILL.md) | ⭐ 117 | `content creation` |
| [Content Guidelines](content-creation/360-content-guidelines_b4f2960a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/references/content-guidelines.md) | ⭐ 117 | `content creation` |
| [Shiny Formatting](content-creation/361-shiny-formatting_30dd05e6/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/references/shiny-formatting.md) | ⭐ 117 | `content creation` |
| [Conversion Blogdown](content-creation/362-conversion-blogdown_038f094b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-blogdown.md) | ⭐ 117 | `data` |
| [Conversion Distill](content-creation/363-conversion-distill_69039bbc/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-distill.md) | ⭐ 117 | `content creation` |
| [Features](content-creation/109-features_40be1130/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/docs/FEATURES.md) | ⭐ 829 | `content creation` |
| [Skill](content-creation/049-name-skill_aa320b97/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/scientific-writing/SKILL.md) | ⭐ 829 | `content creation` |
| [V1.81.12](content-creation/351-v18112_b76af1d8/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/release_notes/v1.81.12.md) | 🔥 36.4k | `content creation` |
| [Skill](content-creation/049-name-skill_232ae8b1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/work-backlog-item/SKILL.md) | ⭐ 18 | `content creation` |
| [Whats New](content-creation/354-whats-new_6a65e564/) | [NTCoding/claude-skillz](https://raw.githubusercontent.com/NTCoding/claude-skillz/main/claude-code-updates/commands/whats-new.md) | ⭐ 239 | `content creation` |
| [Skill](content-creation/049-name-skill_a0bce514/) | [jim60105/copilot-prompt](https://raw.githubusercontent.com/jim60105/copilot-prompt/master/skills/create-blog-post/SKILL.md) | ⭐ 17 | `content creation` |
| [Writing Guidelines](content-creation/355-writing-guidelines_49868b1b/) | [jim60105/copilot-prompt](https://raw.githubusercontent.com/jim60105/copilot-prompt/master/skills/create-blog-post/references/writing-guidelines.md) | ⭐ 17 | `content creation` |
| [P5 Guide](content-creation/353-p5-guide_f110bc04/) | [dmccreary/claude-skills](https://raw.githubusercontent.com/dmccreary/claude-skills/main/skills/microsim-generator/references/p5-guide.md) | ⭐ 49 | `content creation` |
| [02 Scraping](content-creation/354-02-scraping_eb8fe18d/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/02-scraping.md) | 🔥 9.7k | `content creation` |
| [Usage](content-creation/042-usage_cd9a1109/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/archive/legacy/USAGE.md) | 🔥 9.7k | `content creation` |
| [Index](content-creation/index_a0e2ee39/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/index.md) | ⭐ 36 | `content creation` |
| [Toolsets](content-creation/toolsets_5e962330/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/toolsets.md) | ⭐ 36 | `content creation` |
| [Console Toolset](content-creation/console-toolset_a2a7b2fa/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/console-toolset.md) | ⭐ 36 | `content creation` |
| [Backends](content-creation/backends_8631f33a/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/backends.md) | ⭐ 36 | `content creation` |
| [Docker](content-creation/docker_94863f0f/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/docker.md) | ⭐ 36 | `content creation` |
| [Permissions](content-creation/permissions_f36bcca1/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/permissions.md) | ⭐ 36 | `content creation` |
| [Types](content-creation/types_57d018af/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/types.md) | ⭐ 36 | `content creation` |

### Daily Assistant (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Package Configs](daily-assistant/294-package-configs_66bf674a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/package-configs.md) | ⭐ 102 | `daily assistant` |
| [Crewai](daily-assistant/295-crewai_c93ce999/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/crewai.md) | ⭐ 3.3k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_165a91fb/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/create-release-checklist/SKILL.md) | ⭐ 117 | `daily assistant` |
| [Claude](daily-assistant/037-claude_60cb35f0/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/CLAUDE.md) | ⭐ 1.4k | `daily assistant` |
| [Codex Models](daily-assistant/270-codex-models_4cf972bb/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-orchestrator/references/codex-models.md) | ⭐ 10 | `daily assistant` |

### Data Analysis (54 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](data-analysis/226-name-skill_2414c4db/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-core/skills/drawio-logical-diagrams/SKILL.md) | ⭐ 102 | `data analysis` |
| [Raw Typescript Lambda](data-analysis/498-raw-typescript-lambda_c2fffb80/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/raw-typescript-lambda.md) | ⭐ 102 | `data analysis` |
| [Performance Architecture](data-analysis/499-performance-architecture_f977d60e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/performance-architecture.md) | ⭐ 3.3k | `data analysis` |
| [Dcr](data-analysis/500-dcr_14b0e9e4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/dcr.md) | ⭐ 3.3k | `data analysis` |
| [Postgres Upgrade Process](data-analysis/501-postgres-upgrade-process_60a48da6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/postgres-upgrade-process.md) | ⭐ 3.3k | `data analysis` |
| [007 Pluggable Cache Backend](data-analysis/502-007-pluggable-cache-backend_6ee2c5f9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/007-pluggable-cache-backend.md) | ⭐ 3.3k | `data analysis` |
| [Index](data-analysis/113-index_974fad2a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/index.md) | ⭐ 3.3k | `data analysis` |
| [Internal Observability](data-analysis/503-internal-observability_328ec0d8/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/internal-observability.md) | ⭐ 3.3k | `data analysis` |
| [Plugins](data-analysis/396-plugins_5f80238a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/plugins.md) | ⭐ 3.3k | `data analysis` |
| [Csv Pandas Chat Server](data-analysis/504-csv-pandas-chat-server_454d219f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/csv-pandas-chat-server.md) | ⭐ 3.3k | `data analysis` |
| [Data Analysis Server](data-analysis/505-data-analysis-server_cef0d00a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/data-analysis-server.md) | ⭐ 3.3k | `data analysis` |
| [Docx Server](data-analysis/506-docx-server_a79fb0e0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/docx-server.md) | ⭐ 3.3k | `data analysis` |
| [Graphviz Server](data-analysis/507-graphviz-server_b48cf167/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/graphviz-server.md) | ⭐ 3.3k | `data analysis` |
| [Plotly Server](data-analysis/508-plotly-server_3921f69a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/plotly-server.md) | ⭐ 3.3k | `data analysis` |
| [Pptx Server](data-analysis/509-pptx-server_9c6d336a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/pptx-server.md) | ⭐ 3.3k | `data analysis` |
| [Python Sandbox Server](data-analysis/510-python-sandbox-server_563cf18c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/python-sandbox-server.md) | ⭐ 3.3k | `data analysis` |
| [Brand Yml In R](data-analysis/479-brand-yml-in-r_780fb3bb/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/brand-yml-in-r.md) | ⭐ 117 | `data analysis` |
| [Quarto](data-analysis/480-quarto_81427cd4/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/quarto.md) | ⭐ 117 | `data analysis` |
| [Shiny R](data-analysis/481-shiny-r_0109581d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/shiny-r.md) | ⭐ 117 | `data analysis` |
| [Skill](data-analysis/226-name-skill_533d0b3b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/SKILL.md) | ⭐ 117 | `data analysis` |
| [Callouts](data-analysis/482-callouts_e4d4e5e1/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/callouts.md) | ⭐ 117 | `data analysis` |
| [Conversion Bookdown](data-analysis/483-conversion-bookdown_d994073b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-bookdown.md) | ⭐ 117 | `data analysis` |
| [Conversion Rmarkdown](data-analysis/484-conversion-rmarkdown_aee6105e/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-rmarkdown.md) | ⭐ 117 | `data analysis` |
| [Cross References](data-analysis/485-cross-references_202cb6fe/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/cross-references.md) | ⭐ 117 | `data analysis` |
| [Diagrams](data-analysis/442-diagrams_c3d9252d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/diagrams.md) | ⭐ 117 | `data analysis` |
| [Shortcodes](data-analysis/486-shortcodes_3dd7acd7/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/shortcodes.md) | ⭐ 117 | `data analysis` |
| [Tables](data-analysis/487-tables_15e26f8d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/tables.md) | ⭐ 117 | `data analysis` |
| [Accordions](data-analysis/488-accordions_c86bfde7/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/accordions.md) | ⭐ 117 | `data analysis` |
| [Navigation](data-analysis/489-navigation_d4984c12/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/navigation.md) | ⭐ 117 | `data analysis` |
| [Page Layouts](data-analysis/490-page-layouts_7aa85757/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/page-layouts.md) | ⭐ 117 | `data analysis` |
| [Theming](data-analysis/491-theming_77a63bbc/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/theming.md) | ⭐ 117 | `data analysis` |
| [Value Boxes](data-analysis/492-value-boxes_3b530b15/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/value-boxes.md) | ⭐ 117 | `data analysis` |
| [Claude](data-analysis/036-claude_8dccc45b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/CLAUDE.md) | ⭐ 18 | `data analysis` |
| [Orchestrator Discipline Patterns](data-analysis/486-orchestrator-discipline-patterns_95fa0f32/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/codebase/orchestrator-discipline-patterns.md) | ⭐ 18 | `data analysis` |
| [Research Scripts Refs](data-analysis/487-research-scripts-refs_9880003d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/plan/the-rewrite-room/research-scripts-refs.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_3854a340/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/research-curator/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_0637d133/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/swarm-primitives/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_a73903d5/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agentskill-kaizen/skills/kaizen-improvement/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_4667e05f/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/dasel/skills/dasel-reference/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_0f557dc7/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/audit-skill-lifecycle/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_5d22706c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/plugin-creator/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_7bfb55e1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/summarizer/skills/file-summarization/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_7326fbc9/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/summarizer/skills/image-summarization/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_58a16c4c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/development-harness/skills/workflows/execution/SKILL.md) | ⭐ 18 | `data analysis` |
| [Prompt Optimization](data-analysis/488-prompt-optimization_f32f3366/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/the-rewrite-room/skills/the-rewrite-room/workflows/prompt-optimization.md) | ⭐ 18 | `data analysis` |
| [Research Utilities](data-analysis/489-research-utilities_a32adef1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/the-rewrite-room/skills/the-rewrite-room/workflows/research-utilities.md) | ⭐ 18 | `data analysis` |
| [Contextual Ai Documentation Optimizer Adapter](data-analysis/490-contextual-ai-documentation-optimizer-adapter_394973bd/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/the-rewrite-room/skills/the-rewrite-room/workflows/adapters/contextual-ai-documentation-optimizer-adapter.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_bf9acbc8/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/mempool-analyzer/skills/analyzing-mempool/SKILL.md) | ⭐ 1.4k | `data analysis` |
| [Skill](data-analysis/226-name-skill_fecbac2d/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/nft-rarity-analyzer/skills/analyzing-nft-rarity/SKILL.md) | ⭐ 1.4k | `data analysis` |
| [Skill](data-analysis/226-name-skill_3ce7b1ff/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/core-development/claude-usage/SKILL.md) | ⭐ 10 | `data analysis` |
| [014 Security Headers Cors Middleware](data-analysis/481-014-security-headers-cors-middleware_f2c11a7a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/014-security-headers-cors-middleware.md) | ⭐ 3.3k | `data analysis` |
| [Skill](data-analysis/226-name-skill_95cf7107/) | [bowenliang123/md_exporter](https://raw.githubusercontent.com/bowenliang123/md_exporter/main/SKILL.md) | ⭐ 182 | `data analysis` |
| [Docker](data-analysis/docker_d9e87c50/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/docker.md) | ⭐ 36 | `data analysis` |
| [Docker Sandbox](data-analysis/docker-sandbox_0981cd42/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/examples/docker-sandbox.md) | ⭐ 36 | `data analysis` |

### Development (104 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Multitenancy](development/1203-multitenancy_d26020b1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/multitenancy.md) | ⭐ 3.3k | `development` |
| [Plugins](development/698-plugins_b81593c7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins.md) | ⭐ 3.3k | `openai` `moderation` `content-safety` |
| [Selecting An Mcp Gateway](development/955-selecting-an-mcp-gateway_cec67c55/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/selecting-an-mcp-gateway.md) | ⭐ 3.3k | `development` |
| [Postgresql Schema Configuration](development/2977-postgresql-schema-configuration_9f6cd11d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/postgresql-schema-configuration.md) | ⭐ 3.3k | `development` |
| [Mcp Developer Guide Json Rpc](development/2978-mcp-developer-guide-json-rpc_2bfea8c1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/mcp-developer-guide-json-rpc.md) | ⭐ 3.3k | `development` |
| [Profiling](development/2979-profiling_d8b52edd/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/profiling.md) | ⭐ 3.3k | `development` |
| [Sso Entra Role Mapping](development/2980-sso-entra-role-mapping_c6294cc1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-entra-role-mapping.md) | ⭐ 3.3k | `development` |
| [Sso Github Tutorial](development/2981-sso-github-tutorial_68c131fe/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-github-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Ibm Tutorial](development/2982-sso-ibm-tutorial_61b0c795/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-ibm-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Microsoft Entra Id Tutorial](development/2983-sso-microsoft-entra-id-tutorial_e579028c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-microsoft-entra-id-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Okta Tutorial](development/2984-sso-okta-tutorial_76fd8f99/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-okta-tutorial.md) | ⭐ 3.3k | `development` |
| [Teams](development/2985-teams_c57a4406/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/teams.md) | ⭐ 3.3k | `development` |
| [003 Expose Multi Transport Endpoints](development/2986-003-expose-multi-transport-endpoints_9637d88e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/003-expose-multi-transport-endpoints.md) | ⭐ 3.3k | `development` |
| [005 Vscode Devcontainer Support](development/2987-005-vscode-devcontainer-support_6e77ea8b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/005-vscode-devcontainer-support.md) | ⭐ 3.3k | `development` |
| [008 Federation Discovery](development/2988-008-federation-discovery_753ace5b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/008-federation-discovery.md) | ⭐ 3.3k | `development` |
| [010 Observability Prometheus](development/2989-010-observability-prometheus_fcec712d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/010-observability-prometheus.md) | ⭐ 3.3k | `development` |
| [011 Tool Federation](development/2990-011-tool-federation_24908a26/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/011-tool-federation.md) | ⭐ 3.3k | `development` |
| [016 Plugin Framework Ai Middleware](development/2991-016-plugin-framework-ai-middleware_f40a7f2c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/016-plugin-framework-ai-middleware.md) | ⭐ 3.3k | `development` |
| [017 Adopt Orjson Json Serialization](development/2992-017-adopt-orjson-json-serialization_7a66024e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/017-adopt-orjson-json-serialization.md) | ⭐ 3.3k | `development` |
| [022 Elicitation Passthrough Implementation](development/2993-022-elicitation-passthrough-implementation_288c5524/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/022-elicitation-passthrough-implementation.md) | ⭐ 3.3k | `development` |
| [024 Uvicorn Standard Extras](development/2994-024-uvicorn-standard-extras_6d24b0a8/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/024-uvicorn-standard-extras.md) | ⭐ 3.3k | `development` |
| [040 Flexible Admin Ui Sections](development/2995-040-flexible-admin-ui-sections_361375d2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/040-flexible-admin-ui-sections.md) | ⭐ 3.3k | `development` |
| [Gateway Hooks](development/2996-gateway-hooks_9aa359bb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins/gateway-hooks.md) | ⭐ 3.3k | `development` |
| [Security Hooks](development/2997-security-hooks_aa01ca2a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins/security-hooks.md) | ⭐ 3.3k | `development` |
| [A2A](development/2998-a2a_9a2ca7a0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/a2a.md) | ⭐ 3.3k | `development` |
| [Cline](development/2999-cline_8713195d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/cline.md) | ⭐ 3.3k | `development` |
| [Llm Chat](development/3000-llm-chat_b4fec4c4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/llm-chat.md) | ⭐ 3.3k | `development` |
| [Index](development/468-index_79b6cf6d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/index.md) | ⭐ 3.3k | `custom` `filter` |
| [Lifecycle](development/2173-lifecycle_17a38751/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/lifecycle.md) | ⭐ 3.3k | `plugin` |
| [Rust Plugins](development/3001-rust-plugins_b6ecb04f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/rust-plugins.md) | ⭐ 3.3k | `development` |
| [Index](development/468-index_817e60d4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/index.md) | ⭐ 3.3k | `development` |
| [Calculator Server](development/3002-calculator-server_84a343ea/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/calculator-server.md) | ⭐ 3.3k | `development` |
| [Box](development/3003-box_1e5cdfee/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/box/box.md) | ⭐ 3.3k | `development` |
| [Monday Mcp](development/3004-monday-mcp_91796fb7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/monday/monday-mcp.md) | ⭐ 3.3k | `development` |
| [Skill](development/1178-name-skill_17e6beab/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/SKILL.md) | ⭐ 117 | `development` |
| [Brand Yml Spec](development/2903-brand-yml-spec_4c1c8800/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/brand-yml-spec.md) | ⭐ 117 | `development` |
| [Shiny Python](development/2904-shiny-python_6d034f3b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/shiny-python.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_d8e79aed/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/posit-dev/critical-code-reviewer/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_e5da82b8/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/posit-dev/describe-design/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_06ff523f/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/posit-dev/pr-create/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_b817f91f/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_faade469/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cran-extrachecks/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_1a9af13a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/lifecycle/SKILL.md) | ⭐ 117 | `development` |
| [Tidyverse Formatting](development/2905-tidyverse-formatting_6f9d8aae/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/references/tidyverse-formatting.md) | ⭐ 117 | `package-name` `category` |
| [Code Cells](development/2906-code-cells_4883683d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/code-cells.md) | ⭐ 117 | `development` |
| [Conditional Content](development/2907-conditional-content_3c229e8a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conditional-content.md) | ⭐ 117 | `development` |
| [Conversion Xaringan](development/2908-conversion-xaringan_39d57dc1/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-xaringan.md) | ⭐ 117 | `development` |
| [Extensions](development/2909-extensions_7560b3bd/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/extensions.md) | ⭐ 117 | `development` |
| [Figures](development/2910-figures_c0e7f57e/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/figures.md) | ⭐ 117 | `development` |
| [Markdown Linting](development/2911-markdown-linting_0dff392d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/markdown-linting.md) | ⭐ 117 | `development` |
| [Ansi Operations](development/2912-ansi-operations_3b86a218/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/ansi-operations.md) | ⭐ 117 | `development` |
| [Inline Markup](development/2913-inline-markup_1ed5cb03/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/inline-markup.md) | ⭐ 117 | `development` |
| [Progress](development/714-progress_d7998082/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/progress.md) | ⭐ 117 | `development` |
| [Lifecycle Stages](development/2914-lifecycle-stages_4f024f57/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/lifecycle/references/lifecycle-stages.md) | ⭐ 117 | `development` |
| [Migration](development/2204-migration_55087718/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/migration.md) | ⭐ 117 | `development` |
| [Feature Context Validate Orchestrator Discipline](development/2977-feature-context-validate-orchestrator-discipline_1964559c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/feature-context-validate-orchestrator-discipline.md) | ⭐ 18 | `development` |
| [Task File Format](development/2797-task_file_format_c1ef2af8/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/docs/TASK_FILE_FORMAT.md) | ⭐ 18 | `development` |
| [Claude](development/140-claude_76000a30/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/CLAUDE.md) | ⭐ 18 | `development` |
| [Research Agents](development/2978-research-agents_1d3d3c9d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/plan/the-rewrite-room/research-agents.md) | ⭐ 18 | `development` |
| [Research Architecture Pattern](development/2979-research-architecture-pattern_edc4bb85/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/plan/the-rewrite-room/research-architecture-pattern.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_3bc11637/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/agent-creator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_771f59ba/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/create-merge-request-changelog/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_4a8cf475/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/external-pattern-integrator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_45901b1b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/orchestrating-swarms/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_013658d1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/scientific-thinking/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_d1022c7b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/swarm-operations/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_d132f8a1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/swarm-patterns/SKILL.md) | ⭐ 18 | `development` |
| [Claude](development/140-claude_f69e49c5/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/scripts/CLAUDE.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_798aebc4/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agent-orchestration/skills/how-to-delegate/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_fde04101/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/bash-development/skills/bash-53-features/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_bbf386d9/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/bash-development/skills/bash-lint/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_8e1c6987/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/bash-development/skills/bash-logging/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_e6f15e2d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/commitlint/skills/commitlint/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_67bd8e87/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/holistic-linting/skills/holistic-linting-orchestrator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_fe04c1e3/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/holistic-linting/skills/holistic-linting/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_6e18624e/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/litellm/skills/litellm/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_4fd6044a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/llamafile/skills/llamafile/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_ff9fe8fb/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/agent-creator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_c856efe0/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/audit-skill-completeness/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_1640461c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/lint/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_5f0e587a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/memory-and-rules/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_becad5b3/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/token-launch-tracker/skills/tracking-token-launches/SKILL.md) | ⭐ 1.4k | `development` |
| [Response Api](development/1220-response_api_7468ad4b/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/response_api.md) | 🔥 36.4k | `development` |
| [Auto Routing](development/2892-auto_routing_e7ca64fd/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/auto_routing.md) | 🔥 36.4k | `development` |
| [Grooming 2026 02 21](development/2869-grooming-2026-02-21_f58716b3/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/grooming-reports/grooming-2026-02-21.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_7063d1fd/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/gh/SKILL.md) | ⭐ 18 | `development` |
| [Cto](development/2917-cto_087484d0/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-cto/agents/cto.md) | ⭐ 10 | `development` |
| [Subagent Patterns](development/2918-subagent-patterns_4b0a81e2/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-orchestrator/references/subagent-patterns.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_4a33313e/) | [artwist-polyakov/polyakov-claude-skills](https://raw.githubusercontent.com/artwist-polyakov/polyakov-claude-skills/main/plugins/codex-review/skills/codex-review/SKILL.md) | ⭐ 38 | `development` |
| [Documentation Updates Summary](development/2874-documentation_updates_summary_665ee62f/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/DOCUMENTATION_UPDATES_SUMMARY.md) | 🔥 9.7k | `development` |
| [Multi Source](development/2875-multi-source_9f256d80/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/advanced/multi-source.md) | 🔥 9.7k | `development` |
| [Unified Scraping](development/1108-unified_scraping_1d019192/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/features/UNIFIED_SCRAPING.md) | 🔥 9.7k | `development` |
| [01 Installation](development/2876-01-installation_39803a31/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/01-installation.md) | 🔥 9.7k | `development` |
| [03 Your First Skill](development/2877-03-your-first-skill_c9fe1c7e/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/03-your-first-skill.md) | 🔥 9.7k | `development` |
| [Cli Reference](development/1593-cli_reference_ba8e75e8/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/CLI_REFERENCE.md) | 🔥 9.7k | `development` |
| [Config Format](development/2878-config_format_e39c2327/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/CONFIG_FORMAT.md) | 🔥 9.7k | `development` |
| [06 Troubleshooting](development/2879-06-troubleshooting_78ddfebd/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/06-troubleshooting.md) | 🔥 9.7k | `development` |
| [Cli Reference](development/1593-cli_reference_2e61f4f8/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/zh-CN/reference/CLI_REFERENCE.md) | 🔥 9.7k | `development` |
| [Config Format](development/2878-config_format_af066045/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/zh-CN/reference/CONFIG_FORMAT.md) | 🔥 9.7k | `development` |
| [Claude](development/claude_8e8685ee/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/CLAUDE.md) | ⭐ 36 | `development` |
| [Abbreviations](development/abbreviations_cdb8d3d2/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/includes/abbreviations.md) | ⭐ 36 | `development` |
| [Installation](development/installation_b151f30f/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/installation.md) | ⭐ 36 | `development` |
| [Cli Agent](development/cli-agent_8017030b/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/examples/cli-agent.md) | ⭐ 36 | `development` |
| [Local Backend](development/local-backend_752c9262/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/examples/local-backend.md) | ⭐ 36 | `development` |

### Development/Devops (102 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agent Deployment Guide](development/devops/229-agent-deployment-guide_fc15347e/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-deploy/references/agent-deployment-guide.md) | ⭐ 80 | `development` |
| [Guide Skills Frontend](development/devops/376-guide-skills-frontend_2f63a088/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/docs/guide-skills-frontend.md) | ⭐ 102 | `development` |
| [Guide Skills Monorepo](development/devops/377-guide-skills-monorepo_83c4b53f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/docs/guide-skills-monorepo.md) | ⭐ 102 | `development` |
| [Skill](development/devops/014-name-skill_1ec7fc18/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `java` |
| [Skill](development/devops/014-name-skill_9bc6e00e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `php` |
| [Skill](development/devops/014-name-skill_3c0620fc/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `python` |
| [Skill](development/devops/014-name-skill_2658dda5/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `typescript` |
| [Skill](development/devops/014-name-skill_62d193fa/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/SKILL.md) | ⭐ 102 | `nextjs` `next.js` `deployment` |
| [Skill](development/devops/014-name-skill_d0e8129e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/SKILL.md) | ⭐ 102 | `nx` `monorepo` `typescript` |
| [Skill](development/devops/014-name-skill_4133613a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/SKILL.md) | ⭐ 102 | `development` |
| [Raw Java Lambda](development/devops/378-raw-java-lambda_bc77ac25/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/raw-java-lambda.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_7a94e770/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_0963c52b/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Chalice Lambda](development/devops/380-chalice-lambda_73723400/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/chalice-lambda.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_64158df7/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_5ae8bdfe/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Docker Patterns](development/devops/381-docker-patterns_6cc95fae/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/references/docker-patterns.md) | ⭐ 102 | `development` |
| [Basics](development/devops/382-basics_94f6201e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/basics.md) | ⭐ 102 | `development` |
| [Mcp Server Python](development/devops/383-mcp-server-python_5f654215/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/mcp-server-python.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_ccaab12c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/index.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_3c572c12/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/index.md) | ⭐ 3.3k | `development` |
| [Roadmap](development/devops/097-roadmap_c9aa8183/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/roadmap.md) | ⭐ 3.3k | `development` |
| [Security Features](development/devops/384-security-features_012c8460/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/security-features.md) | ⭐ 3.3k | `development` |
| [Input Validation](development/devops/228-input-validation_01d87216/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/input-validation.md) | ⭐ 3.3k | `development` |
| [Mcp Architecture Patterns](development/devops/385-mcp-architecture-patterns_230d8859/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/mcp-architecture-patterns.md) | ⭐ 3.3k | `development` |
| [Argocd](development/devops/386-argocd_31a49405/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/argocd.md) | ⭐ 3.3k | `development` |
| [Cforge Gateway](development/devops/387-cforge-gateway_ac5a8cdb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/cforge-gateway.md) | ⭐ 3.3k | `security` `filter` |
| [Compose](development/devops/388-compose_366cfc28/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/compose.md) | ⭐ 3.3k | `development` |
| [Fly Io](development/devops/389-fly-io_e5ccc7e2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/fly-io.md) | ⭐ 3.3k | `development` |
| [Google Cloud Run](development/devops/390-google-cloud-run_8be8613f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/google-cloud-run.md) | ⭐ 3.3k | `development` |
| [Helm](development/devops/391-helm_2945bf9d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/helm.md) | ⭐ 3.3k | `development` |
| [Ibm Code Engine](development/devops/392-ibm-code-engine_fe3cb3dc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/ibm-code-engine.md) | ⭐ 3.3k | `development` |
| [Local](development/devops/393-local_93c0d987/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/local.md) | ⭐ 3.3k | `development` |
| [Minikube](development/devops/394-minikube_0c8a4b03/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/minikube.md) | ⭐ 3.3k | `development` |
| [Proxy Auth](development/devops/395-proxy-auth_7d43d3de/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/proxy-auth.md) | ⭐ 3.3k | `development` |
| [Tls Configuration](development/devops/201-tls-configuration_54ff4560/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/tls-configuration.md) | ⭐ 3.3k | `development` |
| [Developer Onboarding](development/devops/396-developer-onboarding_61358b2f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/developer-onboarding.md) | ⭐ 3.3k | `development` |
| [Github](development/devops/397-github_487e80b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/github.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_a20a1f7e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/index.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_0010150a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/faq/index.md) | ⭐ 3.3k | `development` |
| [Admin Ui Customization](development/devops/398-admin-ui-customization_b4697c47/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/admin-ui-customization.md) | ⭐ 3.3k | `development` |
| [Api Usage](development/devops/399-api-usage_930cc30b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/api-usage.md) | ⭐ 3.3k | `development` |
| [Configuration](development/devops/009-configuration_ed70f974/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/configuration.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_7c9efd7f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/index.md) | ⭐ 3.3k | `development` |
| [Logging](development/devops/400-logging_7601ca0a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/logging.md) | ⭐ 3.3k | `development` |
| [Proxy](development/devops/401-proxy_5a83e05a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/proxy.md) | ⭐ 3.3k | `development` |
| [Rbac](development/devops/010-rbac_46494f27/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/rbac.md) | ⭐ 3.3k | `development` |
| [Scale](development/devops/402-scale_5375b413/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/scale.md) | ⭐ 3.3k | `development` |
| [Securing](development/devops/403-securing_264f96d4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/securing.md) | ⭐ 3.3k | `development` |
| [Self Signed Certificates](development/devops/404-self-signed-certificates_c29006e4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/self-signed-certificates.md) | ⭐ 3.3k | `development` |
| [Sso Generic Oidc Tutorial](development/devops/405-sso-generic-oidc-tutorial_1cf7119b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-generic-oidc-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Google Tutorial](development/devops/406-sso-google-tutorial_ec0b87b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-google-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Keycloak Tutorial](development/devops/407-sso-keycloak-tutorial_bc429daa/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-keycloak-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso](development/devops/408-sso_3ff6196e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso.md) | ⭐ 3.3k | `development` |
| [Tuning](development/devops/409-tuning_b13077bc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/tuning.md) | ⭐ 3.3k | `development` |
| [Upgrade](development/devops/410-upgrade_b3150af6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/upgrade.md) | ⭐ 3.3k | `development` |
| [Well Known Uris](development/devops/411-well-known-uris_b7a6cdf4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/well-known-uris.md) | ⭐ 3.3k | `development` |
| [Config Validation](development/devops/412-config-validation_dc2663b4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/operations/config-validation.md) | ⭐ 3.3k | `development` |
| [Cpu Spin Loop Mitigation](development/devops/413-cpu-spin-loop-mitigation_fdb3b575/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/operations/cpu-spin-loop-mitigation.md) | ⭐ 3.3k | `development` |
| [Features](development/devops/361-features_4c3fbb11/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/features.md) | ⭐ 3.3k | `development` |
| [Admin UI](development/devops/414-ui_223c1cd2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/ui.md) | ⭐ 3.3k | `development` |
| [Argocd Helm Deployment Ibm Cloud Iks](development/devops/415-argocd-helm-deployment-ibm-cloud-iks_8a45b06f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/argocd-helm-deployment-ibm-cloud-iks.md) | ⭐ 3.3k | `development` |
| [Openwebui Tutorial](development/devops/416-openwebui-tutorial_ace5026d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/openwebui-tutorial.md) | ⭐ 3.3k | `development` |
| [Mcpgateway Translate](development/devops/417-mcpgateway-translate_199775f5/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/mcpgateway-translate.md) | ⭐ 3.3k | `development` |
| [Multi Auth Headers](development/devops/418-multi-auth-headers_d4462e80/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/multi-auth-headers.md) | ⭐ 3.3k | `development` |
| [Reverse Proxy](development/devops/419-reverse-proxy_a49bea67/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/reverse-proxy.md) | ⭐ 3.3k | `development` |
| [014 Security Headers Cors Middleware](development/devops/420-014-security-headers-cors-middleware_009499c3/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/014-security-headers-cors-middleware.md) | ⭐ 3.3k | `development` |
| [015 Well Known Uri Handler](development/devops/421-015-well-known-uri-handler_198ce9fb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/015-well-known-uri-handler.md) | ⭐ 3.3k | `development` |
| [018 Built In Response Compression](development/devops/422-018-built-in-response-compression_161e1ccb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/018-built-in-response-compression.md) | ⭐ 3.3k | `development` |
| [019 Modular Architecture Split](development/devops/423-019-modular-architecture-split_606af208/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/019-modular-architecture-split.md) | ⭐ 3.3k | `development` |
| [025 Granian Http Server](development/devops/424-025-granian-http-server_7849918e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/025-granian-http-server.md) | ⭐ 3.3k | `development` |
| [032 Mcp Session Pool](development/devops/425-032-mcp-session-pool_fd39527c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/032-mcp-session-pool.md) | ⭐ 3.3k | `development` |
| [036 Bootstrap Custom Roles](development/devops/426-036-bootstrap-custom-roles_98fdf334/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/036-bootstrap-custom-roles.md) | ⭐ 3.3k | `development` |
| [038 Multi Worker Session Affinity](development/devops/193-038-multi-worker-session-affinity_de1547c3/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/038-multi-worker-session-affinity.md) | ⭐ 3.3k | `development` |
| [Observability](development/devops/427-observability_358deba7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/observability.md) | ⭐ 3.3k | `development` |
| [Phoenix](development/devops/428-phoenix_a6b52052/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/phoenix.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_49786343/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/social/index.md) | ⭐ 3.3k | `development` |
| [Bee](development/devops/429-bee_bd521c8a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/bee.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_6da5357e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/index.md) | ⭐ 3.3k | `development` |
| [Claude Desktop](development/devops/430-claude-desktop_33f0496c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/claude-desktop.md) | ⭐ 3.3k | `development` |
| [Continue](development/devops/431-continue_4f371801/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/continue.md) | ⭐ 3.3k | `development` |
| [Grpc Transport](development/devops/226-grpc-transport_c50d7487/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/grpc-transport.md) | ⭐ 3.3k | `development` |
| [Mtls](development/devops/432-mtls_54eac01a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/mtls.md) | ⭐ 3.3k | `development` |
| [Terraform](development/devops/433-terraform_86d5b2a9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/hashicorp/terraform.md) | ⭐ 3.3k | `development` |
| [Langflow Server](development/devops/434-langflow-server_4739c90c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/third-party/langflow-server.md) | ⭐ 3.3k | `development` |
| [Github](development/devops/397-github_cda2891c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/microsoft/github.md) | ⭐ 3.3k | `development` |
| [Rollback](development/devops/364-rollback_4f482aad/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/troubleshoot/rollback.md) | 🔥 36.4k | `development` |
| [Cli Commands](development/devops/319-cli-commands_4a42a248/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/cli-commands.md) | ⭐ 61 | `development` |
| [Multi Tenancy](development/devops/361-multi-tenancy_677a145e/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/multi-tenancy.md) | ⭐ 61 | ``tenant:${tenantId}`` `role:support` |
| [Skill](development/devops/014-name-skill_7c48cfe9/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-orchestrator/SKILL.md) | ⭐ 10 | `development` |
| [Readme Es](development/devops/435-readme_es_1474cf17/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_es.md) | ⭐ 820 | `development` |
| [Readme Ja](development/devops/436-readme_ja_a026d9c8/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_ja.md) | ⭐ 820 | `development` |
| [Readme Pt Br](development/devops/437-readme_pt-br_698660d2/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_pt-BR.md) | ⭐ 820 | `development` |
| [Readme Zh Cn](development/devops/438-readme_zh-cn_b5fce0bf/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_zh-CN.md) | ⭐ 820 | `development` |
| [Mcp](development/devops/025-mcp_b3a0b417/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/docs/docs/concepts/mcp.md) | ⭐ 820 | `development` |
| [Securing](development/devops/361-securing_0ab63085/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/securing.md) | ⭐ 3.3k | `development` |
| [Agents](development/devops/053-agents_05228376/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/AGENTS.md) | 🔥 9.7k | `development` |
| [Environment Variables](development/devops/362-environment_variables_588fc08a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/ENVIRONMENT_VARIABLES.md) | 🔥 9.7k | `development` |
| [05 Workflows](development/devops/363-05-workflows_75f7af4a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/05-workflows.md) | 🔥 9.7k | `development` |
| [Quick Reference](development/devops/364-quick_reference_07ac890f/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/archive/legacy/QUICK_REFERENCE.md) | 🔥 9.7k | `development` |
| [05 Workflows](development/devops/363-05-workflows_f0c92b2a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/zh-CN/user-guide/05-workflows.md) | 🔥 9.7k | `development` |
| [Index](development/devops/index_e539cf55/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/examples/index.md) | ⭐ 36 | `development` |

### Development/Testing (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/testing/002-name-skill_187b30c5/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nestjs-drizzle-crud-generator/SKILL.md) | ⭐ 102 | `development` |
| [Nestjs Config](development/testing/089-nestjs-config_ab544f6a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/nestjs-config.md) | ⭐ 102 | `development` |
| [Db Performance](development/testing/090-db-performance_9aa0f078/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/db-performance.md) | ⭐ 3.3k | `development` |
| [Backends](development/testing/backends_3471597a/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/backends.md) | ⭐ 36 | `development` |
| [Index](development/testing/index_9516d09b/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/index.md) | ⭐ 36 | `development` |

### Development/Tools (55 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Typescript](development/tools/349-typescript_d85d7b93/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/typescript.md) | ⭐ 102 | `development` |
| [Skill](development/tools/002-name-skill_034a479c/) | [netresearch/jira-skill](https://raw.githubusercontent.com/netresearch/jira-skill/main/skills/jira-communication/SKILL.md) | ⭐ 24 | `development` |
| [Developing](development/tools/350-developing_4eef14b0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/DEVELOPING.md) | ⭐ 3.3k | `development` |
| [Enable Payload Logging](development/tools/351-enable_payload_logging_c75da92d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/enable_payload_logging.md) | ⭐ 3.3k | `development` |
| [Mcpgateway](development/tools/010-mcpgateway_f02247c0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/mcpgateway.md) | ⭐ 3.3k | `development` |
| [Plugins Llms](development/tools/352-plugins-llms_8fca01fe/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/plugins-llms.md) | ⭐ 3.3k | `development` |
| [Oauth Authorization Code Ui Design](development/tools/011-oauth-authorization-code-ui-design_ec055b8d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/oauth-authorization-code-ui-design.md) | ⭐ 3.3k | `development` |
| [Oauth Design](development/tools/012-oauth-design_dd0d2a41/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/oauth-design.md) | ⭐ 3.3k | `development` |
| [Bulk Import](development/tools/353-bulk-import_8f7144ba/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/bulk-import.md) | ⭐ 3.3k | `development` |
| [Export Import Reference](development/tools/354-export-import-reference_cd3b304c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import-reference.md) | ⭐ 3.3k | `development` |
| [Metadata Tracking](development/tools/355-metadata-tracking_93c516b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/metadata-tracking.md) | ⭐ 3.3k | `development` |
| [Oauth](development/tools/356-oauth_2ff3a5f9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/oauth.md) | ⭐ 3.3k | `development` |
| [Tags](development/tools/357-tags_0898f285/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/tags.md) | ⭐ 3.3k | `development` |
| [006 Gateway Tool Rate Limiting](development/tools/358-006-gateway-tool-rate-limiting_bb786b5a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/006-gateway-tool-rate-limiting.md) | ⭐ 3.3k | `development` |
| [021 Built In Proxy Vs Service Mesh](development/tools/359-021-built-in-proxy-vs-service-mesh_87238005/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/021-built-in-proxy-vs-service-mesh.md) | ⭐ 3.3k | `development` |
| [023 One Time Authentication Servers](development/tools/360-023-one-time-authentication-servers_6f4caffb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/023-one-time-authentication-servers.md) | ⭐ 3.3k | `development` |
| [027 Migrate Psycopg3](development/tools/361-027-migrate-psycopg3_ae99362f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/027-migrate-psycopg3.md) | ⭐ 3.3k | `development` |
| [035 Query Parameter Authentication](development/tools/362-035-query-parameter-authentication_7da5eb7c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/035-query-parameter-authentication.md) | ⭐ 3.3k | `development` |
| [Llamaindex](development/tools/363-llamaindex_332803db/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/llamaindex.md) | ⭐ 3.3k | `development` |
| [Http Auth Hooks](development/tools/364-http-auth-hooks_18dede4f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/http-auth-hooks.md) | ⭐ 3.3k | `development` |
| [Fast Time Server](development/tools/365-fast-time-server_c35f9972/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/fast-time-server.md) | ⭐ 3.3k | `development` |
| [Instana](development/tools/366-instana_c81a89af/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/ibm/instana.md) | ⭐ 3.3k | `development` |
| [Mermaid Server](development/tools/367-mermaid-server_f0be2d55/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/mermaid-server.md) | ⭐ 3.3k | `development` |
| [Index](development/tools/062-index_ec82d631/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/open/index.md) | ⭐ 3.3k | `development` |
| [Skill](development/tools/002-name-skill_fe5d0c93/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib-theming/SKILL.md) | ⭐ 117 | `development` |
| [Conditions](development/tools/329-conditions_fbf1013a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/conditions.md) | ⭐ 117 | `development` |
| [Sass And Css Variables](development/tools/330-sass-and-css-variables_db885c80/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib-theming/references/sass-and-css-variables.md) | ⭐ 117 | `development` |
| [Backlog](development/tools/336-backlog_5f579c4c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/BACKLOG.md) | ⭐ 18 | `development` |
| [Copilot Instructions](development/tools/337-copilot-instructions_727ce275/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.github/copilot-instructions.md) | ⭐ 18 | `development` |
| [Architect Validate Orchestrator Discipline](development/tools/338-architect-validate-orchestrator-discipline_2adb21c5/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/architect-validate-orchestrator-discipline.md) | ⭐ 18 | `development` |
| [Tasks 4 Validate Orchestrator Discipline](development/tools/339-tasks-4-validate-orchestrator-discipline_55f8f536/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/tasks-4-validate-orchestrator-discipline.md) | ⭐ 18 | `development` |
| [Orchestrator Discipline Grooming 2026 02 20](development/tools/340-orchestrator-discipline-grooming-2026-02-20_115fd850/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/grooming-reports/orchestrator-discipline-grooming-2026-02-20.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_f8e7b13a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/commit-staged/SKILL.md) | ⭐ 18 | `development` |
| [2026 02 20 Duckdb Lock Scope Flag](development/tools/341-2026-02-20-duckdb-lock-scope-flag_18db1838/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agentskill-kaizen/docs/plans/2026-02-20-duckdb-lock-scope-flag.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_0fffb33b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/clang-format/skills/clang-format/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_98d88f7e/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/conventional-commits/skills/conventional-commits/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_97365666/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/development-harness/skills/implementation-manager/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_b65d9a7d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/fastmcp-creator/skills/fastmcp-creator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_900d15ae/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/agentskills/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_da240963/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_ab8de90a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/hatchling/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_7c7e3416/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/implementation-manager/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_2be1a33c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/python3-development/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_d5053e84/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/toml-python/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_6f3596ac/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/workshops/.claude/skills/embedded-debug-tools/SKILL.md) | ⭐ 18 | `development` |
| [Fact Checker](development/tools/329-fact-checker_7699dca6/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/agents/fact-checker.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_f380706f/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/fact-check/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_07dab008/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-cto/SKILL.md) | ⭐ 10 | `development` |
| [Agents](development/tools/015-agents_d9ec1552/) | [bowenliang123/md_exporter](https://raw.githubusercontent.com/bowenliang123/md_exporter/main/AGENTS.md) | ⭐ 182 | `development` |
| [Claude](development/tools/017-claude_69bd09db/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/CLAUDE.md) | 🔥 9.7k | `development` |
| [Architecture](development/tools/051-architecture_97bc174a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/ARCHITECTURE.md) | 🔥 9.7k | `development` |
| [Mcp Server](development/tools/330-mcp-server_11035d67/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/advanced/mcp-server.md) | 🔥 9.7k | `development` |
| [Mcp Reference](development/tools/331-mcp_reference_9b2d06aa/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/MCP_REFERENCE.md) | 🔥 9.7k | `development` |
| [Getting Help](development/tools/getting-help_75a6ba4c/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/getting-help.md) | ⭐ 36 | `development` |
| [Index](development/tools/index_8f2ce9d0/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/index.md) | ⭐ 36 | `development` |

### Investment (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](investment/021-name-skill_65b23c64/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/whale-alert-monitor/skills/monitoring-whale-activity/SKILL.md) | ⭐ 1.4k | `investment` |

### Productivity (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Parallel Session Cleanup](productivity/179-parallel-session-cleanup_d83f616a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/parallel-session-cleanup.md) | ⭐ 3.3k | `productivity` |
| [Dark Mode](productivity/176-dark-mode_77760605/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib-theming/references/dark-mode.md) | ⭐ 117 | `productivity` |
| [Tooltips Popovers](productivity/177-tooltips-popovers_f123d4bd/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/tooltips-popovers.md) | ⭐ 117 | `productivity` |

### Research (23 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Pandoc Server](research/265-pandoc-server_a65d4d68/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/pandoc-server.md) | ⭐ 3.3k | `research` |
| [Latex Server](research/266-latex-server_e6f7093b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/latex-server.md) | ⭐ 3.3k | `research` |
| [Citations](research/261-citations_00e0d6bf/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/citations.md) | ⭐ 117 | `research` |
| [Divs And Spans](research/262-divs-and-spans_a11251da/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/divs-and-spans.md) | ⭐ 117 | `research` |
| [Layout](research/263-layout_f4e230b1/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/layout.md) | ⭐ 117 | `research` |
| [Yaml Front Matter](research/264-yaml-front-matter_8c124133/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/yaml-front-matter.md) | ⭐ 117 | `research` |
| [Claude](research/015-claude_b0ed45bf/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/CLAUDE.md) | ⭐ 829 | `research` |
| [Writer](research/144-writer_1664b873/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/.claude/WRITER.md) | ⭐ 829 | `research` |
| [Api](research/146-api_a79b5e1a/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/docs/API.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_83163c0c/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/citation-management/SKILL.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_995bd408/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/SKILL.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_ed983e74/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/research-lookup/SKILL.md) | ⭐ 829 | `research` |
| [Citation Validation](research/014-citation_validation_5901b5ff/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/citation-management/references/citation_validation.md) | ⭐ 829 | `research` |
| [Api Reference](research/007-api_reference_bd3cf235/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/api_reference.md) | ⭐ 829 | `research` |
| [Deep Research Guide](research/258-deep_research_guide_c3d92b15/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/deep_research_guide.md) | ⭐ 829 | `research` |
| [Extraction Patterns](research/259-extraction_patterns_49ce48cf/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/extraction_patterns.md) | ⭐ 829 | `research` |
| [Search Best Practices](research/260-search_best_practices_8ba3e7fe/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/search_best_practices.md) | ⭐ 829 | `research` |
| [Workflow Recipes](research/261-workflow_recipes_0676ea23/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/workflow_recipes.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_41583fad/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/cross-chain-bridge-monitor/skills/monitoring-cross-chain-bridges/SKILL.md) | ⭐ 1.4k | `research` |
| [Fleet Config](research/257-fleet-config_ee28fcec/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/fleet-config.md) | ⭐ 61 | `]                            # Optional: key:value strings for filtering
    memory_blocks: [` |
| [Introduction](research/267-introduction_9d11f1af/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/docs/docs/introduction.md) | ⭐ 820 | `research` |
| [02 Quick Start](research/258-02-quick-start_aa356bca/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/02-quick-start.md) | 🔥 9.7k | `research` |
| [01 Core Concepts](research/259-01-core-concepts_8eaa4c6e/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/01-core-concepts.md) | 🔥 9.7k | `research` |

## How Skills Are Organized

Skills are automatically categorized based on their purpose:

- **Development**: Coding, debugging, testing, and developer tools
- **Daily Assistant**: Task management, scheduling, and reminders
- **Content Creation**: Writing, editing, and content generation
- **Data Analysis**: Visualization, statistics, and data processing
- **Automation**: Workflows, scripts, and task automation
- **Research**: Academic tools, citations, and literature
- **Communication**: Email, messaging, and collaboration
- **Productivity**: Efficiency tools and optimization
- **Commercial**: E-commerce and business tools
- **Investment**: Trading, stocks, and financial analysis

## Manual Usage

These skills can also be used directly without installing patches:

1. Browse the category folders to find relevant skills
2. Navigate to a skill's subdirectory
3. Read the skill's README.md for metadata and description
4. Use the skill's .md file content with Claude Code or similar AI assistants

## File Naming Convention

Each skill is stored in a subdirectory named: `source_name_hashprefix/`

- `source_name`: The original filename (sanitized)
- `hashprefix`: First 8 characters of the content hash (ensures uniqueness)

The hash-based naming ensures that:
- The same skill content always maps to the same directory
- Updated skills automatically replace old versions
- No duplicate directories for the same content

## Skill Index

This repository includes a `.index.json` file that tracks all skills and their locations.
This index enables:
- Incremental updates (only writing changed skills)
- Efficient change detection
- Proper handling of skill updates from source repositories

## Contributing

This repository is automatically maintained by [SkillFlow](https://github.com/tools-only/SkillFlow). Skills are aggregated from open-source repositories.

---

*Last updated: 2026-02-22 09:17:28 UTC*
*Automatically maintained by SkillFlow*
