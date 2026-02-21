# X-Skills

A curated collection of **301 AI-powered skills** organized into 12 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Workflow** (10 skills)
- **Commercial** (23 skills)
- **Communication** (7 skills)
- **Content Creation** (21 skills)
- **Daily Assistant** (2 skills)
- **Data Analysis** (20 skills)
- **Development** (65 skills)
- **Development/Devops** (109 skills)
- **Development/Testing** (7 skills)
- **Development/Tools** (34 skills)
- **Productivity** (1 skill)
- **Research** (2 skills)

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


### Automation/Workflow (10 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Vfunc Sig](automation/workflow/136-vfunc_sig_9dfecfa2/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.serena/memories/vfunc_sig.md) | ⭐ 19 | `automation` |
| [Skill](automation/workflow/002-name-skill_2854fe0b/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-deploy/SKILL.md) | ⭐ 80 | `automation` |
| [Cli Guide](automation/workflow/137-cli-guide_69513c0f/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-ai-agentscript/references/cli-guide.md) | ⭐ 80 | `automation` |
| [Export Import Architecture](automation/workflow/export-import-architecture_1fc68b51/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/export-import-architecture.md) | ⭐ 3.3k | `automation` |
| [Export Import Tutorial](automation/workflow/export-import-tutorial_372778e0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import-tutorial.md) | ⭐ 3.3k | `automation` |
| [Export Import](automation/workflow/export-import_af809418/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import.md) | ⭐ 3.3k | `automation` |
| [Mcp Cli](automation/workflow/mcp-cli_dbea27bc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/mcp-cli.md) | ⭐ 3.3k | `automation` |
| [Libreoffice Server](automation/workflow/libreoffice-server_2c0ac43b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/libreoffice-server.md) | ⭐ 3.3k | `automation` |
| [Index](automation/workflow/index_f2308b5e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/testimonials/index.md) | ⭐ 3.3k | `automation` |
| [Openai Sdk](automation/workflow/openai-sdk_4fb90c4e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/openai-sdk.md) | ⭐ 3.3k | `automation` |

### Commercial (23 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Express Adapter](commercial/express-adapter_6db5888c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/express-adapter.md) | ⭐ 102 | `commercial` |
| [Passkey](commercial/passkey_b16def64/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/passkey.md) | ⭐ 102 | `commercial` |
| [Routing Patterns](commercial/routing-patterns_959a3406/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/routing-patterns.md) | ⭐ 102 | `commercial` |
| [Caching Strategies](commercial/caching-strategies_3decca32/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/caching-strategies.md) | ⭐ 102 | `products` |
| [Server Components](commercial/server-components_249b33ee/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/server-components.md) | ⭐ 102 | `commercial` |
| [Advanced](commercial/advanced_87d29bd8/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/advanced.md) | ⭐ 102 | `commercial` |
| [Migration 0.7.0](commercial/migration-070_0d56b710/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/MIGRATION-0.7.0.md) | ⭐ 3.3k | `commercial` |
| [Kubernetes](commercial/kubernetes_0493cab9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/kubernetes.md) | ⭐ 3.3k | `commercial` |
| [Developer Workstation](commercial/developer-workstation_2f62c709/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/developer-workstation.md) | ⭐ 3.3k | `commercial` |
| [Oauth Troubleshooting](commercial/oauth-troubleshooting_2786d167/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/oauth-troubleshooting.md) | ⭐ 3.3k | `commercial` |
| [Supported Databases](commercial/supported-databases_d8a56a4d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/supported-databases.md) | ⭐ 3.3k | `commercial` |
| [Index](commercial/index_833de47b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/index.md) | ⭐ 3.3k | `commercial` |
| [Grpc Services](commercial/grpc-services_a969a94f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/grpc-services.md) | ⭐ 3.3k | `commercial` |
| [Query Param Auth](commercial/query-param-auth_e000b191/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/query-param-auth.md) | ⭐ 3.3k | `commercial` |
| [Diagram Templates](commercial/diagram-templates_0d34a579/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-core/skills/drawio-logical-diagrams/references/diagram-templates.md) | ⭐ 102 | `commercial` |
| [Caching Strategies](commercial/caching-strategies_d292ad8f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/caching-strategies.md) | ⭐ 102 | `commercial` |
| [Data Fetching](commercial/data-fetching_cd99416e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/data-fetching.md) | ⭐ 102 | `users` |
| [Streaming Suspense](commercial/streaming-suspense_d5acf02e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/streaming-suspense.md) | ⭐ 102 | `commercial` |
| [Ci Cd](commercial/ci-cd_3445fe57/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/ci-cd.md) | ⭐ 102 | `commercial` |
| [React](commercial/react_019731cb/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/react.md) | ⭐ 102 | `commercial` |
| [Ci Cd](commercial/ci-cd_850d9586/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/ci-cd.md) | ⭐ 102 | `commercial` |
| [Aws](commercial/aws_235d4424/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/aws.md) | ⭐ 3.3k | `commercial` |
| [Backup](commercial/backup_2ce0262a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/backup.md) | ⭐ 3.3k | `commercial` |

### Communication (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Schema](communication/schema_737de31c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/SCHEMA.md) | ⭐ 102 | `communication` |
| [Nestjs Setup](communication/nestjs-setup_476bb950/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/nestjs-setup.md) | ⭐ 102 | `communication` |
| [Nextjs Setup](communication/nextjs-setup_c674b3e1/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/nextjs-setup.md) | ⭐ 102 | `communication` |
| [Server Actions](communication/server-actions_4b85acf1/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/server-actions.md) | ⭐ 102 | `communication` |
| [Database Adapter](communication/database-adapter_0b083935/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-authentication/references/database-adapter.md) | ⭐ 102 | `communication` |
| [Rfc9728 Compliance](communication/rfc9728-compliance_0ee721f6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/rfc9728-compliance.md) | ⭐ 3.3k | `communication` |
| [Password Management](communication/password-management_bc5d8451/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/password-management.md) | ⭐ 3.3k | `communication` |

### Content Creation (21 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](content-creation/name-skill_df94dda4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/drizzle-orm-patterns/SKILL.md) | ⭐ 102 | `drizzle` `orm` `database` |
| [Skill](content-creation/name-skill_feaebcdb/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/SKILL.md) | ⭐ 102 | `nextjs` `next.js` `app-router` |
| [Skill](content-creation/name-skill_2cfa63f4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-data-fetching/SKILL.md) | ⭐ 102 | `posts` |
| [App Router Fundamentals](content-creation/app-router-fundamentals_5db2d98d/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/app-router-fundamentals.md) | ⭐ 102 | `content creation` |
| [Metadata Api](content-creation/metadata-api_464de7c9/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/metadata-api.md) | ⭐ 102 | `content creation` |
| [Nextjs16 Migration](content-creation/nextjs16-migration_807e73c6/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/nextjs16-migration.md) | ⭐ 102 | `content creation` |
| [React Query](content-creation/react-query_66bf81c6/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-data-fetching/references/REACT-QUERY.md) | ⭐ 102 | `content creation` |
| [Core Web Vitals](content-creation/core-web-vitals_c726a13b/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/core-web-vitals.md) | ⭐ 102 | `content creation` |
| [Image Optimization](content-creation/image-optimization_8a4da78f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/image-optimization.md) | ⭐ 102 | `content creation` |
| [Metadata Seo](content-creation/metadata-seo_eea42aa7/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/metadata-seo.md) | ⭐ 102 | `content creation` |
| [Nextjs 16 Patterns](content-creation/nextjs-16-patterns_1faa29b3/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/nextjs-16-patterns.md) | ⭐ 102 | `content creation` |
| [Nextjs Config](content-creation/nextjs-config_d4dbc259/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/nextjs-config.md) | ⭐ 102 | `content creation` |
| [Troubleshooting](content-creation/troubleshooting_7c93e51b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/troubleshooting.md) | ⭐ 3.3k | `content creation` |
| [Passthrough](content-creation/passthrough_6dcf2b83/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/passthrough.md) | ⭐ 3.3k | `content creation` |
| [Tool Annotations](content-creation/tool-annotations_9f3db504/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/tool-annotations.md) | ⭐ 3.3k | `content creation` |
| [038 Experimental Rust Transport Backend](content-creation/038-experimental-rust-transport-backend_0eab4f20/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/038-experimental-rust-transport-backend.md) | ⭐ 3.3k | `content creation` |
| [Index](content-creation/index_c4edf73d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/kit/index.md) | ⭐ 3.3k | `content creation` |
| [Index](content-creation/index_1e66c9e7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/press/index.md) | ⭐ 3.3k | `content creation` |
| [Chunker Server](content-creation/chunker-server_334bf908/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/chunker-server.md) | ⭐ 3.3k | `content creation` |
| [Eval Server](content-creation/eval-server_35d409dc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/eval-server.md) | ⭐ 3.3k | `content creation` |
| [Url To Markdown Server](content-creation/url-to-markdown-server_55cb28c7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/url-to-markdown-server.md) | ⭐ 3.3k | `content creation` |

### Daily Assistant (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Package Configs](daily-assistant/package-configs_66bf674a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/package-configs.md) | ⭐ 102 | `daily assistant` |
| [Crewai](daily-assistant/crewai_c93ce999/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/crewai.md) | ⭐ 3.3k | `daily assistant` |

### Data Analysis (20 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](data-analysis/name-skill_2414c4db/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-core/skills/drawio-logical-diagrams/SKILL.md) | ⭐ 102 | `data analysis` |
| [Raw Typescript Lambda](data-analysis/raw-typescript-lambda_c2fffb80/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/raw-typescript-lambda.md) | ⭐ 102 | `data analysis` |
| [Performance Architecture](data-analysis/performance-architecture_f977d60e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/performance-architecture.md) | ⭐ 3.3k | `data analysis` |
| [Dcr](data-analysis/dcr_14b0e9e4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/dcr.md) | ⭐ 3.3k | `data analysis` |
| [Postgres Upgrade Process](data-analysis/postgres-upgrade-process_60a48da6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/postgres-upgrade-process.md) | ⭐ 3.3k | `data analysis` |
| [007 Pluggable Cache Backend](data-analysis/007-pluggable-cache-backend_6ee2c5f9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/007-pluggable-cache-backend.md) | ⭐ 3.3k | `data analysis` |
| [Index](data-analysis/index_974fad2a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/index.md) | ⭐ 3.3k | `data analysis` |
| [Internal Observability](data-analysis/internal-observability_328ec0d8/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/internal-observability.md) | ⭐ 3.3k | `data analysis` |
| [Plugins](data-analysis/plugins_5f80238a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/plugins.md) | ⭐ 3.3k | `data analysis` |
| [Csv Pandas Chat Server](data-analysis/csv-pandas-chat-server_454d219f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/csv-pandas-chat-server.md) | ⭐ 3.3k | `data analysis` |
| [Data Analysis Server](data-analysis/data-analysis-server_cef0d00a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/data-analysis-server.md) | ⭐ 3.3k | `data analysis` |
| [Docx Server](data-analysis/docx-server_a79fb0e0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/docx-server.md) | ⭐ 3.3k | `data analysis` |
| [Graphviz Server](data-analysis/graphviz-server_b48cf167/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/graphviz-server.md) | ⭐ 3.3k | `data analysis` |
| [Plotly Server](data-analysis/plotly-server_3921f69a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/plotly-server.md) | ⭐ 3.3k | `data analysis` |
| [Pptx Server](data-analysis/pptx-server_9c6d336a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/pptx-server.md) | ⭐ 3.3k | `data analysis` |
| [Python Sandbox Server](data-analysis/python-sandbox-server_563cf18c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/python-sandbox-server.md) | ⭐ 3.3k | `data analysis` |
| [Nestjs](data-analysis/nestjs_ae32ef45/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/nestjs.md) | ⭐ 102 | `data analysis` |
| [Charts](data-analysis/charts_99941b48/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/deployment/CHARTS.md) | ⭐ 3.3k | `data analysis` |
| [Performance Strategy](data-analysis/performance_strategy_821a3c9c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/tests/performance/PERFORMANCE_STRATEGY.md) | ⭐ 3.3k | `data analysis` |
| [Xlsx Server](data-analysis/xlsx-server_2e1df286/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/xlsx-server.md) | ⭐ 3.3k | `data analysis` |

### Development (65 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Multitenancy](development/multitenancy_d26020b1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/multitenancy.md) | ⭐ 3.3k | `development` |
| [Plugins](development/plugins_b81593c7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins.md) | ⭐ 3.3k | `openai` `moderation` `content-safety` |
| [Selecting An Mcp Gateway](development/selecting-an-mcp-gateway_cec67c55/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/selecting-an-mcp-gateway.md) | ⭐ 3.3k | `development` |
| [Postgresql Schema Configuration](development/postgresql-schema-configuration_9f6cd11d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/postgresql-schema-configuration.md) | ⭐ 3.3k | `development` |
| [Mcp Developer Guide Json Rpc](development/mcp-developer-guide-json-rpc_2bfea8c1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/mcp-developer-guide-json-rpc.md) | ⭐ 3.3k | `development` |
| [Profiling](development/profiling_d8b52edd/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/profiling.md) | ⭐ 3.3k | `development` |
| [Sso Entra Role Mapping](development/sso-entra-role-mapping_c6294cc1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-entra-role-mapping.md) | ⭐ 3.3k | `development` |
| [Sso Github Tutorial](development/sso-github-tutorial_68c131fe/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-github-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Ibm Tutorial](development/sso-ibm-tutorial_61b0c795/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-ibm-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Microsoft Entra Id Tutorial](development/sso-microsoft-entra-id-tutorial_e579028c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-microsoft-entra-id-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Okta Tutorial](development/sso-okta-tutorial_76fd8f99/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-okta-tutorial.md) | ⭐ 3.3k | `development` |
| [Teams](development/teams_c57a4406/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/teams.md) | ⭐ 3.3k | `development` |
| [003 Expose Multi Transport Endpoints](development/003-expose-multi-transport-endpoints_9637d88e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/003-expose-multi-transport-endpoints.md) | ⭐ 3.3k | `development` |
| [005 Vscode Devcontainer Support](development/005-vscode-devcontainer-support_6e77ea8b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/005-vscode-devcontainer-support.md) | ⭐ 3.3k | `development` |
| [008 Federation Discovery](development/008-federation-discovery_753ace5b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/008-federation-discovery.md) | ⭐ 3.3k | `development` |
| [010 Observability Prometheus](development/010-observability-prometheus_fcec712d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/010-observability-prometheus.md) | ⭐ 3.3k | `development` |
| [011 Tool Federation](development/011-tool-federation_24908a26/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/011-tool-federation.md) | ⭐ 3.3k | `development` |
| [016 Plugin Framework Ai Middleware](development/016-plugin-framework-ai-middleware_f40a7f2c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/016-plugin-framework-ai-middleware.md) | ⭐ 3.3k | `development` |
| [017 Adopt Orjson Json Serialization](development/017-adopt-orjson-json-serialization_7a66024e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/017-adopt-orjson-json-serialization.md) | ⭐ 3.3k | `development` |
| [022 Elicitation Passthrough Implementation](development/022-elicitation-passthrough-implementation_288c5524/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/022-elicitation-passthrough-implementation.md) | ⭐ 3.3k | `development` |
| [024 Uvicorn Standard Extras](development/024-uvicorn-standard-extras_6d24b0a8/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/024-uvicorn-standard-extras.md) | ⭐ 3.3k | `development` |
| [040 Flexible Admin Ui Sections](development/040-flexible-admin-ui-sections_361375d2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/040-flexible-admin-ui-sections.md) | ⭐ 3.3k | `development` |
| [Gateway Hooks](development/gateway-hooks_9aa359bb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins/gateway-hooks.md) | ⭐ 3.3k | `development` |
| [Security Hooks](development/security-hooks_aa01ca2a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins/security-hooks.md) | ⭐ 3.3k | `development` |
| [A2A](development/a2a_9a2ca7a0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/a2a.md) | ⭐ 3.3k | `development` |
| [Cline](development/cline_8713195d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/cline.md) | ⭐ 3.3k | `development` |
| [Llm Chat](development/llm-chat_b4fec4c4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/llm-chat.md) | ⭐ 3.3k | `development` |
| [Index](development/index_79b6cf6d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/index.md) | ⭐ 3.3k | `custom` `filter` |
| [Lifecycle](development/lifecycle_17a38751/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/lifecycle.md) | ⭐ 3.3k | `plugin` |
| [Rust Plugins](development/rust-plugins_b6ecb04f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/rust-plugins.md) | ⭐ 3.3k | `development` |
| [Index](development/index_817e60d4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/index.md) | ⭐ 3.3k | `development` |
| [Calculator Server](development/calculator-server_84a343ea/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/calculator-server.md) | ⭐ 3.3k | `development` |
| [Box](development/box_1e5cdfee/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/box/box.md) | ⭐ 3.3k | `development` |
| [Monday Mcp](development/monday-mcp_91796fb7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/monday/monday-mcp.md) | ⭐ 3.3k | `development` |
| [Installation](development/installation_38060e27/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-core/docs/installation.md) | ⭐ 102 | `development` |
| [Guide Skills Architecture](development/guide-skills-architecture_ffe890a0/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/docs/guide-skills-architecture.md) | ⭐ 102 | `development` |
| [Skill](development/name-skill_e5618400/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/SKILL.md) | ⭐ 102 | `authentication` `better-auth` `nestjs` |
| [Skill](development/name-skill_0aa4f73f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-authentication/SKILL.md) | ⭐ 102 | `development` |
| [Skill](development/name-skill_4ad4e6d6/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/SKILL.md) | ⭐ 102 | `products` |
| [Shape Styles](development/shape-styles_d578b142/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-core/skills/drawio-logical-diagrams/references/shape-styles.md) | ⭐ 102 | `development` |
| [Micronaut Lambda](development/micronaut-lambda_0d942b2e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/micronaut-lambda.md) | ⭐ 102 | `development` |
| [Bref Lambda](development/bref-lambda_abe50d59/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/references/bref-lambda.md) | ⭐ 102 | `development` |
| [Raw Php Lambda](development/raw-php-lambda_af10d762/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/references/raw-php-lambda.md) | ⭐ 102 | `development` |
| [Testing Lambda](development/testing-lambda_d82074a4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/references/testing-lambda.md) | ⭐ 102 | `development` |
| [Raw Python Lambda](development/raw-python-lambda_15f76685/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/raw-python-lambda.md) | ⭐ 102 | `development` |
| [Fastify Adapter](development/fastify-adapter_b2dd6c13/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/fastify-adapter.md) | ⭐ 102 | `development` |
| [Nestjs Lambda](development/nestjs-lambda_932b0d76/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/nestjs-lambda.md) | ⭐ 102 | `development` |
| [Testing](development/testing_7e93cab4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/testing.md) | ⭐ 102 | `development` |
| [Plugins](development/plugins_908cf704/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/PLUGINS.md) | ⭐ 102 | `development` |
| [Mfa 2Fa](development/mfa-2fa_9fe02da9/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/mfa-2fa.md) | ⭐ 102 | `development` |
| [Social Providers](development/social-providers_9d4193d7/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/social-providers.md) | ⭐ 102 | `development` |
| [Authjs Setup](development/authjs-setup_ea51025c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-authentication/references/authjs-setup.md) | ⭐ 102 | `development` |
| [Oauth Providers](development/oauth-providers_dd5610b4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-authentication/references/oauth-providers.md) | ⭐ 102 | `development` |
| [Deployment Platforms](development/deployment-platforms_113fe037/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/references/deployment-platforms.md) | ⭐ 102 | `development` |
| [Monitoring](development/monitoring_ef29a8d8/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/references/monitoring.md) | ⭐ 102 | `development` |
| [Api Routes](development/api-routes_ad607ee3/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/api-routes.md) | ⭐ 102 | `development` |
| [Bundle Optimization](development/bundle-optimization_bf03068e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/bundle-optimization.md) | ⭐ 102 | `development` |
| [Font Optimization](development/font-optimization_2269aa12/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/font-optimization.md) | ⭐ 102 | `development` |
| [Agents](development/agents_5970f48f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/AGENTS.md) | ⭐ 3.3k | `development` |
| [Api](development/api_08409a71/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/api.md) | ⭐ 3.3k | `development` |
| [Audit Db Transaction Management](development/audit-db-transaction-management_caf0c54e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/audit-db-transaction-management.md) | ⭐ 3.3k | `development` |
| [Doctest Coverage](development/doctest-coverage_83608663/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/doctest-coverage.md) | ⭐ 3.3k | `development` |
| [Acceptance](development/acceptance_cba125ac/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/testing/acceptance.md) | ⭐ 3.3k | `development` |
| [Semantic Kernel](development/semantic-kernel_2983535d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/semantic-kernel.md) | ⭐ 3.3k | `development` |
| [Code Splitter Server](development/code-splitter-server_0c8d6e48/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/code-splitter-server.md) | ⭐ 3.3k | `development` |

### Development/Devops (109 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agent Deployment Guide](development/devops/229-agent-deployment-guide_fc15347e/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-deploy/references/agent-deployment-guide.md) | ⭐ 80 | `development` |
| [Guide Skills Frontend](development/devops/guide-skills-frontend_2f63a088/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/docs/guide-skills-frontend.md) | ⭐ 102 | `development` |
| [Guide Skills Monorepo](development/devops/guide-skills-monorepo_83c4b53f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/docs/guide-skills-monorepo.md) | ⭐ 102 | `development` |
| [Skill](development/devops/name-skill_1ec7fc18/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `java` |
| [Skill](development/devops/name-skill_9bc6e00e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `php` |
| [Skill](development/devops/name-skill_3c0620fc/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `python` |
| [Skill](development/devops/name-skill_2658dda5/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `typescript` |
| [Skill](development/devops/name-skill_62d193fa/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/SKILL.md) | ⭐ 102 | `nextjs` `next.js` `deployment` |
| [Skill](development/devops/name-skill_d0e8129e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/SKILL.md) | ⭐ 102 | `nx` `monorepo` `typescript` |
| [Skill](development/devops/name-skill_4133613a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/SKILL.md) | ⭐ 102 | `development` |
| [Raw Java Lambda](development/devops/raw-java-lambda_bc77ac25/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/raw-java-lambda.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/serverless-deployment_7a94e770/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/serverless-deployment_0963c52b/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Chalice Lambda](development/devops/chalice-lambda_73723400/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/chalice-lambda.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/serverless-deployment_64158df7/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/serverless-deployment_5ae8bdfe/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Docker Patterns](development/devops/docker-patterns_6cc95fae/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/references/docker-patterns.md) | ⭐ 102 | `development` |
| [Basics](development/devops/basics_94f6201e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/basics.md) | ⭐ 102 | `development` |
| [Mcp Server Python](development/devops/mcp-server-python_5f654215/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/mcp-server-python.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_ccaab12c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/index.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_3c572c12/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/index.md) | ⭐ 3.3k | `development` |
| [Roadmap](development/devops/roadmap_c9aa8183/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/roadmap.md) | ⭐ 3.3k | `development` |
| [Security Features](development/devops/security-features_012c8460/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/security-features.md) | ⭐ 3.3k | `development` |
| [Input Validation](development/devops/input-validation_01d87216/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/input-validation.md) | ⭐ 3.3k | `development` |
| [Mcp Architecture Patterns](development/devops/mcp-architecture-patterns_230d8859/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/mcp-architecture-patterns.md) | ⭐ 3.3k | `development` |
| [Argocd](development/devops/argocd_31a49405/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/argocd.md) | ⭐ 3.3k | `development` |
| [Cforge Gateway](development/devops/cforge-gateway_ac5a8cdb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/cforge-gateway.md) | ⭐ 3.3k | `security` `filter` |
| [Compose](development/devops/compose_366cfc28/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/compose.md) | ⭐ 3.3k | `development` |
| [Fly Io](development/devops/fly-io_e5ccc7e2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/fly-io.md) | ⭐ 3.3k | `development` |
| [Google Cloud Run](development/devops/google-cloud-run_8be8613f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/google-cloud-run.md) | ⭐ 3.3k | `development` |
| [Helm](development/devops/helm_2945bf9d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/helm.md) | ⭐ 3.3k | `development` |
| [Ibm Code Engine](development/devops/ibm-code-engine_fe3cb3dc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/ibm-code-engine.md) | ⭐ 3.3k | `development` |
| [Local](development/devops/local_93c0d987/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/local.md) | ⭐ 3.3k | `development` |
| [Minikube](development/devops/minikube_0c8a4b03/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/minikube.md) | ⭐ 3.3k | `development` |
| [Proxy Auth](development/devops/proxy-auth_7d43d3de/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/proxy-auth.md) | ⭐ 3.3k | `development` |
| [Tls Configuration](development/devops/tls-configuration_54ff4560/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/tls-configuration.md) | ⭐ 3.3k | `development` |
| [Developer Onboarding](development/devops/developer-onboarding_61358b2f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/developer-onboarding.md) | ⭐ 3.3k | `development` |
| [Github](development/devops/github_487e80b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/github.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_a20a1f7e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/index.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_0010150a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/faq/index.md) | ⭐ 3.3k | `development` |
| [Admin Ui Customization](development/devops/admin-ui-customization_b4697c47/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/admin-ui-customization.md) | ⭐ 3.3k | `development` |
| [Api Usage](development/devops/api-usage_930cc30b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/api-usage.md) | ⭐ 3.3k | `development` |
| [Configuration](development/devops/configuration_ed70f974/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/configuration.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_7c9efd7f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/index.md) | ⭐ 3.3k | `development` |
| [Logging](development/devops/logging_7601ca0a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/logging.md) | ⭐ 3.3k | `development` |
| [Proxy](development/devops/proxy_5a83e05a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/proxy.md) | ⭐ 3.3k | `development` |
| [Rbac](development/devops/rbac_46494f27/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/rbac.md) | ⭐ 3.3k | `development` |
| [Scale](development/devops/scale_5375b413/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/scale.md) | ⭐ 3.3k | `development` |
| [Securing](development/devops/securing_264f96d4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/securing.md) | ⭐ 3.3k | `development` |
| [Self Signed Certificates](development/devops/self-signed-certificates_c29006e4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/self-signed-certificates.md) | ⭐ 3.3k | `development` |
| [Sso Generic Oidc Tutorial](development/devops/sso-generic-oidc-tutorial_1cf7119b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-generic-oidc-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Google Tutorial](development/devops/sso-google-tutorial_ec0b87b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-google-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Keycloak Tutorial](development/devops/sso-keycloak-tutorial_bc429daa/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-keycloak-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso](development/devops/sso_3ff6196e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso.md) | ⭐ 3.3k | `development` |
| [Tuning](development/devops/tuning_b13077bc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/tuning.md) | ⭐ 3.3k | `development` |
| [Upgrade](development/devops/upgrade_b3150af6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/upgrade.md) | ⭐ 3.3k | `development` |
| [Well Known Uris](development/devops/well-known-uris_b7a6cdf4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/well-known-uris.md) | ⭐ 3.3k | `development` |
| [Config Validation](development/devops/config-validation_dc2663b4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/operations/config-validation.md) | ⭐ 3.3k | `development` |
| [Cpu Spin Loop Mitigation](development/devops/cpu-spin-loop-mitigation_fdb3b575/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/operations/cpu-spin-loop-mitigation.md) | ⭐ 3.3k | `development` |
| [Features](development/devops/features_4c3fbb11/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/features.md) | ⭐ 3.3k | `development` |
| [Admin UI](development/devops/ui_223c1cd2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/ui.md) | ⭐ 3.3k | `development` |
| [Argocd Helm Deployment Ibm Cloud Iks](development/devops/argocd-helm-deployment-ibm-cloud-iks_8a45b06f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/argocd-helm-deployment-ibm-cloud-iks.md) | ⭐ 3.3k | `development` |
| [Openwebui Tutorial](development/devops/openwebui-tutorial_ace5026d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/openwebui-tutorial.md) | ⭐ 3.3k | `development` |
| [Mcpgateway Translate](development/devops/mcpgateway-translate_199775f5/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/mcpgateway-translate.md) | ⭐ 3.3k | `development` |
| [Multi Auth Headers](development/devops/multi-auth-headers_d4462e80/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/multi-auth-headers.md) | ⭐ 3.3k | `development` |
| [Reverse Proxy](development/devops/reverse-proxy_a49bea67/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/reverse-proxy.md) | ⭐ 3.3k | `development` |
| [014 Security Headers Cors Middleware](development/devops/014-security-headers-cors-middleware_009499c3/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/014-security-headers-cors-middleware.md) | ⭐ 3.3k | `development` |
| [015 Well Known Uri Handler](development/devops/015-well-known-uri-handler_198ce9fb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/015-well-known-uri-handler.md) | ⭐ 3.3k | `development` |
| [018 Built In Response Compression](development/devops/018-built-in-response-compression_161e1ccb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/018-built-in-response-compression.md) | ⭐ 3.3k | `development` |
| [019 Modular Architecture Split](development/devops/019-modular-architecture-split_606af208/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/019-modular-architecture-split.md) | ⭐ 3.3k | `development` |
| [025 Granian Http Server](development/devops/025-granian-http-server_7849918e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/025-granian-http-server.md) | ⭐ 3.3k | `development` |
| [032 Mcp Session Pool](development/devops/032-mcp-session-pool_fd39527c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/032-mcp-session-pool.md) | ⭐ 3.3k | `development` |
| [036 Bootstrap Custom Roles](development/devops/036-bootstrap-custom-roles_98fdf334/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/036-bootstrap-custom-roles.md) | ⭐ 3.3k | `development` |
| [038 Multi Worker Session Affinity](development/devops/038-multi-worker-session-affinity_de1547c3/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/038-multi-worker-session-affinity.md) | ⭐ 3.3k | `development` |
| [Observability](development/devops/observability_358deba7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/observability.md) | ⭐ 3.3k | `development` |
| [Phoenix](development/devops/phoenix_a6b52052/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/phoenix.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_49786343/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/social/index.md) | ⭐ 3.3k | `development` |
| [Bee](development/devops/bee_bd521c8a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/bee.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_6da5357e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/index.md) | ⭐ 3.3k | `development` |
| [Claude Desktop](development/devops/claude-desktop_33f0496c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/claude-desktop.md) | ⭐ 3.3k | `development` |
| [Continue](development/devops/continue_4f371801/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/continue.md) | ⭐ 3.3k | `development` |
| [Grpc Transport](development/devops/grpc-transport_c50d7487/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/grpc-transport.md) | ⭐ 3.3k | `development` |
| [Mtls](development/devops/mtls_54eac01a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/mtls.md) | ⭐ 3.3k | `development` |
| [Terraform](development/devops/terraform_86d5b2a9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/hashicorp/terraform.md) | ⭐ 3.3k | `development` |
| [Langflow Server](development/devops/langflow-server_4739c90c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/third-party/langflow-server.md) | ⭐ 3.3k | `development` |
| [Github](development/devops/github_cda2891c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/microsoft/github.md) | ⭐ 3.3k | `development` |
| [Testing Lambda](development/devops/testing-lambda_cb2bfa70/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/testing-lambda.md) | ⭐ 102 | `development` |
| [Testing Lambda](development/devops/testing-lambda_67a7202d/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/testing-lambda.md) | ⭐ 102 | `development` |
| [Serverless Config](development/devops/serverless-config_abb2a10a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/serverless-config.md) | ⭐ 102 | `development` |
| [Github Actions](development/devops/github-actions_c3f7d0df/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/references/github-actions.md) | ⭐ 102 | `v*` |
| [Testing Config](development/devops/testing-config_f6011718/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/testing-config.md) | ⭐ 102 | `development` |
| [Quickstart](development/devops/quickstart_4b1a224b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/plugins_rust/QUICKSTART.md) | ⭐ 3.3k | `development` |
| [Manual Testing](development/devops/manual_testing_ce3c70ad/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/tests/performance/MANUAL_TESTING.md) | ⭐ 3.3k | `development` |
| [Azure](development/devops/azure_0c39cc0f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/azure.md) | ⭐ 3.3k | `development` |
| [Container](development/devops/container_18dafad2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/container.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_27e0ccd7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/index.md) | ⭐ 3.3k | `development` |
| [Openshift](development/devops/openshift_e3f6d529/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/openshift.md) | ⭐ 3.3k | `development` |
| [Packaging](development/devops/packaging_61cc17c6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/packaging.md) | ⭐ 3.3k | `development` |
| [Logging Examples](development/devops/logging-examples_25339325/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/logging-examples.md) | ⭐ 3.3k | `development` |
| [Observability](development/devops/observability_e6a7dc7d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability.md) | ⭐ 3.3k | `development` |
| [Quick Start](development/devops/quick_start_885d4bd6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/quick_start.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_f25773c2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/testing/index.md) | ⭐ 3.3k | `development` |
| [Performance](development/devops/performance_d50c481e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/testing/performance.md) | ⭐ 3.3k | `development` |
| [Dcr Hyprmcp](development/devops/dcr-hyprmcp_7ceee4b1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/dcr-hyprmcp.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/index_a595636a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/index.md) | ⭐ 3.3k | `development` |
| [009 Built In Health Checks](development/devops/009-built-in-health-checks_05433232/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/009-built-in-health-checks.md) | ⭐ 3.3k | `development` |
| [026 Hiredis Redis Parser](development/devops/026-hiredis-redis-parser_7c184216/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/026-hiredis-redis-parser.md) | ⭐ 3.3k | `development` |
| [Copilot](development/devops/copilot_ebc99473/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/copilot.md) | ⭐ 3.3k | `development` |
| [Openwebui](development/devops/openwebui_fd32a81d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/openwebui.md) | ⭐ 3.3k | `development` |

### Development/Testing (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/testing/name-skill_187b30c5/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nestjs-drizzle-crud-generator/SKILL.md) | ⭐ 102 | `development` |
| [Nestjs Config](development/testing/nestjs-config_ab544f6a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/nestjs-config.md) | ⭐ 102 | `development` |
| [Db Performance](development/testing/db-performance_9aa0f078/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/db-performance.md) | ⭐ 3.3k | `development` |
| [Testing Patterns](development/testing/testing-patterns_f02f2d5c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-authentication/references/testing-patterns.md) | ⭐ 102 | `development` |
| [Generators](development/testing/generators_7452a2df/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/generators.md) | ⭐ 102 | `development` |
| [Testing](development/testing/testing_6ef3e854/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/TESTING.md) | ⭐ 3.3k | `development` |
| [Testing](development/testing/testing_b7b74390/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/plugins/webhook_notification/TESTING.md) | ⭐ 3.3k | `development` |

### Development/Tools (34 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Typescript](development/tools/typescript_d85d7b93/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/typescript.md) | ⭐ 102 | `development` |
| [Skill](development/tools/name-skill_034a479c/) | [netresearch/jira-skill](https://raw.githubusercontent.com/netresearch/jira-skill/main/skills/jira-communication/SKILL.md) | ⭐ 24 | `development` |
| [Developing](development/tools/developing_4eef14b0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/DEVELOPING.md) | ⭐ 3.3k | `development` |
| [Enable Payload Logging](development/tools/enable_payload_logging_c75da92d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/enable_payload_logging.md) | ⭐ 3.3k | `development` |
| [Mcpgateway](development/tools/mcpgateway_f02247c0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/mcpgateway.md) | ⭐ 3.3k | `development` |
| [Plugins Llms](development/tools/plugins-llms_8fca01fe/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/plugins-llms.md) | ⭐ 3.3k | `development` |
| [Oauth Authorization Code Ui Design](development/tools/oauth-authorization-code-ui-design_ec055b8d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/oauth-authorization-code-ui-design.md) | ⭐ 3.3k | `development` |
| [Oauth Design](development/tools/oauth-design_dd0d2a41/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/oauth-design.md) | ⭐ 3.3k | `development` |
| [Bulk Import](development/tools/bulk-import_8f7144ba/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/bulk-import.md) | ⭐ 3.3k | `development` |
| [Export Import Reference](development/tools/export-import-reference_cd3b304c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import-reference.md) | ⭐ 3.3k | `development` |
| [Metadata Tracking](development/tools/metadata-tracking_93c516b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/metadata-tracking.md) | ⭐ 3.3k | `development` |
| [Oauth](development/tools/oauth_2ff3a5f9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/oauth.md) | ⭐ 3.3k | `development` |
| [Tags](development/tools/tags_0898f285/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/tags.md) | ⭐ 3.3k | `development` |
| [006 Gateway Tool Rate Limiting](development/tools/006-gateway-tool-rate-limiting_bb786b5a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/006-gateway-tool-rate-limiting.md) | ⭐ 3.3k | `development` |
| [021 Built In Proxy Vs Service Mesh](development/tools/021-built-in-proxy-vs-service-mesh_87238005/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/021-built-in-proxy-vs-service-mesh.md) | ⭐ 3.3k | `development` |
| [023 One Time Authentication Servers](development/tools/023-one-time-authentication-servers_6f4caffb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/023-one-time-authentication-servers.md) | ⭐ 3.3k | `development` |
| [027 Migrate Psycopg3](development/tools/027-migrate-psycopg3_ae99362f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/027-migrate-psycopg3.md) | ⭐ 3.3k | `development` |
| [035 Query Parameter Authentication](development/tools/035-query-parameter-authentication_7da5eb7c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/035-query-parameter-authentication.md) | ⭐ 3.3k | `development` |
| [Llamaindex](development/tools/llamaindex_332803db/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/llamaindex.md) | ⭐ 3.3k | `development` |
| [Http Auth Hooks](development/tools/http-auth-hooks_18dede4f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/http-auth-hooks.md) | ⭐ 3.3k | `development` |
| [Fast Time Server](development/tools/fast-time-server_c35f9972/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/fast-time-server.md) | ⭐ 3.3k | `development` |
| [Instana](development/tools/instana_c81a89af/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/ibm/instana.md) | ⭐ 3.3k | `development` |
| [Mermaid Server](development/tools/mermaid-server_f0be2d55/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/mermaid-server.md) | ⭐ 3.3k | `development` |
| [Index](development/tools/index_ec82d631/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/open/index.md) | ⭐ 3.3k | `development` |
| [Basic](development/tools/basic_5cb3231f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/testing/basic.md) | ⭐ 3.3k | `development` |
| [Fuzzing](development/tools/fuzzing_cb54bd1f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/testing/fuzzing.md) | ⭐ 3.3k | `development` |
| [Unittest](development/tools/unittest_dee6b743/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/testing/unittest.md) | ⭐ 3.3k | `development` |
| [Index](development/tools/index_fb8ba861/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/index.md) | ⭐ 3.3k | `development` |
| [Testing](development/tools/testing_def21240/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/mcp-servers/python/output_schema_test_server/TESTING.md) | ⭐ 3.3k | `development` |
| [001 Adopt Fastapi Pydantic](development/tools/001-adopt-fastapi-pydantic_af0cde87/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/001-adopt-fastapi-pydantic.md) | ⭐ 3.3k | `development` |
| [Autogen](development/tools/autogen_354daee9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/autogen.md) | ⭐ 3.3k | `development` |
| [Langchain](development/tools/langchain_ec9efd8f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/langchain.md) | ⭐ 3.3k | `development` |
| [Langgraph](development/tools/langgraph_bf2b363a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/langgraph.md) | ⭐ 3.3k | `development` |
| [Index](development/tools/index_0af3c028/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/index.md) | ⭐ 3.3k | `development` |

### Productivity (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Parallel Session Cleanup](productivity/parallel-session-cleanup_d83f616a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/parallel-session-cleanup.md) | ⭐ 3.3k | `productivity` |

### Research (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Pandoc Server](research/pandoc-server_a65d4d68/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/pandoc-server.md) | ⭐ 3.3k | `research` |
| [Latex Server](research/latex-server_e6f7093b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/latex-server.md) | ⭐ 3.3k | `research` |

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

*Last updated: 2026-02-21 19:40:24 UTC*
*Automatically maintained by SkillFlow*
