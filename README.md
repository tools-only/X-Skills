# X-Skills

A curated collection of **234 AI-powered skills** organized into 14 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (26 skills)
- **Automation/Workflow** (22 skills)
- **Commercial** (17 skills)
- **Communication** (3 skills)
- **Content Creation** (29 skills)
- **Daily Assistant** (19 skills)
- **Data Analysis** (13 skills)
- **Development** (46 skills)
- **Development/Devops** (24 skills)
- **Development/Testing** (1 skill)
- **Development/Tools** (13 skills)
- **Investment** (7 skills)
- **Productivity** (8 skills)
- **Research** (6 skills)

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


### Automation/Scripting (26 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/scripting/name-skill_c94cab98/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-speech-to-text-rest-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_18050c5d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-queue-ts/SKILL.md) | 🔥 15.4k | `automation` |
| [Web Spec](automation/scripting/web-spec_e7518f36/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/web-spec.md) | ⭐ 95 | `automation` |
| [Reference](automation/scripting/reference_704f68d2/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/cf-edit/reference.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_e622e2f6/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/mxl-info/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_e110d174/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/skd-info/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_aa22b354/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-textanalytics-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_834b27de/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-translation-document-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_212882c8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-compute-batch-java/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_70dd7219/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-data-tables-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_db3eb0ef/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-eventgrid-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_3023e815/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-eventhub-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_fe6cf97f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-ingestion-java/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_1afe3a01/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-query-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_54c897f4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-search-documents-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_34944c09/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-servicebus-dotnet/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_7f6d2f34/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/database-migrations-sql-migrations/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/name-skill_395e980f/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/cf-init/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_846f3bc5/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/cfe-borrow/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_8a6bda74/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/cfe-validate/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_b90e6dd1/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/img-grid/SKILL.md) | ⭐ 95 | `automation` |
| [Child Operations](automation/scripting/child-operations_cfe2e2f1/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/meta-edit/child-operations.md) | ⭐ 95 | `automation` |
| [Properties Reference](automation/scripting/properties-reference_fbd8941c/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/meta-edit/properties-reference.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_bf90d9db/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/subsystem-validate/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_2b2c8262/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/web-info/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/name-skill_f34b65e7/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/web-stop/SKILL.md) | ⭐ 95 | `automation` |

### Automation/Workflow (22 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/name-skill_1b8e6a33/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ai-agent-development/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_0d6abb80/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ai-ml/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_46db224b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/aws-cost-optimizer/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_4472d1b2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-ml-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_1d0e6f8c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-messaging-webpubsubservice-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_ec3c62a4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-fabric-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_fd49fef2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/business-analyst/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_274d2321/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/database/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_8d459ba9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/e2e-testing/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_0847c16f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hybrid-cloud-architect/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_58f43b98/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/incident-responder/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_e11ae46a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/m365-agents-dotnet/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_5b803d51/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/office-productivity/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_e301bd77/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/postgresql-optimization/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_78bec55c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security-audit/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_7f586da1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/terraform-infrastructure/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_ed36a8dc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/web-security-testing/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_3fede59c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/libreoffice/base/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_b3725c56/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/libreoffice/draw/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_0e17d0c7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/libreoffice/writer/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/workflow/name-skill_7c1f4e2d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security/aws-compliance-checker/SKILL.md) | 🔥 15.4k | `aws` `compliance` `audit` |
| [Skill](automation/workflow/name-skill_42f96e8a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security/aws-iam-best-practices/SKILL.md) | 🔥 15.4k | `aws` `iam` `security` |

### Commercial (17 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/andruia-skill_bc3d6ec8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/20-andruia-niche-intelligence/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_76bcd7f7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-appconfiguration-java/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_f9d2fd81/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-apimanagement-py/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_1a0ca218/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-sql-dotnet/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_c41137ec/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/customs-trade-compliance/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_cdd40689/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/market-sizing-analysis/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_1204d25e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/returns-reverse-logistics/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_c147e6cf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/shopify-development/SKILL.md) | 🔥 15.4k | `commercial` |
| [Skill](commercial/name-skill_d63cb7e3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-business-analyst-financial-projections/SKILL.md) | 🔥 15.4k | `commercial` |
| [Communication Templates](commercial/communication-templates_438a48cf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/carrier-relationship-management/references/communication-templates.md) | 🔥 15.4k | `commercial` |
| [Communication Templates](commercial/communication-templates_df567df3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/customs-trade-compliance/references/communication-templates.md) | 🔥 15.4k | `commercial` |
| [Edge Cases](commercial/edge-cases_d67b9803/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/customs-trade-compliance/references/edge-cases.md) | 🔥 15.4k | `commercial` |
| [Decision Frameworks](commercial/decision-frameworks_46db9d31/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/logistics-exception-management/references/decision-frameworks.md) | 🔥 15.4k | `commercial` |
| [Communication Templates](commercial/communication-templates_807c2d2f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/production-scheduling/references/communication-templates.md) | 🔥 15.4k | `commercial` |
| [Edge Cases](commercial/edge-cases_ba2efa50/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/production-scheduling/references/edge-cases.md) | 🔥 15.4k | `commercial` |
| [Communication Templates](commercial/communication-templates_199f3df6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/returns-reverse-logistics/references/communication-templates.md) | 🔥 15.4k | `commercial` |
| [Decision Frameworks](commercial/decision-frameworks_e707229f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/returns-reverse-logistics/references/decision-frameworks.md) | 🔥 15.4k | `commercial` |

### Communication (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](communication/name-skill_2b647276/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-java/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/name-skill_45420053/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-arizeaiobservabilityeval-dotnet/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/name-skill_373d16f0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-botservice-py/SKILL.md) | 🔥 15.4k | `communication` |

### Content Creation (29 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](content-creation/name-skill_3e14f061/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-contentunderstanding-py/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_5df3e2ea/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-file-datalake-py/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_49bd2701/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-file-share-ts/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_ec5b8b01/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/copywriting/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_8239255c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-technologies/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_fe6b5383/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/page-cro/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_5e603001/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/programmatic-seo/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_9fbb8949/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-content-refresher/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_9566e909/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-fundamentals/SKILL.md) | 🔥 15.4k | `content creation` |
| [Skill](content-creation/name-skill_383f8e6a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-structure-architect/SKILL.md) | 🔥 15.4k | `content creation` |
| [Edge Cases](content-creation/edge-cases_1ae015db/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/quality-nonconformance/references/edge-cases.md) | 🔥 15.4k | `content creation` |
| [Iteration Guide](content-creation/iteration-guide_b1a688bc/) | [tripleyak/SkillForge](https://raw.githubusercontent.com/tripleyak/SkillForge/main/references/iteration-guide.md) | ⭐ 516 | `content creation` |
| [1C Extension Spec](content-creation/1c-extension-spec_65b84edc/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/1c-extension-spec.md) | ⭐ 95 | `content creation` |
| [Cf Guide](content-creation/cf-guide_9240bd02/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/cf-guide.md) | ⭐ 95 | `content creation` |
| [Form Dsl Spec](content-creation/form-dsl-spec_8193cb94/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/form-dsl-spec.md) | ⭐ 95 | `content creation` |
| [Form Patterns](content-creation/form-patterns_62560bb2/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/form-patterns.md) | ⭐ 95 | `content creation` |
| [Role Dsl Spec](content-creation/role-dsl-spec_6bada2b0/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/role-dsl-spec.md) | ⭐ 95 | `content creation` |
| [Subsystem Guide](content-creation/subsystem-guide_884e2548/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/subsystem-guide.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_dfc4168f/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/epf-add-form/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_f969c4ad/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/form-compile/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_7944b047/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/form-patterns/SKILL.md) | ⭐ 95 | `content creation` |
| [Reference](content-creation/reference_852f6a58/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/interface-edit/reference.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_f0e6fadc/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/mxl-compile/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_d93afdee/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/role-compile/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_0f35115f/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/role-info/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_fe226ce7/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/skd-compile/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_9f965a62/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/subsystem-compile/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_d66b2b98/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/subsystem-edit/SKILL.md) | ⭐ 95 | `content creation` |
| [Skill](content-creation/name-skill_2417e66d/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/template-add/SKILL.md) | ⭐ 95 | `content creation` |

### Daily Assistant (19 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Edge Cases](daily-assistant/288-edge-cases_50a38ae5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/inventory-demand-planning/references/edge-cases.md) | 🔥 15.4k | `daily assistant` |
| [Decision Frameworks](daily-assistant/289-decision-frameworks_90e5dd3e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/production-scheduling/references/decision-frameworks.md) | 🔥 15.4k | `daily assistant` |
| [Decision Frameworks](daily-assistant/289-decision-frameworks_e50a63c7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/quality-nonconformance/references/decision-frameworks.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_76644efe/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/legal-advisor/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_01c2e76e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/risk-manager/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_25fc0d0a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-authority-builder/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_ac125d57/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-meta-optimizer/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_27ccbbd6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/track-management/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Decision Frameworks](daily-assistant/decision-frameworks_750ee441/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/inventory-demand-planning/references/decision-frameworks.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_dc345072/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/c4-component/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_2e2a1f5f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/c4-context/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_ca657afb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/competitive-landscape/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_8b65100d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conductor-validator/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_0a25874f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/debugger/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_82bbace0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/linux-troubleshooting/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Communication Templates](daily-assistant/communication-templates_b9a597dd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/quality-nonconformance/references/communication-templates.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/name-skill_3ab9a0b1/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/cfe-init/SKILL.md) | ⭐ 95 | `daily assistant` |
| [Skill](daily-assistant/name-skill_da832524/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/epf-bsp-init/SKILL.md) | ⭐ 95 | `daily assistant` |
| [Skill](daily-assistant/name-skill_2b9e204b/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/meta-remove/SKILL.md) | ⭐ 95 | `daily assistant` |

### Data Analysis (13 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](data-analysis/name-skill_3159036d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/android-jetpack-compose-expert/SKILL.md) | 🔥 15.4k | `data analysis` |
| [Skill](data-analysis/name-skill_a17f8188/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-controls/SKILL.md) | 🔥 15.4k | `data analysis` |
| [Skill](data-analysis/name-skill_7edb395a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mermaid-expert/SKILL.md) | 🔥 15.4k | `data analysis` |
| [Skill](data-analysis/name-skill_a5a64d2d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-content-auditor/SKILL.md) | 🔥 15.4k | `data analysis` |
| [Advanced Attacks](data-analysis/advanced-attacks_56418a81/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/web-app/public/skills/active-directory-attacks/references/advanced-attacks.md) | 🔥 15.4k | `data analysis` |
| [Skill](data-analysis/name-skill_7c8eb33b/) | [tripleyak/SkillForge](https://raw.githubusercontent.com/tripleyak/SkillForge/main/SKILL.md) | ⭐ 516 | `data analysis` |
| [Regression Questions](data-analysis/regression-questions_b934bd95/) | [tripleyak/SkillForge](https://raw.githubusercontent.com/tripleyak/SkillForge/main/references/regression-questions.md) | ⭐ 516 | `data analysis` |
| [Script Patterns Catalog](data-analysis/script-patterns-catalog_ac036d19/) | [tripleyak/SkillForge](https://raw.githubusercontent.com/tripleyak/SkillForge/main/references/script-patterns-catalog.md) | ⭐ 516 | `data analysis` |
| [1C Configuration Spec](data-analysis/1c-configuration-spec_9983109b/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/1c-configuration-spec.md) | ⭐ 95 | `data analysis` |
| [1C Specs Index](data-analysis/1c-specs-index_b7b8a8d3/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/1c-specs-index.md) | ⭐ 95 | `data analysis` |
| [1C Subsystem Spec](data-analysis/1c-subsystem-spec_985a6bc0/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/1c-subsystem-spec.md) | ⭐ 95 | `data analysis` |
| [Mxl Guide](data-analysis/mxl-guide_ec86d6cd/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/mxl-guide.md) | ⭐ 95 | `data analysis` |
| [Skd Dsl Spec](data-analysis/skd-dsl-spec_b9e1becb/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/skd-dsl-spec.md) | ⭐ 95 | `data analysis` |

### Development (46 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Catalog](development/catalog_25603c93/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/CATALOG.md) | 🔥 15.4k | `development` |
| [Bundles](development/bundles_e1d442f9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/docs/BUNDLES.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_5e1679cc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/api-documenter/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_fdbc231b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/aws-cost-cleanup/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_f7d17c1f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-projects-dotnet/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_ab353286/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-vision-imageanalysis-py/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_c6c48c4c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-rust/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_66f7496c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-identity-dotnet/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_1f0230f0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-identity-rust/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_19d6d394/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-keyvault-keys-rust/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_eb72f545/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-keyvault-secrets-rust/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_e17fed7f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-opentelemetry-exporter-java/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_99290b4a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-opentelemetry-py/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_de2e670d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-cosmosdb-dotnet/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_6de00132/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-blob-rust/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_cd5f8cda/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/backend-security-coder/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_6643a8eb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bash-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_141a2291/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cdk-patterns/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_2b6ed6ca/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cloudflare-workers-expert/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_75d666fe/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/csharp-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_0f2a6378/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/docs-architect/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_0739037a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/flutter-expert/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_cbb16e54/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/frontend-developer/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_f3b482c0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/godot-4-migration/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_a5f29da5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-foundations/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_88018179/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ios-developer/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_c3a5efb1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/javascript-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_3e7883fd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/kotlin-coroutines-expert/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_4b2dbc05/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/minecraft-bukkit-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_61a74054/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mobile-security-coder/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_4689f81b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/nerdzao-elite-gemini-high/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_047ac83f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/python-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_85965eaf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/rust-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_ab4d5702/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/scala-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_1d9f066f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security-scanning-security-sast/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_edf52845/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/typescript-expert/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_4fabba11/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ui-ux-designer/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_34209777/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/unity-developer/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_589fcb1c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wordpress-plugin-development/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_6b6f277a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wordpress-theme-development/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/name-skill_1a5c1860/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wordpress/SKILL.md) | 🔥 15.4k | `development` |
| [Advanced Aws Pentesting](development/advanced-aws-pentesting_11d414ad/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/web-app/public/skills/aws-penetration-testing/references/advanced-aws-pentesting.md) | 🔥 15.4k | `development` |
| [Advanced Cloud Scripts](development/advanced-cloud-scripts_7eca868d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/web-app/public/skills/cloud-penetration-testing/references/advanced-cloud-scripts.md) | 🔥 15.4k | `development` |
| [Synthesis Protocol](development/synthesis-protocol_6141a6e8/) | [tripleyak/SkillForge](https://raw.githubusercontent.com/tripleyak/SkillForge/main/references/synthesis-protocol.md) | ⭐ 516 | `development` |
| [Db Guide](development/db-guide_56e641aa/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/db-guide.md) | ⭐ 95 | `development` |
| [Meta Guide](development/meta-guide_f098eb62/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/meta-guide.md) | ⭐ 95 | `development` |

### Development/Devops (24 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Workflow  Bundlesreadme](development/devops/workflow-bundlesreadme_74af78cc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/workflow-%20bundlesREADME.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_045c4e45/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-apicenter-py/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_1b6a88ff/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cloud-architect/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_4ccc74ea/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-driven-development/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_bbc10b0f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/database-admin/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_bb70a7d1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/deployment-engineer/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_76d0f17c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/devops-troubleshooter/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_a461824a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/dotnet-architect/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_34ffe725/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/fastapi-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_24e6fa8c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/golang-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_1267ddd9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/grpc-golang/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_0adc0eb8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/kubernetes-deployment/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_7b6a4559/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/m365-agents-ts/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_16fc635e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mlops-engineer/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_7278c29a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/network-engineer/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_d4daeb69/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/temporal-python-pro/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/devops/name-skill_975bca83/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/test-automator/SKILL.md) | 🔥 15.4k | `development` |
| [1C Epf Spec](development/devops/1c-epf-spec_5adac7bb/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/1c-epf-spec.md) | ⭐ 95 | `development` |
| [Skill](development/devops/name-skill_3264772e/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/db-dump-cf/SKILL.md) | ⭐ 95 | `development` |
| [Skill](development/devops/name-skill_ee4b6b18/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/db-list/SKILL.md) | ⭐ 95 | `development` |
| [Skill](development/devops/name-skill_fcd643bf/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/db-load-git/SKILL.md) | ⭐ 95 | `development` |
| [Skill](development/devops/name-skill_e92700cd/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/db-run/SKILL.md) | ⭐ 95 | `development` |
| [Skill](development/devops/name-skill_410555d2/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/erf-build/SKILL.md) | ⭐ 95 | `development` |
| [Skill](development/devops/name-skill_698119fe/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/erf-dump/SKILL.md) | ⭐ 95 | `development` |

### Development/Testing (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/testing/name-skill_0af17ce7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-orchestrator/SKILL.md) | 🔥 15.4k | `development` |

### Development/Tools (13 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/tools/name-skill_c8db859c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/appdeploy/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/tools/name-skill_50ec7b20/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-agents-persistent-java/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/tools/name-skill_e84e5829/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-voicelive-java/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/tools/name-skill_6f3ad8e1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-layout/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/tools/name-skill_319a5cdc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-search/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/tools/name-skill_4ce1b2d9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-inputs/SKILL.md) | 🔥 15.4k | `development` |
| [Skill](development/tools/name-skill_f965679f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-nextjs-development/SKILL.md) | 🔥 15.4k | `development` |
| [Implementation Playbook](development/tools/implementation-playbook_7db0d79c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/grpc-golang/resources/implementation-playbook.md) | 🔥 15.4k | `development` |
| [Skill](development/tools/name-skill_15d20563/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security/aws-security-audit/SKILL.md) | 🔥 15.4k | `aws` `security` `audit` |
| [1C Help Spec](development/tools/1c-help-spec_6cebe8d3/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/1c-help-spec.md) | ⭐ 95 | `development` |
| [Skill](development/tools/name-skill_24950a8f/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/epf-dump/SKILL.md) | ⭐ 95 | `development` |
| [Skill](development/tools/name-skill_49bbefb7/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/form-edit/SKILL.md) | ⭐ 95 | `development` |
| [Skill](development/tools/name-skill_4ba6b448/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/meta-compile/SKILL.md) | ⭐ 95 | `development` |

### Investment (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](investment/name-skill_3e0f379b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-scientist/SKILL.md) | 🔥 15.4k | `investment` |
| [Skill](investment/name-skill_0db988ba/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/SKILL.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/decision-frameworks_95ad10a0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/carrier-relationship-management/references/decision-frameworks.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/decision-frameworks_30ceb8b9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/references/decision-frameworks.md) | 🔥 15.4k | `investment` |
| [Skill](investment/name-skill_4f5c9321/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/quant-analyst/SKILL.md) | 🔥 15.4k | `investment` |
| [Skill](investment/name-skill_605596b3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-financial-modeling/SKILL.md) | 🔥 15.4k | `investment` |
| [Communication Templates](investment/communication-templates_b7c10e59/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/references/communication-templates.md) | 🔥 15.4k | `investment` |

### Productivity (8 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](productivity/name-skill_541ef5b9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/arm-cortex-expert/SKILL.md) | 🔥 15.4k | `productivity` |
| [Skill](productivity/name-skill_cfa4e81d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-mysql-dotnet/SKILL.md) | 🔥 15.4k | `productivity` |
| [Skill](productivity/name-skill_c527e215/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-postgresql-dotnet/SKILL.md) | 🔥 15.4k | `productivity` |
| [Skill](productivity/name-skill_c5bbef87/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bevy-ecs-expert/SKILL.md) | 🔥 15.4k | `productivity` |
| [Skill](productivity/name-skill_a88e6107/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-system/SKILL.md) | 🔥 15.4k | `productivity` |
| [Skill](productivity/name-skill_2eafb96c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-platforms/SKILL.md) | 🔥 15.4k | `productivity` |
| [Skill](productivity/name-skill_ea52a75f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/php-pro/SKILL.md) | 🔥 15.4k | `productivity` |
| [Skill](productivity/name-skill_1df978d7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/web-app/public/skills/bevy-ecs-expert/SKILL.md) | 🔥 15.4k | `productivity` |

### Research (6 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](research/139-name-skill_40343acb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-business-analyst-market-opportunity/SKILL.md) | 🔥 15.4k | `research` |
| [Build Spec](research/265-build-spec_f3072f59/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/build-spec.md) | ⭐ 95 | `research` |
| [Skill](research/name-skill_b6e33b5c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-analyst/SKILL.md) | 🔥 15.4k | `research` |
| [Skill](research/name-skill_ed811839/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/epf-validate/SKILL.md) | ⭐ 95 | `research` |
| [Skill](research/name-skill_5f9c53d8/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/erf-validate/SKILL.md) | ⭐ 95 | `research` |
| [Skill](research/name-skill_8063f048/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/form-validate/SKILL.md) | ⭐ 95 | `research` |

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

*Last updated: 2026-02-26 18:44:15 UTC*
*Automatically maintained by SkillFlow*
