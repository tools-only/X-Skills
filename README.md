# X-Skills

A curated collection of **33 AI-powered skills** organized into 10 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Workflow** (6 skills)
- **Commercial** (5 skills)
- **Communication** (1 skill)
- **Content Creation** (5 skills)
- **Daily Assistant** (2 skills)
- **Data Analysis** (1 skill)
- **Development** (5 skills)
- **Development/Devops** (1 skill)
- **Development/Tools** (2 skills)
- **Research** (5 skills)

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


### Automation/Workflow (6 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/name-skill_ce4bb9a9/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-deploy/SKILL.md) | ⭐ 65 | `automation` |
| [Cli Guide](automation/workflow/cli-guide_49550e7a/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/docs/cli-guide.md) | ⭐ 65 | `automation` |
| [Skill](automation/workflow/name-skill_09aea88b/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/human-checkpoint/SKILL.md) | ⭐ 197 | `automation` |
| [Skill](automation/workflow/name-skill_4052fdf4/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/pipeline-router/SKILL.md) | ⭐ 197 | `automation` |
| [Skill](automation/workflow/name-skill_09555a86/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/screening-manager/SKILL.md) | ⭐ 197 | `automation` |
| [Test Spec Reference](automation/workflow/test-spec-reference_86af974b/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentforce-testing/resources/test-spec-reference.md) | ⭐ 65 | `automation` |

### Commercial (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/name-skill_10221f10/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/manuscript-ingest/SKILL.md) | ⭐ 197 | `commercial` |
| [Security Policy](commercial/371-security_policy_f8dea503/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/SECURITY_POLICY.md) | ⭐ 18 | `commercial` |
| [Customer Success](commercial/372-customer_success_85216f92/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/customer_success.md) | ⭐ 18 | `commercial` |
| [Executive](commercial/373-executive_e9668604/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/executive.md) | ⭐ 18 | `commercial` |
| [Finance](commercial/374-finance_8fdccf44/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/finance.md) | ⭐ 18 | `commercial` |

### Communication (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Syntax Reference](communication/syntax-reference_a16eba1d/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/syntax-reference.md) | ⭐ 65 | `communication` |

### Content Creation (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill Index](content-creation/skill_index_09523f5a/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/SKILL_INDEX.md) | ⭐ 197 | `content creation` |
| [Idea Finder.Pipeline](content-creation/idea-finderpipeline_1e97a927/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/pipelines/idea-finder.pipeline.md) | ⭐ 197 | `content creation` |
| [Lit Snapshot.Pipeline](content-creation/lit-snapshotpipeline_cfa446bb/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/pipelines/lit-snapshot.pipeline.md) | ⭐ 197 | `content creation` |
| [Tutorial.Pipeline](content-creation/tutorialpipeline_a1a30b37/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/pipelines/tutorial.pipeline.md) | ⭐ 197 | `content creation` |
| [Skill](content-creation/name-skill_e03f8008/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/idea-shortlist-curator/SKILL.md) | ⭐ 197 | `content creation` |

### Daily Assistant (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Peer Review.Pipeline](daily-assistant/peer-reviewpipeline_47e029c1/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/pipelines/peer-review.pipeline.md) | ⭐ 197 | `daily assistant` |
| [Skill](daily-assistant/name-skill_a81b8190/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/deliverable-selfloop/SKILL.md) | ⭐ 197 | `daily assistant` |

### Data Analysis (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](data-analysis/name-skill_d0682987/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/idea-pool-expander/SKILL.md) | ⭐ 197 | `data analysis` |

### Development (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Credits](development/credits_7d3b11da/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/CREDITS.md) | ⭐ 65 | `development` |
| [Actions Reference](development/actions-reference_733b0bdd/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/actions-reference.md) | ⭐ 65 | `development` |
| [Migration Guide](development/migration-guide_6f7f2a5a/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/migration-guide.md) | ⭐ 65 | `development` |
| [Skill](development/name-skill_590ac30d/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentforce-testing/SKILL.md) | ⭐ 65 | `development` |
| [Test Spec Guide](development/test-spec-guide_17704202/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentforce-testing/docs/test-spec-guide.md) | ⭐ 65 | `development` |

### Development/Devops (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Known Issues](development/devops/known-issues_97ef6daa/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/known-issues.md) | ⭐ 65 | `development` |

### Development/Tools (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/tools/name-skill_727accdb/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/SKILL.md) | ⭐ 65 | `development` |
| [Official Sources](development/tools/official-sources_aab46d33/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/official-sources.md) | ⭐ 65 | `development` |

### Research (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skills Standard](research/skills_standard_dea82938/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/SKILLS_STANDARD.md) | ⭐ 197 | `research` |
| [Systematic Review.Pipeline](research/systematic-reviewpipeline_28159c96/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/pipelines/systematic-review.pipeline.md) | ⭐ 197 | `research` |
| [Skill](research/name-skill_c6fee295/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/idea-brief/SKILL.md) | ⭐ 197 | `research` |
| [Skill](research/name-skill_d7b34b7c/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/protocol-writer/SKILL.md) | ⭐ 197 | `research` |
| [Skill](research/name-skill_70d2c27b/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/research-pipeline-runner/SKILL.md) | ⭐ 197 | `research` |

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

*Last updated: 2026-02-12 20:23:46 UTC*
*Automatically maintained by SkillFlow*
