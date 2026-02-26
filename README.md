# X-Skills

A curated collection of **40 AI-powered skills** organized into 9 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (8 skills)
- **Communication** (3 skills)
- **Content Creation** (7 skills)
- **Daily Assistant** (11 skills)
- **Development** (1 skill)
- **Development/Testing** (2 skills)
- **Development/Tools** (1 skill)
- **Investment** (4 skills)
- **Research** (3 skills)

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


### Automation/Scripting (8 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/scripting/003-name-skill_c94cab98/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-speech-to-text-rest-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/003-name-skill_18050c5d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-queue-ts/SKILL.md) | 🔥 15.4k | `automation` |
| [Web Spec](automation/scripting/086-web-spec_e7518f36/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/web-spec.md) | ⭐ 95 | `automation` |
| [Reference](automation/scripting/087-reference_704f68d2/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/cf-edit/reference.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/003-name-skill_e622e2f6/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/mxl-info/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/003-name-skill_e110d174/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/skd-info/SKILL.md) | ⭐ 95 | `automation` |
| [Interview Openspec](automation/scripting/interview-openspec_ccdcaf6f/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/docs/skills/interview-openspec.md) | ⭐ 10 | `automation` |
| [Interview Principles](automation/scripting/interview_principles_544f0960/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/resources/INTERVIEW_PRINCIPLES.md) | ⭐ 10 | `automation` |

### Communication (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](communication/127-name-skill_2b647276/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-java/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/127-name-skill_45420053/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-arizeaiobservabilityeval-dotnet/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/127-name-skill_373d16f0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-botservice-py/SKILL.md) | 🔥 15.4k | `communication` |

### Content Creation (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Build Program](content-creation/366-build_program_6794f052/) | [julianghadially/CodeEvolver](https://raw.githubusercontent.com/julianghadially/CodeEvolver/main/specs/analysis/build_program.md) | ⭐ 12 | `content creation` |
| [Claude](content-creation/claude_bd86cf0c/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/CLAUDE.md) | ⭐ 10 | `content creation` |
| [Readme Cn](content-creation/readme_cn_e8c6ee87/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/README_CN.md) | ⭐ 10 | `content creation` |
| [Index](content-creation/index_4fd68aa3/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/docs/skills/index.md) | ⭐ 10 | `content creation` |
| [Index](content-creation/index_28dd84e3/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/docs/zh/skills/index.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/name-skill_f2f65ed3/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/SKILL.md) | ⭐ 10 | `content creation` |
| [Examples](content-creation/examples_35d2e9f5/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/resources/EXAMPLES.md) | ⭐ 10 | `content creation` |

### Daily Assistant (11 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Edge Cases](daily-assistant/288-edge-cases_50a38ae5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/inventory-demand-planning/references/edge-cases.md) | 🔥 15.4k | `daily assistant` |
| [Decision Frameworks](daily-assistant/289-decision-frameworks_90e5dd3e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/production-scheduling/references/decision-frameworks.md) | 🔥 15.4k | `daily assistant` |
| [Decision Frameworks](daily-assistant/289-decision-frameworks_e50a63c7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/quality-nonconformance/references/decision-frameworks.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_76644efe/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/legal-advisor/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_01c2e76e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/risk-manager/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_25fc0d0a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-authority-builder/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_ac125d57/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-meta-optimizer/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_27ccbbd6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/track-management/SKILL.md) | 🔥 15.4k | `daily assistant` |
| [Decision Frameworks](daily-assistant/288-decision-frameworks_750ee441/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/inventory-demand-planning/references/decision-frameworks.md) | 🔥 15.4k | `daily assistant` |
| [Interview Openspec](daily-assistant/interview-openspec_da8dc2d0/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/docs/zh/skills/interview-openspec.md) | ⭐ 10 | `daily assistant` |
| [Artifact Templates](daily-assistant/artifact_templates_6c8c9b85/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/resources/ARTIFACT_TEMPLATES.md) | ⭐ 10 | `daily assistant` |

### Development (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Openspec Interview Dimensions](development/openspec_interview_dimensions_693c5b0e/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/resources/OPENSPEC_INTERVIEW_DIMENSIONS.md) | ⭐ 10 | `development` |

### Development/Testing (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/testing/002-name-skill_0af17ce7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-orchestrator/SKILL.md) | 🔥 15.4k | `development` |
| [Risk Protection](development/testing/risk_protection_ff4e951b/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/resources/RISK_PROTECTION.md) | ⭐ 10 | `development` |

### Development/Tools (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agents](development/tools/agents_15ebd473/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/AGENTS.md) | ⭐ 10 | `development` |

### Investment (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](investment/021-name-skill_3e0f379b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-scientist/SKILL.md) | 🔥 15.4k | `investment` |
| [Skill](investment/021-name-skill_0db988ba/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/SKILL.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/049-decision-frameworks_95ad10a0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/carrier-relationship-management/references/decision-frameworks.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/049-decision-frameworks_30ceb8b9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/references/decision-frameworks.md) | 🔥 15.4k | `investment` |

### Research (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](research/139-name-skill_40343acb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-business-analyst-market-opportunity/SKILL.md) | 🔥 15.4k | `research` |
| [Build Spec](research/265-build-spec_f3072f59/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/build-spec.md) | ⭐ 95 | `research` |
| [Skill](research/139-name-skill_b6e33b5c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-analyst/SKILL.md) | 🔥 15.4k | `research` |

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

*Last updated: 2026-02-26 20:48:02 UTC*
*Automatically maintained by SkillFlow*
