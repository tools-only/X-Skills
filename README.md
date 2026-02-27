# X-Skills

A curated collection of **924 AI-powered skills** organized into 15 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (61 skills)
- **Automation/Workflow** (176 skills)
- **Commercial** (54 skills)
- **Communication** (28 skills)
- **Content Creation** (76 skills)
- **Daily Assistant** (50 skills)
- **Data Analysis** (47 skills)
- **Development** (224 skills)
- **Development/Devops** (73 skills)
- **Development/Testing** (15 skills)
- **Development/Tools** (51 skills)
- **Investment** (19 skills)
- **Other** (2 skills)
- **Productivity** (21 skills)
- **Research** (27 skills)

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


### Automation/Scripting (61 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/scripting/003-name-skill_c94cab98/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-speech-to-text-rest-py/SKILL.md) | 🔥 15.4k | `automation` |
| [Skill](automation/scripting/003-name-skill_18050c5d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-queue-ts/SKILL.md) | 🔥 15.4k | `automation` |
| [Web Spec](automation/scripting/086-web-spec_e7518f36/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/web-spec.md) | ⭐ 95 | `automation` |
| [Reference](automation/scripting/087-reference_704f68d2/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/cf-edit/reference.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/003-name-skill_e622e2f6/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/mxl-info/SKILL.md) | ⭐ 95 | `automation` |
| [Skill](automation/scripting/003-name-skill_e110d174/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/.claude/skills/skd-info/SKILL.md) | ⭐ 95 | `automation` |
| [Interview Openspec](automation/scripting/090-interview-openspec_ccdcaf6f/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/docs/skills/interview-openspec.md) | ⭐ 10 | `automation` |
| [Phase 2 Compiler Analysis](automation/scripting/090-phase-2-compiler-analysis_9ac1b095/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/workflows/phase-2-compiler-analysis.md) | ⭐ 2.9k | `automation` |
| [Slurm Directives](automation/scripting/090-slurm_directives_ef58f60f/) | [HeshamFS/materials-simulation-skills](https://raw.githubusercontent.com/HeshamFS/materials-simulation-skills/main/skills/hpc-deployment/slurm-job-script-generator/references/slurm_directives.md) | ⭐ 20 | `automation` |
| [Skill](automation/scripting/003-name-skill_67621ec1/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/skills/system/get-skill-info/SKILL.md) | ⭐ 763 | `automation` |
| [Backlog](automation/scripting/088-backlog_193d131d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/BACKLOG.md) | ⭐ 20 | `automation` |
| [P0 Reduce Session Start Context Load Via Rules Path Scoping And](automation/scripting/089-p0-reduce-session-start-context-load-via-rules-path-scoping-and_a0cd61f8/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p0-reduce-session-start-context-load-via-rules-path-scoping-and.md) | ⭐ 20 | `automation` |
| [P1 P1 Plugin Validator Pre Commit Output Is Too Noisy](automation/scripting/090-p1-p1-plugin-validator-pre-commit-output-is-too-noisy_bd37edca/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-p1-plugin-validator-pre-commit-output-is-too-noisy.md) | ⭐ 20 | `automation` |
| [P1 Sam Human Escalation Criteria](automation/scripting/091-p1-sam-human-escalation-criteria_22c63258/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p1-sam-human-escalation-criteria.md) | ⭐ 20 | `automation` |
| [P2 Conventional Commits Fix Skillmd Frontmatter Broken Links An](automation/scripting/092-p2-conventional-commits-fix-skillmd-frontmatter-broken-links-an_a2fe1a4b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-conventional-commits-fix-skillmd-frontmatter-broken-links-an.md) | ⭐ 20 | `automation` |
| [P2 Plan Artifact Diverges From Implementation Without Update Me](automation/scripting/093-p2-plan-artifact-diverges-from-implementation-without-update-me_e309f768/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/backlog/p2-plan-artifact-diverges-from-implementation-without-update-me.md) | ⭐ 20 | `automation` |
| [Session](automation/scripting/091-session_17d91954/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/commands/session.md) | ⭐ 543 | `automation` |
| [Json Output](automation/scripting/094-json-output_9f8d7838/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/json-output.md) | ⭐ 643 | `automation` |
| [Orc](automation/scripting/095-orc_a5d95af8/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/agent/orc.md) | ⭐ 643 | `automation` |
| [Skill](automation/scripting/003-name-skill_582450bc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ai-wrapper-product/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_d3aaf2e2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/algolia-search/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_e2e66acf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/aws-serverless/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_d1364535/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-anomalydetector-java/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_de6e1d28/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-contentsafety-java/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_b754d714/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-translation-document-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_8a1b2928/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-translation-ts/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_4228a2a9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-communication-sms-java/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_bc9140de/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-compute-batch-java/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_2a7ee078/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-data-tables-java/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_95396749/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-data-tables-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_54acb0c7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-eventgrid-dotnet/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_4b93371d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-eventhub-dotnet/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_0c02ed93/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-eventhub-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_f58d80fb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-ingestion-java/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_c3649151/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-opentelemetry-ts/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_63519d8c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-query-java/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_145e2d3c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-query-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_dabb3d89/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-search-documents-dotnet/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_66e50d19/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-search-documents-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_41b30733/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-search-documents-ts/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_98d0d783/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-servicebus-dotnet/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_9a081cad/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-speech-to-text-rest-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_944ad65a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-queue-ts/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_9880f848/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/claude-win11-speckit-update-skill/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_a663a5b2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cost-optimization/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_5dde8ed0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/event-store-design/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_3ee0dbb1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/executing-plans/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_572c001e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/finishing-a-development-branch/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_14356b40/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/firecrawl-scraper/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_30d87ea2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/game-development/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_5bdee658/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hubspot-integration/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_14ec4b9a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/microservices-patterns/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_e37403ed/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/multi-agent-patterns/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_bd9843b3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/nft-standards/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_b711c3da/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/segment-cdp/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_2aec8798/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/skill-seekers/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_f60cb467/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/varlock-claude-skill/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_47931fb1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/game-development/pc-games/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_43af504e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/game-development/web-games/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/scripting/003-name-skill_a471a895/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/libreoffice/calc/SKILL.md) | 🔥 13.9k | `automation` |
| [Uniswap Driver](automation/scripting/uniswap-driver_6ed6ed96/) | [Uniswap/uniswap-ai](https://raw.githubusercontent.com/Uniswap/uniswap-ai/main/docs/plugins/uniswap-driver.md) | ⭐ 130 | `automation` |

### Automation/Workflow (176 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Cross File Resolution](automation/workflow/144-cross-file-resolution_c15a46bb/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/cross-file-resolution.md) | ⭐ 2.9k | `automation` |
| [Macos Arm64E Workaround](automation/workflow/145-macos-arm64e-workaround_2b4e18ae/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/codeql/references/macos-arm64e-workaround.md) | ⭐ 2.9k | `automation` |
| [Task](automation/workflow/146-task_dd739ef5/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/prompts/task.md) | ⭐ 2.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_247714ad/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/slides-generator/SKILL.md) | ⭐ 11 | `automation` |
| [Lp Executor Guide](automation/workflow/139-lp_executor_guide_0040aeee/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/lp-agent/references/lp_executor_guide.md) | ⭐ 11 | `automation` |
| [Skill](automation/workflow/002-name-skill_4c881e3b/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/skills/system/create-plan/SKILL.md) | ⭐ 763 | `automation` |
| [Skill](automation/workflow/002-name-skill_0e789cbb/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/.github/skills/classify-temp-entries-to-section/SKILL.md) | ⭐ 392 | `automation` |
| [Skill](automation/workflow/133-andruia-skill_d3638bbf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/00-andruia-consultant/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_2557f868/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/3d-web-experience/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_55e15cb6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/active-directory-attacks/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_48ddc190/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/activecampaign-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_fb1eb36b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/agent-evaluation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_d076502c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/agent-memory-systems/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_aa0748a6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/agent-orchestration-multi-agent-optimize/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_ef32d8c4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ai-ml/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_18237c41/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/airtable-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_78ae17e1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/amplitude-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b8b568eb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/antigravity-workflows/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_4fd06e45/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/architecture/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_2cf3896a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/asana-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_27fc8d4a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/audio-transcriber/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_738e02a0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/avalonia-viewmodels-zafiro/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_8b9844a6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-ml-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_302501b2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-appconfiguration-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_ff169171/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-eventhub-ts/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_8d30b32d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-keyvault-keys-ts/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_a5112cd9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-keyvault-secrets-ts/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_991192e7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-messaging-webpubsubservice-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_d4235c61/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-fabric-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_6e8cf33a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-file-share-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_cdc1dea6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-queue-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_54aed437/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bamboohr-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_43faa950/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/basecamp-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_f5133ce2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bitbucket-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_cf403a77/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/box-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_5d727724/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/brevo-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_1d4f6542/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bullmq-specialist/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_dfabaa47/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/business-analyst/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_56ac0b37/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cal-com-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_76f3ad22/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/calendly-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_59c84d04/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/canva-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_ce13ee05/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cc-skill-backend-patterns/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_c1dc9577/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/changelog-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_fef8088f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/circleci-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_47cd1e7d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/clerk-auth/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_f64f42d6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/clickup-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_57a74bf1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/close-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_bee54def/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/coda-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_4d498659/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-review-ai-ai-review/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_1b521bab/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/codebase-cleanup-deps-audit/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_70870fa1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/comprehensive-review-full-review/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b584225c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conductor-manage/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_43d9c078/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/confluence-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_59ff14bb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-window-management/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_1a1b4543/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conversation-memory/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_c900e2c7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/convertkit-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_115ff7d3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/crypto-bd-agent/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_8a35613a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-engineering-data-driven-feature/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_56ad5bcf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/database-architect/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_a6f87377/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/database-optimizer/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_892a6eb7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/datadog-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_9d9f5a77/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/dependency-management-deps-audit/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_20dc62ec/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/deployment-validation-config-validate/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_7e17ca21/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/design-orchestration/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_58225486/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/discord-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_f2fe94a3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/docusign-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_cd999390/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/email-sequence/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_42b188f9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/figma-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_89391b5e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/file-uploads/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b08b3872/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/framework-migration-code-migrate/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_4343a237/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/framework-migration-deps-upgrade/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_be131533/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/freshdesk-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_fbc27572/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/freshservice-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_143f8016/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/git-pr-workflows-pr-enhance/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_43910758/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/git-pushing/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_8bf692bc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/github-workflow-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_ccb69342/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/gitlab-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_2368aa70/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/gitops-workflow/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_0b766a91/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/gmail-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_01cb3e25/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/google-analytics-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_109401cc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/google-calendar-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_e16a552c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/google-drive-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_bb980b65/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/googlesheets-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_e5862b60/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/helpdesk-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_112c2371/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hubspot-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_cc22c364/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hugging-face-cli/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_8c4b707a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hugging-face-jobs/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_1c18cf8b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hybrid-cloud-architect/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_810d5e47/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/incident-responder/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_7519a71b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/incident-response-smart-fix/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_ee48d35a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/inngest/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_0e5bce9f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/intercom-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_3e3c4ff0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/javascript-typescript-typescript-scaffold/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_a076bc6c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/jira-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b967839e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/k8s-security-policies/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_f4c84259/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/klaviyo-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_874a319c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/linear-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_a4473581/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/linkedin-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_31ffba4e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/linkedin-cli/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b56b195b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/linux-privilege-escalation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_09dd8a18/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/linux-shell-scripting/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b8427991/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/logistics-exception-management/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_d7e9e388/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/m365-agents-dotnet/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_c49e0521/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/m365-agents-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_ebc81458/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mailchimp-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_1b5004d3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/make-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_6f2b5704/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/miro-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_9f1c8a7f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mixpanel-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_4809eb75/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ml-engineer/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_be877626/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ml-pipeline-workflow/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_d5ef3d12/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/monday-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b46d7891/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/monorepo-architect/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_66e77a9f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/neon-postgres/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_7c1f3a43/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/nextjs-supabase-auth/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_939b5249/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/notion-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_0006a906/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/observe-whatsapp/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_787f2107/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/office-productivity/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_7b4d33bb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/one-drive-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_962bea24/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/pagerduty-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_80e5057d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/performance-profiling/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_c973b620/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/pipedrive-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_703b7f74/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/plaid-fintech/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_89325155/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/posthog-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_3d51350b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/postmark-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_0a93fea0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/prompt-caching/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_dd9c1629/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/prompt-engineering-patterns/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_01b4f66c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/pydantic-models-py/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_a8abb260/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/rag-engineer/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_d2a84854/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-best-practices/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_e2b26a4b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/red-team-tools/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_985965f3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/reddit-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_fd1b96c7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/render-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_a443ad11/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/salesforce-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_1bb3bf02/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/salesforce-development/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_65639e01/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/scanning-tools/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_8d74ce70/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security-audit/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b7fc34b3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security-compliance-compliance-check/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_891e843c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security-scanning-security-hardening/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_720fd150/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/segment-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_39ae8340/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/sendgrid-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_cd509d43/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/senior-architect/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_8c9ca936/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/senior-fullstack/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_55a656af/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/sentry-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_c6f99a7f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/service-mesh-expert/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_db7a625b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/skill-creator/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_15111732/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/slack-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_f3bf2786/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/square-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_13492db4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/supabase-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_e3586171/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/telegram-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_9c601d0f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/terraform-specialist/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_40998f5f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/theme-factory/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_3cd6580b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/todoist-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_7f6abd65/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/trello-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_6821b12b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/trigger-dev/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_6693dbe9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/twilio-communications/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_4a9d20d9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/twitter-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_da95b7dc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ui-visual-validator/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_f9b50702/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/vercel-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_46052def/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/voice-agents/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_2c0bbd98/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/vulnerability-scanner/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_bafc0132/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wiki-vitepress/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_43600599/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wrike-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_73526343/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/youtube-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_937e4b07/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/zendesk-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_fde66cf3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/zoho-crm-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_547e8663/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/zoom-automation/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_893cdd0d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/game-development/game-design/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_73214b76/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/game-development/mobile-games/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_e0b29456/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/game-development/multiplayer/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b7ff54ed/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/libreoffice/draw/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_a28a8196/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/libreoffice/impress/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_b7bfad60/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/libreoffice/writer/SKILL.md) | 🔥 13.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_58a9e1a9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security/aws-iam-best-practices/SKILL.md) | 🔥 13.9k | `automation` |
| [Asset Bundles](automation/workflow/147-asset-bundles_e0cde027/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks/asset-bundles.md) | ⭐ 13 | `automation` |
| [Data Exploration](automation/workflow/148-data-exploration_68b4b3e4/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks/data-exploration.md) | ⭐ 13 | `automation` |
| [Databricks Cli Auth](automation/workflow/149-databricks-cli-auth_d30d5288/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks/databricks-cli-auth.md) | ⭐ 13 | `automation` |

### Commercial (54 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/210-name-skill_5fcd131e/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/SKILL.md) | ⭐ 2.9k | `commercial` |
| [Ir Analysis](commercial/372-ir-analysis_473b57e4/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/references/ir-analysis.md) | ⭐ 2.9k | `commercial` |
| [17 The Star](commercial/373-17-the-star_0cc63419/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/17-the-star.md) | ⭐ 2.9k | `commercial` |
| [King Of Pentacles](commercial/374-king-of-pentacles_ee2168eb/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/pentacles/king-of-pentacles.md) | ⭐ 2.9k | `commercial` |
| [Claude](commercial/036-claude_10641fac/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/CLAUDE.md) | ⭐ 543 | `commercial` |
| [Skill](commercial/210-name-skill_e1fcc7c4/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/shopify/skills/shopify-setup/SKILL.md) | ⭐ 543 | `commercial` |
| [Skill](commercial/210-name-skill_1b1d59b7/) | [timescale/pg-aiguide](https://raw.githubusercontent.com/timescale/pg-aiguide/main/skills/design-postgis-tables/SKILL.md) | ⭐ 1.6k | `commercial` |
| [Skill](commercial/368-andruia-skill_d3ca8342/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/20-andruia-niche-intelligence/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_508ed84f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/angular-best-practices/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_d3cf4f16/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/angular-state-management/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_73c04116/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-document-intelligence-dotnet/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_f7b83489/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-document-intelligence-ts/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_e98c1ae6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-formrecognizer-java/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_cc601eb7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-appconfiguration-java/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_aaa8a4a1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-appconfiguration-ts/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_4eef4cf8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-ts/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_d108988d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-eventhub-java/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_52b07764/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-apimanagement-py/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_fc1b020c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-postgres-ts/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_8ff58956/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-sql-dotnet/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_9cc24f7c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-security-keyvault-keys-dotnet/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_6bae6591/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-security-keyvault-keys-java/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_e1b393a5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-security-keyvault-secrets-java/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_4c9ca815/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-optimization/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_f2086a9f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/customs-trade-compliance/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_afae905b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ddd-context-mapping/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_e041ea43/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/domain-driven-design/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_801ff84f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/event-sourcing-architect/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_0fc4180b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/expo-deployment/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_a432f5ee/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/laravel-security-audit/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_63d7e96c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/market-sizing-analysis/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_b4dc8ff7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/memory-systems/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_fe778798/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/network-101/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_3c65ffc4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/payment-integration/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_c12d3b0c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/paywall-upgrade-cro/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_327d0197/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/pci-compliance/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_4908afc4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/postgresql/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_1fc2168b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-flow-architect/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_8a453392/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-patterns/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_7223c024/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/red-team-tactics/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_078360e9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/returns-reverse-logistics/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_73562642/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/saga-orchestration/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_3d52c44e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/shopify-apps/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_8f7e36a5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/shopify-automation/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_e1741ea7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/shopify-development/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_e0879773/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/signup-flow-cro/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_f7a84577/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-business-analyst-financial-projections/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_2b8514ab/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/stitch-ui-design/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_ac4b6a48/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/stripe-automation/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_970b0c5d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/stripe-integration/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_bbc1d0c4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/telegram-bot-builder/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_dd15e392/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/terraform-aws-modules/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_c25c4f67/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/workflow-orchestration-patterns/SKILL.md) | 🔥 13.9k | `commercial` |
| [Skill](commercial/210-name-skill_9ff101ae/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/game-development/game-art/SKILL.md) | 🔥 13.9k | `commercial` |

### Communication (28 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](communication/127-name-skill_2b647276/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-java/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/127-name-skill_45420053/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-arizeaiobservabilityeval-dotnet/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/127-name-skill_373d16f0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-botservice-py/SKILL.md) | 🔥 15.4k | `communication` |
| [Foundations](communication/252-foundations_9e0eb90b/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/foundations.md) | ⭐ 2.9k | `communication` |
| [Authentication](communication/016-authentication_432b0de4/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/user-guide/authentication.md) | ⭐ 177 | `communication` |
| [Readme Cn](communication/252-readme_cn_f4fb8d9b/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/README_CN.md) | ⭐ 763 | `communication` |
| [Skill](communication/127-name-skill_bb401c33/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/skills/system/get-chat-history/SKILL.md) | ⭐ 763 | `communication` |
| [Skill](communication/127-name-skill_facee7f5/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/writing/skills/aussie-business-english/SKILL.md) | ⭐ 543 | `communication` |
| [Recipes](communication/253-recipes_07b045a2/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/frontend/skills/shadcn-ui/references/recipes.md) | ⭐ 543 | `communication` |
| [Skill](communication/127-name-skill_bb4e4dfd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/angular-ui-patterns/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_156d7047/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-communication-chat-java/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_a829ab5d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-java/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_7dcc6e1f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-messaging-webpubsub-java/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_849ce4fa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-arizeaiobservabilityeval-dotnet/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_cd492dff/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-botservice-dotnet/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_b74f700d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-botservice-py/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_26477cd2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-weightsandbiases-dotnet/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_2e737a03/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-web-pubsub-ts/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_af005b8f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/commit/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_362f39d8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/customer-support/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_6e4fedb5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/email-systems/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_13e9cdf1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/form-cro/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_52bc8b8a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/launch-strategy/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_571aba9a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/n8n-node-configuration/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_59eb1273/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/popup-cro/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_660df10a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/slack-bot-builder/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_3fe12461/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/slack-gif-creator/SKILL.md) | 🔥 13.9k | `communication` |
| [Skill](communication/127-name-skill_e2326fed/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/upstash-qstash/SKILL.md) | 🔥 13.9k | `communication` |

### Content Creation (76 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Build Program](content-creation/366-build_program_6794f052/) | [julianghadially/CodeEvolver](https://raw.githubusercontent.com/julianghadially/CodeEvolver/main/specs/analysis/build_program.md) | ⭐ 12 | `content creation` |
| [Claude](content-creation/007-claude_bd86cf0c/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/CLAUDE.md) | ⭐ 10 | `content creation` |
| [Readme Cn](content-creation/366-readme_cn_e8c6ee87/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/README_CN.md) | ⭐ 10 | `content creation` |
| [Index](content-creation/019-index_4fd68aa3/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/docs/skills/index.md) | ⭐ 10 | `content creation` |
| [Index](content-creation/019-index_28dd84e3/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/docs/zh/skills/index.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/049-name-skill_f2f65ed3/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/SKILL.md) | ⭐ 10 | `content creation` |
| [4 Report Assembler](content-creation/384-4-report-assembler_b99ad74e/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/4-report-assembler.md) | ⭐ 2.9k | `content creation` |
| [5 Poc Generator](content-creation/385-5-poc-generator_18b5c964/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/5-poc-generator.md) | ⭐ 2.9k | `content creation` |
| [Skill](content-creation/049-name-skill_a4c5b6f1/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/codeql/SKILL.md) | ⭐ 2.9k | `content creation` |
| [Action Profiles](content-creation/386-action-profiles_ea755730/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/action-profiles.md) | ⭐ 2.9k | `content creation` |
| [Vector G Eval Of Ai Output](content-creation/387-vector-g-eval-of-ai-output_b9a47a41/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/vector-g-eval-of-ai-output.md) | ⭐ 2.9k | `content creation` |
| [Vector I Wildcard Allowlists](content-creation/388-vector-i-wildcard-allowlists_68a4cd77/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/vector-i-wildcard-allowlists.md) | ⭐ 2.9k | `content creation` |
| [21 The World](content-creation/389-21-the-world_ccef15a0/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/21-the-world.md) | ⭐ 2.9k | `content creation` |
| [Config Settings](content-creation/387-config_settings_8affec2b/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/config_settings.md) | 🔥 36.9k | `cache_hit` `cache_key` `proxy_base_url` |
| [Models Research](content-creation/384-models_research_1bcb2c64/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/section/models_research.md) | ⭐ 392 | `content creation` |
| [Johnw](content-creation/388-johnw_0c635fd9/) | [jwiegley/claude-prompts](https://raw.githubusercontent.com/jwiegley/claude-prompts/main/commands/johnw.md) | ⭐ 12 | `content creation` |
| [Skill](content-creation/049-name-skill_90a79ecb/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/process-siren/skills/improve-processes/SKILL.md) | ⭐ 20 | `content creation` |
| [Api](content-creation/072-api_a71b2d70/) | [timescale/pg-aiguide](https://raw.githubusercontent.com/timescale/pg-aiguide/main/API.md) | ⭐ 1.6k | `content creation` |
| [Skill](content-creation/049-name-skill_23bd51e3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ai-engineer/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_09fe7c87/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/algorithmic-art/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_e177214e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/architecture-decision-records/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_87e357ef/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-contentsafety-py/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_80fc976f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-contentsafety-ts/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_2415f7ea/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-contentunderstanding-py/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_9b6e60de/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-communication-callautomation-java/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_c375f893/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-containerregistry-py/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_af34007e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-blob-java/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_b012f396/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-blob-py/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_011d370e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-file-datalake-py/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_ccdc1044/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-file-share-ts/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_37bc1734/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bazel-build-optimization/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_e8d0421d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/canvas-design/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_ef6cc525/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/competitor-alternatives/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_91a54c0a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conductor-new-track/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_850b267e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conductor-revert/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_2c3ed8c8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/content-creator/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_a629db3d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/content-marketer/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_f4408fa8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-degradation/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_bb94b383/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/copy-editing/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_117f78bc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/copywriting/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_c8060c6d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/create-pr/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_e5585ed5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/crewai/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_4e7fd83e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/doc-coauthoring/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_d4b2e5af/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/docx-official/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_25d6b28f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/dropbox-automation/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_8caf6e98/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/file-organizer/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_dff5b5d5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/gemini-api-dev/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_0950c3bb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/git-advanced-workflows/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_edb8513c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-technologies/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_61445369/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/instagram-automation/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_ef746650/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/interactive-portfolio/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_3bf39775/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/internal-comms-anthropic/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_4ad85d4a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/internal-comms-community/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_0204a076/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/last30days/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_9fd97ba9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/lint-and-validate/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_b97cb9dc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/loki-mode/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_6125afcd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/marketing-psychology/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_a99d444e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/page-cro/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_781f3f82/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/paid-ads/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_bfd0a845/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/pdf-official/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_937b5d3e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/planning-with-files/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_0ac92a8d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/postmortem-writing/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_4bdaabd2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/programmatic-seo/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_51512227/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-flow-node-ts/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_e7f94bbc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/referral-program/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_279f1419/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/schema-markup/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_5ac416d9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-content-refresher/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_c6f715e4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-fundamentals/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_a2114334/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/social-content/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_bfa88b4f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tailwind-patterns/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_07ae9c3a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tiktok-automation/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_c6884e66/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/webflow-automation/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_845f1dc5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/whatsapp-automation/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_c68eb579/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/xlsx-official/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/049-name-skill_73ad2192/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/xss-html-injection/SKILL.md) | 🔥 13.9k | `content creation` |
| [Skill](content-creation/name-skill_0b6d8cb4/) | [marimo-team/skills](https://raw.githubusercontent.com/marimo-team/skills/main/skills/add-molab-badge/SKILL.md) | ⭐ 52 | `content creation` |

### Daily Assistant (50 skills)

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
| [Rust Zeroization Patterns](daily-assistant/270-rust-zeroization-patterns_5a86f5ac/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/references/rust-zeroization-patterns.md) | ⭐ 2.9k | `daily assistant` |
| [Phase 1 Source Analysis](daily-assistant/271-phase-1-source-analysis_23c0b46e/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/workflows/phase-1-source-analysis.md) | ⭐ 2.9k | `daily assistant` |
| [Phase 3 Interim Report](daily-assistant/272-phase-3-interim-report_f898af29/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/workflows/phase-3-interim-report.md) | ⭐ 2.9k | `daily assistant` |
| [Phase 4 Poc Generation](daily-assistant/273-phase-4-poc-generation_2f7a0daf/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/workflows/phase-4-poc-generation.md) | ⭐ 2.9k | `daily assistant` |
| [Phase 6 Final Report](daily-assistant/274-phase-6-final-report_9a318c4d/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/workflows/phase-6-final-report.md) | ⭐ 2.9k | `daily assistant` |
| [Ten Of Wands](daily-assistant/275-ten-of-wands_85a48ee3/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/wands/ten-of-wands.md) | ⭐ 2.9k | `daily assistant` |
| [Systemd Setup](daily-assistant/266-systemd-setup_f4c67252/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/systemd-setup.md) | ⭐ 177 | `daily assistant` |
| [Report](daily-assistant/289-report_a65e3fc4/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/agent/report.md) | ⭐ 643 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_8477c2b7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ai-agents-architect/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_0e2f8752/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/app-store-optimization/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_62bccbc2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/attack-tree-construction/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_bb8d9eab/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/behavioral-modes/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_44be1f6d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cloudformation-best-practices/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_a587e177/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/concise-planning/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_6a59d853/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conductor-status/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_67ce5626/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ddd-strategic-design/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_170f7d07/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/frontend-mobile-development-component-scaffold/SKILL.md) | 🔥 13.9k | `autodocs` |
| [Skill](daily-assistant/032-name-skill_d568dcf0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/go-concurrency-patterns/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_de98e3a5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hr-pro/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_9011e8d7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/incident-response-incident-response/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_ab86a6e6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/inventory-demand-planning/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_0a515ab0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/k8s-manifest-generator/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_437766ff/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/legal-advisor/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_d7ad3cb1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/microsoft-teams-automation/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_f25512fb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/on-call-handoff-patterns/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_8e8f8ba0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/onboarding-cro/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_464a133d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/outlook-automation/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_441ba572/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/outlook-calendar-automation/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_37c2dc8b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/production-scheduling/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_76929eb3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/risk-manager/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_520c0a2a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-authority-builder/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_94be79aa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-cannibalization-detector/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_bc9ccb43/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-meta-optimizer/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_0d12c35f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/subagent-driven-development/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_40beab76/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/team-collaboration-standup-notes/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_53d1142f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/threat-mitigation-mapping/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_7b8aa6ec/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/track-management/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_986e42f5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/using-superpowers/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_299e3396/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/verification-before-completion/SKILL.md) | 🔥 13.9k | `daily assistant` |
| [Index](daily-assistant/index_464fd71b/) | [Uniswap/uniswap-ai](https://raw.githubusercontent.com/Uniswap/uniswap-ai/main/docs/skills/index.md) | ⭐ 130 | `daily assistant` |
| [Claude](daily-assistant/claude_967a7d58/) | [Uniswap/uniswap-ai](https://raw.githubusercontent.com/Uniswap/uniswap-ai/main/packages/plugins/uniswap-driver/CLAUDE.md) | ⭐ 130 | `daily assistant` |

### Data Analysis (47 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [2 Source Analyzer](data-analysis/495-2-source-analyzer_ca0407ae/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/2-source-analyzer.md) | ⭐ 2.9k | `data analysis` |
| [3 Tu Compiler Analyzer](data-analysis/496-3-tu-compiler-analyzer_a875558c/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/3-tu-compiler-analyzer.md) | ⭐ 2.9k | `data analysis` |
| [3B Rust Compiler Analyzer](data-analysis/497-3b-rust-compiler-analyzer_f2a7fd92/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/3b-rust-compiler-analyzer.md) | ⭐ 2.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_369b998e/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/SKILL.md) | ⭐ 2.9k | `data analysis` |
| [System](data-analysis/498-system_4e70fd23/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/prompts/system.md) | ⭐ 2.9k | `data analysis` |
| [Compile Commands](data-analysis/499-compile-commands_ee40f46f/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/references/compile-commands.md) | ⭐ 2.9k | `data analysis` |
| [Detection Strategy](data-analysis/500-detection-strategy_acdf4242/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/references/detection-strategy.md) | ⭐ 2.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_d516ab63/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/lp-agent/SKILL.md) | ⭐ 11 | `data analysis` |
| [Memory Architecture](data-analysis/377-memory_architecture_c2a7a6b8/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/docs/memory_architecture.md) | ⭐ 763 | `data analysis` |
| [Nft Monitoring](data-analysis/495-nft-monitoring_be7eb631/) | [mensfeld/code-on-incus](https://raw.githubusercontent.com/mensfeld/code-on-incus/master/docs/NFT-MONITORING.md) | ⭐ 268 | `data analysis` |
| [Ui Ux Pro Max Skill](data-analysis/501-ui-ux-pro-max-skill_b408db6a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/ai-design-tools/ui-ux-pro-max-skill.md) | ⭐ 20 | `data analysis` |
| [Skill](data-analysis/226-name-skill_a7129ee9/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/integrations/skills/google-chat-messages/SKILL.md) | ⭐ 543 | `data analysis` |
| [Nfsd](data-analysis/491-nfsd_0b74b7f0/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/nfsd.md) | ⭐ 643 | `data analysis` |
| [Skill](data-analysis/226-name-skill_d8365148/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/android-jetpack-compose-expert/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_aad78aa3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-maps-search-dotnet/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_ed631af1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cc-skill-clickhouse-io/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_5e89a28b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/claude-d3js-skill/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_0a13ff30/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-management-context-save/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_e33a859e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-storytelling/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_43a9489e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/dependency-upgrade/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_0fb2410c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/design-md/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_ab0bbff8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/embedding-strategies/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_3e91f158/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/find-bugs/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_535ae6d5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/frontend-dev-guidelines/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_40de832a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/frontend-slides/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_ea881379/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/grafana-dashboards/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_8efae08f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-controls/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_07bd2ba1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-patterns/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_a74ffd4a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/kpi-dashboard-design/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_090982c0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/n8n-code-python/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_b54fd026/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/nextjs-best-practices/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_2a497fce/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/pptx-official/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_371dd4d6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-ui-patterns/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_40332b2f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/remotion-best-practices/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_ed230aa6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/screenshots/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_af36d77e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/service-mesh-observability/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_4fb4936f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/shodan-reconnaissance/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_2bb5a6b1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ui-ux-pro-max/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_1cd75c08/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/using-git-worktrees/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Skill](data-analysis/226-name-skill_00e4d228/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wireshark-analysis/SKILL.md) | 🔥 13.9k | `data analysis` |
| [Frontend](data-analysis/502-frontend_2a5cfb09/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks-apps/references/appkit/frontend.md) | ⭐ 13 | `data analysis` |
| [Overview](data-analysis/503-overview_0974fcfd/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks-apps/references/appkit/overview.md) | ⭐ 13 | `data analysis` |
| [Trpc](data-analysis/504-trpc_e55724ce/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks-apps/references/appkit/trpc.md) | ⭐ 13 | `data analysis` |
| [Phase4 Detection](data-analysis/phase4-detection_7dee7d81/) | [FrancyJGLisboa/agent-skill-creator](https://raw.githubusercontent.com/FrancyJGLisboa/agent-skill-creator/main/references/phase4-detection.md) | ⭐ 293 | `data analysis` |
| [Pipeline Phases](data-analysis/pipeline-phases_01041310/) | [FrancyJGLisboa/agent-skill-creator](https://raw.githubusercontent.com/FrancyJGLisboa/agent-skill-creator/main/references/pipeline-phases.md) | ⭐ 293 | `data analysis` |
| [Quality Standards](data-analysis/quality-standards_bf8c1ab8/) | [FrancyJGLisboa/agent-skill-creator](https://raw.githubusercontent.com/FrancyJGLisboa/agent-skill-creator/main/references/quality-standards.md) | ⭐ 293 | `data analysis` |
| [Pyodide Packages](data-analysis/pyodide-packages_bbb0f687/) | [marimo-team/skills](https://raw.githubusercontent.com/marimo-team/skills/main/skills/wasm-compatibility/references/pyodide-packages.md) | ⭐ 52 | `data analysis` |

### Development (224 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Openspec Interview Dimensions](development/2943-openspec_interview_dimensions_693c5b0e/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/content/skills/workflow-skills/interview-openspec/resources/OPENSPEC_INTERVIEW_DIMENSIONS.md) | ⭐ 10 | `development` |
| [Readme.Ja](development/1201-readmeja_35452de7/) | [japan1988/multi-agent-mediation](https://raw.githubusercontent.com/japan1988/multi-agent-mediation/main/README.ja.md) | ⭐ 29 | `development` |
| [0 Preflight](development/2943-0-preflight_b2407984/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/0-preflight.md) | ⭐ 2.9k | `development` |
| [5C Poc Verifier](development/2944-5c-poc-verifier_5b7725c8/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/5c-poc-verifier.md) | ⭐ 2.9k | `development` |
| [Skill](development/1178-name-skill_4d81229c/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/supply-chain-risk-auditor/skills/supply-chain-risk-auditor/SKILL.md) | ⭐ 2.9k | `development` |
| [Vector B Direct Expression Injection](development/2945-vector-b-direct-expression-injection_d074c95e/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/vector-b-direct-expression-injection.md) | ⭐ 2.9k | `development` |
| [Vector D Pr Target Checkout](development/2946-vector-d-pr-target-checkout_fc029037/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/vector-d-pr-target-checkout.md) | ⭐ 2.9k | `development` |
| [Extension Yaml Format](development/2947-extension-yaml-format_a263fa41/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/codeql/references/extension-yaml-format.md) | ⭐ 2.9k | `development` |
| [Quality Assessment](development/2948-quality-assessment_1dabaec5/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/codeql/references/quality-assessment.md) | ⭐ 2.9k | `development` |
| [Run Analysis](development/2831-run-analysis_8da60238/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/codeql/workflows/run-analysis.md) | ⭐ 2.9k | `development` |
| [Mcp Analysis](development/2949-mcp-analysis_704a9c6b/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/references/mcp-analysis.md) | ⭐ 2.9k | `development` |
| [Poc Generation](development/2950-poc-generation_61fb008c/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/references/poc-generation.md) | ⭐ 2.9k | `development` |
| [Phase 0 Preflight](development/2951-phase-0-preflight_cdf85ce4/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/workflows/phase-0-preflight.md) | ⭐ 2.9k | `development` |
| [Phase 5 Poc Validation](development/2952-phase-5-poc-validation_0cd58c2b/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/workflows/phase-5-poc-validation.md) | ⭐ 2.9k | `development` |
| [Queen Of Swords](development/2953-queen-of-swords_0b757ae9/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/swords/queen-of-swords.md) | ⭐ 2.9k | `development` |
| [Configuration](development/191-configuration_ede081fc/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/CONFIGURATION.md) | ⭐ 4.0k | `development` |
| [Conventions](development/2930-conventions_b0e7f580/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/CONVENTIONS.md) | ⭐ 177 | `development` |
| [Plugin Authoring](development/2931-plugin_authoring_766dc887/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/PLUGIN_AUTHORING.md) | ⭐ 177 | `development` |
| [Plugin System Documentation](development/2932-plugin_system_documentation_f214d0af/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/PLUGIN_SYSTEM_DOCUMENTATION.md) | ⭐ 177 | `development` |
| [Documentation Guide](development/2933-documentation-guide_b3c0b700/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/documentation-guide.md) | ⭐ 177 | `development` |
| [Metrics Api](development/2934-metrics-api_954803c6/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/metrics-api.md) | ⭐ 177 | `development` |
| [Debugging With Proxy](development/2935-debugging-with-proxy_e4c3deb0/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/development/debugging-with-proxy.md) | ⭐ 177 | `development` |
| [Api Usage](development/047-api-usage_0fa5b854/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/user-guide/api-usage.md) | ⭐ 177 | `development` |
| [Claude Sdk Compatibility](development/2936-claude-sdk-compatibility_2d74e6d3/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/user-guide/claude-sdk-compatibility.md) | ⭐ 177 | `development` |
| [Codex Api](development/2937-codex-api_5ad2064a/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/user-guide/codex-api.md) | ⭐ 177 | `development` |
| [Applications](development/2944-applications_c92cf5e8/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/section/applications.md) | ⭐ 392 | `development` |
| [Azure](development/084-azure_e856bbe7/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/section/azure.md) | ⭐ 392 | `development` |
| [Best Practices](development/103-best_practices_998cf88f/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/section/best_practices.md) | ⭐ 392 | `development` |
| [Tools Extra](development/2945-tools_extra_469d0147/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/section/tools_extra.md) | ⭐ 392 | `development` |
| [X Popular Papers](development/2946-x_popular_papers_f6579aef/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/section/x_popular_papers.md) | ⭐ 392 | `development` |
| [Fix](development/2942-fix_4334049b/) | [jwiegley/claude-prompts](https://raw.githubusercontent.com/jwiegley/claude-prompts/main/commands/fix.md) | ⭐ 12 | `development` |
| [Azure](development/084-azure_bb4e38b3/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/section/azure.md) | ⭐ 392 | `development` |
| [Research Agent Assessment](development/2959-research-agent-assessment_c6d81a5d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/audits/research-agent-assessment.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_9d69f716/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/1178-name-skill_07bba6a0/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/vite-flare-starter/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/1178-name-skill_1891cee6/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/gemini-guide/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/1178-name-skill_6c7e72be/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/github-release/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/1178-name-skill_b0e0b759/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/skill-creator/SKILL.md) | ⭐ 543 | `development` |
| [Deployment](development/272-deployment_28c3c43e/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/references/deployment.md) | ⭐ 543 | `development` |
| [Argument Parsing](development/2978-argument-parsing_f46ba92e/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/argument-parsing.md) | ⭐ 643 | `development` |
| [Coding Style](development/2979-coding-style_1338db47/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/coding-style.md) | ⭐ 643 | `development` |
| [Common Bugs](development/2980-common-bugs_e52fe6c4/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/common-bugs.md) | ⭐ 643 | `development` |
| [Kernel Compat](development/2981-kernel-compat_2db6869f/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/kernel-compat.md) | ⭐ 643 | `development` |
| [Netlink](development/2982-netlink_b4900310/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/netlink.md) | ⭐ 643 | `development` |
| [Patch Submission](development/2983-patch-submission_32be650f/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/patch-submission.md) | ⭐ 643 | `development` |
| [Review Core](development/2984-review-core_89eff333/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/review-core.md) | ⭐ 643 | `development` |
| [Technical Patterns](development/2985-technical-patterns_0e516b8b/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/technical-patterns.md) | ⭐ 643 | `development` |
| [Iproute Verify](development/2986-iproute-verify_fdf36dc1/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/slash-commands/iproute-verify.md) | ⭐ 643 | `development` |
| [Mm Folio](development/2987-mm-folio_2fa17c29/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/mm-folio.md) | ⭐ 643 | `development` |
| [Subsystem](development/2988-subsystem_424e4d1c/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/subsystem/subsystem.md) | ⭐ 643 | `development` |
| [Catalog](development/126-catalog_6a96f5ba/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/CATALOG.md) | 🔥 13.9k | `development` |
| [Skill Anatomy](development/984-skill_anatomy_0e01cdaf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/docs/SKILL_ANATOMY.md) | 🔥 13.9k | `react` `typescript` |
| [Skill](development/1178-name-skill_a9f1a082/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/address-github-comments/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9dbfc671/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/agent-orchestration-improve-agent/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f4ffc185/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/agent-tool-builder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_4e419747/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/agentfolio/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c339f994/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ai-product/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9f56f6a5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/angular-migration/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_bf04cfe0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/api-documentation-generator/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_dc68be9b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/api-documentation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f2bef58c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/api-documenter/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_4c8901e5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/api-fuzzing-bug-bounty/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_86c5334c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/api-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5d797c24/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/api-security-best-practices/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_ea7b1db8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/automate-whatsapp/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_3413f991/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/autonomous-agent-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f64b74ce/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/aws-cost-cleanup/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_83098164/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azd-deployment/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_edf6b56b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-agents-persistent-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c3fcd976/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-openai-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_b246c8fa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-projects-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_b5ea8ef8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-projects-java/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_b4c4874a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-projects-ts/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_0c1522af/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-vision-imageanalysis-java/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_20012eb4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-vision-imageanalysis-py/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_bb097fb7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-voicelive-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_dbb2e139/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-voicelive-py/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_000e8ddd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-voicelive-ts/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_a1e9b4f2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-rust/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_09381c9e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-functions/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_870befaa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-identity-java/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_72d8739b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-keyvault-certificates-rust/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_cecadb52/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-keyvault-secrets-rust/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d4bf7b6e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-apimanagement-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_b04dec20/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-applicationinsights-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f0a94ebb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-fabric-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_b7d74878/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-monitor-opentelemetry-exporter-java/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_babf686c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-cosmosdb-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5bf7cc76/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-durabletask-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d15424ed/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-playwright-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_bc6406b2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-redis-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_13fceb94/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-storage-blob-ts/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_1dba608e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/backend-dev-guidelines/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_7c1bd7c7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/backend-development-feature-development/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_1d6732a9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/backend-security-coder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_985a50e1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bash-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_1d9db141/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/blockchain-developer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_07f04bea/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/brainstorming/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d036c78c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/broken-authentication/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_3b1776ad/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/browser-extension-builder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_fb1cb725/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bun-development/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d2752909/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/c4-architecture-c4-architecture/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_cdce2c63/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cc-skill-coding-standards/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_cce5ee0b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cc-skill-frontend-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_1fb4423a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/clean-code/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_1cbe0a46/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cloudflare-workers-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_a6bee09c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-documentation-code-explain/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_3a6a7155/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-refactoring-refactor-clean/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c6770b4c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-refactoring-tech-debt/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_49277182/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-review-checklist/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_473455d6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-review-excellence/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_09296810/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-reviewer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_37bd8f3b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/codebase-cleanup-tech-debt/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_8850345f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/codex-review/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5ae223f8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conductor-implement/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_afd3bef5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/conductor-setup/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_369221f3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-compression/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d4814a18/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-fundamentals/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_8f93adcf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-manager/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9aa9bed1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/copilot-sdk/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_8431741b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-structure-protocol/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f2586412/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/dbos-python/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_ef69b5d9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ddd-tactical-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_cf02ed3e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/debugging-toolkit-smart-debug/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9175a877/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/discord-bot-architect/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d4057c47/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/docs-architect/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_29a5c57c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/documentation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_32cd3ff4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/dotnet-backend/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c7f2cbe2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/environment-setup-guide/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_39151f3a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/error-debugging-multi-agent-review/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_563f1df5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/error-diagnostics-smart-debug/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_0a85b8f2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/file-path-traversal/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_232265a8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/firebase/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_a6e1b348/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/firmware-analyst/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9f47ef79/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/fix-review/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_326752d0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/flutter-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_da6fa41e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/fp-ts-react/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_30b5ebbe/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/framework-migration-legacy-modernize/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_6569e2ee/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/free-tool-strategy/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_4809bfe3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/frontend-design/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_036b2dcf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/frontend-developer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_cfa04aef/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/frontend-security-coder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_22ffede3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/full-stack-orchestration-full-stack-feature/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_84496c7b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/git-pr-workflows-git-workflow/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_851c46ae/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/godot-4-migration/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_660eff3f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/graphql-architect/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_914875b7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/graphql/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_0b0bbb22/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-dialogs/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_a4699768/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-status/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_89fe563e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-project-context/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_bb3a7f31/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ios-developer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_85219543/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/iterate-pr/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_1b457144/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/java-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_0d5b3556/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/javascript-mastery/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d4bc4852/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/kaizen/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_6800921c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/langgraph/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9aac3505/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/laravel-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_6cc96eb1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/legacy-modernizer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_d2e826f0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/llm-app-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_b5bfb36b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/llm-application-dev-langchain-agent/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_08c3aa80/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/malware-analyst/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_65267a13/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/memory-forensics/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_4168bd08/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/metasploit-framework/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_2740aa53/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/microsoft-azure-webjobs-extensions-authentication-events-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5657deb1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/minecraft-bukkit-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f87d24b6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mobile-design/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f00f532f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mobile-developer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_44985d91/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mobile-security-coder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_43824f5f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/moodle-external-api-development/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f63ffff0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/multi-platform-apps-multi-platform/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_db3272aa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/paypal-integration/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5b3d7ab4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/personal-tool-builder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c9616ee3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/podcast-generation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_709f17a0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/posix-shell-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_65ec4f21/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/prompt-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_86cd5a8e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/prompt-engineering/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5df40654/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/prompt-library/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5c07186e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/python-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_f6781a2c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/radix-ui-design-system/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9caf8e39/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-modernization/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_433c235f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/receiving-code-review/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_99146faa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/reference-builder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5383bdb8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/requesting-code-review/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_106acb3a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/research-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_48a3fd0a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/reverse-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_0ca4bfc2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/sast-configuration/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_6f8fde3b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/scala-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c9e03685/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security-auditor/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_a5adcc6d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security-scanning-security-sast/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_db3e7cf2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/sharp-edges/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_b87a4cd8/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/shellcheck-configuration/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_69db7e48/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/skill-creator-ms/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_3aee5180/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/skill-developer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_04eed10c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/skill-rails-upgrade/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c695fdcb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/software-architecture/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_23085fea/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/spark-optimization/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_c7fad5db/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/swiftui-expert-skill/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_e8e96c7f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/systematic-debugging/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_61b250c4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/systems-programming-rust-project/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_42614ddb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-workflow/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_17369ad0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-workflows-tdd-refactor/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_5c047297/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/team-composition-analysis/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_bad27e11/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/telegram-mini-app/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_0777edd5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tool-design/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_12ba576e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/top-web-vulnerabilities/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_610e2b12/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tutorial-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_a56d5e86/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/typescript-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_64b3c507/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ui-ux-designer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_0834ae01/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/unity-developer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_fc131b59/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/upgrading-expo/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_99b83d7b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/voice-ai-development/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_9c35a326/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/voice-ai-engine-development/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_658a3ca4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wordpress/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_80f90720/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/workflow-automation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_2f24891a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/writing-plans/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_8934aa35/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/youtube-summarizer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_7f15798a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/zapier-make-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_fdae94f6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security/aws-secrets-rotation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/1178-name-skill_a95e8f64/) | [tsaol/awesome-claude](https://raw.githubusercontent.com/tsaol/awesome-claude/main/skills/git-commit/SKILL.md) | ⭐ 41 | `development` |
| [Skill](development/name-skill_d3573840/) | [FrancyJGLisboa/agent-skill-creator](https://raw.githubusercontent.com/FrancyJGLisboa/agent-skill-creator/main/SKILL.md) | ⭐ 293 | `development` |
| [Agent Skill Creator Full Brief](development/agent-skill-creator-full-brief_70aab80a/) | [FrancyJGLisboa/agent-skill-creator](https://raw.githubusercontent.com/FrancyJGLisboa/agent-skill-creator/main/Dynamous/Content-Ideation/agent-skill-creator-full-brief.md) | ⭐ 293 | `development` |
| [Overview](development/overview_5c2e86e3/) | [Uniswap/uniswap-ai](https://raw.githubusercontent.com/Uniswap/uniswap-ai/main/docs/OVERVIEW.md) | ⭐ 130 | `development` |
| [Skill](development/name-skill_8ea64245/) | [Uniswap/uniswap-ai](https://raw.githubusercontent.com/Uniswap/uniswap-ai/main/packages/plugins/uniswap-trading/skills/swap-integration/SKILL.md) | ⭐ 130 | `development` |
| [Sql](development/sql_5c2a2f0a/) | [marimo-team/skills](https://raw.githubusercontent.com/marimo-team/skills/main/skills/marimo-notebook/references/SQL.md) | ⭐ 52 | `development` |

### Development/Devops (73 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Scan Workflow](development/devops/374-scan-workflow_de63106a/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/semgrep/workflows/scan-workflow.md) | ⭐ 2.9k | `development` |
| [Claude](development/devops/205-claude_c917be58/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/CLAUDE.md) | ⭐ 11 | `development` |
| [Release Notes Generation Instructions](development/devops/020-release_notes_generation_instructions_65a0ecf2/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/cookbook/misc/RELEASE_NOTES_GENERATION_INSTRUCTIONS.md) | 🔥 36.9k | `development` |
| [Configuration](development/devops/009-configuration_39c227e8/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/getting-started/configuration.md) | ⭐ 177 | `development` |
| [Mcp Integration](development/devops/172-mcp-integration_d8ae2e35/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/user-guide/mcp-integration.md) | ⭐ 177 | `development` |
| [Mcp Debug Research](development/devops/387-mcp-debug-research_fa166c9a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.planning/mcp-debug-research.md) | ⭐ 20 | `development` |
| [Micro Agent](development/devops/388-micro-agent_19ed3c9a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/agent-frameworks/micro-agent.md) | ⭐ 20 | `development` |
| [Ra Aid](development/devops/111-ra-aid_ea31434d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/agent-frameworks/ra-aid.md) | ⭐ 20 | `development` |
| [Motia](development/devops/389-motia_636fe247/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/api-frameworks/motia.md) | ⭐ 20 | `development` |
| [Pocketbase](development/devops/390-pocketbase_ed64764b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/api-frameworks/pocketbase.md) | ⭐ 20 | `development` |
| [Openhands](development/devops/112-openhands_3db9bc92/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/coding-agents/openhands.md) | ⭐ 20 | `development` |
| [Microsoft Graphrag](development/devops/113-microsoft-graphrag_72e6a096/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/context-management/microsoft-graphrag.md) | ⭐ 20 | `development` |
| [Claude Quickstarts](development/devops/391-claude-quickstarts_b1f6f524/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/developer-tools/claude-quickstarts.md) | ⭐ 20 | `development` |
| [Copier Astral](development/devops/114-copier-astral_c79630e2/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/developer-tools/copier-astral.md) | ⭐ 20 | `development` |
| [Devenv](development/devops/392-devenv_1889434b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/developer-tools/devenv.md) | ⭐ 20 | `development` |
| [Everything Claude Code](development/devops/393-everything-claude-code_0d9c4c88/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/developer-tools/everything-claude-code.md) | ⭐ 20 | `development` |
| [Bifrost](development/devops/394-bifrost_81c5f947/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/llm-infrastructure/bifrost.md) | ⭐ 20 | `development` |
| [Docs Mcp Server](development/devops/117-docs-mcp-server_79dbf0ef/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/mcp-ecosystem/docs-mcp-server.md) | ⭐ 20 | `development` |
| [Mcpjam](development/devops/118-mcpjam_628da6fa/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/mcp-ecosystem/mcpjam.md) | ⭐ 20 | `development` |
| [Skill](development/devops/014-name-skill_f2fb6f03/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/project-health/SKILL.md) | ⭐ 543 | `development` |
| [Architecture](development/devops/032-architecture_8a5a2c86/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/references/architecture.md) | ⭐ 543 | `development` |
| [Customization Guide](development/devops/375-customization-guide_69b157f0/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/vite-flare-starter/references/customization-guide.md) | ⭐ 543 | `development` |
| [Devto V060](development/devops/376-devto-v060_d158348a/) | [doramirdor/NadirClaw](https://raw.githubusercontent.com/doramirdor/NadirClaw/main/docs/devto-v060.md) | ⭐ 236 | `development` |
| [Skill](development/devops/014-name-skill_c061bf2e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/architect-review/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_fbf3ef60/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-apicenter-dotnet/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_f78f007f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-apicenter-py/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_0f3a25c0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/backend-architect/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_97e53614/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cc-skill-security-review/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_71988ee2/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/cloud-architect/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_ea6d824b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/computer-use-agents/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_04398e50/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-driven-development/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_4ce73560/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_405793ae/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-engineering-data-pipeline/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_ed2662a7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/database-admin/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_106bf3ce/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/deployment-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_7098fdfc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/deployment-pipeline-design/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_061aa032/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/deployment-procedures/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_dd68e389/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/development/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_79342d79/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/devops-troubleshooter/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_1a536efe/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/distributed-tracing/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_2926002a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/docker-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_1c99c612/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/dotnet-architect/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_25451b40/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/fastapi-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_30273e6b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/gcp-cloud-run/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_7f9e1bfb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/git-pr-workflows-onboard/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_6ff88568/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/github-automation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_5fc9acdd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/gitlab-ci-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_2b6a738e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/golang-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_bb43b2dd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/grpc-golang/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_c75119b7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/istio-traffic-management/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_5e753cb9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/kubernetes-architect/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_aafbcfa0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/llm-evaluation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_ab1ca318/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/m365-agents-ts/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_6ca2a02b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/machine-learning-ops-ml-pipeline/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_d72f8e60/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mcp-builder-ms/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_03101a4b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mlops-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_14ba1db4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/multi-cloud-architecture/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_89573200/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/n8n-mcp-tools-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_aa9d558a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/network-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_9ea0fcf6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/nodejs-best-practices/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_7eff94bb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/parallel-agents/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_70458c9a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/performance-engineer/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_cd3bccd5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/prisma-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_2b178567/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/production-code-audit/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_915a91f6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/prometheus-configuration/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_d388fe57/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/readme/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_87496fec/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/secrets-management/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_11abfede/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/server-management/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_0e1f1460/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/sql-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_02b51c7a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/temporal-python-pro/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_71341d0f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/terraform-skill/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_a2d598f6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/vercel-deployment/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/devops/014-name-skill_73fec275/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks-apps/SKILL.md) | ⭐ 13 | `development` |

### Development/Testing (15 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/testing/002-name-skill_0af17ce7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-orchestrator/SKILL.md) | 🔥 15.4k | `development` |
| [Oauth Plugin Architecture](development/testing/084-oauth_plugin_architecture_8d4366d2/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/OAUTH_PLUGIN_ARCHITECTURE.md) | ⭐ 177 | `development` |
| [Auth Providers](development/testing/085-auth-providers_346cccf7/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/development/auth-providers.md) | ⭐ 177 | `development` |
| [Skill](development/testing/002-name-skill_21bdc8b7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/angular/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_921fbe82/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/anti-reversing-techniques/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_5a2a0c0b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-db-py/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_060af7e9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/browser-automation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_a836fb2a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/core-components/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_0f6a6b1d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/go-rod-master/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_36f22f62/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/langchain-architecture/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_577d770f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/plan-writing/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_295a00a1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/python-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_23c7b055/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-orchestrator/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_608ca966/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-workflows-tdd-cycle/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/testing/002-name-skill_3667e087/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-workflows-tdd-red/SKILL.md) | 🔥 13.9k | `development` |

### Development/Tools (51 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agents](development/tools/015-agents_15ebd473/) | [bahayonghang/my-claude-code-settings](https://raw.githubusercontent.com/bahayonghang/my-claude-code-settings/master/AGENTS.md) | ⭐ 10 | `development` |
| [2B Rust Source Analyzer](development/tools/368-2b-rust-source-analyzer_220f9c1f/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/2b-rust-source-analyzer.md) | ⭐ 2.9k | `development` |
| [Skill](development/tools/002-name-skill_0764b487/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/seatbelt-sandboxer/skills/seatbelt-sandboxer/SKILL.md) | ⭐ 2.9k | `development` |
| [Skill](development/tools/002-name-skill_4d485cc0/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/second-opinion/skills/second-opinion/SKILL.md) | ⭐ 2.9k | `development` |
| [Skill](development/tools/002-name-skill_b12e2130/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/semgrep/SKILL.md) | ⭐ 2.9k | `development` |
| [Vector A Env Var Intermediary](development/tools/369-vector-a-env-var-intermediary_7903992f/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/vector-a-env-var-intermediary.md) | ⭐ 2.9k | `development` |
| [Vector C Cli Data Fetch](development/tools/370-vector-c-cli-data-fetch_f11d2333/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/vector-c-cli-data-fetch.md) | ⭐ 2.9k | `development` |
| [Gemini Invocation](development/tools/312-gemini-invocation_e9fd259e/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/second-opinion/skills/second-opinion/references/gemini-invocation.md) | ⭐ 2.9k | `development` |
| [Create Data Extensions](development/tools/371-create-data-extensions_2210ce14/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/codeql/workflows/create-data-extensions.md) | ⭐ 2.9k | `development` |
| [Code Index](development/tools/282-code_index_d1c06dc8/) | [plexe-ai/plexe](https://raw.githubusercontent.com/plexe-ai/plexe/main/plexe/CODE_INDEX.md) | ⭐ 2.5k | `development` |
| [Skill](development/tools/002-name-skill_ad5c103c/) | [HeshamFS/materials-simulation-skills](https://raw.githubusercontent.com/HeshamFS/materials-simulation-skills/main/skills/ontology/ontology-explorer/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/tools/002-name-skill_b244adfc/) | [HeshamFS/materials-simulation-skills](https://raw.githubusercontent.com/HeshamFS/materials-simulation-skills/main/skills/ontology/ontology-validator/SKILL.md) | ⭐ 20 | `development` |
| [Vertex Realtime](development/tools/366-vertex_realtime_fdde2fa3/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/providers/vertex_realtime.md) | 🔥 36.9k | `development` |
| [Agents](development/tools/015-agents_647ea66b/) | [gptme/gptme](https://raw.githubusercontent.com/gptme/gptme/master/AGENTS.md) | ⭐ 4.2k | `development` |
| [Claude Code Options](development/tools/367-claude-code-options_3b4b3a1f/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/user-guide/claude-code-options.md) | ⭐ 177 | `development` |
| [Helpers](development/tools/368-helpers_24dafc5c/) | [anthropics/anthropic-sdk-python](https://raw.githubusercontent.com/anthropics/anthropic-sdk-python/main/helpers.md) | ⭐ 2.8k | `development` |
| [Skill](development/tools/002-name-skill_d043a7e0/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/gemini-peer-review/SKILL.md) | ⭐ 543 | `development` |
| [Iproute Review](development/tools/342-iproute-review_943594e2/) | [masoncl/review-prompts](https://raw.githubusercontent.com/masoncl/review-prompts/main/iproute/slash-commands/iproute-review.md) | ⭐ 643 | `development` |
| [Skill](development/tools/002-name-skill_3d187d78/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/agents-v2-py/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_d1270001/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/analytics-tracking/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_daae890f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/appdeploy/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_763107a0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/autonomous-agents/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_db86a1a4/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-ai-voicelive-java/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_6b2ca0d5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-py/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_5379e338/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/database-migration/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_475153d3/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/ethical-hacking-methodology/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_8c3d602d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/go-playwright/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_8eaa9bd1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-layout/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_c57c7b3d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-menus/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_2904aeab/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-search/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_ff337a33/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-inputs/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_70cb86ec/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hosted-agents-v2-py/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_2be94c36/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/langfuse/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_17cf7a40/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/linear-claude-skill/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_02587bbb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/mcp-builder/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_7d75334d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/nestjs-expert/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_3c992361/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/oss-hunter/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_c70a6be1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/playwright-skill/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_add8b685/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/privilege-escalation-methods/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_e70ed24a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/python-development-python-scaffold/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_1920babd/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/react-state-management/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_b2a44d9f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/scroll-experience/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_129ef2d9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/threejs-skills/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_0d7ca31d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/windows-privilege-escalation/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_0bac9fc0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/workflow-patterns/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_bbe62bfa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/writing-skills/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_aba7d708/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/security/aws-security-audit/SKILL.md) | 🔥 13.9k | `development` |
| [Skill](development/tools/002-name-skill_fabb3077/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks/SKILL.md) | ⭐ 13 | `development` |
| [Databricks Cli Install](development/tools/372-databricks-cli-install_ad904585/) | [databricks/databricks-agent-skills](https://raw.githubusercontent.com/databricks/databricks-agent-skills/main/skills/databricks/databricks-cli-install.md) | ⭐ 13 | `development` |
| [Code Index](development/tools/282-code_index_a8e950fb/) | [plexe-ai/plexe](https://raw.githubusercontent.com/plexe-ai/plexe/main/plexe/CODE_INDEX.md) | ⭐ 2.5k | `development` |
| [Skill](development/tools/name-skill_522e2ffc/) | [marimo-team/skills](https://raw.githubusercontent.com/marimo-team/skills/main/skills/wasm-compatibility/SKILL.md) | ⭐ 52 | `development` |

### Investment (19 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](investment/021-name-skill_3e0f379b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-scientist/SKILL.md) | 🔥 15.4k | `investment` |
| [Skill](investment/021-name-skill_0db988ba/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/SKILL.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/049-decision-frameworks_95ad10a0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/carrier-relationship-management/references/decision-frameworks.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/049-decision-frameworks_30ceb8b9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/references/decision-frameworks.md) | 🔥 15.4k | `investment` |
| [Seven Of Cups](investment/052-seven-of-cups_852b35c0/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/cups/seven-of-cups.md) | ⭐ 2.9k | `investment` |
| [Nine Of Pentacles](investment/053-nine-of-pentacles_272b240f/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/pentacles/nine-of-pentacles.md) | ⭐ 2.9k | `investment` |
| [Skill](investment/021-name-skill_ea8c4627/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/connectors-available/SKILL.md) | ⭐ 11 | `investment` |
| [Server Functions](investment/054-server-functions_535f27be/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/references/server-functions.md) | ⭐ 543 | `investment` |
| [Skill](investment/021-name-skill_a4c1aa27/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/carrier-relationship-management/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_897e53ad/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-scientist/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_3a78ddaa/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_da9ae21c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/marketing-ideas/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_dc977f7c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/pricing-strategy/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_1102cc5f/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/product-manager-toolkit/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_8301cf1d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/quality-nonconformance/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_54dac2de/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/risk-metrics-calculation/SKILL.md) | 🔥 13.9k | `investment` |
| [Skill](investment/021-name-skill_d7f692ed/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-business-analyst-business-case/SKILL.md) | 🔥 13.9k | `investment` |
| [Multi Agent Guide](investment/multi-agent-guide_db50cf80/) | [FrancyJGLisboa/agent-skill-creator](https://raw.githubusercontent.com/FrancyJGLisboa/agent-skill-creator/main/references/multi-agent-guide.md) | ⭐ 293 | `investment` |
| [Skill](investment/name-skill_ccd2df2b/) | [FrancyJGLisboa/agent-skill-creator](https://raw.githubusercontent.com/FrancyJGLisboa/agent-skill-creator/main/references/examples/stock-analyzer/SKILL.md) | ⭐ 293 | `investment` |

### Other (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [04 The Emperor](other/036-04-the-emperor_1b48cba1/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/04-the-emperor.md) | ⭐ 2.9k | `other` |
| [Eight Of Swords](other/037-eight-of-swords_24e30c87/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/swords/eight-of-swords.md) | ⭐ 2.9k | `other` |

### Productivity (21 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Four Of Cups](productivity/174-four-of-cups_e7c80a02/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/cups/four-of-cups.md) | ⭐ 2.9k | `productivity` |
| [07 The Chariot](productivity/175-07-the-chariot_fd9f65cd/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/07-the-chariot.md) | ⭐ 2.9k | `productivity` |
| [Two Of Pentacles](productivity/176-two-of-pentacles_95edc988/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/pentacles/two-of-pentacles.md) | ⭐ 2.9k | `productivity` |
| [Knight Of Wands](productivity/177-knight-of-wands_4ae17136/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/wands/knight-of-wands.md) | ⭐ 2.9k | `productivity` |
| [Page Of Wands](productivity/178-page-of-wands_6da5cc15/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/wands/page-of-wands.md) | ⭐ 2.9k | `productivity` |
| [Lp Rebalancer Guide](productivity/174-lp_rebalancer_guide_63b97f02/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/lp-agent/references/lp_rebalancer_guide.md) | ⭐ 11 | `productivity` |
| [Skill](productivity/093-name-skill_35b8645b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/arm-cortex-expert/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_b4255575/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-mysql-dotnet/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_75206fac/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-resource-manager-postgresql-dotnet/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_6d4ae24a/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/bevy-ecs-expert/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_bf1af702/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-refactoring-context-restore/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_4d31fa1b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context-management-context-restore/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_1d1fb220/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/geo-fundamentals/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_e3d48d2b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/hig-components-system/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_a2d6bb43/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/micro-saas-launcher/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_a804d53e/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/multi-agent-brainstorming/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_f9f8f703/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/nosql-expert/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_5f22a935/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/seo-audit/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_a7855f23/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/threat-modeling-expert/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_5a2cd2a5/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/viral-generator-builder/SKILL.md) | 🔥 13.9k | `productivity` |
| [Skill](productivity/093-name-skill_8f90834d/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wiki-changelog/SKILL.md) | 🔥 13.9k | `productivity` |

### Research (27 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](research/139-name-skill_40343acb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-business-analyst-market-opportunity/SKILL.md) | 🔥 15.4k | `research` |
| [Build Spec](research/265-build-spec_f3072f59/) | [Nikolay-Shirokov/cc-1c-skills](https://raw.githubusercontent.com/Nikolay-Shirokov/cc-1c-skills/main/docs/build-spec.md) | ⭐ 95 | `research` |
| [Skill](research/139-name-skill_b6e33b5c/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-analyst/SKILL.md) | 🔥 15.4k | `research` |
| [Skill](research/139-name-skill_12e80dcd/) | [fivetaku/deep-research-kit](https://raw.githubusercontent.com/fivetaku/deep-research-kit/main/skills/deep-research-main/SKILL.md) | ⭐ 15 | `research` |
| [Skill](research/139-name-skill_767f4ca9/) | [fivetaku/deep-research-kit](https://raw.githubusercontent.com/fivetaku/deep-research-kit/main/skills/deep-research-query/SKILL.md) | ⭐ 15 | `research` |
| [Citation Rules](research/257-citation_rules_cb28ab52/) | [fivetaku/deep-research-kit](https://raw.githubusercontent.com/fivetaku/deep-research-kit/main/skills/deep-research-main/references/citation_rules.md) | ⭐ 15 | `research` |
| [Quality Rubric](research/258-quality_rubric_5c8b32c0/) | [fivetaku/deep-research-kit](https://raw.githubusercontent.com/fivetaku/deep-research-kit/main/skills/deep-research-main/references/quality_rubric.md) | ⭐ 15 | `research` |
| [Query Generator](research/259-query_generator_0192eac0/) | [fivetaku/deep-research-kit](https://raw.githubusercontent.com/fivetaku/deep-research-kit/main/skills/deep-research-main/references/query_generator.md) | ⭐ 15 | `research` |
| [Tool Strategy](research/260-tool_strategy_14a1b963/) | [fivetaku/deep-research-kit](https://raw.githubusercontent.com/fivetaku/deep-research-kit/main/skills/deep-research-main/references/tool_strategy.md) | ⭐ 15 | `research` |
| [09 The Hermit](research/261-09-the-hermit_9d747c35/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/09-the-hermit.md) | ⭐ 2.9k | `research` |
| [Two Of Swords](research/262-two-of-swords_a633f81b/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/swords/two-of-swords.md) | ⭐ 2.9k | `research` |
| [Asmo Guide](research/258-asmo_guide_0a4c01c8/) | [HeshamFS/materials-simulation-skills](https://raw.githubusercontent.com/HeshamFS/materials-simulation-skills/main/skills/ontology/ontology-explorer/references/asmo_guide.md) | ⭐ 20 | `research` |
| [Skill](research/139-name-skill_e113c43e/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/.github/skills/add-new-entry/SKILL.md) | ⭐ 392 | `research` |
| [Skill](research/139-name-skill_1ee806ea/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/.github/skills/update-cite-count/SKILL.md) | ⭐ 392 | `research` |
| [Completeness Report Plugin Creator](research/263-completeness-report-plugin-creator_bf91f9a2/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/audits/completeness-report-plugin-creator.md) | ⭐ 20 | `research` |
| [Skill](research/139-name-skill_3edd3a77/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-mongodbatlas-dotnet/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_19691250/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/blockrun/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_9dc09518/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/context7-auto-research/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_358060bc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/daily-news-report/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_569b0269/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/evaluation/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_0eb19edf/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/exa-search/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_703ba365/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/notebooklm/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_092aec56/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-analyst/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_11a300c1/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/startup-business-analyst-market-opportunity/SKILL.md) | 🔥 13.9k | `research` |
| [Skill](research/139-name-skill_34edf9bc/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/wiki-onboarding/SKILL.md) | 🔥 13.9k | `research` |
| [Liquidity Planner](research/liquidity-planner_39a0a3d3/) | [Uniswap/uniswap-ai](https://raw.githubusercontent.com/Uniswap/uniswap-ai/main/docs/skills/liquidity-planner.md) | ⭐ 130 | `research` |
| [Swap Planner](research/swap-planner_5e8b8727/) | [Uniswap/uniswap-ai](https://raw.githubusercontent.com/Uniswap/uniswap-ai/main/docs/skills/swap-planner.md) | ⭐ 130 | `research` |

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

*Last updated: 2026-02-27 16:14:42 UTC*
*Automatically maintained by SkillFlow*
