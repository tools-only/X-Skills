# X-Skills

A curated collection of **78 AI-powered skills** organized into 8 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Workflow** (2 skills)
- **Communication** (12 skills)
- **Content Creation** (17 skills)
- **Development** (31 skills)
- **Development/Devops** (2 skills)
- **Development/Testing** (1 skill)
- **Development/Tools** (12 skills)
- **Research** (1 skill)

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


### Automation/Workflow (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Executive Pitch](automation/workflow/067-executive-pitch_61fed1da/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/executive-pitch.md) | ⭐ 60 | `automation` |
| [Detecting Llm Hallucinations In Ci](automation/workflow/138-detecting-llm-hallucinations-in-ci_471c775c/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/guides/detecting-llm-hallucinations-in-ci.md) | ⭐ 43 | `automation` |

### Communication (12 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Help](communication/251-help_dcae5407/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/help.md) | 🔥 40.5k | `communication` |
| [Model Warnings](communication/252-model-warnings_ba625717/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/model-warnings.md) | 🔥 40.5k | `communication` |
| [Multi Line](communication/253-multi-line_5264f142/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/multi-line.md) | 🔥 40.5k | `communication` |
| [Git](communication/254-git_1469e351/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/git.md) | 🔥 40.5k | `communication` |
| [Docker](communication/255-docker_24d9f59d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/install/docker.md) | 🔥 40.5k | `communication` |
| [Analytics](communication/256-analytics_073a41ba/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/analytics.md) | 🔥 40.5k | `communication` |
| [Edit Errors](communication/257-edit-errors_29cca0f5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/edit-errors.md) | 🔥 40.5k | `communication` |
| [Caching](communication/258-caching_50e5e286/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/caching.md) | 🔥 40.5k | `communication` |
| [Notifications](communication/144-notifications_38a3472f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/notifications.md) | 🔥 40.5k | `communication` |
| [Cost Tracking](communication/259-cost_tracking_b8c7b063/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/COST_TRACKING.md) | ⭐ 43 | `communication` |
| [Debugging](communication/260-debugging_6f2c48f9/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DEBUGGING.md) | ⭐ 43 | `communication` |
| [Trace Spec](communication/261-trace_spec_24289c84/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TRACE_SPEC.md) | ⭐ 43 | `communication` |

### Content Creation (17 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [No Heredoc.Instructions](content-creation/357-no-heredocinstructions_a1bdea6d/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/no-heredoc.instructions.md) | ⭐ 60 | `content creation` |
| [Skill](content-creation/049-name-skill_83d53e39/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/docs-writer/SKILL.md) | ⭐ 60 | `content creation` |
| [2023 10 22 Repomap](content-creation/358-2023-10-22-repomap_89a9b7a9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-10-22-repomap.md) | 🔥 40.5k | `content creation` |
| [2024 05 02 Browser](content-creation/359-2024-05-02-browser_bd617fca/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-02-browser.md) | 🔥 40.5k | `content creation` |
| [2024 05 24 Self Assembly](content-creation/360-2024-05-24-self-assembly_37b50370/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-24-self-assembly.md) | 🔥 40.5k | `content creation` |
| [2024 07 01 Sonnet Not Lazy](content-creation/361-2024-07-01-sonnet-not-lazy_630e8d66/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-07-01-sonnet-not-lazy.md) | 🔥 40.5k | `content creation` |
| [2024 11 21 Quantization](content-creation/362-2024-11-21-quantization_7b1bdc60/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-11-21-quantization.md) | 🔥 40.5k | `content creation` |
| [2024 12 03 Qwq](content-creation/363-2024-12-03-qwq_223f2745/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-12-03-qwq.md) | 🔥 40.5k | `content creation` |
| [2025 01 28 Deepseek Down](content-creation/364-2025-01-28-deepseek-down_b46d65b3/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-01-28-deepseek-down.md) | 🔥 40.5k | `content creation` |
| [2025 05 08 Qwen3](content-creation/365-2025-05-08-qwen3_b8c3e75c/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-05-08-qwen3.md) | 🔥 40.5k | `content creation` |
| [Benchmarks 0125](content-creation/366-benchmarks-0125_c20592e5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/benchmarks-0125.md) | 🔥 40.5k | `content creation` |
| [Ctags](content-creation/367-ctags_9dd05206/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/ctags.md) | 🔥 40.5k | `content creation` |
| [Faq](content-creation/368-faq_af3fbd08/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/faq.md) | 🔥 40.5k | `content creation` |
| [Unified Diffs](content-creation/369-unified-diffs_16b627c7/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/unified-diffs.md) | 🔥 40.5k | `content creation` |
| [Index](content-creation/019-index_2e12c12a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/index.md) | 🔥 40.5k | `content creation` |
| [Browser](content-creation/370-browser_27e56b0f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/browser.md) | 🔥 40.5k | `content creation` |
| [Modes](content-creation/371-modes_77c3d350/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/modes.md) | 🔥 40.5k | `content creation` |

### Development (31 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Scripting](development/2875-scripting_163d09dc/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/scripting.md) | 🔥 40.5k | `development` |
| [Usage](development/2876-usage_8a5d1b8d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage.md) | 🔥 40.5k | `development` |
| [Index](development/468-index_3f9d5a0f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/share/index.md) | 🔥 40.5k | `development` |
| [Adv Model Settings](development/2877-adv-model-settings_45ec2acd/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/adv-model-settings.md) | 🔥 40.5k | `development` |
| [Reasoning](development/2878-reasoning_a05be33d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/reasoning.md) | 🔥 40.5k | `development` |
| [Contrib](development/2879-contrib_40379351/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/contrib.md) | 🔥 40.5k | `development` |
| [Edit](development/2880-edit_92b57853/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/edit.md) | 🔥 40.5k | `development` |
| [Notes](development/2881-notes_aebac3fc/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/notes.md) | 🔥 40.5k | `development` |
| [Refactor](development/827-refactor_1a6d950f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/refactor.md) | 🔥 40.5k | `development` |
| [Privacy](development/2882-privacy_a4277601/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/legal/privacy.md) | 🔥 40.5k | `development` |
| [Anthropic](development/2185-anthropic_4e802ab8/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/anthropic.md) | 🔥 40.5k | `development` |
| [Gemini](development/1371-gemini_48668ea1/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/gemini.md) | 🔥 40.5k | `development` |
| [Github](development/2317-github_bc027382/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/github.md) | 🔥 40.5k | `development` |
| [Ollama](development/1338-ollama_22770bc8/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/ollama.md) | 🔥 40.5k | `development` |
| [Edit Formats](development/2883-edit-formats_0d3c5081/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/edit-formats.md) | 🔥 40.5k | `development` |
| [Auto Accept Architect](development/2884-auto-accept-architect_4d97fe45/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/auto-accept-architect.md) | 🔥 40.5k | `development` |
| [Dont Drop Original Read Files](development/2885-dont-drop-original-read-files_e4e2aebf/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/dont-drop-original-read-files.md) | 🔥 40.5k | `development` |
| [Model Accepts Settings](development/2886-model-accepts-settings_4bf52d6e/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/model-accepts-settings.md) | 🔥 40.5k | `development` |
| [Tree Sitter Language Pack](development/2887-tree-sitter-language-pack_d7fc4c53/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/tree-sitter-language-pack.md) | 🔥 40.5k | `development` |
| [Support](development/2888-support_2a9ff792/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/support.md) | 🔥 40.5k | `development` |
| [Token Limits](development/2889-token-limits_f51ac904/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/token-limits.md) | 🔥 40.5k | `development` |
| [Copypaste](development/2890-copypaste_845f243a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/copypaste.md) | 🔥 40.5k | `development` |
| [Tips](development/2891-tips_b26427c6/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/tips.md) | 🔥 40.5k | `development` |
| [Tutorials](development/1100-tutorials_64a073b5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/tutorials.md) | 🔥 40.5k | `development` |
| [Voice](development/2892-voice_06d5fba4/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/voice.md) | 🔥 40.5k | `development` |
| [Watch](development/2893-watch_395175e2/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/watch.md) | 🔥 40.5k | `development` |
| [Adapters](development/2894-adapters_beef3be8/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/ADAPTERS.md) | ⭐ 43 | `development` |
| [Backend Requirements](development/2895-backend_requirements_51bd5367/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/BACKEND_REQUIREMENTS.md) | ⭐ 43 | `development` |
| [Database Setup](development/2896-database_setup_061e7377/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DATABASE_SETUP.md) | ⭐ 43 | `development` |
| [Framework Support](development/2897-framework_support_cd166e68/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FRAMEWORK_SUPPORT.md) | ⭐ 43 | `development` |
| [Generic Guardrail Api](development/generic_guardrail_api_60da4bb2/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/adding_provider/generic_guardrail_api.md) | 🔥 35.8k | `development` |

### Development/Devops (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Configs](development/devops/configs_40652549/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/configs.md) | 🔥 35.8k | `development` |
| [Prod](development/devops/prod_27e95b4b/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/prod.md) | 🔥 35.8k | `development` |

### Development/Testing (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [2024 04 09 Gpt 4 Turbo](development/testing/082-2024-04-09-gpt-4-turbo_9e2cf0f0/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-04-09-gpt-4-turbo.md) | 🔥 40.5k | `development` |

### Development/Tools (12 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Docs.Instructions](development/tools/229-docsinstructions_c7327425/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/docs.instructions.md) | ⭐ 60 | `development` |
| [2025 01 15 Uv](development/tools/319-2025-01-15-uv_954a0942/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-01-15-uv.md) | 🔥 40.5k | `development` |
| [Aider Conf](development/tools/320-aider_conf_49df9413/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/aider_conf.md) | 🔥 40.5k | `development` |
| [Dotenv](development/tools/321-dotenv_0088c147/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/dotenv.md) | 🔥 40.5k | `development` |
| [Options](development/tools/322-options_6bd20c61/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/options.md) | 🔥 40.5k | `development` |
| [Commands](development/tools/252-commands_57549608/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/commands.md) | 🔥 40.5k | `development` |
| [Agents](development/tools/015-agents_1bbd0b4d/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/AGENTS.md) | ⭐ 43 | `development` |
| [Evaluation Metrics](development/tools/323-evaluation_metrics_0d6d74bb/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/EVALUATION_METRICS.md) | ⭐ 43 | `development` |
| [Faq](development/tools/324-faq_1839ca4b/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FAQ.md) | ⭐ 43 | `development` |
| [Mcp Contracts](development/tools/325-mcp_contracts_53288f40/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/MCP_CONTRACTS.md) | ⭐ 43 | `development` |
| [Tool Categories](development/tools/326-tool_categories_6c5054ba/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TOOL_CATEGORIES.md) | ⭐ 43 | `development` |
| [Troubleshooting](development/tools/205-troubleshooting_271bcac1/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TROUBLESHOOTING.md) | ⭐ 43 | `development` |

### Research (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Infinite Output](research/259-infinite-output_a05bef43/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/infinite-output.md) | 🔥 40.5k | `research` |

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

*Last updated: 2026-02-12 18:21:07 UTC*
*Automatically maintained by SkillFlow*
