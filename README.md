# X-Skills

A curated collection of **196 AI-powered skills** organized into 14 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (4 skills)
- **Automation/Workflow** (3 skills)
- **Commercial** (4 skills)
- **Communication** (22 skills)
- **Content Creation** (25 skills)
- **Daily Assistant** (2 skills)
- **Data Analysis** (10 skills)
- **Development** (80 skills)
- **Development/Devops** (7 skills)
- **Development/Testing** (6 skills)
- **Development/Tools** (18 skills)
- **Other** (12 skills)
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


### Automation/Scripting (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Recording](automation/scripting/recording_762ea469/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/recording.md) | 🔥 40.5k | `automation` |
| [Config](automation/scripting/config_18092253/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config.md) | 🔥 40.5k | `automation` |
| [Troubleshooting](automation/scripting/troubleshooting_51816008/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting.md) | 🔥 40.5k | `automation` |
| [Aider Not Found](automation/scripting/aider-not-found_23d5fb29/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/aider-not-found.md) | 🔥 40.5k | `automation` |

### Automation/Workflow (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Executive Pitch](automation/workflow/executive-pitch_61fed1da/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/executive-pitch.md) | ⭐ 60 | `automation` |
| [Detecting Llm Hallucinations In Ci](automation/workflow/detecting-llm-hallucinations-in-ci_471c775c/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/guides/detecting-llm-hallucinations-in-ci.md) | ⭐ 43 | `automation` |
| [Agent Testing Plan.Prompt](automation/workflow/agent-testing-planprompt_8740a039/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/prompts/agent-testing-plan.prompt.md) | ⭐ 60 | `automation` |

### Commercial (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Terraform Roadmap](commercial/terraform-roadmap_c6b6b671/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/terraform-roadmap.md) | ⭐ 60 | `commercial` |
| [Keys](commercial/keys_b2c6966b/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/keys.md) | 🔥 40.5k | `commercial` |
| [Api Keys](commercial/api-keys_b8e51df6/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/api-keys.md) | 🔥 40.5k | `commercial` |
| [Contributor Agreement](commercial/contributor-agreement_a721a840/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/legal/contributor-agreement.md) | 🔥 40.5k | `commercial` |

### Communication (22 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Help](communication/help_dcae5407/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/help.md) | 🔥 40.5k | `communication` |
| [Model Warnings](communication/model-warnings_ba625717/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/model-warnings.md) | 🔥 40.5k | `communication` |
| [Multi Line](communication/multi-line_5264f142/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/multi-line.md) | 🔥 40.5k | `communication` |
| [Git](communication/git_1469e351/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/git.md) | 🔥 40.5k | `communication` |
| [Docker](communication/docker_24d9f59d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/install/docker.md) | 🔥 40.5k | `communication` |
| [Analytics](communication/analytics_073a41ba/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/analytics.md) | 🔥 40.5k | `communication` |
| [Edit Errors](communication/edit-errors_29cca0f5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/edit-errors.md) | 🔥 40.5k | `communication` |
| [Caching](communication/caching_50e5e286/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/caching.md) | 🔥 40.5k | `communication` |
| [Notifications](communication/notifications_38a3472f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/notifications.md) | 🔥 40.5k | `communication` |
| [Cost Tracking](communication/cost_tracking_b8c7b063/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/COST_TRACKING.md) | ⭐ 43 | `communication` |
| [Debugging](communication/debugging_6f2c48f9/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DEBUGGING.md) | ⭐ 43 | `communication` |
| [Trace Spec](communication/trace_spec_24289c84/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TRACE_SPEC.md) | ⭐ 43 | `communication` |
| [Replit Pipx](communication/replit-pipx_f4ad3a13/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/replit-pipx.md) | 🔥 40.5k | `communication` |
| [Works Best](communication/works-best_6d95348a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/works-best.md) | 🔥 40.5k | `communication` |
| [Chat Transcript Css](communication/chat-transcript-css_737c75a9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/chat-transcript-css.md) | 🔥 40.5k | `communication` |
| [Pong](communication/pong_3e22e4bd/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/pong.md) | 🔥 40.5k | `communication` |
| [Models And Keys](communication/models-and-keys_fa5239de/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/models-and-keys.md) | 🔥 40.5k | `communication` |
| [Images Urls](communication/images-urls_86ce6a9e/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/images-urls.md) | 🔥 40.5k | `communication` |
| [Lint Test](communication/lint-test_b7ac6076/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/lint-test.md) | 🔥 40.5k | `communication` |
| [Langgraph Cloud](communication/langgraph_cloud_acf6f5cb/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/LANGGRAPH_CLOUD.md) | ⭐ 43 | `communication` |
| [Suite Types](communication/suite_types_2f4a15a8/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/SUITE_TYPES.md) | ⭐ 43 | `communication` |
| [Backend Implementations](communication/backend-implementations_e9fd1e24/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/examples/backend-implementations.md) | ⭐ 43 | `communication` |

### Content Creation (25 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [No Heredoc.Instructions](content-creation/no-heredocinstructions_a1bdea6d/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/no-heredoc.instructions.md) | ⭐ 60 | `content creation` |
| [Skill](content-creation/name-skill_83d53e39/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/docs-writer/SKILL.md) | ⭐ 60 | `content creation` |
| [2023 10 22 Repomap](content-creation/2023-10-22-repomap_89a9b7a9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-10-22-repomap.md) | 🔥 40.5k | `content creation` |
| [2024 05 02 Browser](content-creation/2024-05-02-browser_bd617fca/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-02-browser.md) | 🔥 40.5k | `content creation` |
| [2024 05 24 Self Assembly](content-creation/2024-05-24-self-assembly_37b50370/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-24-self-assembly.md) | 🔥 40.5k | `content creation` |
| [2024 07 01 Sonnet Not Lazy](content-creation/2024-07-01-sonnet-not-lazy_630e8d66/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-07-01-sonnet-not-lazy.md) | 🔥 40.5k | `content creation` |
| [2024 11 21 Quantization](content-creation/2024-11-21-quantization_7b1bdc60/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-11-21-quantization.md) | 🔥 40.5k | `content creation` |
| [2024 12 03 Qwq](content-creation/2024-12-03-qwq_223f2745/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-12-03-qwq.md) | 🔥 40.5k | `content creation` |
| [2025 01 28 Deepseek Down](content-creation/2025-01-28-deepseek-down_b46d65b3/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-01-28-deepseek-down.md) | 🔥 40.5k | `content creation` |
| [2025 05 08 Qwen3](content-creation/2025-05-08-qwen3_b8c3e75c/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-05-08-qwen3.md) | 🔥 40.5k | `content creation` |
| [Benchmarks 0125](content-creation/benchmarks-0125_c20592e5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/benchmarks-0125.md) | 🔥 40.5k | `content creation` |
| [Ctags](content-creation/ctags_9dd05206/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/ctags.md) | 🔥 40.5k | `content creation` |
| [Faq](content-creation/faq_af3fbd08/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/faq.md) | 🔥 40.5k | `content creation` |
| [Unified Diffs](content-creation/unified-diffs_16b627c7/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/unified-diffs.md) | 🔥 40.5k | `content creation` |
| [Index](content-creation/index_2e12c12a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/index.md) | 🔥 40.5k | `content creation` |
| [Browser](content-creation/browser_27e56b0f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/browser.md) | 🔥 40.5k | `content creation` |
| [Modes](content-creation/modes_77c3d350/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/modes.md) | 🔥 40.5k | `content creation` |
| [Visual Elements Guide](content-creation/visual-elements-guide_f53bfc70/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/visual-elements-guide.md) | ⭐ 60 | `content creation` |
| [2024 05 22 Draft](content-creation/2024-05-22-draft_adee17ea/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-22-draft.md) | 🔥 40.5k | `content creation` |
| [Repomap](content-creation/repomap_ac41c884/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/repomap.md) | 🔥 40.5k | `content creation` |
| [Census](content-creation/census_7cf5c0e1/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/census.md) | 🔥 40.5k | `content creation` |
| [Css Exercises](content-creation/css-exercises_aac743bc/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/css-exercises.md) | 🔥 40.5k | `content creation` |
| [Hello World Flask](content-creation/hello-world-flask_d80df16b/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/hello-world-flask.md) | 🔥 40.5k | `content creation` |
| [Hello](content-creation/hello_df5b61c4/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/hello.md) | 🔥 40.5k | `content creation` |
| [Skill](content-creation/name-skill_9a2c44cf/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/examples/agent-test/SKILL.md) | ⭐ 43 | `content creation` |

### Daily Assistant (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agents Definitions.Instructions](daily-assistant/agents-definitionsinstructions_4e39deb3/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/agents-definitions.instructions.md) | ⭐ 60 | `daily assistant` |
| [Freshness Checklist](daily-assistant/freshness-checklist_fe681f3c/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/docs-writer/references/freshness-checklist.md) | ⭐ 60 | `daily assistant` |

### Data Analysis (10 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Workflow](data-analysis/workflow_9313af94/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/workflow.md) | ⭐ 60 | `data analysis` |
| [Shell.Instructions](data-analysis/shellinstructions_81e494f0/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/shell.instructions.md) | ⭐ 60 | `data analysis` |
| [2024 05 13 Models Over Time](data-analysis/2024-05-13-models-over-time_a4091a0a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-13-models-over-time.md) | 🔥 40.5k | `data analysis` |
| [2024 08 26 Sonnet Seems Fine](data-analysis/2024-08-26-sonnet-seems-fine_85e22a86/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-08-26-sonnet-seems-fine.md) | 🔥 40.5k | `data analysis` |
| [2024 09 12 O1](data-analysis/2024-09-12-o1_8336c0a4/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-09-12-o1.md) | 🔥 40.5k | `data analysis` |
| [2025 05 07 Gemini Cost](data-analysis/2025-05-07-gemini-cost_4e409d97/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-05-07-gemini-cost.md) | 🔥 40.5k | `data analysis` |
| [Chat Mode](data-analysis/chat_mode_de4dc6e1/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/CHAT_MODE.md) | ⭐ 43 | `data analysis` |
| [Cli Reference](data-analysis/cli_reference_813f9afa/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/CLI_REFERENCE.md) | ⭐ 43 | `data analysis` |
| [Getting Started](data-analysis/getting_started_f07a70f1/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/GETTING_STARTED.md) | ⭐ 43 | `data analysis` |
| [Bug Report](data-analysis/bug_report_31a16204/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/.github/ISSUE_TEMPLATE/bug_report.md) | ⭐ 43 | `data analysis` |

### Development (80 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Scripting](development/scripting_163d09dc/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/scripting.md) | 🔥 40.5k | `development` |
| [Usage](development/usage_8a5d1b8d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage.md) | 🔥 40.5k | `development` |
| [Index](development/index_3f9d5a0f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/share/index.md) | 🔥 40.5k | `development` |
| [Adv Model Settings](development/adv-model-settings_45ec2acd/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/adv-model-settings.md) | 🔥 40.5k | `development` |
| [Reasoning](development/reasoning_a05be33d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/reasoning.md) | 🔥 40.5k | `development` |
| [Contrib](development/contrib_40379351/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/contrib.md) | 🔥 40.5k | `development` |
| [Edit](development/edit_92b57853/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/edit.md) | 🔥 40.5k | `development` |
| [Notes](development/notes_aebac3fc/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/notes.md) | 🔥 40.5k | `development` |
| [Refactor](development/refactor_1a6d950f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/refactor.md) | 🔥 40.5k | `development` |
| [Privacy](development/privacy_a4277601/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/legal/privacy.md) | 🔥 40.5k | `development` |
| [Anthropic](development/anthropic_4e802ab8/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/anthropic.md) | 🔥 40.5k | `development` |
| [Gemini](development/gemini_48668ea1/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/gemini.md) | 🔥 40.5k | `development` |
| [Github](development/github_bc027382/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/github.md) | 🔥 40.5k | `development` |
| [Ollama](development/ollama_22770bc8/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/ollama.md) | 🔥 40.5k | `development` |
| [Edit Formats](development/edit-formats_0d3c5081/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/edit-formats.md) | 🔥 40.5k | `development` |
| [Auto Accept Architect](development/auto-accept-architect_4d97fe45/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/auto-accept-architect.md) | 🔥 40.5k | `development` |
| [Dont Drop Original Read Files](development/dont-drop-original-read-files_e4e2aebf/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/dont-drop-original-read-files.md) | 🔥 40.5k | `development` |
| [Model Accepts Settings](development/model-accepts-settings_4bf52d6e/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/model-accepts-settings.md) | 🔥 40.5k | `development` |
| [Tree Sitter Language Pack](development/tree-sitter-language-pack_d7fc4c53/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/tree-sitter-language-pack.md) | 🔥 40.5k | `development` |
| [Support](development/support_2a9ff792/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/support.md) | 🔥 40.5k | `development` |
| [Token Limits](development/token-limits_f51ac904/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/token-limits.md) | 🔥 40.5k | `development` |
| [Copypaste](development/copypaste_845f243a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/copypaste.md) | 🔥 40.5k | `development` |
| [Tips](development/tips_b26427c6/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/tips.md) | 🔥 40.5k | `development` |
| [Tutorials](development/tutorials_64a073b5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/tutorials.md) | 🔥 40.5k | `development` |
| [Voice](development/voice_06d5fba4/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/voice.md) | 🔥 40.5k | `development` |
| [Watch](development/watch_395175e2/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/watch.md) | 🔥 40.5k | `development` |
| [Adapters](development/adapters_beef3be8/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/ADAPTERS.md) | ⭐ 43 | `development` |
| [Backend Requirements](development/backend_requirements_51bd5367/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/BACKEND_REQUIREMENTS.md) | ⭐ 43 | `development` |
| [Database Setup](development/database_setup_061e7377/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DATABASE_SETUP.md) | ⭐ 43 | `development` |
| [Framework Support](development/framework_support_cd166e68/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FRAMEWORK_SUPPORT.md) | ⭐ 43 | `development` |
| [Code Review.Instructions](development/code-reviewinstructions_0df9a6cd/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/code-review.instructions.md) | ⭐ 60 | `development` |
| [Self Explanatory Code Commenting.Instructions](development/self-explanatory-code-commentinginstructions_9983f412/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/self-explanatory-code-commenting.instructions.md) | ⭐ 60 | `development` |
| [Update Docs On Code Change.Instructions](development/update-docs-on-code-changeinstructions_49383808/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/update-docs-on-code-change.instructions.md) | ⭐ 60 | `development` |
| [Pilot Success Checklist](development/pilot-success-checklist_aa27242f/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/pilot-success-checklist.md) | ⭐ 60 | `development` |
| [Doc Standards](development/doc-standards_a89fc630/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/docs-writer/references/doc-standards.md) | ⭐ 60 | `development` |
| [Chat History](development/chat-history_a9a89dfd/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/tests/fixtures/chat-history.md) | 🔥 40.5k | `development` |
| [Blame](development/blame_adf37d12/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/blame.md) | 🔥 40.5k | `development` |
| [Get Started](development/get-started_23121008/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/get-started.md) | 🔥 40.5k | `development` |
| [2024 03 08 Claude 3](development/2024-03-08-claude-3_d37f7881/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-03-08-claude-3.md) | 🔥 40.5k | `development` |
| [2024 05 22 Linting](development/2024-05-22-linting_3762dc83/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-22-linting.md) | 🔥 40.5k | `development` |
| [2024 05 22 Swe Bench Lite](development/2024-05-22-swe-bench-lite_3d48ceef/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-22-swe-bench-lite.md) | 🔥 40.5k | `development` |
| [2024 06 02 Main Swe Bench](development/2024-06-02-main-swe-bench_924f7b49/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-06-02-main-swe-bench.md) | 🔥 40.5k | `development` |
| [2024 07 25 New Models](development/2024-07-25-new-models_aab5b30e/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-07-25-new-models.md) | 🔥 40.5k | `development` |
| [2024 08 14 Code In Json](development/2024-08-14-code-in-json_1d0bf2c9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-08-14-code-in-json.md) | 🔥 40.5k | `development` |
| [2024 09 26 Architect](development/2024-09-26-architect_39d04275/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-09-26-architect.md) | 🔥 40.5k | `development` |
| [2024 12 21 Polyglot](development/2024-12-21-polyglot_0ecb6801/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-12-21-polyglot.md) | 🔥 40.5k | `development` |
| [2025 01 24 R1 Sonnet](development/2025-01-24-r1-sonnet_f7fc241a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-01-24-r1-sonnet.md) | 🔥 40.5k | `development` |
| [Benchmarks 1106](development/benchmarks-1106_75ed7cba/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/benchmarks-1106.md) | 🔥 40.5k | `development` |
| [Benchmarks Speed 1106](development/benchmarks-speed-1106_376fccf9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/benchmarks-speed-1106.md) | 🔥 40.5k | `development` |
| [Benchmarks](development/benchmarks_9a51e720/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/benchmarks.md) | 🔥 40.5k | `development` |
| [Index](development/index_274a87d1/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/index.md) | 🔥 40.5k | `development` |
| [Install](development/install_e82f091f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/install.md) | 🔥 40.5k | `development` |
| [Languages](development/languages_c0aef960/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/languages.md) | 🔥 40.5k | `development` |
| [Llms](development/llms_7ead152f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms.md) | 🔥 40.5k | `development` |
| [2048 Game](development/2048-game_fd4c9956/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/2048-game.md) | 🔥 40.5k | `development` |
| [Add Test](development/add-test_db69cc4b/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/add-test.md) | 🔥 40.5k | `development` |
| [Asciinema](development/asciinema_eb47bb6d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/asciinema.md) | 🔥 40.5k | `development` |
| [Complex Change](development/complex-change_4b6c008b/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/complex-change.md) | 🔥 40.5k | `development` |
| [No Color](development/no-color_f4aaecc2/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/no-color.md) | 🔥 40.5k | `development` |
| [Semantic Search Replace](development/semantic-search-replace_7afde258/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/semantic-search-replace.md) | 🔥 40.5k | `development` |
| [Update Docs](development/update-docs_ba400159/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/examples/update-docs.md) | 🔥 40.5k | `development` |
| [Editor](development/editor_4a12477d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/editor.md) | 🔥 40.5k | `development` |
| [Codespaces](development/codespaces_6ba394b6/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/install/codespaces.md) | 🔥 40.5k | `development` |
| [Optional](development/optional_b52b739f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/install/optional.md) | 🔥 40.5k | `development` |
| [By Release Date](development/by-release-date_0ca74c85/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/by-release-date.md) | 🔥 40.5k | `development` |
| [Azure](development/azure_c62b37bf/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/azure.md) | 🔥 40.5k | `development` |
| [Cohere](development/cohere_6e3d63ef/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/cohere.md) | 🔥 40.5k | `development` |
| [Deepseek](development/deepseek_53424aa6/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/deepseek.md) | 🔥 40.5k | `development` |
| [Groq](development/groq_0db51e8d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/groq.md) | 🔥 40.5k | `development` |
| [Lm Studio](development/lm-studio_43952268/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/lm-studio.md) | 🔥 40.5k | `development` |
| [Openai Compat](development/openai-compat_ed0c73ac/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/openai-compat.md) | 🔥 40.5k | `development` |
| [Openai](development/openai_7e1580c8/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/openai.md) | 🔥 40.5k | `development` |
| [Openrouter](development/openrouter_8fa09faf/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/openrouter.md) | 🔥 40.5k | `development` |
| [Other](development/other_307adeb4/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/other.md) | 🔥 40.5k | `development` |
| [Vertex](development/vertex_228723a9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/vertex.md) | 🔥 40.5k | `development` |
| [Xai](development/xai_d15274d9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/xai.md) | 🔥 40.5k | `development` |
| [Index](development/index_21ede540/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/index.md) | 🔥 40.5k | `development` |
| [Conventions](development/conventions_b584901e/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/conventions.md) | 🔥 40.5k | `development` |
| [Quickstart Huggingface](development/quickstart_huggingface_5542d967/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/QUICKSTART_HUGGINGFACE.md) | ⭐ 43 | `development` |
| [Skill](development/name-skill_057fe041/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/examples/skills/test-skill/SKILL.md) | ⭐ 43 | `development` |

### Development/Devops (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Quickstart](development/devops/quickstart_16c8c0f8/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/quickstart.md) | ⭐ 60 | `development` |
| [Troubleshooting](development/devops/troubleshooting_3a0c559b/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/troubleshooting.md) | ⭐ 60 | `development` |
| [Repo Architecture](development/devops/repo-architecture_a17fa352/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/docs-writer/references/repo-architecture.md) | ⭐ 60 | `development` |
| [Not Code](development/devops/not-code_b3322beb/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/not-code.md) | 🔥 40.5k | `development` |
| [Ci Cd](development/devops/ci_cd_be43cfca/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/CI_CD.md) | ⭐ 43 | `development` |
| [Golden Traces](development/devops/golden_traces_36206255/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/GOLDEN_TRACES.md) | ⭐ 43 | `development` |
| [Quickstart Langgraph](development/devops/quickstart_langgraph_91ad62a4/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/QUICKSTART_LANGGRAPH.md) | ⭐ 43 | `development` |

### Development/Testing (6 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [2024 04 09 Gpt 4 Turbo](development/testing/2024-04-09-gpt-4-turbo_9e2cf0f0/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-04-09-gpt-4-turbo.md) | 🔥 40.5k | `development` |
| [Model Aliases](development/testing/model-aliases_f8bc573a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/model-aliases.md) | 🔥 40.5k | `development` |
| [Imports](development/testing/imports_eb889c49/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/imports.md) | 🔥 40.5k | `development` |
| [Pytest For Ai Agents Langgraph Ci](development/testing/pytest-for-ai-agents-langgraph-ci_5b466c0c/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/guides/pytest-for-ai-agents-langgraph-ci.md) | ⭐ 43 | `development` |
| [Adapter Request](development/testing/adapter_request_b9c7014a/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/.github/ISSUE_TEMPLATE/adapter_request.md) | ⭐ 43 | `development` |
| [Feature Request](development/testing/feature_request_b63a05b6/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/.github/ISSUE_TEMPLATE/feature_request.md) | ⭐ 43 | `development` |

### Development/Tools (18 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Docs.Instructions](development/tools/docsinstructions_c7327425/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/docs.instructions.md) | ⭐ 60 | `development` |
| [2025 01 15 Uv](development/tools/2025-01-15-uv_954a0942/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-01-15-uv.md) | 🔥 40.5k | `development` |
| [Aider Conf](development/tools/aider_conf_49df9413/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/aider_conf.md) | 🔥 40.5k | `development` |
| [Dotenv](development/tools/dotenv_0088c147/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/dotenv.md) | 🔥 40.5k | `development` |
| [Options](development/tools/options_6bd20c61/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/options.md) | 🔥 40.5k | `development` |
| [Commands](development/tools/commands_57549608/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/commands.md) | 🔥 40.5k | `development` |
| [Agents](development/tools/agents_1bbd0b4d/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/AGENTS.md) | ⭐ 43 | `development` |
| [Evaluation Metrics](development/tools/evaluation_metrics_0d6d74bb/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/EVALUATION_METRICS.md) | ⭐ 43 | `development` |
| [Faq](development/tools/faq_1839ca4b/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FAQ.md) | ⭐ 43 | `development` |
| [Mcp Contracts](development/tools/mcp_contracts_53288f40/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/MCP_CONTRACTS.md) | ⭐ 43 | `development` |
| [Tool Categories](development/tools/tool_categories_6c5054ba/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TOOL_CATEGORIES.md) | ⭐ 43 | `development` |
| [Troubleshooting](development/tools/troubleshooting_271bcac1/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TROUBLESHOOTING.md) | ⭐ 43 | `development` |
| [Bedrock](development/tools/bedrock_982a37fd/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/bedrock.md) | 🔥 40.5k | `development` |
| [Behavior Coverage](development/tools/behavior_coverage_518953c0/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/BEHAVIOR_COVERAGE.md) | ⭐ 43 | `development` |
| [Setup Langgraph Example](development/tools/setup_langgraph_example_f1952c46/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/SETUP_LANGGRAPH_EXAMPLE.md) | ⭐ 43 | `development` |
| [Skills Testing](development/tools/skills_testing_55b36405/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/SKILLS_TESTING.md) | ⭐ 43 | `development` |
| [Statistical Mode](development/tools/statistical_mode_5e5af8fc/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/STATISTICAL_MODE.md) | ⭐ 43 | `development` |
| [Test Generation](development/tools/test_generation_4cbb5efd/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TEST_GENERATION.md) | ⭐ 43 | `development` |

### Other (12 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Help Tip](other/help-tip_629546a2/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/help-tip.md) | 🔥 40.5k | `other` |
| [Install](other/install_53852552/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/install.md) | 🔥 40.5k | `other` |
| [Python M Aider](other/python-m-aider_586906d0/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/python-m-aider.md) | 🔥 40.5k | `other` |
| [2023 05 25 Ctags](other/2023-05-25-ctags_1790f6d2/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-05-25-ctags.md) | 🔥 40.5k | `other` |
| [2023 07 02 Benchmarks](other/2023-07-02-benchmarks_872110b7/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-07-02-benchmarks.md) | 🔥 40.5k | `other` |
| [2023 11 06 Benchmarks 1106](other/2023-11-06-benchmarks-1106_1c99353d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-11-06-benchmarks-1106.md) | 🔥 40.5k | `other` |
| [2023 12 21 Unified Diffs](other/2023-12-21-unified-diffs_824af77a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-12-21-unified-diffs.md) | 🔥 40.5k | `other` |
| [2024 01 25 Benchmarks 0125](other/2024-01-25-benchmarks-0125_7706f742/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-01-25-benchmarks-0125.md) | 🔥 40.5k | `other` |
| [More Info](other/more-info_acd7a814/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more-info.md) | 🔥 40.5k | `other` |
| [Replit](other/replit_6323a4b2/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/install/replit.md) | 🔥 40.5k | `other` |
| [Warnings](other/warnings_15586043/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/warnings.md) | 🔥 40.5k | `other` |
| [Warnings](other/warnings_7877f069/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/warnings.md) | 🔥 40.5k | `other` |

### Productivity (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [2023 11 06 Benchmarks Speed 1106](productivity/2023-11-06-benchmarks-speed-1106_8e820ec7/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-11-06-benchmarks-speed-1106.md) | 🔥 40.5k | `productivity` |

### Research (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Infinite Output](research/infinite-output_a05bef43/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/infinite-output.md) | 🔥 40.5k | `research` |
| [Yaml Schema](research/yaml_schema_e7e77e14/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/YAML_SCHEMA.md) | ⭐ 43 | `research` |

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

*Last updated: 2026-02-12 17:21:22 UTC*
*Automatically maintained by SkillFlow*
