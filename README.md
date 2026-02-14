# X-Skills

A curated collection of **58 AI-powered skills** organized into 10 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Workflow** (5 skills)
- **Commercial** (3 skills)
- **Communication** (3 skills)
- **Content Creation** (3 skills)
- **Daily Assistant** (1 skill)
- **Data Analysis** (4 skills)
- **Development** (23 skills)
- **Development/Devops** (11 skills)
- **Development/Testing** (2 skills)
- **Development/Tools** (3 skills)

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


### Automation/Workflow (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/002-name-skill_768717c3/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Created](automation/workflow/140-readme_flat_skills_created_518f2683/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_CREATED.md) | 🔥 23.5k | `automation` |
| [Skill](automation/workflow/002-name-skill_82d9d209/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui-integrator/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Az](automation/workflow/136-readme_flat_skills_az_1dd4094b/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_AZ.md) | 🔥 23.5k | `automation` |
| [Overview](automation/workflow/016-overview_cf2b1aeb/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/docs/ai-agents/overview.md) | ⭐ 17 | `automation` |

### Commercial (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/210-name-skill_f3e89c56/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/session-handoff/SKILL.md) | ⭐ 15 | `commercial` |
| [Interface Design](commercial/374-interface-design_764c5ff0/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/interface-design.md) | ⭐ 15 | `commercial` |
| [Dev Endpoint](commercial/375-dev-endpoint_ad8342d7/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/agent-telemetry/references/dev-endpoint.md) | ⭐ 15 | `commercial` |

### Communication (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Troubleshoot](communication/221-troubleshoot_7a8329e7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/troubleshoot.md) | 🔥 35.9k | `communication` |
| [Web Search](communication/206-web_search_6a939ba7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/completion/web_search.md) | 🔥 35.9k | `communication` |
| [Claude](communication/024-claude_8fbab90f/) | [pocketpaw/pocketpaw](https://raw.githubusercontent.com/pocketpaw/pocketpaw/main/docs/CLAUDE.md) | ⭐ 31 | `communication` |

### Content Creation (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Openai](content-creation/359-openai_4ef0cd71/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/providers/openai.md) | 🔥 35.9k | `content creation` |
| [Skill](content-creation/049-name-skill_a9cca04f/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/daiv/automation/agent/skills/skill-creator/SKILL.md) | ⭐ 17 | `content creation` |
| [Target User Agentic Harness V0](content-creation/375-target_user_agentic_harness_v0_a84a827c/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/prds/target_user_agentic_harness_v0.md) | ⭐ 51 | `content creation` |

### Daily Assistant (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Index](daily-assistant/052-index_21c18601/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/index.md) | ⭐ 10 | `daily assistant` |

### Data Analysis (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Architecture](data-analysis/009-architecture_f62cb9de/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/architecture.md) | ⭐ 12 | `data analysis` |
| [Configuration](data-analysis/046-configuration_cda8812a/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/configuration.md) | ⭐ 12 | `data analysis` |
| [Transformation](data-analysis/481-transformation_39e9ae41/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/input_data/transformation.md) | ⭐ 12 | `data analysis` |
| [Index](data-analysis/113-index_152f012b/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/concepts/index.md) | ⭐ 10 | `data analysis` |

### Development (23 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/1178-name-skill_1d0a1439/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/SKILL.md) | ⭐ 15 | `development` |
| [Claude Code Beta Headers](development/2894-claude_code_beta_headers_8fd17b7e/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/tutorials/claude_code_beta_headers.md) | 🔥 35.9k | `development` |
| [Refactoring](development/2896-refactoring_54fced22/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/refactoring.md) | ⭐ 15 | `development` |
| [Responses Api](development/2895-responses_api_ff9997ac/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/providers/openai/responses_api.md) | 🔥 35.9k | `development` |
| [Readme Flat Commands Az](development/785-readme_flat_commands_az_b6610c40/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Skills Releases](development/799-readme_flat_skills_releases_209e0922/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat Claude Md Az](development/777-readme_flat_claude-md_az_3245cf3a/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Claude Md Created](development/778-readme_flat_claude-md_created_68a43aec/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Claude Md Updated](development/780-readme_flat_claude-md_updated_755d5596/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Commands Created](development/786-readme_flat_commands_created_906a7634/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Commands Updated](development/788-readme_flat_commands_updated_62eafa02/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_UPDATED.md) | 🔥 23.5k | `development` |
| [Skills](development/1285-agent-skill_cbde8500/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/docs/ai-agents/skills.md) | ⭐ 17 | `development` |
| [Env Config](development/1286-env-config_86fc6fe4/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/docs/configuration/env-config.md) | ⭐ 17 | `development` |
| [Yaml Config](development/1287-yaml-config_e44080e8/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/docs/configuration/yaml-config.md) | ⭐ 17 | `development` |
| [Skill](development/1178-name-skill_a4acc08a/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/daiv/automation/agent/skills/init/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1178-name-skill_32981586/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/daiv/automation/agent/skills/security-audit/SKILL.md) | ⭐ 17 | `development` |
| [System Prompt Debug Tool B01Af939.Plan](development/2922-system_prompt_debug_tool_b01af939plan_e409ef2e/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/.cursor/plans/system_prompt_debug_tool_b01af939.plan.md) | ⭐ 12 | `development` |
| [Claude](development/140-claude_2b80b6b3/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/CLAUDE.md) | ⭐ 10 | `development` |
| [Getting Help](development/2897-getting-help_29533856/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/getting-help.md) | ⭐ 10 | `development` |
| [Installation](development/474-installation_58370c32/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/installation.md) | ⭐ 10 | `development` |
| [Conditional Middleware](development/2898-conditional-middleware_5b175205/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/conditional-middleware.md) | ⭐ 10 | `development` |
| [Pipeline Spec](development/2899-pipeline-spec_a1346c91/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/pipeline-spec.md) | ⭐ 10 | `development` |
| [Install](development/install_41fb8cee/) | [frmoretto/hardstop](https://raw.githubusercontent.com/frmoretto/hardstop/main/install.md) | ⭐ 15 | `development` |

### Development/Devops (11 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Readme Awesome](development/devops/154-readme_awesome_6e232cc0/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_AWESOME.md) | 🔥 23.5k | `development` |
| [Readme Classic](development/devops/155-readme_classic_ecb74b82/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_CLASSIC.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Az](development/devops/161-readme_flat_tooling_az_ed0b8064/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Releases](development/devops/163-readme_flat_tooling_releases_8847650e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat All Az](development/devops/157-readme_flat_all_az_26e4ead3/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat All Created](development/devops/158-readme_flat_all_created_2abaf022/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat All Releases](development/devops/159-readme_flat_all_releases_91433d1e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat All Updated](development/devops/160-readme_flat_all_updated_1b3f5678/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Created](development/devops/162-readme_flat_tooling_created_c17d7a5f/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Updated](development/devops/164-readme_flat_tooling_updated_91ee5233/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_UPDATED.md) | 🔥 23.5k | `development` |
| [Pipeline Spec](development/devops/367-pipeline_spec_60683b0e/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/api/pipeline_spec.md) | ⭐ 10 | `development` |

### Development/Testing (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Middleware Chains](development/testing/085-middleware-chains_c239a6f7/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/middleware-chains.md) | ⭐ 10 | `development` |
| [Parallel Execution](development/testing/086-parallel-execution_5d3062a6/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/parallel-execution.md) | ⭐ 10 | `development` |

### Development/Tools (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Mcp Troubleshoot](development/tools/324-mcp_troubleshoot_6a5a0c41/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/mcp_troubleshoot.md) | 🔥 35.9k | `development` |
| [Desktop](development/tools/328-desktop_f1b0efd6/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/desktop.md) | ⭐ 51 | `development` |
| [Chain](development/tools/325-chain_02e6ae60/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/api/chain.md) | ⭐ 10 | `development` |

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

*Last updated: 2026-02-14 01:23:20 UTC*
*Automatically maintained by SkillFlow*
