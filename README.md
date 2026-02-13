# X-Skills

A curated collection of **55 AI-powered skills** organized into 10 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Workflow** (5 skills)
- **Commercial** (4 skills)
- **Communication** (2 skills)
- **Content Creation** (2 skills)
- **Development** (11 skills)
- **Development/Devops** (11 skills)
- **Development/Testing** (1 skill)
- **Development/Tools** (17 skills)
- **Other** (1 skill)
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


### Automation/Workflow (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/name-skill_768717c3/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Created](automation/workflow/readme_flat_skills_created_518f2683/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_CREATED.md) | 🔥 23.5k | `automation` |
| [Skill](automation/workflow/002-name-skill_82d9d209/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui-integrator/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Az](automation/workflow/136-readme_flat_skills_az_1dd4094b/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_AZ.md) | 🔥 23.5k | `automation` |
| [Readme Flat Skills Updated](automation/workflow/readme_flat_skills_updated_60c4bdd3/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_UPDATED.md) | 🔥 23.5k | `automation` |

### Commercial (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/name-skill_f3e89c56/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/session-handoff/SKILL.md) | ⭐ 15 | `commercial` |
| [Interface Design](commercial/374-interface-design_764c5ff0/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/interface-design.md) | ⭐ 15 | `commercial` |
| [Dev Endpoint](commercial/dev-endpoint_ad8342d7/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/agent-telemetry/references/dev-endpoint.md) | ⭐ 15 | `commercial` |
| [Tests](commercial/tests_e12182f5/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/tests.md) | ⭐ 15 | `commercial` |

### Communication (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Troubleshoot](communication/221-troubleshoot_7a8329e7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/troubleshoot.md) | 🔥 35.9k | `communication` |
| [Web Search](communication/206-web_search_6a939ba7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/completion/web_search.md) | 🔥 35.9k | `communication` |

### Content Creation (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Openai](content-creation/359-openai_4ef0cd71/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/providers/openai.md) | 🔥 35.9k | `content creation` |
| [Components Catalog](content-creation/components-catalog_ea008930/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui-integrator/references/components-catalog.md) | ⭐ 15 | `content creation` |

### Development (11 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/1178-name-skill_1d0a1439/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/SKILL.md) | ⭐ 15 | `development` |
| [Claude Code Beta Headers](development/2894-claude_code_beta_headers_8fd17b7e/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/tutorials/claude_code_beta_headers.md) | 🔥 35.9k | `development` |
| [Refactoring](development/2896-refactoring_54fced22/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/refactoring.md) | ⭐ 15 | `development` |
| [Responses Api](development/2895-responses_api_ff9997ac/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/providers/openai/responses_api.md) | 🔥 35.9k | `development` |
| [Readme Flat Commands Az](development/785-readme_flat_commands_az_b6610c40/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Skills Releases](development/799-readme_flat_skills_releases_209e0922/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat Claude Md Az](development/readme_flat_claude-md_az_3245cf3a/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Claude Md Created](development/readme_flat_claude-md_created_68a43aec/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Claude Md Updated](development/readme_flat_claude-md_updated_755d5596/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Commands Created](development/readme_flat_commands_created_906a7634/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Commands Updated](development/readme_flat_commands_updated_62eafa02/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_UPDATED.md) | 🔥 23.5k | `development` |

### Development/Devops (11 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Readme Awesome](development/devops/readme_awesome_6e232cc0/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_AWESOME.md) | 🔥 23.5k | `development` |
| [Readme Classic](development/devops/readme_classic_ecb74b82/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_CLASSIC.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Az](development/devops/161-readme_flat_tooling_az_ed0b8064/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Releases](development/devops/163-readme_flat_tooling_releases_8847650e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat All Az](development/devops/readme_flat_all_az_26e4ead3/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat All Created](development/devops/readme_flat_all_created_2abaf022/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat All Releases](development/devops/readme_flat_all_releases_91433d1e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat All Updated](development/devops/readme_flat_all_updated_1b3f5678/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Created](development/devops/readme_flat_tooling_created_c17d7a5f/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Tooling Updated](development/devops/readme_flat_tooling_updated_91ee5233/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Extra](development/devops/readme_extra_e82069ae/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_EXTRA.md) | 🔥 23.5k | `development` |

### Development/Testing (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Mocking](development/testing/mocking_3ceb807f/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/mocking.md) | ⭐ 15 | `development` |

### Development/Tools (17 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Mcp Troubleshoot](development/tools/324-mcp_troubleshoot_6a5a0c41/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/mcp_troubleshoot.md) | 🔥 35.9k | `development` |
| [Skill](development/tools/name-skill_f3f18b25/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/agent-telemetry/SKILL.md) | ⭐ 15 | `development` |
| [Readme Flat Claude Md Releases](development/tools/readme_flat_claude-md_releases_6e122444/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat Clients Az](development/tools/readme_flat_clients_az_074d13be/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLIENTS_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Clients Created](development/tools/readme_flat_clients_created_2ff63d52/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLIENTS_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Clients Releases](development/tools/readme_flat_clients_releases_185dc808/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLIENTS_RELEASES.md) | 🔥 23.5k | `development` |
| [Readme Flat Clients Updated](development/tools/readme_flat_clients_updated_911ccf3f/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLIENTS_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Docs Created](development/tools/readme_flat_docs_created_f8a108ff/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_DOCS_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Docs Updated](development/tools/readme_flat_docs_updated_92ed30d1/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_DOCS_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Hooks Created](development/tools/readme_flat_hooks_created_7d842c39/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_HOOKS_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Hooks Updated](development/tools/readme_flat_hooks_updated_fd3e65bc/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_HOOKS_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Statusline Created](development/tools/readme_flat_statusline_created_b8bf5a38/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STATUSLINE_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Statusline Updated](development/tools/readme_flat_statusline_updated_147b3674/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STATUSLINE_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Styles Created](development/tools/readme_flat_styles_created_9443e240/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STYLES_CREATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Styles Updated](development/tools/readme_flat_styles_updated_e4da8065/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STYLES_UPDATED.md) | 🔥 23.5k | `development` |
| [Readme Flat Workflows Az](development/tools/readme_flat_workflows_az_d2ecb9e3/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_WORKFLOWS_AZ.md) | 🔥 23.5k | `development` |
| [Readme Flat Workflows Releases](development/tools/readme_flat_workflows_releases_313c0f9e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_WORKFLOWS_RELEASES.md) | 🔥 23.5k | `development` |

### Other (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Deep Modules](other/deep-modules_f2123700/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/deep-modules.md) | ⭐ 15 | `other` |

### Research (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Recipes](research/recipes_358c34c7/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui-integrator/references/recipes.md) | ⭐ 15 | `research` |

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

*Last updated: 2026-02-13 23:21:38 UTC*
*Automatically maintained by SkillFlow*
