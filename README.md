# X-Skills

A curated collection of **181 AI-powered skills** organized into 10 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Workflow** (9 skills)
- **Commercial** (8 skills)
- **Communication** (5 skills)
- **Content Creation** (14 skills)
- **Daily Assistant** (4 skills)
- **Data Analysis** (7 skills)
- **Development** (98 skills)
- **Development/Devops** (25 skills)
- **Development/Testing** (3 skills)
- **Development/Tools** (8 skills)

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


### Automation/Workflow (9 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/002-name-skill_768717c3/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Created](automation/workflow/140-readme_flat_skills_created_518f2683/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_CREATED.md) | 🔥 23.5k | `automation` |
| [Skill](automation/workflow/002-name-skill_82d9d209/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui-integrator/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Az](automation/workflow/136-readme_flat_skills_az_1dd4094b/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_AZ.md) | 🔥 23.5k | `automation` |
| [Overview](automation/workflow/016-overview_cf2b1aeb/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/docs/ai-agents/overview.md) | ⭐ 17 | `automation` |
| [Security Scorecard](automation/workflow/035-security_scorecard_a19e5df8/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/.github/SECURITY_SCORECARD.md) | ⭐ 4.0k | `automation` |
| [Knowledge Architecture V5.2](automation/workflow/knowledge-architecture-v52_d2959d9f/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/KNOWLEDGE-ARCHITECTURE-v5.2.md) | ⭐ 151 | `automation` |
| [Control Set 07 Logging](automation/workflow/control-set-07-logging_4409d6b1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-07-logging.md) | ⭐ 151 | `automation` |
| [Scan Headers](automation/workflow/scan-headers_7108340d/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/security/penetration-tester/commands/scan-headers.md) | ⭐ 1.3k | `automation` |

### Commercial (8 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/210-name-skill_f3e89c56/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/session-handoff/SKILL.md) | ⭐ 15 | `commercial` |
| [Interface Design](commercial/374-interface-design_764c5ff0/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/interface-design.md) | ⭐ 15 | `commercial` |
| [Dev Endpoint](commercial/375-dev-endpoint_ad8342d7/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/agent-telemetry/references/dev-endpoint.md) | ⭐ 15 | `commercial` |
| [Control Set 04 Output Encoding](commercial/control-set-04-output-encoding_d2c0a426/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-04-output-encoding.md) | ⭐ 151 | `commercial` |
| [Control Set 10 Data Protection](commercial/control-set-10-data-protection_cbab8bd5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-10-data-protection.md) | ⭐ 151 | `commercial` |
| [Control Set Ext 01 02 Auth Patterns](commercial/control-set-ext-01_02-auth-patterns_20aa6001/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-01_02-auth-patterns.md) | ⭐ 151 | `commercial` |
| [Reference Set 05 Third Party Javascript](commercial/reference-set-05-third-party-javascript_31a81d77/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-third-party-javascript.md) | ⭐ 151 | `commercial` |
| [Reference Set 09 Rest Security](commercial/reference-set-09-rest-security_57c8774a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-rest-security.md) | ⭐ 151 | `commercial` |

### Communication (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Troubleshoot](communication/221-troubleshoot_7a8329e7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/troubleshoot.md) | 🔥 35.9k | `communication` |
| [Web Search](communication/206-web_search_6a939ba7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/completion/web_search.md) | 🔥 35.9k | `communication` |
| [Claude](communication/024-claude_8fbab90f/) | [pocketpaw/pocketpaw](https://raw.githubusercontent.com/pocketpaw/pocketpaw/main/docs/CLAUDE.md) | ⭐ 31 | `communication` |
| [Control Set 01 Authentication](communication/control-set-01-authentication_315e4b3c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-01-authentication.md) | ⭐ 151 | `communication` |
| [Reference Set 02 Idor Prevention](communication/reference-set-02-idor-prevention_344c8152/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-idor-prevention.md) | ⭐ 151 | `communication` |

### Content Creation (14 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Openai](content-creation/359-openai_4ef0cd71/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/providers/openai.md) | 🔥 35.9k | `content creation` |
| [Skill](content-creation/049-name-skill_a9cca04f/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/daiv/automation/agent/skills/skill-creator/SKILL.md) | ⭐ 17 | `content creation` |
| [Target User Agentic Harness V0](content-creation/375-target_user_agentic_harness_v0_a84a827c/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/prds/target_user_agentic_harness_v0.md) | ⭐ 51 | `content creation` |
| [Claude](content-creation/007-claude_a0274342/) | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/CLAUDE.md) | 🔥 20.1k | `content creation` |
| [Reference Set 01 Credential Stuffing Prevention](content-creation/reference-set-01-credential-stuffing-prevention_06dc56b3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-credential-stuffing-prevention.md) | ⭐ 151 | `content creation` |
| [Reference Set 01 Saml Security](content-creation/reference-set-01-saml-security_5bf46ff0/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-saml-security.md) | ⭐ 151 | `content creation` |
| [Reference Set 01 Security Questions](content-creation/reference-set-01-security-questions_b5b13f0b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-security-questions.md) | ⭐ 151 | `content creation` |
| [Reference Set 03 Bean Validation](content-creation/reference-set-03-bean-validation_cb7cb648/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-bean-validation.md) | ⭐ 151 | `content creation` |
| [Reference Set 05 Csp](content-creation/reference-set-05-csp_a185fb6c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-csp.md) | ⭐ 151 | `content creation` |
| [Reference Set 05 Xss Filter Evasion](content-creation/reference-set-05-xss-filter-evasion_ba27d20a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-xss-filter-evasion.md) | ⭐ 151 | `content creation` |
| [Reference Set 07 Logging](content-creation/reference-set-07-logging_b8ad9e2b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-07-logging.md) | ⭐ 151 | `content creation` |
| [Reference Set Ext 11 Docker Security](content-creation/reference-set-ext-11-docker-security_0ecaf542/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-11-docker-security.md) | ⭐ 151 | `content creation` |
| [Reference Set Ext 13 Ai Agent Security](content-creation/reference-set-ext-13-ai-agent-security_43ee4a42/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-13-ai-agent-security.md) | ⭐ 151 | `content creation` |
| [Reference Set Ext 11 Iac Security](content-creation/reference-set-ext-11-iac-security_dbcc6cbf/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-11-iac-security.md) | ⭐ 151 | `content creation` |

### Daily Assistant (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Index](daily-assistant/052-index_21c18601/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/index.md) | ⭐ 10 | `daily assistant` |
| [Workflow](daily-assistant/workflow_097c5b9a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/WORKFLOW.md) | ⭐ 151 | `daily assistant` |
| [P3 Trust Boundary](daily-assistant/p3-trust-boundary_69377f71/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P3-TRUST-BOUNDARY.md) | ⭐ 151 | `daily assistant` |
| [P5 Stride Analysis](daily-assistant/p5-stride-analysis_0687bde5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P5-STRIDE-ANALYSIS.md) | ⭐ 151 | `daily assistant` |

### Data Analysis (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Architecture](data-analysis/009-architecture_f62cb9de/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/architecture.md) | ⭐ 12 | `data analysis` |
| [Configuration](data-analysis/046-configuration_cda8812a/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/configuration.md) | ⭐ 12 | `data analysis` |
| [Transformation](data-analysis/481-transformation_39e9ae41/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/input_data/transformation.md) | ⭐ 12 | `data analysis` |
| [Index](data-analysis/113-index_152f012b/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/concepts/index.md) | ⭐ 10 | `data analysis` |
| [P2 Dfd Analysis](data-analysis/p2-dfd-analysis_3d231074/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P2-DFD-ANALYSIS.md) | ⭐ 151 | `data analysis` |
| [P8 Report Generation](data-analysis/p8-report-generation_bd398302/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P8-REPORT-GENERATION.md) | ⭐ 151 | `data analysis` |
| [Control Set Ext 13 Ai Llm](data-analysis/control-set-ext-13-ai-llm_9ed85b8d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-13-ai-llm.md) | ⭐ 151 | `data analysis` |

### Development (98 skills)

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
| [Claude](development/140-claude_01ceef5e/) | [plastic-labs/honcho](https://raw.githubusercontent.com/plastic-labs/honcho/main/CLAUDE.md) | ⭐ 349 | `development` |
| [Skill Architecture Design Cn](development/skill-architecture-design-cn_1e77c8e5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/SKILL-ARCHITECTURE-DESIGN-cn.md) | ⭐ 151 | `development` |
| [Skill Architecture Design](development/skill-architecture-design_caff0b9a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/SKILL-ARCHITECTURE-DESIGN.md) | ⭐ 151 | `development` |
| [Skillset Threat Modeling Tour Cn V5](development/skillset-threat-modeling-tour-cn-v5_cadc4f29/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/Skillset-threat-modeling-tour-cn-v5.md) | ⭐ 151 | `development` |
| [Skillset Threat Modeling Tour V5](development/skillset-threat-modeling-tour-v5_3a0190a3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/Skillset-threat-modeling-tour-v5.md) | ⭐ 151 | `development` |
| [P1 Project Understanding](development/p1-project-understanding_3463d1f5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P1-PROJECT-UNDERSTANDING.md) | ⭐ 151 | `development` |
| [P6 Risk Validation](development/p6-risk-validation_68e1c938/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P6-RISK-VALIDATION.md) | ⭐ 151 | `development` |
| [P7 Mitigation Planning](development/p7-mitigation-planning_f47272a6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P7-MITIGATION-PLANNING.md) | ⭐ 151 | `development` |
| [P8R Detailed Report](development/p8r-detailed-report_079f4af8/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P8R-DETAILED-REPORT.md) | ⭐ 151 | `development` |
| [Control Set 02 Authorization](development/control-set-02-authorization_d1c8fa3a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-02-authorization.md) | ⭐ 151 | `development` |
| [Control Set 03 Input Validation](development/control-set-03-input-validation_8f6d69d1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-03-input-validation.md) | ⭐ 151 | `development` |
| [Control Set 05 Client Side](development/control-set-05-client-side_44068fc6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-05-client-side.md) | ⭐ 151 | `development` |
| [Control Set 06 Cryptography](development/control-set-06-cryptography_ebb7424d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-06-cryptography.md) | ⭐ 151 | `development` |
| [Control Set 08 Error Handling](development/control-set-08-error-handling_8daca94b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-08-error-handling.md) | ⭐ 151 | `development` |
| [Control Set Ext 10 Hardcoded Credentials](development/control-set-ext-10-hardcoded-credentials_5b8d87dc/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-10-hardcoded-credentials.md) | ⭐ 151 | `development` |
| [Control Set Ext 14 Mobile](development/control-set-ext-14-mobile_56ce0614/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-14-mobile.md) | ⭐ 151 | `development` |
| [Control Set Ext 16 Agentic](development/control-set-ext-16-agentic_6f05df72/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-16-agentic.md) | ⭐ 151 | `development` |
| [Reference Set 01 Authentication](development/reference-set-01-authentication_75581456/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-authentication.md) | ⭐ 151 | `development` |
| [Reference Set 01 Cookie Theft Mitigation](development/reference-set-01-cookie-theft-mitigation_017a4462/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-cookie-theft-mitigation.md) | ⭐ 151 | `development` |
| [Reference Set 01 Forgot Password](development/reference-set-01-forgot-password_ac32c6b6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-forgot-password.md) | ⭐ 151 | `development` |
| [Reference Set 01 Jaas](development/reference-set-01-jaas_2a444833/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-jaas.md) | ⭐ 151 | `development` |
| [Reference Set 01 Jwt Java](development/reference-set-01-jwt-java_e46983c1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-jwt-java.md) | ⭐ 151 | `development` |
| [Reference Set 01 Multifactor Authentication](development/reference-set-01-multifactor-authentication_196f4fda/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-multifactor-authentication.md) | ⭐ 151 | `development` |
| [Reference Set 01 Password Storage](development/reference-set-01-password-storage_bba43c08/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-password-storage.md) | ⭐ 151 | `development` |
| [Reference Set 01 Session Management](development/reference-set-01-session-management_7aaafaee/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-session-management.md) | ⭐ 151 | `development` |
| [Reference Set 02 Access Control](development/reference-set-02-access-control_8dd17c5f/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-access-control.md) | ⭐ 151 | `development` |
| [Reference Set 02 Authorization](development/reference-set-02-authorization_dde4a0c3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-authorization.md) | ⭐ 151 | `development` |
| [Reference Set 02 Transaction Authorization](development/reference-set-02-transaction-authorization_350752a3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-transaction-authorization.md) | ⭐ 151 | `development` |
| [Reference Set 03 Deserialization](development/reference-set-03-deserialization_a0c3c75c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-deserialization.md) | ⭐ 151 | `development` |
| [Reference Set 03 File Upload](development/reference-set-03-file-upload_dde7e823/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-file-upload.md) | ⭐ 151 | `development` |
| [Reference Set 03 Injection Prevention Java](development/reference-set-03-injection-prevention-java_0247ca9c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-injection-prevention-java.md) | ⭐ 151 | `development` |
| [Reference Set 03 Injection Prevention](development/reference-set-03-injection-prevention_9a702874/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-injection-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Input Validation](development/reference-set-03-input-validation_b602fff0/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-input-validation.md) | ⭐ 151 | `development` |
| [Reference Set 03 Ldap Injection Prevention](development/reference-set-03-ldap-injection-prevention_1ae008ac/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-ldap-injection-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Mass Assignment](development/reference-set-03-mass-assignment_350e41f4/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-mass-assignment.md) | ⭐ 151 | `development` |
| [Reference Set 03 Os Command Injection Defense](development/reference-set-03-os-command-injection-defense_1542d73c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-os-command-injection-defense.md) | ⭐ 151 | `development` |
| [Reference Set 03 Query Parameterization](development/reference-set-03-query-parameterization_b2eba6f1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-query-parameterization.md) | ⭐ 151 | `development` |
| [Reference Set 03 Sql Injection Prevention](development/reference-set-03-sql-injection-prevention_80e8fc20/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-sql-injection-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Ssrf Prevention](development/reference-set-03-ssrf-prevention_2816fc60/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-ssrf-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Unvalidated Redirects](development/reference-set-03-unvalidated-redirects_aab33a7d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-unvalidated-redirects.md) | ⭐ 151 | `development` |
| [Reference Set 05 Clickjacking Defense](development/reference-set-05-clickjacking-defense_7353ee0b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-clickjacking-defense.md) | ⭐ 151 | `development` |
| [Reference Set 05 Csrf Prevention](development/reference-set-05-csrf-prevention_d0e31b35/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-csrf-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Css Security](development/reference-set-05-css-security_ae726937/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-css-security.md) | ⭐ 151 | `development` |
| [Reference Set 05 Dom Clobbering Prevention](development/reference-set-05-dom-clobbering-prevention_5f3a6922/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-dom-clobbering-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Dom Xss Prevention](development/reference-set-05-dom-xss-prevention_0d83c37e/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-dom-xss-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Html5 Security](development/reference-set-05-html5-security_2e726920/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-html5-security.md) | ⭐ 151 | `development` |
| [Reference Set 05 Http Headers](development/reference-set-05-http-headers_031e23d6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-http-headers.md) | ⭐ 151 | `development` |
| [Reference Set 05 Prototype Pollution Prevention](development/reference-set-05-prototype-pollution-prevention_91655ed3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-prototype-pollution-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Xss Prevention](development/reference-set-05-xss-prevention_fb45017c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-xss-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 06 Certificate Pinning](development/reference-set-06-certificate-pinning_f0af49c6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-certificate-pinning.md) | ⭐ 151 | `development` |
| [Reference Set 06 Cryptographic Storage](development/reference-set-06-cryptographic-storage_e3562829/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-cryptographic-storage.md) | ⭐ 151 | `development` |
| [Reference Set 06 Key Management](development/reference-set-06-key-management_33d0fc03/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-key-management.md) | ⭐ 151 | `development` |
| [Reference Set 06 Tls Cipher String](development/reference-set-06-tls-cipher-string_18d2f4a7/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-tls-cipher-string.md) | ⭐ 151 | `development` |
| [Reference Set 06 Tls](development/reference-set-06-tls_0321e669/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-tls.md) | ⭐ 151 | `development` |
| [Reference Set 06 Transport Layer Protection](development/reference-set-06-transport-layer-protection_60ac3b0b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-transport-layer-protection.md) | ⭐ 151 | `development` |
| [Reference Set 07 Logging Vocabulary](development/reference-set-07-logging-vocabulary_19ba9628/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-07-logging-vocabulary.md) | ⭐ 151 | `development` |
| [Reference Set 08 Error Handling](development/reference-set-08-error-handling_b2c98807/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-08-error-handling.md) | ⭐ 151 | `development` |
| [Reference Set 08 Xs Leaks](development/reference-set-08-xs-leaks_7525a3da/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-08-xs-leaks.md) | ⭐ 151 | `development` |
| [Reference Set 09 Ajax Security](development/reference-set-09-ajax-security_101c69a6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-ajax-security.md) | ⭐ 151 | `development` |
| [Reference Set 09 Dos Prevention](development/reference-set-09-dos-prevention_4021d4d8/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-dos-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 09 Graphql](development/reference-set-09-graphql_a917b3e3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-graphql.md) | ⭐ 151 | `development` |
| [Reference Set 09 Microservices Security](development/reference-set-09-microservices-security_a2daef28/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-microservices-security.md) | ⭐ 151 | `development` |
| [Reference Set 09 Rest Assessment](development/reference-set-09-rest-assessment_28d67ae6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-rest-assessment.md) | ⭐ 151 | `development` |
| [Reference Set 09 Web Service Security](development/reference-set-09-web-service-security_635ed73a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-web-service-security.md) | ⭐ 151 | `development` |
| [Reference Set 10 Database Security](development/reference-set-10-database-security_e6160141/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-10-database-security.md) | ⭐ 151 | `development` |
| [Reference Set 10 Secrets Management](development/reference-set-10-secrets-management_f17ad0b1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-10-secrets-management.md) | ⭐ 151 | `development` |
| [Reference Set Ext 11 Nodejs Docker](development/reference-set-ext-11-nodejs-docker_28321053/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-11-nodejs-docker.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Dependency Management](development/reference-set-ext-12-dependency-management_b7983d5d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-dependency-management.md) | ⭐ 151 | `development` |
| [Reference Set Ext 14 Automotive Security](development/reference-set-ext-14-automotive-security_6b2d4499/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-14-automotive-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 14 Mobile Security](development/reference-set-ext-14-mobile-security_714b27e2/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-14-mobile-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 15 Cloud Architecture](development/reference-set-ext-15-cloud-architecture_965b5461/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-15-cloud-architecture.md) | ⭐ 151 | `development` |
| [Reference Set Ext 16 Agentic Security](development/reference-set-ext-16-agentic-security_55a161fa/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-16-agentic-security.md) | ⭐ 151 | `development` |
| [Pentest](development/pentest_61dd5ef4/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/security/penetration-tester/commands/pentest.md) | ⭐ 1.3k | `development` |
| [Owasp Top 10](development/owasp_top_10_36e676ab/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/security/penetration-tester/skills/performing-penetration-testing/references/OWASP_TOP_10.md) | ⭐ 1.3k | `development` |
| [Remediation Playbook](development/remediation_playbook_7c3cbc1b/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/security/penetration-tester/skills/performing-penetration-testing/references/REMEDIATION_PLAYBOOK.md) | ⭐ 1.3k | `development` |
| [Security Headers](development/security_headers_88330111/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/security/penetration-tester/skills/performing-penetration-testing/references/SECURITY_HEADERS.md) | ⭐ 1.3k | `development` |

### Development/Devops (25 skills)

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
| [Config Settings](development/devops/033-config_settings_a882cf5b/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/config_settings.md) | 🔥 35.8k | `cache_hit` `cache_key` `proxy_base_url` |
| [Readme Cn](development/devops/readme-cn_162e7ea6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/README-cn.md) | ⭐ 151 | `development` |
| [Sync Anthropic Beta Headers](development/devops/376-sync_anthropic_beta_headers_cc87512d/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/sync_anthropic_beta_headers.md) | 🔥 35.8k | `development` |
| [Architecture Workflow Guide Cn](development/devops/architecture-workflow-guide-cn_1a497d91/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/ARCHITECTURE-WORKFLOW-GUIDE-cn.md) | ⭐ 151 | `development` |
| [Architecture Workflow Guide](development/devops/architecture-workflow-guide_51ef58a8/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/ARCHITECTURE-WORKFLOW-GUIDE.md) | ⭐ 151 | `development` |
| [P4 Security Design Review](development/devops/p4-security-design-review_3ff4329b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P4-SECURITY-DESIGN-REVIEW.md) | ⭐ 151 | `development` |
| [Control Set Ext 11 Infrastructure](development/devops/control-set-ext-11-infrastructure_c39fd27c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-11-infrastructure.md) | ⭐ 151 | `development` |
| [Control Set Ext 12 Supply Chain](development/devops/control-set-ext-12-supply-chain_e1e7576d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-12-supply-chain.md) | ⭐ 151 | `development` |
| [Control Set Ext 15 Cloud](development/devops/control-set-ext-15-cloud_f1136c8f/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-15-cloud.md) | ⭐ 151 | `development` |
| [Reference Set Ext 11 Kubernetes Security](development/devops/reference-set-ext-11-kubernetes-security_f9eb20b5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-11-kubernetes-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Cicd Security](development/devops/reference-set-ext-12-cicd-security_0521ba12/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-cicd-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Npm Security](development/devops/reference-set-ext-12-npm-security_3c3eb33b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-npm-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Supply Chain Security](development/devops/reference-set-ext-12-supply-chain-security_90024891/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-supply-chain-security.md) | ⭐ 151 | `development` |
| [Skill](development/devops/name-skill_9b800c44/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/security/penetration-tester/skills/performing-penetration-testing/SKILL.md) | ⭐ 1.3k | `development` |

### Development/Testing (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Middleware Chains](development/testing/085-middleware-chains_c239a6f7/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/middleware-chains.md) | ⭐ 10 | `development` |
| [Parallel Execution](development/testing/086-parallel-execution_5d3062a6/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/parallel-execution.md) | ⭐ 10 | `development` |
| [Skill](development/testing/threat-skill_0623a026/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/SKILL.md) | ⭐ 151 | `development` |

### Development/Tools (8 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Mcp Troubleshoot](development/tools/324-mcp_troubleshoot_6a5a0c41/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/mcp_troubleshoot.md) | 🔥 35.9k | `development` |
| [Desktop](development/tools/328-desktop_f1b0efd6/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/desktop.md) | ⭐ 51 | `development` |
| [Chain](development/tools/325-chain_02e6ae60/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/api/chain.md) | ⭐ 10 | `development` |
| [Knowledge Architecture V5.2 Cn](development/tools/knowledge-architecture-v52-cn_0d76937b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/KNOWLEDGE-ARCHITECTURE-v5.2-cn.md) | ⭐ 151 | `development` |
| [Report Design](development/tools/report-design_874a71e1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/REPORT-DESIGN.md) | ⭐ 151 | `development` |
| [Control Set 09 Api Security](development/tools/control-set-09-api-security_c77f3363/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-09-api-security.md) | ⭐ 151 | `development` |
| [Reference Set 05 Hsts](development/tools/reference-set-05-hsts_f54544af/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-hsts.md) | ⭐ 151 | `development` |
| [Reference Set 10 User Privacy](development/tools/reference-set-10-user-privacy_724c9ff6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-10-user-privacy.md) | ⭐ 151 | `development` |

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

*Last updated: 2026-02-14 04:21:53 UTC*
*Automatically maintained by SkillFlow*
