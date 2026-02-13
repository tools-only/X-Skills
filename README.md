# X-Skills

A curated collection of **345 AI-powered skills** organized into 14 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (9 skills)
- **Automation/Workflow** (18 skills)
- **Commercial** (50 skills)
- **Communication** (7 skills)
- **Content Creation** (35 skills)
- **Daily Assistant** (17 skills)
- **Data Analysis** (32 skills)
- **Development** (110 skills)
- **Development/Devops** (19 skills)
- **Development/Testing** (5 skills)
- **Development/Tools** (31 skills)
- **Investment** (4 skills)
- **Productivity** (3 skills)
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


### Automation/Scripting (9 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Param Forge Reference](automation/scripting/100-param_forge_reference_35994bfd/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/param_forge_reference.md) | ⭐ 51 | `automation` |
| [Skill](automation/scripting/003-name-skill_7d2e3a7b/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/pdf-processing-pro/SKILL.md) | ⭐ 21 | `automation` |
| [Chrome Installation](automation/scripting/094-chrome-installation_cb622122/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/devtools/references/chrome-installation.md) | ⭐ 21 | `automation` |
| [Workflows](automation/scripting/069-workflows_1f8e9003/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/pdf-processing-pro/references/workflows.md) | ⭐ 21 | `automation` |
| [Api Storage](automation/scripting/095-api-storage_512c2691/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/api-storage.md) | ⭐ 21 | `automation` |
| [Skill](automation/scripting/003-name-skill_98469709/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_IsPlayerPawn/SKILL.md) | ⭐ 17 | `automation` |
| [Skill](automation/scripting/003-name-skill_2f06a357/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_SetStateChanged/SKILL.md) | ⭐ 17 | `automation` |
| [Skill](automation/scripting/003-name-skill_f645e305/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerPawn_GetEyeAngles/SKILL.md) | ⭐ 17 | `automation` |
| [Skill](automation/scripting/003-name-skill_472f8206/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CSource2Server_Init-AND-CGameEventManager_Init-AND-gameeventmanager-AND-s_GameEventManager/SKILL.md) | ⭐ 17 | `automation` |

### Automation/Workflow (18 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/002-name-skill_ce4bb9a9/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-deploy/SKILL.md) | ⭐ 65 | `automation` |
| [Cli Guide](automation/workflow/137-cli-guide_49550e7a/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/docs/cli-guide.md) | ⭐ 65 | `automation` |
| [Skill](automation/workflow/002-name-skill_fc87db77/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/gene-regulatory-networks/coexpression-networks/SKILL.md) | ⭐ 233 | `automation` |
| [Skill](automation/workflow/002-name-skill_007690c0/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/temporal-genomics/circadian-rhythms/SKILL.md) | ⭐ 233 | `automation` |
| [Skill](automation/workflow/002-name-skill_645242d5/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/workflows/edna-pipeline/SKILL.md) | ⭐ 233 | `automation` |
| [Usage Guide](automation/workflow/031-usage-guide_d5b8c5d0/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/workflows/edna-pipeline/usage-guide.md) | ⭐ 233 | `automation` |
| [Usage Guide](automation/workflow/031-usage-guide_c438f645/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/workflows/timecourse-pipeline/usage-guide.md) | ⭐ 233 | `automation` |
| [Readme.Ja](automation/workflow/040-readmeja_e5d1a15a/) | [japan1988/multi-agent-mediation](https://raw.githubusercontent.com/japan1988/multi-agent-mediation/main/README.ja.md) | ⭐ 27 | `automation` |
| [Agents](automation/workflow/073-agents_8aee7225/) | [opendatahub-io/ai-helpers](https://raw.githubusercontent.com/opendatahub-io/ai-helpers/main/AGENTS.md) | ⭐ 13 | `automation` |
| [Skill](automation/workflow/002-name-skill_a9a208cd/) | [adenhq/hive](https://raw.githubusercontent.com/adenhq/hive/main/.claude/skills/hive-create/SKILL.md) | 🔥 7.1k | `automation` |
| [Skill](automation/workflow/002-name-skill_1b2668f1/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-develop-web-game/skills/openai-develop-web-game/SKILL.md) | ⭐ 56 | `automation` |
| [Skill](automation/workflow/002-name-skill_7a0692d3/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/devtools/SKILL.md) | ⭐ 21 | `automation` |
| [Mcp Configuration](automation/workflow/136-mcp-configuration_c1914cda/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/devtools/references/mcp-configuration.md) | ⭐ 21 | `automation` |
| [Git Commands](automation/workflow/137-git-commands_b6a53e87/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/git-commit-helper/references/git-commands.md) | ⭐ 21 | `automation` |
| [Troubleshooting](automation/workflow/138-troubleshooting_fe01a5b8/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/pdf-processing-pro/references/troubleshooting.md) | ⭐ 21 | `automation` |
| [Version](automation/workflow/065-version_b186ef4d/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/VERSION.md) | ⭐ 63 | `automation` |
| [Run Conductor.Prompt](automation/workflow/134-run-conductorprompt_c290deb8/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/prompts/run-conductor.prompt.md) | ⭐ 63 | `automation` |
| [Executive Pitch](automation/workflow/067-executive-pitch_ec3ad432/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/executive-pitch.md) | ⭐ 63 | `automation` |

### Commercial (50 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/210-name-skill_10221f10/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/manuscript-ingest/SKILL.md) | ⭐ 197 | `commercial` |
| [Security Policy](commercial/371-security_policy_f8dea503/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/SECURITY_POLICY.md) | ⭐ 18 | `commercial` |
| [Customer Success](commercial/372-customer_success_85216f92/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/customer_success.md) | ⭐ 18 | `commercial` |
| [Executive](commercial/373-executive_e9668604/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/executive.md) | ⭐ 18 | `commercial` |
| [Finance](commercial/374-finance_8fdccf44/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/finance.md) | ⭐ 18 | `commercial` |
| [Marketing](commercial/375-marketing_499062ff/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/marketing.md) | ⭐ 18 | `commercial` |
| [Agent Coverage](commercial/376-agent_coverage_affae21f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AGENT_COVERAGE.md) | ⭐ 18 | `commercial` |
| [Batch 20260211 212328 Aggregate](commercial/377-batch_20260211_212328_aggregate_b193fed8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_212328_aggregate.md) | ⭐ 18 | `commercial` |
| [Failure Analysis Batch 20260211 212328](commercial/378-failure_analysis_batch_20260211_212328_ebcc0636/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_212328.md) | ⭐ 18 | `commercial` |
| [Cs Risk Mitigation Playbook Eval](commercial/379-cs_risk_mitigation_playbook_eval_9849a6ec/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_risk_mitigation_playbook_eval.md) | ⭐ 18 | `commercial` |
| [Devex Rate Limit Advisor Eval](commercial/380-devex_rate_limit_advisor_eval_e60ff981/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_rate_limit_advisor_eval.md) | ⭐ 18 | `commercial` |
| [Mon Commitment Tracker Eval](commercial/381-mon_commitment_tracker_eval_be78232d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_commitment_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Mon Discount Optimizer Eval](commercial/382-mon_discount_optimizer_eval_5ef6bd58/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_discount_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Plg Progressive Disclosure Eval](commercial/383-plg_progressive_disclosure_eval_829833e5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_progressive_disclosure_eval.md) | ⭐ 18 | `commercial` |
| [Revops Commit Accuracy Eval](commercial/384-revops_commit_accuracy_eval_97f0580c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_commit_accuracy_eval.md) | ⭐ 18 | `commercial` |
| [Vcf Value Discovery Eval](commercial/385-vcf_value_discovery_eval_a49e0e39/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/vcf_value_discovery_eval.md) | ⭐ 18 | `commercial` |
| [Migration Guide V1](commercial/251-migration_guide_v1_845a62d9/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/MIGRATION_GUIDE_v1.md) | ⭐ 4.0k | `commercial` |
| [Aesthetic](commercial/386-aesthetic_beed860c/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/aesthetic.md) | ⭐ 51 | `commercial` |
| [Bank Penetration](commercial/371-bank-penetration_c2e5ee98/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/bank-penetration.md) | ⭐ 56 | `commercial` |
| [Logic Flaws](commercial/372-logic-flaws_d3713f41/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/logic-flaws.md) | ⭐ 56 | `commercial` |
| [Logic Flaws Checklist](commercial/373-logic-flaws-checklist_c6b46334/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/logic-flaws-checklist.md) | ⭐ 56 | `commercial` |
| [Agents](commercial/007-agents_3a4fc98b/) | [maragudk/skills](https://raw.githubusercontent.com/maragudk/skills/main/AGENTS.md) | ⭐ 31 | `commercial` |
| [Skill](commercial/210-name-skill_f5d8c52c/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/SKILL.md) | ⭐ 21 | `commercial` |
| [Skill](commercial/210-name-skill_ff987588/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/SKILL.md) | ⭐ 21 | `commercial` |
| [Audiences](commercial/374-audiences_cdae62b0/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/audiences.md) | ⭐ 21 | `commercial` |
| [Custom Dimensions](commercial/375-custom-dimensions_ac9e88c6/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/custom-dimensions.md) | ⭐ 21 | `commercial` |
| [Custom Events](commercial/376-custom-events_442ce8c1/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/custom-events.md) | ⭐ 21 | `commercial` |
| [Data Management](commercial/377-data-management_1f9e80e0/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/data-management.md) | ⭐ 21 | `commercial` |
| [Debugview](commercial/378-debugview_f2b050db/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/debugview.md) | ⭐ 21 | `commercial` |
| [Events Fundamentals](commercial/379-events-fundamentals_acdeb349/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/events-fundamentals.md) | ⭐ 21 | `commercial` |
| [Measurement Protocol](commercial/380-measurement-protocol_a5b644a5/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/measurement-protocol.md) | ⭐ 21 | `commercial` |
| [Recommended Events](commercial/381-recommended-events_0dccd208/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/recommended-events.md) | ⭐ 21 | `commercial` |
| [Setup](commercial/382-setup_26743b72/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/setup.md) | ⭐ 21 | `commercial` |
| [Datalayer](commercial/383-datalayer_a05f8557/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/datalayer.md) | ⭐ 21 | `commercial` |
| [Debugging](commercial/384-debugging_5966ee22/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/debugging.md) | ⭐ 21 | `commercial` |
| [Setup](commercial/382-setup_d641441d/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/setup.md) | ⭐ 21 | `commercial` |
| [Triggers](commercial/385-triggers_ddc0c0ce/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/triggers.md) | ⭐ 21 | `commercial` |
| [Api Admin](commercial/386-api-admin_cd0708a4/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/api-admin.md) | ⭐ 21 | `commercial` |
| [Api Storefront](commercial/387-api-storefront_b0702ce0/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/api-storefront.md) | ⭐ 21 | `commercial` |
| [App Development](commercial/017-app-development_abcde65b/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/app-development.md) | ⭐ 21 | `commercial` |
| [Debugging](commercial/384-debugging_e4f0435a/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/debugging.md) | ⭐ 21 | `commercial` |
| [Functions](commercial/388-functions_81f7b794/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/functions.md) | ⭐ 21 | `VIP` |
| [Hydrogen](commercial/389-hydrogen_083eca7c/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/hydrogen.md) | ⭐ 21 | `commercial` |
| [Liquid Filters](commercial/390-liquid-filters_85c17799/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/liquid-filters.md) | ⭐ 21 | `commercial` |
| [Liquid Objects](commercial/391-liquid-objects_351e92e4/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/liquid-objects.md) | ⭐ 21 | `commercial` |
| [Liquid Syntax](commercial/392-liquid-syntax_70d26c64/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/liquid-syntax.md) | ⭐ 21 | `commercial` |
| [Performance](commercial/393-performance_6a12d640/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/performance.md) | ⭐ 21 | `commercial` |
| [Theme Development](commercial/394-theme-development_dc454c66/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/shopify-developer/skills/shopify-developer/references/theme-development.md) | ⭐ 21 | `commercial` |
| [Api Cookies](commercial/395-api-cookies_cd598abf/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/api-cookies.md) | ⭐ 21 | `commercial` |
| [Plan Requirements.Prompt](commercial/311-plan-requirementsprompt_2651eb49/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/prompts/plan-requirements.prompt.md) | ⭐ 65 | `commercial` |

### Communication (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Syntax Reference](communication/113-syntax-reference_a16eba1d/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/syntax-reference.md) | ⭐ 65 | `communication` |
| [Project Context](communication/255-project_context_5cc91b78/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/param_forge_ref/param_forge/PROJECT_CONTEXT.md) | ⭐ 51 | `communication` |
| [Readme Cn](communication/256-readme_cn_a75d0219/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/README_CN.md) | ⭐ 57 | `communication` |
| [Weak Password Checklist](communication/255-weak-password-checklist_0d3963ca/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/weak-password-checklist.md) | ⭐ 56 | `communication` |
| [Patterns](communication/084-patterns_6c763e12/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-ads-scripts/skills/google-ads-scripts/references/patterns.md) | ⭐ 21 | `communication` |
| [Patterns](communication/084-patterns_73fdd9ee/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-apps-script/skills/google-apps-script/references/patterns.md) | ⭐ 21 | `communication` |
| [Workshop Checklist](communication/209-workshop-checklist_d89cc5e8/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/workshop-checklist.md) | ⭐ 65 | `communication` |

### Content Creation (35 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Quick Wins](content-creation/353-quick_wins_29f2473a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/QUICK_WINS.md) | ⭐ 18 | `content creation` |
| [Skill Chains](content-creation/354-skill_chains_312b4d53/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/SKILL_CHAINS.md) | ⭐ 18 | `content creation` |
| [Skills Catalog](content-creation/355-skills_catalog_b303f24f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/SKILLS_CATALOG.md) | ⭐ 18 | `content creation` |
| [Ai Agent Composability Analysis](content-creation/356-ai_agent_composability_analysis_8ca26161/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AI_AGENT_COMPOSABILITY_ANALYSIS.md) | ⭐ 18 | `content creation` |
| [Skill](content-creation/049-name-skill_42cdd582/) | [cat-xierluo/legal-skills](https://raw.githubusercontent.com/cat-xierluo/legal-skills/main/skills/universal-media-downloader/SKILL.md) | ⭐ 18 | `content creation` |
| [Spawn V0 Local Only](content-creation/357-spawn_v0_local_only_1cb9fb27/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/prds/spawn_v0_local_only.md) | ⭐ 51 | `content creation` |
| [Ethics](content-creation/372-ethics_dbde6459/) | [opendatahub-io/ai-helpers](https://raw.githubusercontent.com/opendatahub-io/ai-helpers/main/ETHICS.md) | ⭐ 13 | `content creation` |
| [Plan](content-creation/353-plan_860e9528/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/commands/plan.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/049-name-skill_82868fe3/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/humanizer/skills/humanizer/SKILL.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/049-name-skill_d1631d58/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/last30days/skills/last30days/SKILL.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/049-name-skill_76a45975/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/skills/planning-with-files/SKILL.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/049-name-skill_8a2e5b56/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/x-research/skills/x-research/SKILL.md) | ⭐ 56 | `content creation` |
| [Principles](content-creation/354-principles_c2c61108/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/skills/planning-with-files/references/principles.md) | ⭐ 56 | `content creation` |
| [Skill Lifecycle](content-creation/355-skill-lifecycle_119157d0/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/skill-extractor/skills/skill-extractor/references/skill-lifecycle.md) | ⭐ 56 | `content creation` |
| [Command Execution](content-creation/356-command-execution_ffdc6b02/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/command-execution.md) | ⭐ 56 | `content creation` |
| [Xss](content-creation/357-xss_d1bc8bd4/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/xss.md) | ⭐ 56 | `content creation` |
| [Command Execution Checklist](content-creation/358-command-execution-checklist_39b874c0/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/command-execution-checklist.md) | ⭐ 56 | `content creation` |
| [Ssrf Checklist](content-creation/359-ssrf-checklist_113f13e6/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/ssrf-checklist.md) | ⭐ 56 | `content creation` |
| [Validation Gates](content-creation/360-validation-gates_a1b8f9f3/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/agents/validation-gates.md) | ⭐ 21 | `content creation` |
| [Setup Skill Hook](content-creation/361-setup-skill-hook_14752a8c/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/commands/setup-skill-hook.md) | ⭐ 21 | `content creation` |
| [Skill](content-creation/049-name-skill_a4befc32/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/docx/SKILL.md) | ⭐ 21 | `content creation` |
| [Skill](content-creation/049-name-skill_aaac38a4/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/humanizer/skills/humanizer/SKILL.md) | ⭐ 21 | `content creation` |
| [Skill](content-creation/049-name-skill_30c11759/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/message/skills/message/SKILL.md) | ⭐ 21 | `content creation` |
| [Skill](content-creation/049-name-skill_82012a54/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/SKILL.md) | ⭐ 21 | `content creation` |
| [Redlining Workflow](content-creation/362-redlining-workflow_16b37008/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/docx/references/redlining-workflow.md) | ⭐ 21 | `content creation` |
| [Openpyxl Patterns](content-creation/363-openpyxl-patterns_2e568fb2/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/xlsx/references/openpyxl-patterns.md) | ⭐ 21 | `content creation` |
| [Advanced Editing](content-creation/364-advanced-editing_e346f129/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ffmpeg/skills/ffmpeg-cli/references/advanced-editing.md) | ⭐ 21 | `content creation` |
| [Asset Generation](content-creation/365-asset-generation_e5564ccd/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ffmpeg/skills/ffmpeg-cli/references/asset-generation.md) | ⭐ 21 | `content creation` |
| [Core Concepts](content-creation/366-core-concepts_2465de6c/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ffmpeg/skills/ffmpeg-cli/references/core-concepts.md) | ⭐ 21 | `content creation` |
| [Encoding And Settings](content-creation/367-encoding-and-settings_f5c747a9/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ffmpeg/skills/ffmpeg-cli/references/encoding-and-settings.md) | ⭐ 21 | `content creation` |
| [Communication Patterns](content-creation/368-communication-patterns_4b16b781/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/humanizer/skills/humanizer/references/communication-patterns.md) | ⭐ 21 | `content creation` |
| [Content Patterns](content-creation/369-content-patterns_9df20e5f/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/humanizer/skills/humanizer/references/content-patterns.md) | ⭐ 21 | `content creation` |
| [Language Patterns](content-creation/370-language-patterns_3a724dd8/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/humanizer/skills/humanizer/references/language-patterns.md) | ⭐ 21 | `content creation` |
| [Skill](content-creation/049-name-skill_28e46635/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/docs-writer/SKILL.md) | ⭐ 63 | `content creation` |
| [Markdown.Instructions](content-creation/256-markdowninstructions_c76ba598/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/markdown.instructions.md) | ⭐ 65 | `bicep` `iac` `azure` |

### Daily Assistant (17 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Research Notes](daily-assistant/267-research_notes_3883eb9d/) | [taylorsatula/mira-OSS](https://raw.githubusercontent.com/taylorsatula/mira-OSS/main/RESEARCH_NOTES.md) | ⭐ 389 | `daily assistant` |
| [Copilot Instructions](daily-assistant/266-copilot-instructions_1423668e/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/.github/copilot-instructions.md) | ⭐ 51 | `daily assistant` |
| [Brood Aip.Instructions](daily-assistant/267-brood-aipinstructions_be0c3008/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/.github/instructions/brood-aip.instructions.md) | ⭐ 51 | `daily assistant` |
| [Brood Aip](daily-assistant/268-brood-aip_be4ed266/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/.windsurf/rules/brood-aip.md) | ⭐ 51 | `daily assistant` |
| [Status](daily-assistant/268-status_2e342fb5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/commands/status.md) | ⭐ 56 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_47ec12a5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-playwright/skills/openai-playwright/SKILL.md) | ⭐ 56 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_f9a1f323/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-security-threat-model/skills/openai-security-threat-model/SKILL.md) | ⭐ 56 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_92eb7fef/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-sentry/skills/openai-sentry/SKILL.md) | ⭐ 56 | `daily assistant` |
| [Style Patterns](daily-assistant/269-style-patterns_4ce412a5/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/humanizer/skills/humanizer/references/style-patterns.md) | ⭐ 21 | `daily assistant` |
| [Architect.Agent](daily-assistant/204-architectagent_a5b83e17/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/agents/architect.agent.md) | ⭐ 63 | `daily assistant` |
| [Bicep Plan.Agent](daily-assistant/235-bicep-planagent_9e68826d/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/agents/bicep-plan.agent.md) | ⭐ 63 | `Environment` `ManagedBy` `Project` |
| [Design.Agent](daily-assistant/263-designagent_23213bfc/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/agents/design.agent.md) | ⭐ 63 | `daily assistant` |
| [Requirements.Agent](daily-assistant/236-requirementsagent_0c3bbb4a/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/agents/requirements.agent.md) | ⭐ 63 | `daily assistant` |
| [Bicep Code.Agent](daily-assistant/251-bicep-codeagent_9b9ea6c5/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/agents/bicep-code.agent.md) | ⭐ 65 | `daily assistant` |
| [Deploy.Agent](daily-assistant/252-deployagent_a9396531/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/agents/deploy.agent.md) | ⭐ 65 | `daily assistant` |
| [Artifact H2 Reference.Instructions](daily-assistant/207-artifact-h2-referenceinstructions_a393e6c0/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/artifact-h2-reference.instructions.md) | ⭐ 65 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_a1e7e6e0/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/azure-artifacts/SKILL.md) | ⭐ 65 | `daily assistant` |

### Data Analysis (32 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](data-analysis/226-name-skill_d0682987/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/idea-pool-expander/SKILL.md) | ⭐ 197 | `data analysis` |
| [Skill](data-analysis/226-name-skill_5f28dab2/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/differential-expression/timeseries-de/SKILL.md) | ⭐ 233 | `data analysis` |
| [Usage Guide](data-analysis/279-usage-guide_2425cb78/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/differential-expression/timeseries-de/usage-guide.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_577bfb3d/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/biodiversity-metrics/SKILL.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_63c7037d/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/conservation-genetics/SKILL.md) | ⭐ 233 | `data analysis` |
| [Usage Guide](data-analysis/279-usage-guide_3ebb7620/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/conservation-genetics/usage-guide.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_7fabca72/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/species-delimitation/SKILL.md) | ⭐ 233 | `data analysis` |
| [Usage Guide](data-analysis/279-usage-guide_8edb9336/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/species-delimitation/usage-guide.md) | ⭐ 233 | `data analysis` |
| [Usage Guide](data-analysis/279-usage-guide_5a4a9b50/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/gene-regulatory-networks/coexpression-networks/usage-guide.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_74be9c92/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/gene-regulatory-networks/differential-networks/SKILL.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_5bd8514b/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/microbiome/diversity-analysis/SKILL.md) | ⭐ 233 | `data analysis` |
| [Usage Guide](data-analysis/279-usage-guide_82645619/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/microbiome/diversity-analysis/usage-guide.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_e7c52123/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/population-genetics/population-structure/SKILL.md) | ⭐ 233 | `data analysis` |
| [Usage Guide](data-analysis/279-usage-guide_b3403d15/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/population-genetics/population-structure/usage-guide.md) | ⭐ 233 | `data analysis` |
| [Usage Guide](data-analysis/279-usage-guide_423469c6/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/temporal-genomics/circadian-rhythms/usage-guide.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_08444c13/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/temporal-genomics/periodicity-detection/SKILL.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_9a3b98e8/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/temporal-genomics/temporal-grn/SKILL.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_55773a6b/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/temporal-genomics/trajectory-modeling/SKILL.md) | ⭐ 233 | `data analysis` |
| [Skill](data-analysis/226-name-skill_0def0501/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/workflows/timecourse-pipeline/SKILL.md) | ⭐ 233 | `data analysis` |
| [Visual Prompting V0](data-analysis/485-visual_prompting_v0_52572161/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/visual_prompting_v0.md) | ⭐ 51 | `data analysis` |
| [Import Skill](data-analysis/478-import-skill_0f090701/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/.claude/commands/import-skill.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/226-name-skill_5f6df240/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-doc/skills/openai-doc/SKILL.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/226-name-skill_ed5b6f7e/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-pdf/skills/openai-pdf/SKILL.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/226-name-skill_b03aba97/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-security-ownership-map/skills/openai-security-ownership-map/SKILL.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/226-name-skill_94733cad/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-spreadsheet/skills/openai-spreadsheet/SKILL.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/226-name-skill_6eec1f38/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/xlsx/SKILL.md) | ⭐ 21 | `data analysis` |
| [Reporting](data-analysis/479-reporting_37cb7cfd/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/reporting.md) | ⭐ 21 | `data analysis` |
| [Preprocess Func Sig Via Mcp](data-analysis/476-preprocess_func_sig_via_mcp_62e04de7/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.serena/memories/preprocess_func_sig_via_mcp.md) | ⭐ 17 | `data analysis` |
| [Preprocess Gen Func Sig Via Mcp](data-analysis/477-preprocess_gen_func_sig_via_mcp_77fd66b9/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.serena/memories/preprocess_gen_func_sig_via_mcp.md) | ⭐ 17 | `data analysis` |
| [Preprocess Gen Gv Sig Via Mcp](data-analysis/478-preprocess_gen_gv_sig_via_mcp_8a79aa32/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.serena/memories/preprocess_gen_gv_sig_via_mcp.md) | ⭐ 17 | `data analysis` |
| [Skill](data-analysis/226-name-skill_87e653f4/) | [bowenliang123/md_exporter](https://raw.githubusercontent.com/bowenliang123/md_exporter/main/SKILL.md) | ⭐ 182 | `data analysis` |
| [04 Governance Constraints](data-analysis/382-04-governance-constraints_d9ed2874/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/agent-output/static-webapp/04-governance-constraints.md) | ⭐ 65 | `data analysis` |

### Development (110 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Credits](development/239-credits_7d3b11da/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/CREDITS.md) | ⭐ 65 | `development` |
| [Actions Reference](development/010-actions-reference_733b0bdd/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/actions-reference.md) | ⭐ 65 | `development` |
| [Migration Guide](development/584-migration-guide_6f7f2a5a/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/migration-guide.md) | ⭐ 65 | `development` |
| [Api Quickstart](development/1540-api-quickstart_a217f99d/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/api-quickstart.md) | ⭐ 4.0k | `development` |
| [Skill](development/1178-name-skill_33884f18/) | [cat-xierluo/legal-skills](https://raw.githubusercontent.com/cat-xierluo/legal-skills/main/skills/repo-research/SKILL.md) | ⭐ 18 | `development` |
| [Agents](development/028-agents_e77da732/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/AGENTS.md) | ⭐ 51 | `development` |
| [Claude](development/140-claude_84c07acf/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/CLAUDE.md) | ⭐ 51 | `development` |
| [Releasing](development/2885-releasing_f6b07c11/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/RELEASING.md) | ⭐ 51 | `development` |
| [Aip Aws Setup](development/2886-aip_aws_setup_103a72ac/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/aip_aws_setup.md) | ⭐ 51 | `development` |
| [Skill](development/1178-name-skill_d43c3e9c/) | [jwiegley/claude-prompts](https://raw.githubusercontent.com/jwiegley/claude-prompts/main/skills/claude-code/SKILL.md) | ⭐ 10 | `development` |
| [Developer Guide](development/282-developer-guide_ff9975f1/) | [adenhq/hive](https://raw.githubusercontent.com/adenhq/hive/main/docs/developer-guide.md) | 🔥 7.1k | `development` |
| [Metrics](development/575-metrics_32f0d9c1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/METRICS.md) | ⭐ 18 | `development` |
| [Recipe Audit](development/2850-recipe_audit_dd40a0b5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/RECIPE_AUDIT.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_bc6bbf16/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-security-best-practices/skills/openai-security-best-practices/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/1178-name-skill_bbfc4111/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/skill-extractor/skills/skill-extractor/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/1178-name-skill_0d3d1f7c/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/SKILL.md) | ⭐ 56 | `development` |
| [File Upload](development/2874-file-upload_b013ed1e/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/file-upload.md) | ⭐ 56 | `development` |
| [Info Disclosure](development/2875-info-disclosure_d4862aba/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/info-disclosure.md) | ⭐ 56 | `development` |
| [Path Traversal](development/2876-path-traversal_c73eea45/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/path-traversal.md) | ⭐ 56 | `development` |
| [Sql Injection](development/2877-sql-injection_cc9c615c/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/sql-injection.md) | ⭐ 56 | `development` |
| [Telecom Penetration](development/2878-telecom-penetration_5c723184/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/telecom-penetration.md) | ⭐ 56 | `development` |
| [X Api](development/2879-x-api_b97fa218/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/x-research/skills/x-research/references/x-api.md) | ⭐ 56 | `development` |
| [Csrf Checklist](development/2880-csrf-checklist_e821ad41/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/csrf-checklist.md) | ⭐ 56 | `development` |
| [File Upload Checklist](development/2881-file-upload-checklist_95e05da5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/file-upload-checklist.md) | ⭐ 56 | `development` |
| [Info Disclosure Checklist](development/2882-info-disclosure-checklist_aac2c34e/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/info-disclosure-checklist.md) | ⭐ 56 | `development` |
| [Path Traversal Checklist](development/2883-path-traversal-checklist_e420e992/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/path-traversal-checklist.md) | ⭐ 56 | `development` |
| [Rce Checklist](development/2884-rce-checklist_60e3a56a/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/rce-checklist.md) | ⭐ 56 | `development` |
| [Xxe Checklist](development/2885-xxe-checklist_0521eb1b/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/xxe-checklist.md) | ⭐ 56 | `development` |
| [Agents](development/028-agents_a5d5dc88/) | [strands-agents/sdk-python](https://raw.githubusercontent.com/strands-agents/sdk-python/main/AGENTS.md) | 🔥 5.1k | `development` |
| [Style Guide](development/1377-style_guide_a952b097/) | [strands-agents/sdk-python](https://raw.githubusercontent.com/strands-agents/sdk-python/main/docs/STYLE_GUIDE.md) | 🔥 5.1k | `development` |
| [Provider Migration Checklist](development/2852-provider_migration_checklist_97c70964/) | [openlit/openlit](https://raw.githubusercontent.com/openlit/openlit/main/contributors/PROVIDER_MIGRATION_CHECKLIST.md) | ⭐ 2.2k | `development` |
| [Skill](development/1178-name-skill_0cfcbc4b/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/analyzing-data/SKILL.md) | ⭐ 208 | `development` |
| [Migration Patterns](development/1521-migration-patterns_742e6f49/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/migrating-airflow-2-to-3/reference/migration-patterns.md) | ⭐ 208 | `development` |
| [Skill](development/1178-name-skill_7abf2b26/) | [zrt-ai-lab/opencode-skills](https://raw.githubusercontent.com/zrt-ai-lab/opencode-skills/main/log-analyzer/SKILL.md) | ⭐ 42 | `development` |
| [Claude](development/140-claude_9f4a5980/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/CLAUDE.md) | ⭐ 21 | `development` |
| [Reflection](development/2889-reflection_f2683bee/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/claudecode/commands/reflection.md) | ⭐ 21 | `development` |
| [Documentation Manager](development/2890-documentation-manager_987b65fe/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/agents/documentation-manager.md) | ⭐ 21 | `development` |
| [Fullstack Developer](development/2623-fullstack-developer_4efd5c7e/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/agents/fullstack-developer.md) | ⭐ 21 | `{
    type: String` `trim: true` `lowercase: true
  }` |
| [Skill Architect](development/2891-skill-architect_3baae92b/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/agents/skill-architect.md) | ⭐ 21 | `development` |
| [Skill](development/1178-name-skill_f7ea5091/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/git-commit-helper/SKILL.md) | ⭐ 21 | `development` |
| [Skill](development/1178-name-skill_82f514f9/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/git-worktrees/skills/git-worktrees/SKILL.md) | ⭐ 21 | `development` |
| [Skill](development/1178-name-skill_a0f9afc6/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-ads-scripts/skills/google-ads-scripts/SKILL.md) | ⭐ 21 | `development` |
| [Skill](development/1178-name-skill_73dc1dd9/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-apps-script/skills/google-apps-script/SKILL.md) | ⭐ 21 | `development` |
| [Skill](development/1178-name-skill_150a5c43/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/SKILL.md) | ⭐ 21 | `development` |
| [Skill](development/1178-name-skill_cd24e17d/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/react-best-practices/skills/react-best-practices/SKILL.md) | ⭐ 21 | `development` |
| [Financial Model Standards](development/2892-financial-model-standards_058c3646/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/xlsx/references/financial-model-standards.md) | ⭐ 21 | `development` |
| [Ads Api Reference](development/2893-ads-api-reference_6682095f/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-ads-scripts/skills/google-ads-scripts/references/ads-api-reference.md) | ⭐ 21 | `development` |
| [Best Practices](development/102-best-practices_99fa6cef/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-ads-scripts/skills/google-ads-scripts/references/best-practices.md) | ⭐ 21 | `development` |
| [Gtag](development/2894-gtag_bef05560/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/gtag.md) | ⭐ 21 | `development` |
| [Gtm Integration](development/2895-gtm-integration_3c95c7b6/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/gtm-integration.md) | ⭐ 21 | `development` |
| [Privacy](development/2896-privacy_09e2db4c/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/privacy.md) | ⭐ 21 | `development` |
| [Apps Script Api Reference](development/2897-apps-script-api-reference_c61d553d/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-apps-script/skills/google-apps-script/references/apps-script-api-reference.md) | ⭐ 21 | `development` |
| [Best Practices](development/102-best-practices_4fc37192/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-apps-script/skills/google-apps-script/references/best-practices.md) | ⭐ 21 | `development` |
| [Best Practices](development/102-best-practices_e7e2a0be/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/best-practices.md) | ⭐ 21 | `development` |
| [Tags](development/2898-tags_fab07438/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/tags.md) | ⭐ 21 | `development` |
| [Variables](development/2899-variables_4ba0899e/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/variables.md) | ⭐ 21 | `development` |
| [Api Async](development/2900-api-async_acf5fc6a/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/api-async.md) | ⭐ 21 | `development` |
| [Browser Compatibility](development/2901-browser-compatibility_1e5e5ed8/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/browser-compatibility.md) | ⭐ 21 | `development` |
| [Common Pitfalls](development/2902-common-pitfalls_7fcd2c42/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/common-pitfalls.md) | ⭐ 21 | `development` |
| [Debugging](development/266-debugging_de6ce48e/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/debugging.md) | ⭐ 21 | `development` |
| [Header Reference](development/2903-header-reference_f7d994b6/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/header-reference.md) | ⭐ 21 | `development` |
| [Http Requests](development/2904-http-requests_8c188ad3/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/http-requests.md) | ⭐ 21 | `development` |
| [Patterns](development/672-patterns_b6373511/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/patterns.md) | ⭐ 21 | `development` |
| [Security Checklist](development/2905-security-checklist_68a4b99c/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/security-checklist.md) | ⭐ 21 | `development` |
| [Web Requests](development/2906-web-requests_d683cf25/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/web-requests.md) | ⭐ 21 | `development` |
| [Archetypes](development/2907-archetypes_45dcef2e/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/skills/skill-mastery/references/archetypes.md) | ⭐ 21 | `development` |
| [Token Hierarchy](development/2908-token-hierarchy_c583db8a/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/skills/skill-mastery/references/token-hierarchy.md) | ⭐ 21 | `development` |
| [Skill](development/1178-name-skill_dbe46d44/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_StartTouch-AND-CBaseEntity_Touch-AND-CBaseEntity_EndTouch/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1178-name-skill_78278410/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseTrigger_StartTouch-AND-CBaseTrigger_EndTouch/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1178-name-skill_ba137163/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CTriggerPush_Touch/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1178-name-skill_a50a5cdc/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/write-globalvar-as-yaml/SKILL.md) | ⭐ 17 | `development` |
| [Update Docs On Code Change.Instructions](development/2853-update-docs-on-code-changeinstructions_b72cfda7/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/update-docs-on-code-change.instructions.md) | ⭐ 63 | `development` |
| [Diagnose Resources.Prompt](development/2854-diagnose-resourcesprompt_efa3ee1b/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/prompts/diagnose-resources.prompt.md) | ⭐ 63 | `development` |
| [Pilot Success Checklist](development/2125-pilot-success-checklist_d6827b5c/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/pilot-success-checklist.md) | ⭐ 63 | `development` |
| [Skill](development/1178-name-skill_242173f0/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/azure-adr/SKILL.md) | ⭐ 63 | `development` |
| [Infraops Conductor.Agent](development/2701-infraops-conductoragent_97975bdc/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/agents/infraops-conductor.agent.md) | ⭐ 65 | `development` |
| [Assess Architecture.Prompt](development/2879-assess-architectureprompt_df715185/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/prompts/assess-architecture.prompt.md) | ⭐ 65 | `development` |
| [Objection Handling](development/2124-objection-handling_16d1104d/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/objection-handling.md) | ⭐ 65 | `development` |
| [Roi Calculator](development/2126-roi-calculator_8ba8d13d/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/roi-calculator.md) | ⭐ 65 | `development` |
| [Skill](development/1178-name-skill_d0e41534/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/azure-diagrams/SKILL.md) | ⭐ 65 | `development` |
| [Skill](development/1178-name-skill_fb0daf12/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/git-commit/SKILL.md) | ⭐ 65 | `development` |
| [Agents](development/028-agents_15100b81/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/setups/AGENTS.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_be5d1ead/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-python-x402-client/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_1542a458/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-python-x402-facilitator-bazaar/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_2735fccb/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-python-x402-facilitator/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_30a4884f/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-typescript-x402-client/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_b4213d1b/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-typescript-x402-facilitator/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_3f98ecd8/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-typescript-x402-paywall/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_edc38ab1/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/explain-algorand-x402-typescript/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_a8e656dc/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/teach-algorand-x402/SKILL.md) | ⭐ 20 | `development` |
| [Skill](development/1178-name-skill_1c23a658/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/use-typescript-x402-core-avm/SKILL.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_455dec80/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-python-x402-facilitator/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_33ea3f85/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-typescript-x402-client/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_d5f4d996/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-typescript-x402-facilitator/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_b3c00018/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-typescript-x402-nextjs/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_67c246ad/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-typescript-x402-paywall/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_d370bccd/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/explain-algorand-x402-python/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_821d4a5b/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/explain-algorand-x402-typescript/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_5e308754/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/teach-algorand-x402/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_e336e402/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/use-python-x402-core-avm/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Reference](development/828-reference_7d3016f2/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/use-typescript-x402-core-avm/references/REFERENCE.md) | ⭐ 20 | `development` |
| [Linting Root Cause Resolver](development/1699-linting-root-cause-resolver_c81d2991/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/holistic-linting/agents/linting-root-cause-resolver.md) | ⭐ 17 | `development` |
| [Context Refinement](development/1665-context-refinement_b4db7922/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/agents/context-refinement.md) | ⭐ 17 | `development` |
| [Skill](development/1530-description-skill_8071c541/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/commitlint/skills/commitlint/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1530-description-skill_83522e6b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/holistic-linting/skills/holistic-linting-orchestrator/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1530-description-skill_8dc55b7b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/holistic-linting/skills/holistic-linting/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1530-description-skill_1dae8dde/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/litellm/skills/litellm/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1530-description-skill_38fa6f00/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/llamafile/skills/llamafile/SKILL.md) | ⭐ 17 | `development` |
| [Skill](development/1530-description-skill_f704848d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/audit-skill-completeness/SKILL.md) | ⭐ 17 | `development` |
| [Multi Tenant Design](development/multi-tenant-design_4fdd0f01/) | [volcengine/OpenViking](https://raw.githubusercontent.com/volcengine/OpenViking/main/docs/design/multi-tenant-design.md) | ⭐ 1.1k | `development` |

### Development/Devops (19 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Known Issues](development/devops/373-known-issues_97ef6daa/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/known-issues.md) | ⭐ 65 | `development` |
| [Batch Eval Quickstart](development/devops/364-batch_eval_quickstart_19b916c5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/BATCH_EVAL_QUICKSTART.md) | ⭐ 18 | `development` |
| [Troubleshooting](development/devops/093-troubleshooting_a88da6d4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/TROUBLESHOOTING.md) | ⭐ 18 | `development` |
| [Env Configuration](development/devops/091-env_configuration_2db3de1d/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/env_configuration.md) | ⭐ 4.0k | `development` |
| [Troubleshooting Openai Api Key](development/devops/092-troubleshooting-openai-api-key_c455d79c/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/troubleshooting-openai-api-key.md) | ⭐ 4.0k | `development` |
| [Roadmap](development/devops/097-roadmap_aa63c2f2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/roadmap.md) | ⭐ 3.3k | `development` |
| [Skill](development/devops/014-name-skill_ae24568a/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/semgrep/SKILL.md) | ⭐ 2.6k | `development` |
| [Environment Setup](development/devops/200-environment-setup_d0710e36/) | [adenhq/hive](https://raw.githubusercontent.com/adenhq/hive/main/docs/environment-setup.md) | 🔥 7.1k | `development` |
| [Roadmap](development/devops/097-roadmap_62290bac/) | [adenhq/hive](https://raw.githubusercontent.com/adenhq/hive/main/docs/roadmap.md) | 🔥 7.1k | `development` |
| [Wordlists](development/devops/364-wordlists_f64624a5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/ffuf-web-fuzzing/skills/ffuf-web-fuzzing/references/wordlists.md) | ⭐ 56 | `development` |
| [Unauthorized Access](development/devops/365-unauthorized-access_8bef681c/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/unauthorized-access.md) | ⭐ 56 | `development` |
| [Misconfig Checklist](development/devops/366-misconfig-checklist_2aa6973b/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/misconfig-checklist.md) | ⭐ 56 | `development` |
| [Containerize](development/devops/367-containerize_f635c1bc/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/commands/containerize.md) | ⭐ 21 | `development` |
| [Modern Extensions](development/devops/368-modern-extensions_74d49234/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/fifteen-factor-app/references/modern-extensions.md) | ⭐ 21 | `development` |
| [Original Factors](development/devops/369-original-factors_a71bbfe5/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/fifteen-factor-app/references/original-factors.md) | ⭐ 21 | `development` |
| [Api](development/devops/165-api_bb11980e/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-tagmanager/skills/google-tagmanager/references/api.md) | ⭐ 21 | `development` |
| [Advanced Patterns](development/devops/370-advanced-patterns_6c3be033/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/skills/skill-mastery/references/advanced-patterns.md) | ⭐ 21 | `development` |
| [03 Des Cost Estimate](development/devops/210-03-des-cost-estimate_d7a9199c/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/agent-output/pci-dss-gw/03-des-cost-estimate.md) | ⭐ 63 | `development` |
| [Skill](development/devops/014-name-skill_3f42193e/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/azure-defaults/SKILL.md) | ⭐ 63 | `development` |

### Development/Testing (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Batch 20260211 221333 Aggregate](development/testing/085-batch_20260211_221333_aggregate_2012d195/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221333_aggregate.md) | ⭐ 18 | `development` |
| [Batch 20260211 221435 Aggregate](development/testing/086-batch_20260211_221435_aggregate_09b4ca91/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221435_aggregate.md) | ⭐ 18 | `development` |
| [Usage Guide](development/testing/017-usage-guide_324bef25/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/temporal-genomics/periodicity-detection/usage-guide.md) | ⭐ 233 | `development` |
| [Formula Verification](development/testing/085-formula-verification_19feb703/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/documents/skills/xlsx/references/formula-verification.md) | ⭐ 21 | `development` |
| [Version Numbering](development/testing/086-version-numbering_619d38b1/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/version-numbering.md) | ⭐ 21 | `development` |

### Development/Tools (31 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/tools/002-name-skill_727accdb/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/SKILL.md) | ⭐ 65 | `development` |
| [Official Sources](development/tools/320-official-sources_aab46d33/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/official-sources.md) | ⭐ 65 | `development` |
| [Security Analysis](development/tools/321-security_analysis_fdc4ce06/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/security_analysis.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_9f88cf3e/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/authoring-dags/SKILL.md) | ⭐ 208 | `development` |
| [Skill](development/tools/002-name-skill_04e2c484/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/edna-metabarcoding/SKILL.md) | ⭐ 233 | `development` |
| [Usage Guide](development/tools/069-usage-guide_abeec102/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/edna-metabarcoding/usage-guide.md) | ⭐ 233 | `development` |
| [Skill](development/tools/002-name-skill_e212fd0c/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/ecological-genomics/landscape-genomics/SKILL.md) | ⭐ 233 | `development` |
| [Desktop](development/tools/322-desktop_17f69f89/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/desktop.md) | ⭐ 51 | `development` |
| [Release Description](development/tools/322-release_description_d0311e45/) | [michaelbeijer/Supervertaler](https://raw.githubusercontent.com/michaelbeijer/Supervertaler/main/RELEASE_DESCRIPTION.md) | ⭐ 22 | `development` |
| [Contributing Lint Setup](development/tools/321-contributing-lint-setup_2918e77d/) | [adenhq/hive](https://raw.githubusercontent.com/adenhq/hive/main/docs/contributing-lint-setup.md) | 🔥 7.1k | `development` |
| [Skill](development/tools/002-name-skill_d2a89379/) | [adenhq/hive](https://raw.githubusercontent.com/adenhq/hive/main/.claude/skills/hive-debugger/SKILL.md) | 🔥 7.1k | `development` |
| [Claude](development/tools/017-claude_c795e188/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/CLAUDE.md) | ⭐ 56 | `development` |
| [Skill](development/tools/002-name-skill_926f799f/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/ffuf-web-fuzzing/skills/ffuf-web-fuzzing/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/002-name-skill_92250a36/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-jupyter-notebook/skills/openai-jupyter-notebook/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/002-name-skill_df491cf9/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-netlify-deploy/skills/openai-netlify-deploy/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/002-name-skill_69eee0a1/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-screenshot/skills/openai-screenshot/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/002-name-skill_920d7880/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-yeet/skills/openai-yeet/SKILL.md) | ⭐ 56 | `development` |
| [Cli](development/tools/263-cli_a08ee059/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-playwright/skills/openai-playwright/references/cli.md) | ⭐ 56 | `development` |
| [Events Migration Guide](development/tools/319-events_migration_guide_a1d86b29/) | [openlit/openlit](https://raw.githubusercontent.com/openlit/openlit/main/contributors/EVENTS_MIGRATION_GUIDE.md) | ⭐ 2.2k | `development` |
| [Gemini Cli Headless](development/tools/322-gemini-cli-headless_467b26d2/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/gemini-cli-headless/agents/gemini-cli-headless.md) | ⭐ 21 | `development` |
| [Create Skill Ultimate](development/tools/323-create-skill-ultimate_3ac6b128/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/commands/create-skill-ultimate.md) | ⭐ 21 | `development` |
| [Skill](development/tools/002-name-skill_5d89e32d/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/skills/skill-mastery/SKILL.md) | ⭐ 21 | `development` |
| [Troubleshooting](development/tools/205-troubleshooting_16180695/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/developer/skills/devtools/references/troubleshooting.md) | ⭐ 21 | `development` |
| [Bigquery](development/tools/324-bigquery_fadc3efd/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/bigquery.md) | ⭐ 21 | `development` |
| [Api Sync](development/tools/325-api-sync_f3b574db/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/tampermonkey/skills/tampermonkey/references/api-sync.md) | ⭐ 21 | `development` |
| [Api Reference](development/tools/073-api-reference_9899dff9/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/ultimate-skill-creator/skills/skill-mastery/references/api-reference.md) | ⭐ 21 | `development` |
| [Agent Skills.Instructions](development/tools/228-agent-skillsinstructions_95976920/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/agent-skills.instructions.md) | ⭐ 63 | `development` |
| [Docs.Instructions](development/tools/229-docsinstructions_85f2b3f7/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/docs.instructions.md) | ⭐ 63 | `development` |
| [Skill](development/tools/002-name-skill_c88e2aec/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/github-operations/SKILL.md) | ⭐ 63 | `development` |
| [Plan Bicep.Prompt](development/tools/321-plan-bicepprompt_6c8677c3/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/prompts/plan-bicep.prompt.md) | ⭐ 65 | `development` |
| [Reference](development/tools/074-reference_4115642e/) | [algorand-devrel/algorand-agent-skills](https://raw.githubusercontent.com/algorand-devrel/algorand-agent-skills/main/skills/create-python-x402-facilitator-bazaar/references/REFERENCE.md) | ⭐ 20 | `development` |

### Investment (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Ai Agent Summary](investment/049-ai_agent_summary_eb15178c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AI_AGENT_SUMMARY.md) | ⭐ 18 | `investment` |
| [Failure Analysis Batch 20260211 213131](investment/050-failure_analysis_batch_20260211_213131_17df2793/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_213131.md) | ⭐ 18 | `investment` |
| [Ultra Think](investment/049-ultra-think_fa8c9d32/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/claudecode/commands/ultra-think.md) | ⭐ 21 | `investment` |
| [User Tracking](investment/050-user-tracking_5a945df2/) | [henkisdabro/wookstar-claude-plugins](https://raw.githubusercontent.com/henkisdabro/wookstar-claude-plugins/main/plugins/google-analytics/skills/google-analytics/references/user-tracking.md) | ⭐ 21 | `investment` |

### Productivity (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Growth Pm](productivity/173-growth-pm_48f76195/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/growth-pm.md) | ⭐ 18 | `productivity` |
| [Prompts](productivity/174-prompts_5155f536/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/param_forge_ref/param_forge/docs/prompts.md) | ⭐ 51 | `productivity` |
| [Patterns](productivity/160-patterns_01ee8bbb/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/humanizer/skills/humanizer/references/patterns.md) | ⭐ 56 | `productivity` |

### Research (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Showcase](research/258-showcase_cafd0d13/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/SHOWCASE.md) | ⭐ 18 | `research` |
| [Directory](research/259-directory_ffc70b99/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/directory.md) | ⭐ 18 | `research` |
| [Researcher](research/260-researcher_edc4a804/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/researcher.md) | ⭐ 18 | `research` |
| [Batch 20260211 221758 Aggregate](research/261-batch_20260211_221758_aggregate_8cac16a4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221758_aggregate.md) | ⭐ 18 | `research` |
| [Time Savings Evidence](research/225-time-savings-evidence_a6749979/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/time-savings-evidence.md) | ⭐ 63 | `research` |

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

*Last updated: 2026-02-13 10:21:39 UTC*
*Automatically maintained by SkillFlow*
