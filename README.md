# X-Skills

A curated collection of **172 AI-powered skills** organized into 14 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (1 skill)
- **Automation/Workflow** (12 skills)
- **Commercial** (21 skills)
- **Communication** (4 skills)
- **Content Creation** (19 skills)
- **Daily Assistant** (9 skills)
- **Data Analysis** (26 skills)
- **Development** (33 skills)
- **Development/Devops** (13 skills)
- **Development/Testing** (3 skills)
- **Development/Tools** (21 skills)
- **Investment** (2 skills)
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


### Automation/Scripting (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Param Forge Reference](automation/scripting/100-param_forge_reference_35994bfd/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/param_forge_reference.md) | ⭐ 51 | `automation` |

### Automation/Workflow (12 skills)

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
| [Skill](automation/workflow/name-skill_1b2668f1/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-develop-web-game/skills/openai-develop-web-game/SKILL.md) | ⭐ 56 | `automation` |
| [Skill](automation/workflow/name-skill_5ed000bc/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/ghidra-headless/skills/ghidra-headless/SKILL.md) | ⭐ 56 | `automation` |

### Commercial (21 skills)

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
| [Bank Penetration](commercial/bank-penetration_c2e5ee98/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/bank-penetration.md) | ⭐ 56 | `commercial` |
| [Logic Flaws](commercial/logic-flaws_d3713f41/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/logic-flaws.md) | ⭐ 56 | `commercial` |
| [Logic Flaws Checklist](commercial/logic-flaws-checklist_c6b46334/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/logic-flaws-checklist.md) | ⭐ 56 | `commercial` |

### Communication (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Syntax Reference](communication/113-syntax-reference_a16eba1d/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/syntax-reference.md) | ⭐ 65 | `communication` |
| [Project Context](communication/255-project_context_5cc91b78/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/param_forge_ref/param_forge/PROJECT_CONTEXT.md) | ⭐ 51 | `communication` |
| [Readme Cn](communication/256-readme_cn_a75d0219/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/README_CN.md) | ⭐ 57 | `communication` |
| [Weak Password Checklist](communication/weak-password-checklist_0d3963ca/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/weak-password-checklist.md) | ⭐ 56 | `communication` |

### Content Creation (19 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Quick Wins](content-creation/353-quick_wins_29f2473a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/QUICK_WINS.md) | ⭐ 18 | `content creation` |
| [Skill Chains](content-creation/354-skill_chains_312b4d53/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/SKILL_CHAINS.md) | ⭐ 18 | `content creation` |
| [Skills Catalog](content-creation/355-skills_catalog_b303f24f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/SKILLS_CATALOG.md) | ⭐ 18 | `content creation` |
| [Ai Agent Composability Analysis](content-creation/356-ai_agent_composability_analysis_8ca26161/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AI_AGENT_COMPOSABILITY_ANALYSIS.md) | ⭐ 18 | `content creation` |
| [Skill](content-creation/049-name-skill_42cdd582/) | [cat-xierluo/legal-skills](https://raw.githubusercontent.com/cat-xierluo/legal-skills/main/skills/universal-media-downloader/SKILL.md) | ⭐ 18 | `content creation` |
| [Spawn V0 Local Only](content-creation/357-spawn_v0_local_only_1cb9fb27/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/prds/spawn_v0_local_only.md) | ⭐ 51 | `content creation` |
| [Ethics](content-creation/372-ethics_dbde6459/) | [opendatahub-io/ai-helpers](https://raw.githubusercontent.com/opendatahub-io/ai-helpers/main/ETHICS.md) | ⭐ 13 | `content creation` |
| [Plan](content-creation/plan_860e9528/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/commands/plan.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/name-skill_82868fe3/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/humanizer/skills/humanizer/SKILL.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/name-skill_d1631d58/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/last30days/skills/last30days/SKILL.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/name-skill_76a45975/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/skills/planning-with-files/SKILL.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/name-skill_8a2e5b56/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/x-research/skills/x-research/SKILL.md) | ⭐ 56 | `content creation` |
| [Principles](content-creation/principles_c2c61108/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/skills/planning-with-files/references/principles.md) | ⭐ 56 | `content creation` |
| [Skill Lifecycle](content-creation/skill-lifecycle_119157d0/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/skill-extractor/skills/skill-extractor/references/skill-lifecycle.md) | ⭐ 56 | `content creation` |
| [Command Execution](content-creation/command-execution_ffdc6b02/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/command-execution.md) | ⭐ 56 | `content creation` |
| [Xss](content-creation/xss_d1bc8bd4/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/xss.md) | ⭐ 56 | `content creation` |
| [Command Execution Checklist](content-creation/command-execution-checklist_39b874c0/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/command-execution-checklist.md) | ⭐ 56 | `content creation` |
| [Ssrf Checklist](content-creation/ssrf-checklist_113f13e6/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/ssrf-checklist.md) | ⭐ 56 | `content creation` |
| [Skill](content-creation/name-skill_0b19f29a/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-cloudflare-deploy/skills/openai-cloudflare-deploy/SKILL.md) | ⭐ 56 | `content creation` |

### Daily Assistant (9 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Research Notes](daily-assistant/267-research_notes_3883eb9d/) | [taylorsatula/mira-OSS](https://raw.githubusercontent.com/taylorsatula/mira-OSS/main/RESEARCH_NOTES.md) | ⭐ 389 | `daily assistant` |
| [Copilot Instructions](daily-assistant/266-copilot-instructions_1423668e/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/.github/copilot-instructions.md) | ⭐ 51 | `daily assistant` |
| [Brood Aip.Instructions](daily-assistant/267-brood-aipinstructions_be0c3008/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/.github/instructions/brood-aip.instructions.md) | ⭐ 51 | `daily assistant` |
| [Brood Aip](daily-assistant/268-brood-aip_be4ed266/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/.windsurf/rules/brood-aip.md) | ⭐ 51 | `daily assistant` |
| [Status](daily-assistant/status_2e342fb5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/commands/status.md) | ⭐ 56 | `daily assistant` |
| [Skill](daily-assistant/name-skill_47ec12a5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-playwright/skills/openai-playwright/SKILL.md) | ⭐ 56 | `daily assistant` |
| [Skill](daily-assistant/name-skill_f9a1f323/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-security-threat-model/skills/openai-security-threat-model/SKILL.md) | ⭐ 56 | `daily assistant` |
| [Skill](daily-assistant/name-skill_92eb7fef/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-sentry/skills/openai-sentry/SKILL.md) | ⭐ 56 | `daily assistant` |
| [Examples](daily-assistant/examples_c00cd5b6/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/skills/planning-with-files/references/examples.md) | ⭐ 56 | `daily assistant` |

### Data Analysis (26 skills)

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
| [Import Skill](data-analysis/import-skill_0f090701/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/.claude/commands/import-skill.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/name-skill_5f6df240/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-doc/skills/openai-doc/SKILL.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/name-skill_ed5b6f7e/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-pdf/skills/openai-pdf/SKILL.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/name-skill_b03aba97/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-security-ownership-map/skills/openai-security-ownership-map/SKILL.md) | ⭐ 56 | `data analysis` |
| [Skill](data-analysis/name-skill_94733cad/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-spreadsheet/skills/openai-spreadsheet/SKILL.md) | ⭐ 56 | `data analysis` |
| [Review Plugin](data-analysis/review-plugin_4097752f/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/.claude/commands/review-plugin.md) | ⭐ 56 | `data analysis` |

### Development (33 skills)

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
| [Skill](development/name-skill_bc6bbf16/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-security-best-practices/skills/openai-security-best-practices/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/name-skill_bbfc4111/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/skill-extractor/skills/skill-extractor/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/name-skill_0d3d1f7c/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/SKILL.md) | ⭐ 56 | `development` |
| [File Upload](development/file-upload_b013ed1e/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/file-upload.md) | ⭐ 56 | `development` |
| [Info Disclosure](development/info-disclosure_d4862aba/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/info-disclosure.md) | ⭐ 56 | `development` |
| [Path Traversal](development/path-traversal_c73eea45/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/path-traversal.md) | ⭐ 56 | `development` |
| [Sql Injection](development/sql-injection_cc9c615c/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/sql-injection.md) | ⭐ 56 | `development` |
| [Telecom Penetration](development/telecom-penetration_5c723184/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/telecom-penetration.md) | ⭐ 56 | `development` |
| [X Api](development/x-api_b97fa218/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/x-research/skills/x-research/references/x-api.md) | ⭐ 56 | `development` |
| [Csrf Checklist](development/csrf-checklist_e821ad41/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/csrf-checklist.md) | ⭐ 56 | `development` |
| [File Upload Checklist](development/file-upload-checklist_95e05da5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/file-upload-checklist.md) | ⭐ 56 | `development` |
| [Info Disclosure Checklist](development/info-disclosure-checklist_aac2c34e/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/info-disclosure-checklist.md) | ⭐ 56 | `development` |
| [Path Traversal Checklist](development/path-traversal-checklist_e420e992/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/path-traversal-checklist.md) | ⭐ 56 | `development` |
| [Rce Checklist](development/rce-checklist_60e3a56a/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/rce-checklist.md) | ⭐ 56 | `development` |
| [Xxe Checklist](development/xxe-checklist_0521eb1b/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/xxe-checklist.md) | ⭐ 56 | `development` |
| [Notebook Structure](development/notebook-structure_a5fb0429/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-jupyter-notebook/skills/openai-jupyter-notebook/references/notebook-structure.md) | ⭐ 56 | `development` |
| [Templates](development/templates_7a6ccc51/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/planning-with-files/skills/planning-with-files/references/templates.md) | ⭐ 56 | `development` |
| [Quality Guide](development/quality-guide_327ccc8a/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/skill-extractor/skills/skill-extractor/references/quality-guide.md) | ⭐ 56 | `development` |
| [Sql Injection Checklist](development/sql-injection-checklist_d12eeef4/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/sql-injection-checklist.md) | ⭐ 56 | `development` |
| [Xss Checklist](development/xss-checklist_b09898a4/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/xss-checklist.md) | ⭐ 56 | `development` |

### Development/Devops (13 skills)

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
| [Wordlists](development/devops/wordlists_f64624a5/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/ffuf-web-fuzzing/skills/ffuf-web-fuzzing/references/wordlists.md) | ⭐ 56 | `development` |
| [Unauthorized Access](development/devops/unauthorized-access_8bef681c/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/unauthorized-access.md) | ⭐ 56 | `development` |
| [Misconfig Checklist](development/devops/misconfig-checklist_2aa6973b/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/misconfig-checklist.md) | ⭐ 56 | `development` |
| [Unauthorized Access Checklist](development/devops/unauthorized-access-checklist_9f1f7c60/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/wooyun-legacy/skills/wooyun-legacy/references/checklists/unauthorized-access-checklist.md) | ⭐ 56 | `development` |

### Development/Testing (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Batch 20260211 221333 Aggregate](development/testing/085-batch_20260211_221333_aggregate_2012d195/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221333_aggregate.md) | ⭐ 18 | `development` |
| [Batch 20260211 221435 Aggregate](development/testing/086-batch_20260211_221435_aggregate_09b4ca91/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221435_aggregate.md) | ⭐ 18 | `development` |
| [Usage Guide](development/testing/017-usage-guide_324bef25/) | [GPTomics/bioSkills](https://raw.githubusercontent.com/GPTomics/bioSkills/main/temporal-genomics/periodicity-detection/usage-guide.md) | ⭐ 233 | `development` |

### Development/Tools (21 skills)

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
| [Claude](development/tools/claude_c795e188/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/CLAUDE.md) | ⭐ 56 | `development` |
| [Skill](development/tools/name-skill_926f799f/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/ffuf-web-fuzzing/skills/ffuf-web-fuzzing/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/name-skill_92250a36/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-jupyter-notebook/skills/openai-jupyter-notebook/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/name-skill_df491cf9/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-netlify-deploy/skills/openai-netlify-deploy/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/name-skill_69eee0a1/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-screenshot/skills/openai-screenshot/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/name-skill_920d7880/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-yeet/skills/openai-yeet/SKILL.md) | ⭐ 56 | `development` |
| [Cli](development/tools/cli_a08ee059/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-playwright/skills/openai-playwright/references/cli.md) | ⭐ 56 | `development` |
| [Skill](development/tools/name-skill_fcd5e5cf/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-gh-address-comments/skills/openai-gh-address-comments/SKILL.md) | ⭐ 56 | `development` |
| [Skill](development/tools/name-skill_9fd96efd/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/openai-gh-fix-ci/skills/openai-gh-fix-ci/SKILL.md) | ⭐ 56 | `development` |
| [Request Templates](development/tools/request-templates_66fd8938/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/ffuf-web-fuzzing/skills/ffuf-web-fuzzing/references/request-templates.md) | ⭐ 56 | `development` |

### Investment (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Ai Agent Summary](investment/049-ai_agent_summary_eb15178c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AI_AGENT_SUMMARY.md) | ⭐ 18 | `investment` |
| [Failure Analysis Batch 20260211 213131](investment/050-failure_analysis_batch_20260211_213131_17df2793/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_213131.md) | ⭐ 18 | `investment` |

### Productivity (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Growth Pm](productivity/173-growth-pm_48f76195/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/growth-pm.md) | ⭐ 18 | `productivity` |
| [Prompts](productivity/174-prompts_5155f536/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/param_forge_ref/param_forge/docs/prompts.md) | ⭐ 51 | `productivity` |
| [Patterns](productivity/patterns_01ee8bbb/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/humanizer/skills/humanizer/references/patterns.md) | ⭐ 56 | `productivity` |

### Research (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Showcase](research/258-showcase_cafd0d13/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/SHOWCASE.md) | ⭐ 18 | `research` |
| [Directory](research/259-directory_ffc70b99/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/directory.md) | ⭐ 18 | `research` |
| [Researcher](research/260-researcher_edc4a804/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/researcher.md) | ⭐ 18 | `research` |
| [Batch 20260211 221758 Aggregate](research/261-batch_20260211_221758_aggregate_8cac16a4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221758_aggregate.md) | ⭐ 18 | `research` |
| [Skill Template](research/skill-template_7551123b/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/skill-extractor/skills/skill-extractor/references/skill-template.md) | ⭐ 56 | `research` |

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

*Last updated: 2026-02-13 03:22:48 UTC*
*Automatically maintained by SkillFlow*
