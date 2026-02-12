# X-Skills

A curated collection of **361 AI-powered skills** organized into 14 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (36 skills)
- **Automation/Workflow** (10 skills)
- **Commercial** (252 skills)
- **Communication** (1 skill)
- **Content Creation** (7 skills)
- **Daily Assistant** (1 skill)
- **Data Analysis** (11 skills)
- **Development** (22 skills)
- **Development/Devops** (3 skills)
- **Development/Testing** (3 skills)
- **Development/Tools** (4 skills)
- **Investment** (5 skills)
- **Productivity** (1 skill)
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


### Automation/Scripting (36 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Eval Harness Rollout Complete](automation/scripting/eval_harness_rollout_complete_cba15e32/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/EVAL_HARNESS_ROLLOUT_COMPLETE.md) | ⭐ 18 | `automation` |
| [Phase 0 Implementation Summary](automation/scripting/phase_0_implementation_summary_01ba15da/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/PHASE_0_IMPLEMENTATION_SUMMARY.md) | ⭐ 18 | `automation` |
| [Phase 1 Pilot Complete](automation/scripting/phase_1_pilot_complete_ef0fedf8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/PHASE_1_PILOT_COMPLETE.md) | ⭐ 18 | `automation` |
| [Phase 1 Quick Reference](automation/scripting/phase_1_quick_reference_4f5cd4c5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/PHASE_1_QUICK_REFERENCE.md) | ⭐ 18 | `automation` |
| [Batch 20260211 171245 Aggregate](automation/scripting/batch_20260211_171245_aggregate_eb9e62bd/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_171245_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 171548 Aggregate](automation/scripting/batch_20260211_171548_aggregate_12c83d51/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_171548_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 212313 Aggregate](automation/scripting/batch_20260211_212313_aggregate_2c131433/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_212313_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 212318 Aggregate](automation/scripting/batch_20260211_212318_aggregate_183875df/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_212318_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 212326 Aggregate](automation/scripting/batch_20260211_212326_aggregate_fab58b68/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_212326_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 213131 Aggregate](automation/scripting/batch_20260211_213131_aggregate_9b275057/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_213131_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 213340 Aggregate](automation/scripting/batch_20260211_213340_aggregate_97305460/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_213340_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 230434 Aggregate](automation/scripting/batch_20260211_230434_aggregate_7ea48412/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_230434_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 230439 Aggregate](automation/scripting/batch_20260211_230439_aggregate_ff7d0c8e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_230439_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 230628 Aggregate](automation/scripting/batch_20260211_230628_aggregate_d251bd84/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_230628_aggregate.md) | ⭐ 18 | `automation` |
| [Batch 20260211 230749 Aggregate](automation/scripting/batch_20260211_230749_aggregate_e655ecf1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_230749_aggregate.md) | ⭐ 18 | `automation` |
| [Ai Autonomous Outreach Eval](automation/scripting/ai_autonomous_outreach_eval_75266df5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_autonomous_outreach_eval.md) | ⭐ 18 | `automation` |
| [Ai Conversation Intelligence Eval](automation/scripting/ai_conversation_intelligence_eval_311b3c80/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_conversation_intelligence_eval.md) | ⭐ 18 | `automation` |
| [Ai Personalization Engine Eval](automation/scripting/ai_personalization_engine_eval_ce163d4a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_personalization_engine_eval.md) | ⭐ 18 | `automation` |
| [Ai Predictive Lead Scoring Eval](automation/scripting/ai_predictive_lead_scoring_eval_87b4ceab/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_predictive_lead_scoring_eval.md) | ⭐ 18 | `automation` |
| [Cs Churn Prediction Eval](automation/scripting/cs_churn_prediction_eval_88b9b085/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_churn_prediction_eval.md) | ⭐ 18 | `automation` |
| [Cs Expansion Playbook Eval](automation/scripting/cs_expansion_playbook_eval_7825dd7a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_expansion_playbook_eval.md) | ⭐ 18 | `automation` |
| [Cs Health Scoring Eval](automation/scripting/cs_health_scoring_eval_ca7a935a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_health_scoring_eval.md) | ⭐ 18 | `automation` |
| [Cs Nps Followup Eval](automation/scripting/cs_nps_followup_eval_8ad2e9b7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_nps_followup_eval.md) | ⭐ 18 | `automation` |
| [Cs Renewal Orchestration Eval](automation/scripting/cs_renewal_orchestration_eval_22ac9093/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_renewal_orchestration_eval.md) | ⭐ 18 | `automation` |
| [Cs Value Realization Eval](automation/scripting/cs_value_realization_eval_f58a2ee9/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_value_realization_eval.md) | ⭐ 18 | `automation` |
| [Mon Dunning Automation Eval](automation/scripting/mon_dunning_automation_eval_a7644535/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_dunning_automation_eval.md) | ⭐ 18 | `automation` |
| [Mon Limit Notification Eval](automation/scripting/mon_limit_notification_eval_dd5fcc3b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_limit_notification_eval.md) | ⭐ 18 | `automation` |
| [Mon Pricing Optimization Eval](automation/scripting/mon_pricing_optimization_eval_3d4164a2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_pricing_optimization_eval.md) | ⭐ 18 | `automation` |
| [Mon Upgrade Trigger Eval](automation/scripting/mon_upgrade_trigger_eval_20d18ff6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_upgrade_trigger_eval.md) | ⭐ 18 | `automation` |
| [Mon Usage Metering Eval](automation/scripting/mon_usage_metering_eval_394b5e9c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_usage_metering_eval.md) | ⭐ 18 | `automation` |
| [Prodops Feature Adoption Eval](automation/scripting/prodops_feature_adoption_eval_40a96c08/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_feature_adoption_eval.md) | ⭐ 18 | `automation` |
| [Prodops Feedback Synthesis Eval](automation/scripting/prodops_feedback_synthesis_eval_f7976462/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_feedback_synthesis_eval.md) | ⭐ 18 | `automation` |
| [Prodops Roadmap Alignment Eval](automation/scripting/prodops_roadmap_alignment_eval_8b74c070/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_roadmap_alignment_eval.md) | ⭐ 18 | `automation` |
| [Prodops Voc Aggregation Eval](automation/scripting/prodops_voc_aggregation_eval_5547733f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_voc_aggregation_eval.md) | ⭐ 18 | `automation` |
| [Support Macros Optimizer Eval](automation/scripting/support_macros_optimizer_eval_3ded962b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_macros_optimizer_eval.md) | ⭐ 18 | `automation` |
| [Skill](automation/scripting/name-skill_28bd1f73/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/checking-freshness/SKILL.md) | ⭐ 208 | `automation` |

### Automation/Workflow (10 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/002-name-skill_ce4bb9a9/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-deploy/SKILL.md) | ⭐ 65 | `automation` |
| [Cli Guide](automation/workflow/137-cli-guide_49550e7a/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/docs/cli-guide.md) | ⭐ 65 | `automation` |
| [Announcement](automation/workflow/announcement_23270937/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/ANNOUNCEMENT.md) | ⭐ 18 | `automation` |
| [Architecture](automation/workflow/architecture_bef9ba1a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/ARCHITECTURE.md) | ⭐ 18 | `automation` |
| [Phase 2 Business Critical Complete](automation/workflow/phase_2_business_critical_complete_fac6699d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/PHASE_2_BUSINESS_CRITICAL_COMPLETE.md) | ⭐ 18 | `automation` |
| [Legal](automation/workflow/legal_b6e495b7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/legal.md) | ⭐ 18 | `automation` |
| [Cfo](automation/workflow/cfo_6dfc98f1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/cfo.md) | ⭐ 18 | `automation` |
| [Sales Leader](automation/workflow/sales-leader_ae635cbe/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/sales-leader.md) | ⭐ 18 | `automation` |
| [Ai Workflow Automator Eval](automation/workflow/ai_workflow_automator_eval_d8780b5a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_workflow_automator_eval.md) | ⭐ 18 | `automation` |
| [Skill](automation/workflow/name-skill_029a02bb/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/tracing-downstream-lineage/SKILL.md) | ⭐ 208 | `automation` |

### Commercial (252 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/210-name-skill_10221f10/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/manuscript-ingest/SKILL.md) | ⭐ 197 | `commercial` |
| [Security Policy](commercial/security_policy_f8dea503/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/SECURITY_POLICY.md) | ⭐ 18 | `commercial` |
| [Customer Success](commercial/customer_success_85216f92/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/customer_success.md) | ⭐ 18 | `commercial` |
| [Executive](commercial/executive_e9668604/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/executive.md) | ⭐ 18 | `commercial` |
| [Finance](commercial/finance_8fdccf44/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/finance.md) | ⭐ 18 | `commercial` |
| [Marketing](commercial/marketing_499062ff/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/marketing.md) | ⭐ 18 | `commercial` |
| [Agent Coverage](commercial/agent_coverage_affae21f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AGENT_COVERAGE.md) | ⭐ 18 | `commercial` |
| [Batch 20260211 212328 Aggregate](commercial/batch_20260211_212328_aggregate_b193fed8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_212328_aggregate.md) | ⭐ 18 | `commercial` |
| [Failure Analysis Batch 20260211 212328](commercial/failure_analysis_batch_20260211_212328_ebcc0636/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_212328.md) | ⭐ 18 | `commercial` |
| [Cs Risk Mitigation Playbook Eval](commercial/cs_risk_mitigation_playbook_eval_9849a6ec/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_risk_mitigation_playbook_eval.md) | ⭐ 18 | `commercial` |
| [Devex Rate Limit Advisor Eval](commercial/devex_rate_limit_advisor_eval_e60ff981/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_rate_limit_advisor_eval.md) | ⭐ 18 | `commercial` |
| [Mon Commitment Tracker Eval](commercial/mon_commitment_tracker_eval_be78232d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_commitment_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Mon Discount Optimizer Eval](commercial/mon_discount_optimizer_eval_5ef6bd58/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_discount_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Plg Progressive Disclosure Eval](commercial/plg_progressive_disclosure_eval_829833e5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_progressive_disclosure_eval.md) | ⭐ 18 | `commercial` |
| [Revops Commit Accuracy Eval](commercial/revops_commit_accuracy_eval_97f0580c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_commit_accuracy_eval.md) | ⭐ 18 | `commercial` |
| [Vcf Value Discovery Eval](commercial/vcf_value_discovery_eval_a49e0e39/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/vcf_value_discovery_eval.md) | ⭐ 18 | `commercial` |
| [Dedupe Summary](commercial/dedupe_summary_f8258a18/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/dedupe_summary.md) | ⭐ 18 | `commercial` |
| [Product](commercial/product_91710cb6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/product.md) | ⭐ 18 | `commercial` |
| [Batch 20260211 212412 Aggregate](commercial/batch_20260211_212412_aggregate_a26ab2ca/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_212412_aggregate.md) | ⭐ 18 | `commercial` |
| [Batch 20260211 231557 Aggregate](commercial/batch_20260211_231557_aggregate_85f51bf1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_231557_aggregate.md) | ⭐ 18 | `commercial` |
| [Batch 20260212 000513 Aggregate](commercial/batch_20260212_000513_aggregate_07abe36f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260212_000513_aggregate.md) | ⭐ 18 | `commercial` |
| [Failure Analysis Batch 20260211 212412](commercial/failure_analysis_batch_20260211_212412_8a363b35/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_212412.md) | ⭐ 18 | `commercial` |
| [Ai Agent Performance Eval](commercial/ai_agent_performance_eval_4d0a76a9/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_agent_performance_eval.md) | ⭐ 18 | `commercial` |
| [Ai Anomaly Detector Eval](commercial/ai_anomaly_detector_eval_d4c4794a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_anomaly_detector_eval.md) | ⭐ 18 | `commercial` |
| [Ai Chatbot Orchestrator Eval](commercial/ai_chatbot_orchestrator_eval_f545c32d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_chatbot_orchestrator_eval.md) | ⭐ 18 | `commercial` |
| [Ai Content Generator Eval](commercial/ai_content_generator_eval_b1204745/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_content_generator_eval.md) | ⭐ 18 | `commercial` |
| [Ai Data Quality Agent Eval](commercial/ai_data_quality_agent_eval_a966ec67/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_data_quality_agent_eval.md) | ⭐ 18 | `commercial` |
| [Ai Document Parser Eval](commercial/ai_document_parser_eval_27e4b740/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_document_parser_eval.md) | ⭐ 18 | `commercial` |
| [Ai Intent Classifier Eval](commercial/ai_intent_classifier_eval_a749440a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_intent_classifier_eval.md) | ⭐ 18 | `commercial` |
| [Ai Knowledge Synthesizer Eval](commercial/ai_knowledge_synthesizer_eval_907f87fb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_knowledge_synthesizer_eval.md) | ⭐ 18 | `commercial` |
| [Ai Meeting Intelligence Eval](commercial/ai_meeting_intelligence_eval_ea57d034/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_meeting_intelligence_eval.md) | ⭐ 18 | `commercial` |
| [Ai Response Suggester Eval](commercial/ai_response_suggester_eval_9ce4537e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_response_suggester_eval.md) | ⭐ 18 | `commercial` |
| [Ai Summarization Agent Eval](commercial/ai_summarization_agent_eval_af775b26/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_summarization_agent_eval.md) | ⭐ 18 | `commercial` |
| [Ai Ticket Classifier Eval](commercial/ai_ticket_classifier_eval_da569ec9/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_ticket_classifier_eval.md) | ⭐ 18 | `commercial` |
| [Ai Translation Localizer Eval](commercial/ai_translation_localizer_eval_793d5289/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/ai_translation_localizer_eval.md) | ⭐ 18 | `commercial` |
| [Community Ama Orchestrator Eval](commercial/community_ama_orchestrator_eval_1a077bb8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_ama_orchestrator_eval.md) | ⭐ 18 | `commercial` |
| [Community Ambassador Program Eval](commercial/community_ambassador_program_eval_b03de7ce/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_ambassador_program_eval.md) | ⭐ 18 | `commercial` |
| [Community Case Study Finder Eval](commercial/community_case_study_finder_eval_42d530f1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_case_study_finder_eval.md) | ⭐ 18 | `commercial` |
| [Community Champion Identifier Eval](commercial/community_champion_identifier_eval_bd1063bb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_champion_identifier_eval.md) | ⭐ 18 | `commercial` |
| [Community Content Curator Eval](commercial/community_content_curator_eval_0c083f3b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_content_curator_eval.md) | ⭐ 18 | `commercial` |
| [Community Contributor Tracker Eval](commercial/community_contributor_tracker_eval_c76a93b1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_contributor_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Community Event Manager Eval](commercial/community_event_manager_eval_b3b7241a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_event_manager_eval.md) | ⭐ 18 | `commercial` |
| [Community Forum Moderator Eval](commercial/community_forum_moderator_eval_3413b3e3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_forum_moderator_eval.md) | ⭐ 18 | `commercial` |
| [Community Hackathon Manager Eval](commercial/community_hackathon_manager_eval_1e19051c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_hackathon_manager_eval.md) | ⭐ 18 | `commercial` |
| [Community Swag Manager Eval](commercial/community_swag_manager_eval_c1cfc4f1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_swag_manager_eval.md) | ⭐ 18 | `commercial` |
| [Community Ugc Curator Eval](commercial/community_ugc_curator_eval_b8db4f1c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_ugc_curator_eval.md) | ⭐ 18 | `commercial` |
| [Community User Group Manager Eval](commercial/community_user_group_manager_eval_8dac5425/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/community_user_group_manager_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Access Reviewer Eval](commercial/compliance_access_reviewer_eval_e74939c5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_access_reviewer_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Audit Preparer Eval](commercial/compliance_audit_preparer_eval_12177d86/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_audit_preparer_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Breach Response Eval](commercial/compliance_breach_response_eval_b2a852a6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_breach_response_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Consent Manager Eval](commercial/compliance_consent_manager_eval_082c9417/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_consent_manager_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Data Retention Eval](commercial/compliance_data_retention_eval_5e8e0dad/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_data_retention_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Gdpr Manager Eval](commercial/compliance_gdpr_manager_eval_05146667/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_gdpr_manager_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Pii Detector Eval](commercial/compliance_pii_detector_eval_7225b08f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_pii_detector_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Privacy Policy Eval](commercial/compliance_privacy_policy_eval_2c20c668/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_privacy_policy_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Soc2 Tracker Eval](commercial/compliance_soc2_tracker_eval_dbed6123/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_soc2_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Compliance Vendor Risk Eval](commercial/compliance_vendor_risk_eval_8ccf4b15/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/compliance_vendor_risk_eval.md) | ⭐ 18 | `commercial` |
| [Cs Adoption Score Eval](commercial/cs_adoption_score_eval_280f87f6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_adoption_score_eval.md) | ⭐ 18 | `commercial` |
| [Cs Advocacy Identifier Eval](commercial/cs_advocacy_identifier_eval_999db1ec/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_advocacy_identifier_eval.md) | ⭐ 18 | `commercial` |
| [Cs Contract Intelligence Eval](commercial/cs_contract_intelligence_eval_66cd5488/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_contract_intelligence_eval.md) | ⭐ 18 | `commercial` |
| [Cs Csql Generator Eval](commercial/cs_csql_generator_eval_515b83f5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_csql_generator_eval.md) | ⭐ 18 | `commercial` |
| [Cs Customer Journey Orchestrator Eval](commercial/cs_customer_journey_orchestrator_eval_398a7519/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_customer_journey_orchestrator_eval.md) | ⭐ 18 | `commercial` |
| [Cs Customer Maturity Model Eval](commercial/cs_customer_maturity_model_eval_6b9a931d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_customer_maturity_model_eval.md) | ⭐ 18 | `commercial` |
| [Cs Escalation Manager Eval](commercial/cs_escalation_manager_eval_70f5e9d8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_escalation_manager_eval.md) | ⭐ 18 | `commercial` |
| [Cs Executive Alignment Eval](commercial/cs_executive_alignment_eval_0fdbcb18/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_executive_alignment_eval.md) | ⭐ 18 | `commercial` |
| [Cs Multi Product Adoption Eval](commercial/cs_multi_product_adoption_eval_11fbbc59/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_multi_product_adoption_eval.md) | ⭐ 18 | `commercial` |
| [Cs Onboarding Health Eval](commercial/cs_onboarding_health_eval_c10da718/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_onboarding_health_eval.md) | ⭐ 18 | `commercial` |
| [Cs Outcome Tracker Eval](commercial/cs_outcome_tracker_eval_7d48f9f3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_outcome_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Cs Playbook Selector Eval](commercial/cs_playbook_selector_eval_5a63ca4e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_playbook_selector_eval.md) | ⭐ 18 | `commercial` |
| [Cs Product Adoption Guidance Eval](commercial/cs_product_adoption_guidance_eval_a5f3f3ba/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_product_adoption_guidance_eval.md) | ⭐ 18 | `commercial` |
| [Cs Product Gap Reporter Eval](commercial/cs_product_gap_reporter_eval_9099ab11/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_product_gap_reporter_eval.md) | ⭐ 18 | `commercial` |
| [Cs Quarterly Business Review Eval](commercial/cs_quarterly_business_review_eval_2518855a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_quarterly_business_review_eval.md) | ⭐ 18 | `commercial` |
| [Cs Red Flag Detector Eval](commercial/cs_red_flag_detector_eval_527d8161/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_red_flag_detector_eval.md) | ⭐ 18 | `commercial` |
| [Cs Sentiment Analyzer Eval](commercial/cs_sentiment_analyzer_eval_e5a5a2fe/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_sentiment_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Cs Stakeholder Mapper Eval](commercial/cs_stakeholder_mapper_eval_cca2e744/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_stakeholder_mapper_eval.md) | ⭐ 18 | `commercial` |
| [Cs Success Plan Generator Eval](commercial/cs_success_plan_generator_eval_a579f6e6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_success_plan_generator_eval.md) | ⭐ 18 | `commercial` |
| [Cs Time To Impact Eval](commercial/cs_time_to_impact_eval_7d07d7e7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_time_to_impact_eval.md) | ⭐ 18 | `commercial` |
| [Cs Vla Tracker Eval](commercial/cs_vla_tracker_eval_e7500715/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/cs_vla_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Data Anomaly Alerter Eval](commercial/data_anomaly_alerter_eval_3fc60fcd/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_anomaly_alerter_eval.md) | ⭐ 18 | `commercial` |
| [Data Attribution Modeler Eval](commercial/data_attribution_modeler_eval_86a85b5f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_attribution_modeler_eval.md) | ⭐ 18 | `commercial` |
| [Data Cohort Builder Eval](commercial/data_cohort_builder_eval_6cdc7418/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_cohort_builder_eval.md) | ⭐ 18 | `commercial` |
| [Data Dashboard Builder Eval](commercial/data_dashboard_builder_eval_affd5301/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_dashboard_builder_eval.md) | ⭐ 18 | `commercial` |
| [Data Etl Monitor Eval](commercial/data_etl_monitor_eval_c4384dbb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_etl_monitor_eval.md) | ⭐ 18 | `commercial` |
| [Data Event Validator Eval](commercial/data_event_validator_eval_4911ceca/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_event_validator_eval.md) | ⭐ 18 | `commercial` |
| [Data Experiment Analyzer Eval](commercial/data_experiment_analyzer_eval_5b46a191/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_experiment_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Data Funnel Optimizer Eval](commercial/data_funnel_optimizer_eval_9fa51af1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_funnel_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Data Privacy Scanner Eval](commercial/data_privacy_scanner_eval_daedeca3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_privacy_scanner_eval.md) | ⭐ 18 | `commercial` |
| [Data Warehouse Optimizer Eval](commercial/data_warehouse_optimizer_eval_26ee3f0c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/data_warehouse_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Devex Api Onboarding Eval](commercial/devex_api_onboarding_eval_7a3e1942/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_api_onboarding_eval.md) | ⭐ 18 | `commercial` |
| [Devex Changelog Tracker Eval](commercial/devex_changelog_tracker_eval_013f98b6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_changelog_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Devex Code Sample Generator Eval](commercial/devex_code_sample_generator_eval_e7440a10/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_code_sample_generator_eval.md) | ⭐ 18 | `commercial` |
| [Devex Deprecation Notifier Eval](commercial/devex_deprecation_notifier_eval_cec775a3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_deprecation_notifier_eval.md) | ⭐ 18 | `commercial` |
| [Devex Error Explainer Eval](commercial/devex_error_explainer_eval_85b51b3f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_error_explainer_eval.md) | ⭐ 18 | `commercial` |
| [Devex Integration Health Eval](commercial/devex_integration_health_eval_d6236daf/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_integration_health_eval.md) | ⭐ 18 | `commercial` |
| [Devex Latency Advisor Eval](commercial/devex_latency_advisor_eval_2759c019/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_latency_advisor_eval.md) | ⭐ 18 | `commercial` |
| [Devex Migration Assistant Eval](commercial/devex_migration_assistant_eval_0509c864/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_migration_assistant_eval.md) | ⭐ 18 | `commercial` |
| [Devex Oauth Helper Eval](commercial/devex_oauth_helper_eval_250de550/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_oauth_helper_eval.md) | ⭐ 18 | `commercial` |
| [Devex Sandbox Manager Eval](commercial/devex_sandbox_manager_eval_1b027e6c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_sandbox_manager_eval.md) | ⭐ 18 | `commercial` |
| [Devex Sdk Version Monitor Eval](commercial/devex_sdk_version_monitor_eval_8c721ddf/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_sdk_version_monitor_eval.md) | ⭐ 18 | `commercial` |
| [Devex Support Deflector Eval](commercial/devex_support_deflector_eval_966aaa8b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_support_deflector_eval.md) | ⭐ 18 | `commercial` |
| [Devex Webhook Tester Eval](commercial/devex_webhook_tester_eval_6e9c9181/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/devex_webhook_tester_eval.md) | ⭐ 18 | `commercial` |
| [Elg Channel Enablement Eval](commercial/elg_channel_enablement_eval_6abe5ad2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_channel_enablement_eval.md) | ⭐ 18 | `commercial` |
| [Elg Co Sell Trigger Eval](commercial/elg_co_sell_trigger_eval_97bc1b1d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_co_sell_trigger_eval.md) | ⭐ 18 | `commercial` |
| [Elg Deal Registration Eval](commercial/elg_deal_registration_eval_cfe77025/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_deal_registration_eval.md) | ⭐ 18 | `commercial` |
| [Elg Ecosystem Intelligence Eval](commercial/elg_ecosystem_intelligence_eval_37bd319c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_ecosystem_intelligence_eval.md) | ⭐ 18 | `commercial` |
| [Elg Eql Scoring Eval](commercial/elg_eql_scoring_eval_6d858d47/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_eql_scoring_eval.md) | ⭐ 18 | `commercial` |
| [Elg Integration Health Eval](commercial/elg_integration_health_eval_22aa1742/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_integration_health_eval.md) | ⭐ 18 | `commercial` |
| [Elg Joint Marketing Eval](commercial/elg_joint_marketing_eval_3cf8fd60/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_joint_marketing_eval.md) | ⭐ 18 | `commercial` |
| [Elg Marketplace Integration Eval](commercial/elg_marketplace_integration_eval_93777779/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_marketplace_integration_eval.md) | ⭐ 18 | `commercial` |
| [Elg Marketplace Listing Optimizer Eval](commercial/elg_marketplace_listing_optimizer_eval_99114b2f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_marketplace_listing_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Elg Mdf Tracker Eval](commercial/elg_mdf_tracker_eval_6849431b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_mdf_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Elg Nearbound Signal Eval](commercial/elg_nearbound_signal_eval_1975ee29/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_nearbound_signal_eval.md) | ⭐ 18 | `commercial` |
| [Elg Partner Influenced Revenue Eval](commercial/elg_partner_influenced_revenue_eval_978bd04c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_partner_influenced_revenue_eval.md) | ⭐ 18 | `commercial` |
| [Elg Partner Mapping Eval](commercial/elg_partner_mapping_eval_1e2c3c05/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_partner_mapping_eval.md) | ⭐ 18 | `commercial` |
| [Elg Partner Tier Manager Eval](commercial/elg_partner_tier_manager_eval_2ff5a0e7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_partner_tier_manager_eval.md) | ⭐ 18 | `commercial` |
| [Elg Referral Program Eval](commercial/elg_referral_program_eval_9b0a6864/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_referral_program_eval.md) | ⭐ 18 | `commercial` |
| [Elg Tech Partner Finder Eval](commercial/elg_tech_partner_finder_eval_cdcd1cb2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/elg_tech_partner_finder_eval.md) | ⭐ 18 | `commercial` |
| [Finops Arr Waterfall Eval](commercial/finops_arr_waterfall_eval_06982553/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_arr_waterfall_eval.md) | ⭐ 18 | `commercial` |
| [Finops Benchmark Comparator Eval](commercial/finops_benchmark_comparator_eval_a67f3d97/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_benchmark_comparator_eval.md) | ⭐ 18 | `commercial` |
| [Finops Burn Rate Monitor Eval](commercial/finops_burn_rate_monitor_eval_bfbdbb78/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_burn_rate_monitor_eval.md) | ⭐ 18 | `commercial` |
| [Finops Cac Calculator Eval](commercial/finops_cac_calculator_eval_e025f421/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_cac_calculator_eval.md) | ⭐ 18 | `commercial` |
| [Finops Expense Allocator Eval](commercial/finops_expense_allocator_eval_0ef4bc25/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_expense_allocator_eval.md) | ⭐ 18 | `commercial` |
| [Finops Gross Margin Eval](commercial/finops_gross_margin_eval_a6e1a201/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_gross_margin_eval.md) | ⭐ 18 | `commercial` |
| [Finops Investor Metrics Eval](commercial/finops_investor_metrics_eval_7f551918/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_investor_metrics_eval.md) | ⭐ 18 | `commercial` |
| [Finops Magic Number Eval](commercial/finops_magic_number_eval_acef26ff/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_magic_number_eval.md) | ⭐ 18 | `commercial` |
| [Finops Ndr Tracker Eval](commercial/finops_ndr_tracker_eval_bfada2b6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_ndr_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Finops Payback Period Eval](commercial/finops_payback_period_eval_379726c1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_payback_period_eval.md) | ⭐ 18 | `commercial` |
| [Finops Rule Of 40 Eval](commercial/finops_rule_of_40_eval_ada4a21a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_rule_of_40_eval.md) | ⭐ 18 | `commercial` |
| [Finops Scenario Planner Eval](commercial/finops_scenario_planner_eval_213ba6a3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/finops_scenario_planner_eval.md) | ⭐ 18 | `commercial` |
| [Governance Context Sync Eval](commercial/governance_context_sync_eval_489f026a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/governance_context_sync_eval.md) | ⭐ 18 | `commercial` |
| [Mon Cohort Ltv Analyzer Eval](commercial/mon_cohort_ltv_analyzer_eval_21fea3de/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_cohort_ltv_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Mon Consumption Analyzer Eval](commercial/mon_consumption_analyzer_eval_3079333b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_consumption_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Mon Contract Value Tracker Eval](commercial/mon_contract_value_tracker_eval_ac90523d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_contract_value_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Mon Credits Balance Tracker Eval](commercial/mon_credits_balance_tracker_eval_52c95e5f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_credits_balance_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Mon Entitlement Manager Eval](commercial/mon_entitlement_manager_eval_5ed2dc63/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_entitlement_manager_eval.md) | ⭐ 18 | `commercial` |
| [Mon Invoice Explainer Eval](commercial/mon_invoice_explainer_eval_caba6641/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_invoice_explainer_eval.md) | ⭐ 18 | `commercial` |
| [Mon Mrr Movement Tracker Eval](commercial/mon_mrr_movement_tracker_eval_97cd4fcc/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_mrr_movement_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Mon Overage Predictor Eval](commercial/mon_overage_predictor_eval_8b84d7a4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_overage_predictor_eval.md) | ⭐ 18 | `commercial` |
| [Mon Packaging Optimizer Eval](commercial/mon_packaging_optimizer_eval_3fd5569f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_packaging_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Mon Payment Method Optimizer Eval](commercial/mon_payment_method_optimizer_eval_07755e8c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_payment_method_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Mon Price Experimentation Eval](commercial/mon_price_experimentation_eval_171eb119/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_price_experimentation_eval.md) | ⭐ 18 | `commercial` |
| [Mon Revenue Recognition Eval](commercial/mon_revenue_recognition_eval_b76c0e76/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_revenue_recognition_eval.md) | ⭐ 18 | `commercial` |
| [Mon Tax Compliance Eval](commercial/mon_tax_compliance_eval_05165f11/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/mon_tax_compliance_eval.md) | ⭐ 18 | `commercial` |
| [People Candidate Sourcer Eval](commercial/people_candidate_sourcer_eval_ba098af0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_candidate_sourcer_eval.md) | ⭐ 18 | `commercial` |
| [People Comp Benchmarker Eval](commercial/people_comp_benchmarker_eval_2305faa0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_comp_benchmarker_eval.md) | ⭐ 18 | `commercial` |
| [People Dei Tracker Eval](commercial/people_dei_tracker_eval_7c3e4415/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_dei_tracker_eval.md) | ⭐ 18 | `commercial` |
| [People Offboarding Manager Eval](commercial/people_offboarding_manager_eval_4d56224d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_offboarding_manager_eval.md) | ⭐ 18 | `commercial` |
| [People Onboarding Checklist Eval](commercial/people_onboarding_checklist_eval_4fb20a8e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_onboarding_checklist_eval.md) | ⭐ 18 | `commercial` |
| [People Perf Review Generator Eval](commercial/people_perf_review_generator_eval_30bbe8fe/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_perf_review_generator_eval.md) | ⭐ 18 | `commercial` |
| [People Pulse Analyzer Eval](commercial/people_pulse_analyzer_eval_2b867612/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_pulse_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [People Skill Gap Analyzer Eval](commercial/people_skill_gap_analyzer_eval_c11ae693/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/people_skill_gap_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Plg Activation Eval](commercial/plg_activation_eval_8cb5d547/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_activation_eval.md) | ⭐ 18 | `commercial` |
| [Plg Aha Moment Detection Eval](commercial/plg_aha_moment_detection_eval_4c040154/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_aha_moment_detection_eval.md) | ⭐ 18 | `commercial` |
| [Plg Contextual Help Eval](commercial/plg_contextual_help_eval_f250f13d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_contextual_help_eval.md) | ⭐ 18 | `commercial` |
| [Plg Empty State Optimizer Eval](commercial/plg_empty_state_optimizer_eval_04245f84/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_empty_state_optimizer_eval.md) | ⭐ 18 | `commercial` |
| [Plg Feature Request Handler Eval](commercial/plg_feature_request_handler_eval_925c48ef/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_feature_request_handler_eval.md) | ⭐ 18 | `commercial` |
| [Plg Friction Detector Eval](commercial/plg_friction_detector_eval_3da86b05/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_friction_detector_eval.md) | ⭐ 18 | `commercial` |
| [Plg Guided Setup Wizard Eval](commercial/plg_guided_setup_wizard_eval_cd2dae50/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_guided_setup_wizard_eval.md) | ⭐ 18 | `commercial` |
| [Plg Habit Loop Builder Eval](commercial/plg_habit_loop_builder_eval_030385fe/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_habit_loop_builder_eval.md) | ⭐ 18 | `commercial` |
| [Plg Interactive Tour Eval](commercial/plg_interactive_tour_eval_578b92a7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_interactive_tour_eval.md) | ⭐ 18 | `commercial` |
| [Plg Milestone Celebration Eval](commercial/plg_milestone_celebration_eval_ff5ea93f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_milestone_celebration_eval.md) | ⭐ 18 | `commercial` |
| [Plg Network Effect Amplifier Eval](commercial/plg_network_effect_amplifier_eval_0a84f8d4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_network_effect_amplifier_eval.md) | ⭐ 18 | `commercial` |
| [Plg Onboarding Flow Eval](commercial/plg_onboarding_flow_eval_be15c32b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_onboarding_flow_eval.md) | ⭐ 18 | `commercial` |
| [Plg Personalized Checklist Eval](commercial/plg_personalized_checklist_eval_401d8ce0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_personalized_checklist_eval.md) | ⭐ 18 | `commercial` |
| [Plg Pql Scoring Eval](commercial/plg_pql_scoring_eval_203bf062/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_pql_scoring_eval.md) | ⭐ 18 | `commercial` |
| [Plg Quick Win Generator Eval](commercial/plg_quick_win_generator_eval_127c04da/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_quick_win_generator_eval.md) | ⭐ 18 | `commercial` |
| [Plg Reverse Trial Eval](commercial/plg_reverse_trial_eval_5b31a475/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_reverse_trial_eval.md) | ⭐ 18 | `commercial` |
| [Plg Sandbox Environment Eval](commercial/plg_sandbox_environment_eval_b0d8f6c0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_sandbox_environment_eval.md) | ⭐ 18 | `commercial` |
| [Plg Self Serve Expansion Eval](commercial/plg_self_serve_expansion_eval_8cf2446d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_self_serve_expansion_eval.md) | ⭐ 18 | `commercial` |
| [Plg Social Proof Injector Eval](commercial/plg_social_proof_injector_eval_048b6307/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_social_proof_injector_eval.md) | ⭐ 18 | `commercial` |
| [Plg Time To Value Eval](commercial/plg_time_to_value_eval_248fae90/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_time_to_value_eval.md) | ⭐ 18 | `commercial` |
| [Plg Trial Extension Evaluator Eval](commercial/plg_trial_extension_evaluator_eval_5955e746/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_trial_extension_evaluator_eval.md) | ⭐ 18 | `commercial` |
| [Plg Usage Depth Analyzer Eval](commercial/plg_usage_depth_analyzer_eval_3608681a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_usage_depth_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Plg Viral Loop Eval](commercial/plg_viral_loop_eval_a8ab2022/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plg_viral_loop_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Activation Metrics Eval](commercial/plgf_activation_metrics_eval_ee69f103/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_activation_metrics_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Engagement Loops Eval](commercial/plgf_engagement_loops_eval_1d0a7576/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_engagement_loops_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Expansion Revenue Eval](commercial/plgf_expansion_revenue_eval_32c20178/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_expansion_revenue_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Feature Adoption Eval](commercial/plgf_feature_adoption_eval_fb9c4403/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_feature_adoption_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Feature Gating Eval](commercial/plgf_feature_gating_eval_7978800c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_feature_gating_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Free Tool Strategy Eval](commercial/plgf_free_tool_strategy_eval_373c031f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_free_tool_strategy_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Growth Loops Eval](commercial/plgf_growth_loops_eval_29a6f302/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_growth_loops_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Growth Modeling Eval](commercial/plgf_growth_modeling_eval_ff40627e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_growth_modeling_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Ideas Eval](commercial/plgf_ideas_eval_da6fd6ef/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_ideas_eval.md) | ⭐ 18 | `commercial` |
| [Plgf In Product Messaging Eval](commercial/plgf_in_product_messaging_eval_70ee61c8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_in_product_messaging_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Mental Models Eval](commercial/plgf_mental_models_eval_cb7cf3c6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_mental_models_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Metrics Eval](commercial/plgf_metrics_eval_b2dbda10/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_metrics_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Paywall Upgrade Cro Eval](commercial/plgf_paywall_upgrade_cro_eval_deab0da4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_paywall_upgrade_cro_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Pricing Strategy Eval](commercial/plgf_pricing_strategy_eval_b0d164f8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_pricing_strategy_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Product Analytics Eval](commercial/plgf_product_analytics_eval_67df35b0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_product_analytics_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Product Led Sales Eval](commercial/plgf_product_led_sales_eval_422bd000/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_product_led_sales_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Product Onboarding Eval](commercial/plgf_product_onboarding_eval_767ed472/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_product_onboarding_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Referral Program Eval](commercial/plgf_referral_program_eval_fbd66f74/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_referral_program_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Retention Analysis Eval](commercial/plgf_retention_analysis_eval_01cef32c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_retention_analysis_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Self Serve Motion Eval](commercial/plgf_self_serve_motion_eval_3ea329b7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_self_serve_motion_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Signup Flow Cro Eval](commercial/plgf_signup_flow_cro_eval_fca47ffc/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_signup_flow_cro_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Strategy Eval](commercial/plgf_strategy_eval_d1b5c9e4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_strategy_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Trial Optimization Eval](commercial/plgf_trial_optimization_eval_be088442/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_trial_optimization_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Usage Based Pricing Eval](commercial/plgf_usage_based_pricing_eval_6612274f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_usage_based_pricing_eval.md) | ⭐ 18 | `commercial` |
| [Plgf User Segmentation Eval](commercial/plgf_user_segmentation_eval_3635260f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_user_segmentation_eval.md) | ⭐ 18 | `commercial` |
| [Plgf Viral Loops Eval](commercial/plgf_viral_loops_eval_156db78b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/plgf_viral_loops_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Accessibility Auditor Eval](commercial/prodops_accessibility_auditor_eval_537beda1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_accessibility_auditor_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Api Deprecation Eval](commercial/prodops_api_deprecation_eval_ff65f294/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_api_deprecation_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Beta Program Manager Eval](commercial/prodops_beta_program_manager_eval_445776c1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_beta_program_manager_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Competitive Feature Tracker Eval](commercial/prodops_competitive_feature_tracker_eval_c89282bf/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_competitive_feature_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Documentation Health Eval](commercial/prodops_documentation_health_eval_a8f43afb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_documentation_health_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Dogfooding Tracker Eval](commercial/prodops_dogfooding_tracker_eval_6e971c96/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_dogfooding_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Experimentation Eval](commercial/prodops_experimentation_eval_24694948/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_experimentation_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Impact Analyzer Eval](commercial/prodops_impact_analyzer_eval_94e78cab/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_impact_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Incident Analyzer Eval](commercial/prodops_incident_analyzer_eval_074de406/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_incident_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Localization Manager Eval](commercial/prodops_localization_manager_eval_81c4b91b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_localization_manager_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Performance Budget Eval](commercial/prodops_performance_budget_eval_b7369470/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_performance_budget_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Prioritization Engine Eval](commercial/prodops_prioritization_engine_eval_950d8396/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_prioritization_engine_eval.md) | ⭐ 18 | `commercial` |
| [Prodops Release Notes Generator Eval](commercial/prodops_release_notes_generator_eval_a47cf0e1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_release_notes_generator_eval.md) | ⭐ 18 | `commercial` |
| [Prodops User Interview Synthesizer Eval](commercial/prodops_user_interview_synthesizer_eval_d5dfc061/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/prodops_user_interview_synthesizer_eval.md) | ⭐ 18 | `commercial` |
| [Revops Activity Capture Eval](commercial/revops_activity_capture_eval_ca2a468a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_activity_capture_eval.md) | ⭐ 18 | `commercial` |
| [Revops Attribution Tracker Eval](commercial/revops_attribution_tracker_eval_cabe26d8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_attribution_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Revops Competitive Intel Eval](commercial/revops_competitive_intel_eval_221f662c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_competitive_intel_eval.md) | ⭐ 18 | `commercial` |
| [Revops Cpq Quote Generator Eval](commercial/revops_cpq_quote_generator_eval_fb5d3ddb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_cpq_quote_generator_eval.md) | ⭐ 18 | `commercial` |
| [Revops Data Enrichment Eval](commercial/revops_data_enrichment_eval_feaadfa8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_data_enrichment_eval.md) | ⭐ 18 | `commercial` |
| [Revops Deal Inspection Eval](commercial/revops_deal_inspection_eval_2b10d3ac/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_deal_inspection_eval.md) | ⭐ 18 | `commercial` |
| [Revops Deal Velocity Eval](commercial/revops_deal_velocity_eval_560b058a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_deal_velocity_eval.md) | ⭐ 18 | `commercial` |
| [Revops Forecast Intelligence Eval](commercial/revops_forecast_intelligence_eval_01b9b97b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_forecast_intelligence_eval.md) | ⭐ 18 | `commercial` |
| [Revops Handoff Orchestration Eval](commercial/revops_handoff_orchestration_eval_86dd8da2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_handoff_orchestration_eval.md) | ⭐ 18 | `commercial` |
| [Revops Lead Qualification Eval](commercial/revops_lead_qualification_eval_321874cc/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_lead_qualification_eval.md) | ⭐ 18 | `commercial` |
| [Revops Lead Routing Eval](commercial/revops_lead_routing_eval_cf1ae1eb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_lead_routing_eval.md) | ⭐ 18 | `commercial` |
| [Revops Multi Thread Tracker Eval](commercial/revops_multi_thread_tracker_eval_e2dc0333/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_multi_thread_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Revops Next Best Action Eval](commercial/revops_next_best_action_eval_629f400f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_next_best_action_eval.md) | ⭐ 18 | `commercial` |
| [Revops Opportunity Scoring Eval](commercial/revops_opportunity_scoring_eval_cb9974ce/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_opportunity_scoring_eval.md) | ⭐ 18 | `commercial` |
| [Revops Pipeline Health Eval](commercial/revops_pipeline_health_eval_e851e0b3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_pipeline_health_eval.md) | ⭐ 18 | `commercial` |
| [Revops Pricing Guidance Eval](commercial/revops_pricing_guidance_eval_3d236ec8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_pricing_guidance_eval.md) | ⭐ 18 | `commercial` |
| [Revops Quota Setter Eval](commercial/revops_quota_setter_eval_f4a10f2e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_quota_setter_eval.md) | ⭐ 18 | `commercial` |
| [Revops Renewals Handoff Eval](commercial/revops_renewals_handoff_eval_5b49e494/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_renewals_handoff_eval.md) | ⭐ 18 | `commercial` |
| [Revops Sales Coaching Eval](commercial/revops_sales_coaching_eval_0f2c978c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_sales_coaching_eval.md) | ⭐ 18 | `commercial` |
| [Revops Sales Content Recommender Eval](commercial/revops_sales_content_recommender_eval_a35bb63a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_sales_content_recommender_eval.md) | ⭐ 18 | `commercial` |
| [Revops Stage Duration Eval](commercial/revops_stage_duration_eval_80487aca/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_stage_duration_eval.md) | ⭐ 18 | `commercial` |
| [Revops Territory Planner Eval](commercial/revops_territory_planner_eval_20a188ab/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_territory_planner_eval.md) | ⭐ 18 | `commercial` |
| [Revops Win Loss Analyzer Eval](commercial/revops_win_loss_analyzer_eval_0c7d5705/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_win_loss_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Support Bug Linker Eval](commercial/support_bug_linker_eval_5ee899a9/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_bug_linker_eval.md) | ⭐ 18 | `commercial` |
| [Support Csat Analyzer Eval](commercial/support_csat_analyzer_eval_5f340526/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_csat_analyzer_eval.md) | ⭐ 18 | `commercial` |
| [Support Escalation Predictor Eval](commercial/support_escalation_predictor_eval_074f1a90/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_escalation_predictor_eval.md) | ⭐ 18 | `commercial` |
| [Support First Response Eval](commercial/support_first_response_eval_aed8a8b5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_first_response_eval.md) | ⭐ 18 | `commercial` |
| [Support Kb Gap Finder Eval](commercial/support_kb_gap_finder_eval_93f0b49d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_kb_gap_finder_eval.md) | ⭐ 18 | `commercial` |
| [Support Proactive Outreach Eval](commercial/support_proactive_outreach_eval_b7bd93b6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_proactive_outreach_eval.md) | ⭐ 18 | `commercial` |
| [Support Queue Balancer Eval](commercial/support_queue_balancer_eval_958b0224/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_queue_balancer_eval.md) | ⭐ 18 | `commercial` |
| [Support Resolution Suggester Eval](commercial/support_resolution_suggester_eval_666d2826/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_resolution_suggester_eval.md) | ⭐ 18 | `commercial` |
| [Support Sentiment Tracker Eval](commercial/support_sentiment_tracker_eval_0b1d1f55/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_sentiment_tracker_eval.md) | ⭐ 18 | `commercial` |
| [Support Sla Monitor Eval](commercial/support_sla_monitor_eval_6f0eda6e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_sla_monitor_eval.md) | ⭐ 18 | `commercial` |
| [Support Ticket Triage Eval](commercial/support_ticket_triage_eval_5c29b380/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/support_ticket_triage_eval.md) | ⭐ 18 | `commercial` |
| [Vcf Stakeholder Mapping Eval](commercial/vcf_stakeholder_mapping_eval_44191c2a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/vcf_stakeholder_mapping_eval.md) | ⭐ 18 | `commercial` |
| [Vcf Value Metrics Eval](commercial/vcf_value_metrics_eval_c293be34/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/vcf_value_metrics_eval.md) | ⭐ 18 | `commercial` |
| [Code Review Checklist Eval](commercial/code-review-checklist_eval_376bb2b5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/development/code-review-checklist_eval.md) | ⭐ 18 | `commercial` |
| [Codebase Onboarding Eval](commercial/codebase-onboarding_eval_8e35fe31/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/development/codebase-onboarding_eval.md) | ⭐ 18 | `commercial` |

### Communication (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Syntax Reference](communication/113-syntax-reference_a16eba1d/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/syntax-reference.md) | ⭐ 65 | `communication` |

### Content Creation (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Quick Wins](content-creation/quick_wins_29f2473a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/QUICK_WINS.md) | ⭐ 18 | `content creation` |
| [Skill Chains](content-creation/skill_chains_312b4d53/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/SKILL_CHAINS.md) | ⭐ 18 | `content creation` |
| [Skills Catalog](content-creation/skills_catalog_b303f24f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/SKILLS_CATALOG.md) | ⭐ 18 | `content creation` |
| [Ai Agent Composability Analysis](content-creation/ai_agent_composability_analysis_8ca26161/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AI_AGENT_COMPOSABILITY_ANALYSIS.md) | ⭐ 18 | `content creation` |
| [Batch 20260211 214136 Aggregate](content-creation/batch_20260211_214136_aggregate_47b6a9b2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_214136_aggregate.md) | ⭐ 18 | `content creation` |
| [Batch 20260211 214713 Aggregate](content-creation/batch_20260211_214713_aggregate_dfad9010/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_214713_aggregate.md) | ⭐ 18 | `content creation` |
| [Failure Analysis Batch 20260211 214136](content-creation/failure_analysis_batch_20260211_214136_2e9ae572/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_214136.md) | ⭐ 18 | `content creation` |

### Daily Assistant (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Revops Meeting Scheduler Eval](daily-assistant/revops_meeting_scheduler_eval_493debe5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/skills/revops_meeting_scheduler_eval.md) | ⭐ 18 | `daily assistant` |

### Data Analysis (11 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](data-analysis/226-name-skill_d0682987/) | [WILLOSCAR/research-units-pipeline-skills](https://raw.githubusercontent.com/WILLOSCAR/research-units-pipeline-skills/main/.codex/skills/idea-pool-expander/SKILL.md) | ⭐ 197 | `data analysis` |
| [Eval Harness Implementation](data-analysis/eval_harness_implementation_1cfc532d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/EVAL_HARNESS_IMPLEMENTATION.md) | ⭐ 18 | `data analysis` |
| [Eval Harness Phase 0 Complete](data-analysis/eval_harness_phase_0_complete_e2460ca4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/EVAL_HARNESS_PHASE_0_COMPLETE.md) | ⭐ 18 | `data analysis` |
| [Deduplication And Chaining](data-analysis/deduplication_and_chaining_bfadd626/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/DEDUPLICATION_AND_CHAINING.md) | ⭐ 18 | `data analysis` |
| [Example Workflow Diagram](data-analysis/example-workflow-diagram_cf31d473/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/example-workflow-diagram.md) | ⭐ 18 | `data analysis` |
| [Skill Tree](data-analysis/skill-tree_23bf21eb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/skill-tree.md) | ⭐ 18 | `data analysis` |
| [Engineering](data-analysis/engineering_6a9532c0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/engineering.md) | ⭐ 18 | `data analysis` |
| [Sales](data-analysis/sales_cdb3b292/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/sales.md) | ⭐ 18 | `data analysis` |
| [Failure Analysis Batch 20260211 221333](data-analysis/failure_analysis_batch_20260211_221333_bb42ff12/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_221333.md) | ⭐ 18 | `data analysis` |
| [Failure Analysis Batch 20260211 221435](data-analysis/failure_analysis_batch_20260211_221435_8ea4f9f3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_221435.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/name-skill_4a5837ce/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/tracing-upstream-lineage/SKILL.md) | ⭐ 208 | `data analysis` |

### Development (22 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Credits](development/239-credits_7d3b11da/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/CREDITS.md) | ⭐ 65 | `development` |
| [Actions Reference](development/010-actions-reference_733b0bdd/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/actions-reference.md) | ⭐ 65 | `development` |
| [Migration Guide](development/584-migration-guide_6f7f2a5a/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/migration-guide.md) | ⭐ 65 | `development` |
| [Ascii Design System Port](development/ascii_design_system_port_fb98ee4a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/ASCII_DESIGN_SYSTEM_PORT.md) | ⭐ 18 | `development` |
| [E2E Testing Strategy](development/e2e_testing_strategy_a8857663/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/E2E_TESTING_STRATEGY.md) | ⭐ 18 | `development` |
| [Skene Rebrand Summary](development/skene_rebrand_summary_ce6a5e15/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/SKENE_REBRAND_SUMMARY.md) | ⭐ 18 | `development` |
| [Social Announcements](development/social_announcements_2f9cc1d5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/SOCIAL_ANNOUNCEMENTS.md) | ⭐ 18 | `development` |
| [Before After](development/before_after_8e440787/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/BEFORE_AFTER.md) | ⭐ 18 | `development` |
| [Build Your First Skill](development/build_your_first_skill_d77b5552/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/BUILD_YOUR_FIRST_SKILL.md) | ⭐ 18 | `development` |
| [Cli Guide](development/cli_guide_46cf5f6d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/CLI_GUIDE.md) | ⭐ 18 | `development` |
| [Project Summary](development/project_summary_37f81d98/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/PROJECT_SUMMARY.md) | ⭐ 18 | `development` |
| [Value](development/value_95e26ad4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/VALUE.md) | ⭐ 18 | `development` |
| [Welcome Screen](development/welcome_screen_de0a7f4a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/WELCOME_SCREEN.md) | ⭐ 18 | `development` |
| [Data](development/data_eed398f1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/data.md) | ⭐ 18 | `development` |
| [Design](development/design_df0db37d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/design.md) | ⭐ 18 | `development` |
| [Operations](development/operations_e35ec716/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/operations.md) | ⭐ 18 | `development` |
| [Batch 20260211 221835 Aggregate](development/batch_20260211_221835_aggregate_ff3fbdf1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221835_aggregate.md) | ⭐ 18 | `development` |
| [Failure Analysis Batch 20260211 221835](development/failure_analysis_batch_20260211_221835_8382f786/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_221835.md) | ⭐ 18 | `development` |
| [Implementation Summary](development/implementation_summary_988e05dc/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/scripts/ecosystem-generator/IMPLEMENTATION_SUMMARY.md) | ⭐ 18 | `development` |
| [Skill](development/name-skill_47176d9d/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/airflow/SKILL.md) | ⭐ 208 | `development` |
| [Skill](development/name-skill_556763a0/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/debugging-dags/SKILL.md) | ⭐ 208 | `development` |
| [Skill](development/name-skill_7166aaf2/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/testing-dags/SKILL.md) | ⭐ 208 | `development` |

### Development/Devops (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Known Issues](development/devops/373-known-issues_97ef6daa/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/known-issues.md) | ⭐ 65 | `development` |
| [Batch Eval Quickstart](development/devops/batch_eval_quickstart_19b916c5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/BATCH_EVAL_QUICKSTART.md) | ⭐ 18 | `development` |
| [Troubleshooting](development/devops/troubleshooting_a88da6d4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/TROUBLESHOOTING.md) | ⭐ 18 | `development` |

### Development/Testing (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Batch 20260211 221333 Aggregate](development/testing/batch_20260211_221333_aggregate_2012d195/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221333_aggregate.md) | ⭐ 18 | `development` |
| [Batch 20260211 221435 Aggregate](development/testing/batch_20260211_221435_aggregate_09b4ca91/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221435_aggregate.md) | ⭐ 18 | `development` |
| [Hr Skills](development/testing/hr_60505201/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/hr.md) | ⭐ 18 | `development` |

### Development/Tools (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/tools/002-name-skill_727accdb/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/SKILL.md) | ⭐ 65 | `development` |
| [Official Sources](development/tools/320-official-sources_aab46d33/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/sf-ai-agentscript/resources/official-sources.md) | ⭐ 65 | `development` |
| [Security Analysis](development/tools/security_analysis_fdc4ce06/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/security_analysis.md) | ⭐ 18 | `development` |
| [Skill](development/tools/name-skill_9f88cf3e/) | [astronomer/agents](https://raw.githubusercontent.com/astronomer/agents/main/skills/authoring-dags/SKILL.md) | ⭐ 208 | `development` |

### Investment (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Ai Agent Summary](investment/ai_agent_summary_eb15178c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/technical/AI_AGENT_SUMMARY.md) | ⭐ 18 | `investment` |
| [Failure Analysis Batch 20260211 213131](investment/failure_analysis_batch_20260211_213131_17df2793/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_213131.md) | ⭐ 18 | `investment` |
| [Batch 20260211 171252 Aggregate](investment/batch_20260211_171252_aggregate_b701c31c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_171252_aggregate.md) | ⭐ 18 | `investment` |
| [Failure Analysis Batch 20260211 171252](investment/failure_analysis_batch_20260211_171252_cb909f7b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_171252.md) | ⭐ 18 | `investment` |
| [Failure Analysis Batch 20260211 171548](investment/failure_analysis_batch_20260211_171548_acfc1da5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_171548.md) | ⭐ 18 | `investment` |

### Productivity (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Growth Pm](productivity/growth-pm_48f76195/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/growth-pm.md) | ⭐ 18 | `productivity` |

### Research (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Showcase](research/showcase_cafd0d13/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/SHOWCASE.md) | ⭐ 18 | `research` |
| [Directory](research/directory_ffc70b99/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/directory.md) | ⭐ 18 | `research` |
| [Researcher](research/researcher_edc4a804/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/researcher.md) | ⭐ 18 | `research` |
| [Batch 20260211 221758 Aggregate](research/batch_20260211_221758_aggregate_8cac16a4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/batch_20260211_221758_aggregate.md) | ⭐ 18 | `research` |
| [Failure Analysis Batch 20260211 221758](research/failure_analysis_batch_20260211_221758_b09aaf4d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/reports/evals/failure_analysis_batch_20260211_221758.md) | ⭐ 18 | `research` |

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

*Last updated: 2026-02-12 20:23:51 UTC*
*Automatically maintained by SkillFlow*
