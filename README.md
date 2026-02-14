# X-Skills

A curated collection of **365 AI-powered skills** organized into 13 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (4 skills)
- **Automation/Workflow** (19 skills)
- **Commercial** (16 skills)
- **Communication** (5 skills)
- **Content Creation** (26 skills)
- **Daily Assistant** (7 skills)
- **Data Analysis** (58 skills)
- **Development** (132 skills)
- **Development/Devops** (29 skills)
- **Development/Testing** (5 skills)
- **Development/Tools** (14 skills)
- **Productivity** (7 skills)
- **Research** (43 skills)

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
| [Eval Harness Rollout Complete](automation/scripting/097-eval_harness_rollout_complete_f6f97f73/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/internal/EVAL_HARNESS_ROLLOUT_COMPLETE.md) | ⭐ 19 | `automation` |
| [Phase 1 Quick Reference](automation/scripting/098-phase_1_quick_reference_1f8f27ea/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/internal/PHASE_1_QUICK_REFERENCE.md) | ⭐ 19 | `automation` |
| [Instructions](automation/scripting/099-instructions_fd6daff5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/medchem/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/scripting/099-instructions_765c2947/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/rdkit/instructions.md) | ⭐ 19 | `automation` |

### Automation/Workflow (19 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/workflow/002-name-skill_768717c3/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Created](automation/workflow/140-readme_flat_skills_created_518f2683/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_CREATED.md) | 🔥 23.5k | `automation` |
| [Skill](automation/workflow/002-name-skill_82d9d209/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tool-ui-integrator/SKILL.md) | ⭐ 15 | `automation` |
| [Readme Flat Skills Az](automation/workflow/136-readme_flat_skills_az_1dd4094b/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_AZ.md) | 🔥 23.5k | `automation` |
| [Overview](automation/workflow/016-overview_cf2b1aeb/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/docs/ai-agents/overview.md) | ⭐ 17 | `automation` |
| [Security Scorecard](automation/workflow/035-security_scorecard_a19e5df8/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/.github/SECURITY_SCORECARD.md) | ⭐ 4.0k | `automation` |
| [Knowledge Architecture V5.2](automation/workflow/028-knowledge-architecture-v52_d2959d9f/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/KNOWLEDGE-ARCHITECTURE-v5.2.md) | ⭐ 151 | `automation` |
| [Control Set 07 Logging](automation/workflow/029-control-set-07-logging_4409d6b1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-07-logging.md) | ⭐ 151 | `automation` |
| [Asan Validation](automation/workflow/139-asan-validation_e97e1aeb/) | [moonbitlang/moonbit-agent-guide](https://raw.githubusercontent.com/moonbitlang/moonbit-agent-guide/main/moonbit-c-binding/references/asan-validation.md) | ⭐ 45 | `automation` |
| [Cfo](automation/workflow/139-cfo_2e10673f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/cfo.md) | ⭐ 19 | `automation` |
| [Sales Leader](automation/workflow/140-sales-leader_887038e0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/sales-leader.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_cc010b4b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/hypogenic/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_1e33a7ed/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/iso-13485-certification/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_1bff6d69/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/neuropixels-analysis/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_8400bec8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/offer-k-dense-web/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_7da61ecb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pytorch-lightning/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_f584517f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scvi-tools/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_ee39d309/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/transformers/instructions.md) | ⭐ 19 | `automation` |
| [Instructions](automation/workflow/134-instructions_51c596c0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/umap-learn/instructions.md) | ⭐ 19 | `automation` |

### Commercial (16 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/210-name-skill_f3e89c56/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/session-handoff/SKILL.md) | ⭐ 15 | `commercial` |
| [Interface Design](commercial/374-interface-design_764c5ff0/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/tdd/interface-design.md) | ⭐ 15 | `commercial` |
| [Dev Endpoint](commercial/375-dev-endpoint_ad8342d7/) | [petekp/claude-code-setup](https://raw.githubusercontent.com/petekp/claude-code-setup/main/skills/agent-telemetry/references/dev-endpoint.md) | ⭐ 15 | `commercial` |
| [Control Set 04 Output Encoding](commercial/046-control-set-04-output-encoding_d2c0a426/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-04-output-encoding.md) | ⭐ 151 | `commercial` |
| [Control Set 10 Data Protection](commercial/047-control-set-10-data-protection_cbab8bd5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-10-data-protection.md) | ⭐ 151 | `commercial` |
| [Control Set Ext 01 02 Auth Patterns](commercial/048-control-set-ext-01_02-auth-patterns_20aa6001/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-01_02-auth-patterns.md) | ⭐ 151 | `commercial` |
| [Reference Set 05 Third Party Javascript](commercial/168-reference-set-05-third-party-javascript_31a81d77/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-third-party-javascript.md) | ⭐ 151 | `commercial` |
| [Reference Set 09 Rest Security](commercial/169-reference-set-09-rest-security_57c8774a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-rest-security.md) | ⭐ 151 | `commercial` |
| [Requirement Violation](commercial/396-requirement-violation_8c2a25b1/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/requirement-violation.md) | ⭐ 90 | `commercial` |
| [Unchecked Return Values](commercial/397-unchecked-return-values_ae8ce7f8/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unchecked-return-values.md) | ⭐ 90 | `commercial` |
| [Customer Success](commercial/396-customer_success_03440bce/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/customer_success.md) | ⭐ 19 | `commercial` |
| [Executive](commercial/397-executive_48a48df7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/executive.md) | ⭐ 19 | `commercial` |
| [Finance](commercial/398-finance_81b589ce/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/finance.md) | ⭐ 19 | `commercial` |
| [Marketing](commercial/399-marketing_b3e32d9e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/marketing.md) | ⭐ 19 | `commercial` |
| [Instructions](commercial/107-instructions_babffa81/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/simpy/instructions.md) | ⭐ 19 | `commercial` |
| [Goal Setting Frameworks](commercial/095-goal_setting_frameworks_a28a8e14/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/treatment-plans/references/goal_setting_frameworks.md) | ⭐ 19 | `commercial` |

### Communication (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Troubleshoot](communication/221-troubleshoot_7a8329e7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/troubleshoot.md) | 🔥 35.9k | `communication` |
| [Web Search](communication/206-web_search_6a939ba7/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/completion/web_search.md) | 🔥 35.9k | `communication` |
| [Claude](communication/024-claude_8fbab90f/) | [pocketpaw/pocketpaw](https://raw.githubusercontent.com/pocketpaw/pocketpaw/main/docs/CLAUDE.md) | ⭐ 31 | `communication` |
| [Control Set 01 Authentication](communication/037-control-set-01-authentication_315e4b3c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-01-authentication.md) | ⭐ 151 | `communication` |
| [Reference Set 02 Idor Prevention](communication/096-reference-set-02-idor-prevention_344c8152/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-idor-prevention.md) | ⭐ 151 | `communication` |

### Content Creation (26 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Openai](content-creation/359-openai_4ef0cd71/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/providers/openai.md) | 🔥 35.9k | `content creation` |
| [Skill](content-creation/049-name-skill_a9cca04f/) | [srtab/daiv](https://raw.githubusercontent.com/srtab/daiv/main/daiv/automation/agent/skills/skill-creator/SKILL.md) | ⭐ 17 | `content creation` |
| [Target User Agentic Harness V0](content-creation/375-target_user_agentic_harness_v0_a84a827c/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/prds/target_user_agentic_harness_v0.md) | ⭐ 51 | `content creation` |
| [Claude](content-creation/007-claude_a0274342/) | [davila7/claude-code-templates](https://raw.githubusercontent.com/davila7/claude-code-templates/main/CLAUDE.md) | 🔥 20.1k | `content creation` |
| [Reference Set 01 Credential Stuffing Prevention](content-creation/092-reference-set-01-credential-stuffing-prevention_06dc56b3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-credential-stuffing-prevention.md) | ⭐ 151 | `content creation` |
| [Reference Set 01 Saml Security](content-creation/093-reference-set-01-saml-security_5bf46ff0/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-saml-security.md) | ⭐ 151 | `content creation` |
| [Reference Set 01 Security Questions](content-creation/094-reference-set-01-security-questions_b5b13f0b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-security-questions.md) | ⭐ 151 | `content creation` |
| [Reference Set 03 Bean Validation](content-creation/095-reference-set-03-bean-validation_cb7cb648/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-bean-validation.md) | ⭐ 151 | `content creation` |
| [Reference Set 05 Csp](content-creation/096-reference-set-05-csp_a185fb6c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-csp.md) | ⭐ 151 | `content creation` |
| [Reference Set 05 Xss Filter Evasion](content-creation/097-reference-set-05-xss-filter-evasion_ba27d20a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-xss-filter-evasion.md) | ⭐ 151 | `content creation` |
| [Reference Set 07 Logging](content-creation/098-reference-set-07-logging_b8ad9e2b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-07-logging.md) | ⭐ 151 | `content creation` |
| [Reference Set Ext 11 Docker Security](content-creation/099-reference-set-ext-11-docker-security_0ecaf542/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-11-docker-security.md) | ⭐ 151 | `content creation` |
| [Reference Set Ext 13 Ai Agent Security](content-creation/100-reference-set-ext-13-ai-agent-security_43ee4a42/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-13-ai-agent-security.md) | ⭐ 151 | `content creation` |
| [Cheatsheet](content-creation/376-cheatsheet_cd8486fc/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/CHEATSHEET.md) | ⭐ 90 | `content creation` |
| [Moonbit Language Fundamentals.Mbt](content-creation/376-moonbit-language-fundamentalsmbt_12a9c011/) | [moonbitlang/moonbit-agent-guide](https://raw.githubusercontent.com/moonbitlang/moonbit-agent-guide/main/moonbit-agent-guide/references/moonbit-language-fundamentals.mbt.md) | ⭐ 45 | `content creation` |
| [Including C Sources](content-creation/377-including-c-sources_f2b1b21b/) | [moonbitlang/moonbit-agent-guide](https://raw.githubusercontent.com/moonbitlang/moonbit-agent-guide/main/moonbit-c-binding/references/including-c-sources.md) | ⭐ 45 | `content creation` |
| [Skill Chains](content-creation/375-skill_chains_e505df82/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/SKILL_CHAINS.md) | ⭐ 19 | `content creation` |
| [Skills Catalog](content-creation/376-skills_catalog_89484c1f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/SKILLS_CATALOG.md) | ⭐ 19 | `content creation` |
| [Instructions](content-creation/377-instructions_225c135e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pysam/instructions.md) | ⭐ 19 | `content creation` |
| [Case Report Guidelines](content-creation/378-case_report_guidelines_1a1ea7ff/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/case_report_guidelines.md) | ⭐ 19 | `content creation` |
| [Medical Terminology](content-creation/379-medical_terminology_55c7cc7e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/medical_terminology.md) | ⭐ 19 | `content creation` |
| [Broader Impacts](content-creation/005-broader_impacts_7419c446/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/references/broader_impacts.md) | ⭐ 19 | `content creation` |
| [Specific Aims Guide](content-creation/037-specific_aims_guide_d0e2462d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/references/specific_aims_guide.md) | ⭐ 19 | `content creation` |
| [Presentation Structure](content-creation/380-presentation_structure_a9ff0c9a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-slides/references/presentation_structure.md) | ⭐ 19 | `content creation` |
| [Figures Tables](content-creation/381-figures_tables_aeb400fe/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-writing/references/figures_tables.md) | ⭐ 19 | `content creation` |
| [01 About Us](content-creation/350-01-about-us_ae7c9910/) | [volcengine/OpenViking](https://raw.githubusercontent.com/volcengine/OpenViking/main/docs/en/about/01-about-us.md) | ⭐ 1.1k | `content creation` |

### Daily Assistant (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Index](daily-assistant/052-index_21c18601/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/index.md) | ⭐ 10 | `daily assistant` |
| [Workflow](daily-assistant/085-workflow_097c5b9a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/WORKFLOW.md) | ⭐ 151 | `daily assistant` |
| [P3 Trust Boundary](daily-assistant/086-p3-trust-boundary_69377f71/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P3-TRUST-BOUNDARY.md) | ⭐ 151 | `daily assistant` |
| [P5 Stride Analysis](daily-assistant/087-p5-stride-analysis_0687bde5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P5-STRIDE-ANALYSIS.md) | ⭐ 151 | `daily assistant` |
| [Insufficient Gas Griefing](daily-assistant/271-insufficient-gas-griefing_473c0530/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/insufficient-gas-griefing.md) | ⭐ 90 | `daily assistant` |
| [Intervention Guidelines](daily-assistant/271-intervention_guidelines_6bcc583a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/treatment-plans/references/intervention_guidelines.md) | ⭐ 19 | `daily assistant` |
| [Specialty Specific Guidelines](daily-assistant/272-specialty_specific_guidelines_05b82478/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/treatment-plans/references/specialty_specific_guidelines.md) | ⭐ 19 | `daily assistant` |

### Data Analysis (58 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Architecture](data-analysis/009-architecture_f62cb9de/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/architecture.md) | ⭐ 12 | `data analysis` |
| [Configuration](data-analysis/046-configuration_cda8812a/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/configuration.md) | ⭐ 12 | `data analysis` |
| [Transformation](data-analysis/481-transformation_39e9ae41/) | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/docs/guide/input_data/transformation.md) | ⭐ 12 | `data analysis` |
| [Index](data-analysis/113-index_152f012b/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/concepts/index.md) | ⭐ 10 | `data analysis` |
| [P2 Dfd Analysis](data-analysis/251-p2-dfd-analysis_3d231074/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P2-DFD-ANALYSIS.md) | ⭐ 151 | `data analysis` |
| [P8 Report Generation](data-analysis/476-p8-report-generation_bd398302/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P8-REPORT-GENERATION.md) | ⭐ 151 | `data analysis` |
| [Control Set Ext 13 Ai Llm](data-analysis/272-control-set-ext-13-ai-llm_9ed85b8d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-13-ai-llm.md) | ⭐ 151 | `data analysis` |
| [Skill](data-analysis/226-name-skill_248c2e1d/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/react-pdf/skills/react-pdf/SKILL.md) | ⭐ 90 | `data analysis` |
| [Index](data-analysis/113-index_be768ee2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/index.md) | ⭐ 3.3k | `data analysis` |
| [Engineering](data-analysis/481-engineering_b9bc4ec2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/engineering.md) | ⭐ 19 | `data analysis` |
| [Sales](data-analysis/482-sales_4e7ffc3f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/sales.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_ef41724e/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/anndata/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_3148acfa/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/astropy/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_9f404d90/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/cirq/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_58d9df07/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-decision-support/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_54270fc7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_bd804781/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/dask/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_956bead4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/datamol/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_50323d0f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/deepchem/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_c30153e1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/deeptools/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_4314ab70/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/exploratory-data-analysis/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_a5026171/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/flowio/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_72539ad4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/geopandas/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_0231ff33/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/histolab/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_c1da0945/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/imaging-data-commons/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_25226b77/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/lamindb/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_6bd46639/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/market-research-reports/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_507d059f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/markitdown/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_1cba2ce5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/matchms/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_7c72a2c8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/networkx/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_5b964c80/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/neurokit2/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_a816e40d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/plotly/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_170017be/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pydeseq2/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_01e0bf4f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pyhealth/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_6d789ca9/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pymatgen/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_16b2dcb2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pymc/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_844fc5b7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pymoo/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_fa5c6a8a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pyopenms/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_35a0f8f3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/qiskit/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_a3c69014/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-visualization/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_8820b67b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scikit-bio/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_c9a9e4ce/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scikit-learn/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_c35a430b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/seaborn/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_95b84ee0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/shap/instructions.md) | ⭐ 19 | `data analysis` |
| [Instructions](data-analysis/483-instructions_249499ae/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/statsmodels/instructions.md) | ⭐ 19 | `data analysis` |
| [Outcome Analysis](data-analysis/484-outcome_analysis_9d9c20e4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-decision-support/references/outcome_analysis.md) | ⭐ 19 | `data analysis` |
| [Patient Cohort Analysis](data-analysis/485-patient_cohort_analysis_1b4962c1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-decision-support/references/patient_cohort_analysis.md) | ⭐ 19 | `data analysis` |
| [Clinical Trial Reporting](data-analysis/486-clinical_trial_reporting_e6bb4900/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/clinical_trial_reporting.md) | ⭐ 19 | `data analysis` |
| [Data Presentation](data-analysis/064-data_presentation_84af36d4/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/data_presentation.md) | ⭐ 19 | `data analysis` |
| [Diagnostic Reports Standards](data-analysis/487-diagnostic_reports_standards_e3dded66/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/diagnostic_reports_standards.md) | ⭐ 19 | `data analysis` |
| [Poster Content Guide](data-analysis/488-poster_content_guide_6592f544/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/latex-posters/references/poster_content_guide.md) | ⭐ 19 | `data analysis` |
| [Visual Generation Guide](data-analysis/218-visual_generation_guide_a0973a53/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/market-research-reports/references/visual_generation_guide.md) | ⭐ 19 | `data analysis` |
| [Api Reference](data-analysis/008-api_reference_18d13414/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/markitdown/references/api_reference.md) | ⭐ 19 | `data analysis` |
| [File Formats](data-analysis/092-file_formats_051c6450/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/markitdown/references/file_formats.md) | ⭐ 19 | `data analysis` |
| [Best Practices](data-analysis/017-best_practices_7e0d2076/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-schematics/references/best_practices.md) | ⭐ 19 | `data analysis` |
| [Slide Design Principles](data-analysis/489-slide_design_principles_d39c01ab/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-slides/references/slide_design_principles.md) | ⭐ 19 | `data analysis` |
| [Visual Review Workflow](data-analysis/490-visual_review_workflow_579328d5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-slides/references/visual_review_workflow.md) | ⭐ 19 | `data analysis` |
| [Treatment Plan Standards](data-analysis/491-treatment_plan_standards_b6b58a5f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/treatment-plans/references/treatment_plan_standards.md) | ⭐ 19 | `data analysis` |

### Development (132 skills)

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
| [Skill Architecture Design Cn](development/981-skill-architecture-design-cn_1e77c8e5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/SKILL-ARCHITECTURE-DESIGN-cn.md) | ⭐ 151 | `development` |
| [Skill Architecture Design](development/980-skill-architecture-design_caff0b9a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/SKILL-ARCHITECTURE-DESIGN.md) | ⭐ 151 | `development` |
| [Skillset Threat Modeling Tour Cn V5](development/995-skillset-threat-modeling-tour-cn-v5_cadc4f29/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/Skillset-threat-modeling-tour-cn-v5.md) | ⭐ 151 | `development` |
| [Skillset Threat Modeling Tour V5](development/996-skillset-threat-modeling-tour-v5_3a0190a3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/Skillset-threat-modeling-tour-v5.md) | ⭐ 151 | `development` |
| [P1 Project Understanding](development/652-p1-project-understanding_3463d1f5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P1-PROJECT-UNDERSTANDING.md) | ⭐ 151 | `development` |
| [P6 Risk Validation](development/658-p6-risk-validation_68e1c938/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P6-RISK-VALIDATION.md) | ⭐ 151 | `development` |
| [P7 Mitigation Planning](development/659-p7-mitigation-planning_f47272a6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P7-MITIGATION-PLANNING.md) | ⭐ 151 | `development` |
| [P8R Detailed Report](development/2845-p8r-detailed-report_079f4af8/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P8R-DETAILED-REPORT.md) | ⭐ 151 | `development` |
| [Control Set 02 Authorization](development/202-control-set-02-authorization_d1c8fa3a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-02-authorization.md) | ⭐ 151 | `development` |
| [Control Set 03 Input Validation](development/203-control-set-03-input-validation_8f6d69d1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-03-input-validation.md) | ⭐ 151 | `development` |
| [Control Set 05 Client Side](development/204-control-set-05-client-side_44068fc6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-05-client-side.md) | ⭐ 151 | `development` |
| [Control Set 06 Cryptography](development/205-control-set-06-cryptography_ebb7424d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-06-cryptography.md) | ⭐ 151 | `development` |
| [Control Set 08 Error Handling](development/206-control-set-08-error-handling_8daca94b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-08-error-handling.md) | ⭐ 151 | `development` |
| [Control Set Ext 10 Hardcoded Credentials](development/208-control-set-ext-10-hardcoded-credentials_5b8d87dc/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-10-hardcoded-credentials.md) | ⭐ 151 | `development` |
| [Control Set Ext 14 Mobile](development/212-control-set-ext-14-mobile_56ce0614/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-14-mobile.md) | ⭐ 151 | `development` |
| [Control Set Ext 16 Agentic](development/214-control-set-ext-16-agentic_6f05df72/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-16-agentic.md) | ⭐ 151 | `development` |
| [Reference Set 01 Authentication](development/829-reference-set-01-authentication_75581456/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-authentication.md) | ⭐ 151 | `development` |
| [Reference Set 01 Cookie Theft Mitigation](development/830-reference-set-01-cookie-theft-mitigation_017a4462/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-cookie-theft-mitigation.md) | ⭐ 151 | `development` |
| [Reference Set 01 Forgot Password](development/832-reference-set-01-forgot-password_ac32c6b6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-forgot-password.md) | ⭐ 151 | `development` |
| [Reference Set 01 Jaas](development/833-reference-set-01-jaas_2a444833/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-jaas.md) | ⭐ 151 | `development` |
| [Reference Set 01 Jwt Java](development/834-reference-set-01-jwt-java_e46983c1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-jwt-java.md) | ⭐ 151 | `development` |
| [Reference Set 01 Multifactor Authentication](development/835-reference-set-01-multifactor-authentication_196f4fda/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-multifactor-authentication.md) | ⭐ 151 | `development` |
| [Reference Set 01 Password Storage](development/836-reference-set-01-password-storage_bba43c08/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-password-storage.md) | ⭐ 151 | `development` |
| [Reference Set 01 Session Management](development/838-reference-set-01-session-management_7aaafaee/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-01-session-management.md) | ⭐ 151 | `development` |
| [Reference Set 02 Access Control](development/839-reference-set-02-access-control_8dd17c5f/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-access-control.md) | ⭐ 151 | `development` |
| [Reference Set 02 Authorization](development/840-reference-set-02-authorization_dde4a0c3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-authorization.md) | ⭐ 151 | `development` |
| [Reference Set 02 Transaction Authorization](development/841-reference-set-02-transaction-authorization_350752a3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-02-transaction-authorization.md) | ⭐ 151 | `development` |
| [Reference Set 03 Deserialization](development/843-reference-set-03-deserialization_a0c3c75c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-deserialization.md) | ⭐ 151 | `development` |
| [Reference Set 03 File Upload](development/844-reference-set-03-file-upload_dde7e823/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-file-upload.md) | ⭐ 151 | `development` |
| [Reference Set 03 Injection Prevention Java](development/846-reference-set-03-injection-prevention-java_0247ca9c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-injection-prevention-java.md) | ⭐ 151 | `development` |
| [Reference Set 03 Injection Prevention](development/845-reference-set-03-injection-prevention_9a702874/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-injection-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Input Validation](development/847-reference-set-03-input-validation_b602fff0/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-input-validation.md) | ⭐ 151 | `development` |
| [Reference Set 03 Ldap Injection Prevention](development/848-reference-set-03-ldap-injection-prevention_1ae008ac/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-ldap-injection-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Mass Assignment](development/849-reference-set-03-mass-assignment_350e41f4/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-mass-assignment.md) | ⭐ 151 | `development` |
| [Reference Set 03 Os Command Injection Defense](development/850-reference-set-03-os-command-injection-defense_1542d73c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-os-command-injection-defense.md) | ⭐ 151 | `development` |
| [Reference Set 03 Query Parameterization](development/851-reference-set-03-query-parameterization_b2eba6f1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-query-parameterization.md) | ⭐ 151 | `development` |
| [Reference Set 03 Sql Injection Prevention](development/852-reference-set-03-sql-injection-prevention_80e8fc20/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-sql-injection-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Ssrf Prevention](development/853-reference-set-03-ssrf-prevention_2816fc60/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-ssrf-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 03 Unvalidated Redirects](development/854-reference-set-03-unvalidated-redirects_aab33a7d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-03-unvalidated-redirects.md) | ⭐ 151 | `development` |
| [Reference Set 05 Clickjacking Defense](development/855-reference-set-05-clickjacking-defense_7353ee0b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-clickjacking-defense.md) | ⭐ 151 | `development` |
| [Reference Set 05 Csrf Prevention](development/857-reference-set-05-csrf-prevention_d0e31b35/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-csrf-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Css Security](development/858-reference-set-05-css-security_ae726937/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-css-security.md) | ⭐ 151 | `development` |
| [Reference Set 05 Dom Clobbering Prevention](development/859-reference-set-05-dom-clobbering-prevention_5f3a6922/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-dom-clobbering-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Dom Xss Prevention](development/860-reference-set-05-dom-xss-prevention_0d83c37e/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-dom-xss-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Html5 Security](development/862-reference-set-05-html5-security_2e726920/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-html5-security.md) | ⭐ 151 | `development` |
| [Reference Set 05 Http Headers](development/863-reference-set-05-http-headers_031e23d6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-http-headers.md) | ⭐ 151 | `development` |
| [Reference Set 05 Prototype Pollution Prevention](development/864-reference-set-05-prototype-pollution-prevention_91655ed3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-prototype-pollution-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 05 Xss Prevention](development/866-reference-set-05-xss-prevention_fb45017c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-xss-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 06 Certificate Pinning](development/867-reference-set-06-certificate-pinning_f0af49c6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-certificate-pinning.md) | ⭐ 151 | `development` |
| [Reference Set 06 Cryptographic Storage](development/868-reference-set-06-cryptographic-storage_e3562829/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-cryptographic-storage.md) | ⭐ 151 | `development` |
| [Reference Set 06 Key Management](development/869-reference-set-06-key-management_33d0fc03/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-key-management.md) | ⭐ 151 | `development` |
| [Reference Set 06 Tls Cipher String](development/871-reference-set-06-tls-cipher-string_18d2f4a7/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-tls-cipher-string.md) | ⭐ 151 | `development` |
| [Reference Set 06 Tls](development/870-reference-set-06-tls_0321e669/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-tls.md) | ⭐ 151 | `development` |
| [Reference Set 06 Transport Layer Protection](development/872-reference-set-06-transport-layer-protection_60ac3b0b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-06-transport-layer-protection.md) | ⭐ 151 | `development` |
| [Reference Set 07 Logging Vocabulary](development/874-reference-set-07-logging-vocabulary_19ba9628/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-07-logging-vocabulary.md) | ⭐ 151 | `development` |
| [Reference Set 08 Error Handling](development/875-reference-set-08-error-handling_b2c98807/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-08-error-handling.md) | ⭐ 151 | `development` |
| [Reference Set 08 Xs Leaks](development/876-reference-set-08-xs-leaks_7525a3da/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-08-xs-leaks.md) | ⭐ 151 | `development` |
| [Reference Set 09 Ajax Security](development/877-reference-set-09-ajax-security_101c69a6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-ajax-security.md) | ⭐ 151 | `development` |
| [Reference Set 09 Dos Prevention](development/878-reference-set-09-dos-prevention_4021d4d8/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-dos-prevention.md) | ⭐ 151 | `development` |
| [Reference Set 09 Graphql](development/879-reference-set-09-graphql_a917b3e3/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-graphql.md) | ⭐ 151 | `development` |
| [Reference Set 09 Microservices Security](development/880-reference-set-09-microservices-security_a2daef28/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-microservices-security.md) | ⭐ 151 | `development` |
| [Reference Set 09 Rest Assessment](development/881-reference-set-09-rest-assessment_28d67ae6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-rest-assessment.md) | ⭐ 151 | `development` |
| [Reference Set 09 Web Service Security](development/882-reference-set-09-web-service-security_635ed73a/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-09-web-service-security.md) | ⭐ 151 | `development` |
| [Reference Set 10 Database Security](development/883-reference-set-10-database-security_e6160141/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-10-database-security.md) | ⭐ 151 | `development` |
| [Reference Set 10 Secrets Management](development/884-reference-set-10-secrets-management_f17ad0b1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-10-secrets-management.md) | ⭐ 151 | `development` |
| [Reference Set Ext 11 Nodejs Docker](development/888-reference-set-ext-11-nodejs-docker_28321053/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-11-nodejs-docker.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Dependency Management](development/890-reference-set-ext-12-dependency-management_b7983d5d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-dependency-management.md) | ⭐ 151 | `development` |
| [Reference Set Ext 14 Automotive Security](development/894-reference-set-ext-14-automotive-security_6b2d4499/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-14-automotive-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 14 Mobile Security](development/895-reference-set-ext-14-mobile-security_714b27e2/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-14-mobile-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 15 Cloud Architecture](development/896-reference-set-ext-15-cloud-architecture_965b5461/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-15-cloud-architecture.md) | ⭐ 151 | `development` |
| [Reference Set Ext 16 Agentic Security](development/897-reference-set-ext-16-agentic-security_55a161fa/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-16-agentic-security.md) | ⭐ 151 | `development` |
| [Skill](development/1178-name-skill_8e5bcc97/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/SKILL.md) | ⭐ 90 | `development` |
| [Components](development/2938-components_23eca15d/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/react-pdf/skills/react-pdf/references/components.md) | ⭐ 90 | `development` |
| [Asserting Contract From Code Size](development/2939-asserting-contract-from-code-size_53dd37d1/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/asserting-contract-from-code-size.md) | ⭐ 90 | `development` |
| [Dos Gas Limit](development/2940-dos-gas-limit_44ba2d7c/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/dos-gas-limit.md) | ⭐ 90 | `development` |
| [Hash Collision](development/2941-hash-collision_fd105532/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/hash-collision.md) | ⭐ 90 | `development` |
| [Inadherence To Standards](development/2942-inadherence-to-standards_0b295986/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/inadherence-to-standards.md) | ⭐ 90 | `development` |
| [Incorrect Constructor](development/2943-incorrect-constructor_63cf3819/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/incorrect-constructor.md) | ⭐ 90 | `development` |
| [Missing Protection Signature Replay](development/2944-missing-protection-signature-replay_88f560d0/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/missing-protection-signature-replay.md) | ⭐ 90 | `development` |
| [Reentrancy](development/2945-reentrancy_44db9f25/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/reentrancy.md) | ⭐ 90 | `development` |
| [Signature Malleability](development/2946-signature-malleability_2fd4c4c6/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/signature-malleability.md) | ⭐ 90 | `development` |
| [Timestamp Dependence](development/2947-timestamp-dependence_0f401084/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/timestamp-dependence.md) | ⭐ 90 | `development` |
| [Transaction Ordering Dependence](development/2948-transaction-ordering-dependence_33621d21/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/transaction-ordering-dependence.md) | ⭐ 90 | `development` |
| [Unbounded Return Data](development/2949-unbounded-return-data_7da19017/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unbounded-return-data.md) | ⭐ 90 | `development` |
| [Unencrypted Private Data On Chain](development/2950-unencrypted-private-data-on-chain_a43b92aa/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unencrypted-private-data-on-chain.md) | ⭐ 90 | `development` |
| [Unsafe Low Level Call](development/2951-unsafe-low-level-call_176d2a33/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unsafe-low-level-call.md) | ⭐ 90 | `development` |
| [Unsecure Signatures](development/2952-unsecure-signatures_44d435de/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unsecure-signatures.md) | ⭐ 90 | `development` |
| [Unsupported Opcodes](development/2953-unsupported-opcodes_f07f9c1f/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unsupported-opcodes.md) | ⭐ 90 | `development` |
| [Unused Variables](development/2954-unused-variables_15bc83a6/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/unused-variables.md) | ⭐ 90 | `development` |
| [Weak Sources Randomness](development/2955-weak-sources-randomness_34b619d7/) | [trailofbits/skills-curated](https://raw.githubusercontent.com/trailofbits/skills-curated/main/plugins/scv-scan/skills/scv-scan/references/weak-sources-randomness.md) | ⭐ 90 | `development` |
| [Skill](development/1178-name-skill_27d4a744/) | [moonbitlang/moonbit-agent-guide](https://raw.githubusercontent.com/moonbitlang/moonbit-agent-guide/main/moonbit-agent-guide/SKILL.md) | ⭐ 45 | `development` |
| [Skill](development/1178-name-skill_5cb3d28f/) | [moonbitlang/moonbit-agent-guide](https://raw.githubusercontent.com/moonbitlang/moonbit-agent-guide/main/moonbit-c-binding/SKILL.md) | ⭐ 45 | `development` |
| [Skill](development/1178-name-skill_49e95dff/) | [moonbitlang/moonbit-agent-guide](https://raw.githubusercontent.com/moonbitlang/moonbit-agent-guide/main/moonbit-refactoring/SKILL.md) | ⭐ 45 | `development` |
| [Ownership And Memory](development/2924-ownership-and-memory_ae9edde7/) | [moonbitlang/moonbit-agent-guide](https://raw.githubusercontent.com/moonbitlang/moonbit-agent-guide/main/moonbit-c-binding/references/ownership-and-memory.md) | ⭐ 45 | `development` |
| [Ascii Design System Port](development/2935-ascii_design_system_port_a9757625/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/ASCII_DESIGN_SYSTEM_PORT.md) | ⭐ 19 | `development` |
| [Metrics](development/575-metrics_85523ef8/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/METRICS.md) | ⭐ 19 | `development` |
| [Build Your First Skill](development/2936-build_your_first_skill_066b356a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/BUILD_YOUR_FIRST_SKILL.md) | ⭐ 19 | `development` |
| [Value](development/2937-value_647a07db/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/VALUE.md) | ⭐ 19 | `development` |
| [Data](development/2938-data_50773e68/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/data.md) | ⭐ 19 | `development` |
| [Design](development/2939-design_9967d15a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/design.md) | ⭐ 19 | `development` |
| [Operations](development/2310-operations_90fff9a0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/functions/operations.md) | ⭐ 19 | `development` |
| [Instructions](development/1753-instructions_5d420954/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/esm/instructions.md) | ⭐ 19 | `development` |
| [Instructions](development/1753-instructions_a0279196/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pennylane/instructions.md) | ⭐ 19 | `development` |
| [Instructions](development/1753-instructions_b669431a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pylabrobot/instructions.md) | ⭐ 19 | `development` |
| [Instructions](development/1753-instructions_ea2c005a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/zarr-python/instructions.md) | ⭐ 19 | `development` |
| [Services Reference](development/961-services_reference_0109b06a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/bioservices/references/services_reference.md) | ⭐ 19 | `development` |
| [Patient Documentation](development/670-patient_documentation_1e4e5f8c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/patient_documentation.md) | ⭐ 19 | `development` |
| [Bigquery Guide](development/104-bigquery_guide_0ce04067/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/imaging-data-commons/references/bigquery_guide.md) | ⭐ 19 | `development` |
| [Regulatory Compliance](development/901-regulatory_compliance_805f5b21/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/treatment-plans/references/regulatory_compliance.md) | ⭐ 19 | `development` |

### Development/Devops (29 skills)

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
| [Readme Cn](development/devops/069-readme-cn_162e7ea6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/README-cn.md) | ⭐ 151 | `development` |
| [Sync Anthropic Beta Headers](development/devops/376-sync_anthropic_beta_headers_cc87512d/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/sync_anthropic_beta_headers.md) | 🔥 35.8k | `development` |
| [Architecture Workflow Guide Cn](development/devops/070-architecture-workflow-guide-cn_1a497d91/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/ARCHITECTURE-WORKFLOW-GUIDE-cn.md) | ⭐ 151 | `development` |
| [Architecture Workflow Guide](development/devops/071-architecture-workflow-guide_51ef58a8/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/ARCHITECTURE-WORKFLOW-GUIDE.md) | ⭐ 151 | `development` |
| [P4 Security Design Review](development/devops/057-p4-security-design-review_3ff4329b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/phases/P4-SECURITY-DESIGN-REVIEW.md) | ⭐ 151 | `development` |
| [Control Set Ext 11 Infrastructure](development/devops/072-control-set-ext-11-infrastructure_c39fd27c/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-11-infrastructure.md) | ⭐ 151 | `development` |
| [Control Set Ext 12 Supply Chain](development/devops/073-control-set-ext-12-supply-chain_e1e7576d/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-12-supply-chain.md) | ⭐ 151 | `development` |
| [Control Set Ext 15 Cloud](development/devops/074-control-set-ext-15-cloud_f1136c8f/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-ext-15-cloud.md) | ⭐ 151 | `development` |
| [Reference Set Ext 11 Kubernetes Security](development/devops/075-reference-set-ext-11-kubernetes-security_f9eb20b5/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-11-kubernetes-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Cicd Security](development/devops/076-reference-set-ext-12-cicd-security_0521ba12/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-cicd-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Npm Security](development/devops/077-reference-set-ext-12-npm-security_3c3eb33b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-npm-security.md) | ⭐ 151 | `development` |
| [Reference Set Ext 12 Supply Chain Security](development/devops/078-reference-set-ext-12-supply-chain-security_90024891/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-ext-12-supply-chain-security.md) | ⭐ 151 | `development` |
| [Configuration](development/devops/009-configuration_9d0f3514/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/configuration.md) | ⭐ 3.3k | `development` |
| [Batch Eval Quickstart](development/devops/371-batch_eval_quickstart_8a1c7ae5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/BATCH_EVAL_QUICKSTART.md) | ⭐ 19 | `development` |
| [Instructions](development/devops/197-instructions_39d88d4f/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/modal/instructions.md) | ⭐ 19 | `development` |
| [Docker Compose Guide](development/devops/090-docker-compose-guide_72061a3c/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/docker-compose-guide.md) | ⭐ 4.0k | `development` |
| [Agents](development/devops/agents_2f52ed9c/) | [crestalnetwork/intentkit](https://raw.githubusercontent.com/crestalnetwork/intentkit/main/integrations/telegram/AGENTS.md) | 🔥 6.5k | `development` |

### Development/Testing (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Middleware Chains](development/testing/085-middleware-chains_c239a6f7/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/middleware-chains.md) | ⭐ 10 | `development` |
| [Parallel Execution](development/testing/086-parallel-execution_5d3062a6/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/advanced/parallel-execution.md) | ⭐ 10 | `development` |
| [Skill](development/testing/013-threat-skill_0623a026/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/SKILL.md) | ⭐ 151 | `development` |
| [Release Notes](development/testing/release_notes_29151437/) | [crestalnetwork/intentkit](https://raw.githubusercontent.com/crestalnetwork/intentkit/main/release_notes.md) | 🔥 6.5k | `development` |
| [Ops Guide](development/testing/ops_guide_5991f949/) | [crestalnetwork/intentkit](https://raw.githubusercontent.com/crestalnetwork/intentkit/main/agent_docs/ops_guide.md) | 🔥 6.5k | `development` |

### Development/Tools (14 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Mcp Troubleshoot](development/tools/324-mcp_troubleshoot_6a5a0c41/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/mcp_troubleshoot.md) | 🔥 35.9k | `development` |
| [Desktop](development/tools/328-desktop_f1b0efd6/) | [kevinshowkat/brood](https://raw.githubusercontent.com/kevinshowkat/brood/main/docs/desktop.md) | ⭐ 51 | `development` |
| [Chain](development/tools/325-chain_02e6ae60/) | [vstorm-co/pydantic-ai-middleware](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-middleware/main/docs/api/chain.md) | ⭐ 10 | `development` |
| [Knowledge Architecture V5.2 Cn](development/tools/052-knowledge-architecture-v52-cn_0d76937b/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/KNOWLEDGE-ARCHITECTURE-v5.2-cn.md) | ⭐ 151 | `development` |
| [Report Design](development/tools/318-report-design_874a71e1/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/docs/REPORT-DESIGN.md) | ⭐ 151 | `development` |
| [Control Set 09 Api Security](development/tools/053-control-set-09-api-security_c77f3363/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/control-set-09-api-security.md) | ⭐ 151 | `development` |
| [Reference Set 05 Hsts](development/tools/054-reference-set-05-hsts_f54544af/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-05-hsts.md) | ⭐ 151 | `development` |
| [Reference Set 10 User Privacy](development/tools/055-reference-set-10-user-privacy_724c9ff6/) | [fr33d3m0n/threat-modeling](https://raw.githubusercontent.com/fr33d3m0n/threat-modeling/main/knowledge/security-controls/references/reference-set-10-user-privacy.md) | ⭐ 151 | `development` |
| [Instructions](development/tools/326-instructions_67e645e7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/aeon/instructions.md) | ⭐ 19 | `development` |
| [Instructions](development/tools/326-instructions_accc172d/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/geniml/instructions.md) | ⭐ 19 | `development` |
| [Instructions](development/tools/326-instructions_3df34127/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pufferlib/instructions.md) | ⭐ 19 | `development` |
| [Biomarker Classification](development/tools/327-biomarker_classification_32a4dbca/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-decision-support/references/biomarker_classification.md) | ⭐ 19 | `development` |
| [Clinical Decision Algorithms](development/tools/328-clinical_decision_algorithms_638c8f4b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-decision-support/references/clinical_decision_algorithms.md) | ⭐ 19 | `development` |
| [Llm](development/tools/llm_1027c55b/) | [crestalnetwork/intentkit](https://raw.githubusercontent.com/crestalnetwork/intentkit/main/LLM.md) | 🔥 6.5k | `development` |

### Productivity (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Growth Pm](productivity/173-growth-pm_61990fbd/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/growth-pm.md) | ⭐ 19 | `productivity` |
| [Instructions](productivity/142-instructions_8e70d675/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/arboreto/instructions.md) | ⭐ 19 | `productivity` |
| [Instructions](productivity/142-instructions_8423dce7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/get-available-resources/instructions.md) | ⭐ 19 | `productivity` |
| [Instructions](productivity/142-instructions_1cb797ea/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/molfeat/instructions.md) | ⭐ 19 | `productivity` |
| [Instructions](productivity/142-instructions_592d0bff/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/qutip/instructions.md) | ⭐ 19 | `productivity` |
| [Instructions](productivity/142-instructions_34320250/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/stable-baselines3/instructions.md) | ⭐ 19 | `productivity` |
| [Instructions](productivity/142-instructions_a0b1c774/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/vaex/instructions.md) | ⭐ 19 | `productivity` |

### Research (43 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Directory](research/258-directory_a0d73337/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/directory.md) | ⭐ 19 | `research` |
| [Researcher](research/259-researcher_df719c26/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/docs/personas/researcher.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_8e11ef01/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/citation-management/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_4a81772a/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/diffdock/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_bc174929/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/generate-image/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_88dc1368/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/hypothesis-generation/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_edf455e5/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/latex-posters/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_a8a509c6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/literature-review/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_eb5fb503/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/peer-review/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_821ff852/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/perplexity-search/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_2649f6c0/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/pytdc/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_947dbffb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_0b8463b1/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-lookup/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_0c7b82ed/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scholar-evaluation/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_e9bc7902/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-brainstorming/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_1b858180/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-critical-thinking/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_770aa284/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-schematics/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_0f8a0918/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-slides/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_c43a28a7/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-writing/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_9688f3cb/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/statistical-analysis/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_51b72599/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/torch_geometric/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_29e13b31/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/torchdrug/instructions.md) | ⭐ 19 | `research` |
| [Instructions](research/217-instructions_54f0b8e2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/treatment-plans/instructions.md) | ⭐ 19 | `research` |
| [Api Reference](research/007-api_reference_3188878c/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/brenda-database/references/api_reference.md) | ⭐ 19 | `research` |
| [Bibtex Formatting](research/010-bibtex_formatting_72d38876/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/citation-management/references/bibtex_formatting.md) | ⭐ 19 | `research` |
| [Citation Validation](research/014-citation_validation_7b76ecf9/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/citation-management/references/citation_validation.md) | ⭐ 19 | `research` |
| [Google Scholar Search](research/053-google_scholar_search_d741455b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/citation-management/references/google_scholar_search.md) | ⭐ 19 | `research` |
| [Metadata Extraction](research/067-metadata_extraction_71ec43da/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/citation-management/references/metadata_extraction.md) | ⭐ 19 | `research` |
| [Pubmed Search](research/094-pubmed_search_5b0933e6/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/citation-management/references/pubmed_search.md) | ⭐ 19 | `research` |
| [Evidence Synthesis](research/045-evidence_synthesis_8b70f676/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-decision-support/references/evidence_synthesis.md) | ⭐ 19 | `research` |
| [Treatment Recommendations](research/126-treatment_recommendations_f101db52/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-decision-support/references/treatment_recommendations.md) | ⭐ 19 | `research` |
| [Peer Review Standards](research/083-peer_review_standards_87b1180b/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/peer_review_standards.md) | ⭐ 19 | `research` |
| [Regulatory Compliance](research/105-regulatory_compliance_5caa7518/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/clinical-reports/references/regulatory_compliance.md) | ⭐ 19 | `research` |
| [Latex Poster Packages](research/064-latex_poster_packages_564e45c3/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/latex-posters/references/latex_poster_packages.md) | ⭐ 19 | `research` |
| [Poster Design Principles](research/086-poster_design_principles_d09d15ec/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/latex-posters/references/poster_design_principles.md) | ⭐ 19 | `research` |
| [Poster Layout Design](research/087-poster_layout_design_0f4e15b2/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/latex-posters/references/poster_layout_design.md) | ⭐ 19 | `research` |
| [Darpa Guidelines](research/028-darpa_guidelines_547536ad/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/references/darpa_guidelines.md) | ⭐ 19 | `research` |
| [Doe Guidelines](research/040-doe_guidelines_c27d71dc/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/references/doe_guidelines.md) | ⭐ 19 | `research` |
| [Nih Guidelines](research/071-nih_guidelines_7cdbc861/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/references/nih_guidelines.md) | ⭐ 19 | `research` |
| [Nsf Guidelines](research/073-nsf_guidelines_10523438/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/references/nsf_guidelines.md) | ⭐ 19 | `research` |
| [Nstc Guidelines](research/074-nstc_guidelines_0e71e945/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/research-grants/references/nstc_guidelines.md) | ⭐ 19 | `research` |
| [Beamer Guide](research/008-beamer_guide_9084c4b9/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-slides/references/beamer_guide.md) | ⭐ 19 | `research` |
| [Professional Report Formatting](research/093-professional_report_formatting_3dcf76ea/) | [SkeneTechnologies/skene-cookbook](https://raw.githubusercontent.com/SkeneTechnologies/skene-cookbook/main/skills-library/reference/scientific/scientific-writing/references/professional_report_formatting.md) | ⭐ 19 | `research` |

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

*Last updated: 2026-02-14 05:24:07 UTC*
*Automatically maintained by SkillFlow*
