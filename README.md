# X-Skills

A curated collection of **210 AI-powered skills** organized into 15 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (18 skills)
- **Automation/Workflow** (10 skills)
- **Commercial** (7 skills)
- **Communication** (9 skills)
- **Content Creation** (18 skills)
- **Daily Assistant** (17 skills)
- **Data Analysis** (12 skills)
- **Development** (42 skills)
- **Development/Devops** (24 skills)
- **Development/Testing** (3 skills)
- **Development/Tools** (19 skills)
- **Investment** (8 skills)
- **Other** (2 skills)
- **Productivity** (6 skills)
- **Research** (15 skills)

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


### Automation/Scripting (18 skills)

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
| [Session](automation/scripting/session_17d91954/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/commands/session.md) | ⭐ 543 | `automation` |
| [Ux Audit](automation/scripting/ux-audit_63e6bfb2/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/commands/ux-audit.md) | ⭐ 543 | `automation` |

### Automation/Workflow (10 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Cross File Resolution](automation/workflow/144-cross-file-resolution_c15a46bb/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/cross-file-resolution.md) | ⭐ 2.9k | `automation` |
| [Macos Arm64E Workaround](automation/workflow/145-macos-arm64e-workaround_2b4e18ae/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/static-analysis/skills/codeql/references/macos-arm64e-workaround.md) | ⭐ 2.9k | `automation` |
| [Task](automation/workflow/146-task_dd739ef5/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/prompts/task.md) | ⭐ 2.9k | `automation` |
| [Skill](automation/workflow/002-name-skill_247714ad/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/slides-generator/SKILL.md) | ⭐ 11 | `automation` |
| [Lp Executor Guide](automation/workflow/139-lp_executor_guide_0040aeee/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/lp-agent/references/lp_executor_guide.md) | ⭐ 11 | `automation` |
| [Skill](automation/workflow/002-name-skill_4c881e3b/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/skills/system/create-plan/SKILL.md) | ⭐ 763 | `automation` |
| [Skill](automation/workflow/002-name-skill_0e789cbb/) | [kimtth/awesome-azure-openai-llm](https://raw.githubusercontent.com/kimtth/awesome-azure-openai-llm/main/.github/skills/classify-temp-entries-to-section/SKILL.md) | ⭐ 392 | `automation` |
| [Tanstack Start](automation/workflow/tanstack-start_6fca1b22/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/commands/tanstack-start.md) | ⭐ 543 | `automation` |
| [Skill](automation/workflow/name-skill_77d78063/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/design-assets/skills/favicon-gen/SKILL.md) | ⭐ 543 | `automation` |
| [Skill](automation/workflow/name-skill_b30c7212/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/design-assets/skills/image-processing/SKILL.md) | ⭐ 543 | `automation` |

### Commercial (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](commercial/210-name-skill_5fcd131e/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/SKILL.md) | ⭐ 2.9k | `commercial` |
| [Ir Analysis](commercial/372-ir-analysis_473b57e4/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/skills/zeroize-audit/references/ir-analysis.md) | ⭐ 2.9k | `commercial` |
| [17 The Star](commercial/373-17-the-star_0cc63419/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/17-the-star.md) | ⭐ 2.9k | `commercial` |
| [King Of Pentacles](commercial/374-king-of-pentacles_ee2168eb/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/pentacles/king-of-pentacles.md) | ⭐ 2.9k | `commercial` |
| [Claude](commercial/claude_10641fac/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/CLAUDE.md) | ⭐ 543 | `commercial` |
| [Skill](commercial/name-skill_e1fcc7c4/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/shopify/skills/shopify-setup/SKILL.md) | ⭐ 543 | `commercial` |
| [Prompting Guide](commercial/prompting-guide_03d1c91a/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/design-assets/skills/gemini-image-gen/references/prompting-guide.md) | ⭐ 543 | `commercial` |

### Communication (9 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](communication/127-name-skill_2b647276/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-cosmos-java/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/127-name-skill_45420053/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-arizeaiobservabilityeval-dotnet/SKILL.md) | 🔥 15.4k | `communication` |
| [Skill](communication/127-name-skill_373d16f0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/azure-mgmt-botservice-py/SKILL.md) | 🔥 15.4k | `communication` |
| [Foundations](communication/252-foundations_9e0eb90b/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/agentic-actions-auditor/skills/agentic-actions-auditor/references/foundations.md) | ⭐ 2.9k | `communication` |
| [Authentication](communication/016-authentication_432b0de4/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/user-guide/authentication.md) | ⭐ 177 | `communication` |
| [Readme Cn](communication/252-readme_cn_f4fb8d9b/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/README_CN.md) | ⭐ 763 | `communication` |
| [Skill](communication/127-name-skill_bb401c33/) | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/skills/system/get-chat-history/SKILL.md) | ⭐ 763 | `communication` |
| [Skill](communication/name-skill_facee7f5/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/writing/skills/aussie-business-english/SKILL.md) | ⭐ 543 | `communication` |
| [Recipes](communication/recipes_07b045a2/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/frontend/skills/shadcn-ui/references/recipes.md) | ⭐ 543 | `communication` |

### Content Creation (18 skills)

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
| [Skill](content-creation/name-skill_90a79ecb/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/process-siren/skills/improve-processes/SKILL.md) | ⭐ 20 | `content creation` |
| [Skill](content-creation/name-skill_c4d011f1/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/wordpress/skills/wordpress-setup/SKILL.md) | ⭐ 543 | `content creation` |

### Daily Assistant (17 skills)

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
| [Twitter V060](daily-assistant/twitter-v060_00f28c07/) | [doramirdor/NadirClaw](https://raw.githubusercontent.com/doramirdor/NadirClaw/main/docs/twitter-v060.md) | ⭐ 236 | `daily assistant` |

### Data Analysis (12 skills)

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
| [Ui Ux Pro Max Skill](data-analysis/ui-ux-pro-max-skill_b408db6a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/ai-design-tools/ui-ux-pro-max-skill.md) | ⭐ 20 | `data analysis` |
| [Skill](data-analysis/name-skill_a7129ee9/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/integrations/skills/google-chat-messages/SKILL.md) | ⭐ 543 | `data analysis` |

### Development (42 skills)

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
| [Research Agent Assessment](development/research-agent-assessment_c6d81a5d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/audits/research-agent-assessment.md) | ⭐ 20 | `development` |
| [Skill](development/name-skill_9d69f716/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/name-skill_07bba6a0/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/vite-flare-starter/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/name-skill_1891cee6/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/gemini-guide/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/name-skill_6c7e72be/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/github-release/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/name-skill_b0e0b759/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/skill-creator/SKILL.md) | ⭐ 543 | `development` |
| [Deployment](development/deployment_28c3c43e/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/references/deployment.md) | ⭐ 543 | `development` |
| [Skill](development/name-skill_4fd07ecf/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/design-assets/skills/gemini-image-gen/SKILL.md) | ⭐ 543 | `development` |
| [Skill](development/name-skill_b7b9f92f/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/responsiveness-check/SKILL.md) | ⭐ 543 | `development` |
| [Templates](development/templates_c5161a65/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/project-health/references/templates.md) | ⭐ 543 | `development` |

### Development/Devops (24 skills)

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
| [Skill](development/devops/name-skill_f2fb6f03/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/project-health/SKILL.md) | ⭐ 543 | `development` |
| [Architecture](development/devops/architecture_8a5a2c86/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/references/architecture.md) | ⭐ 543 | `development` |
| [Customization Guide](development/devops/customization-guide_69b157f0/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/vite-flare-starter/references/customization-guide.md) | ⭐ 543 | `development` |
| [Devto V060](development/devops/devto-v060_d158348a/) | [doramirdor/NadirClaw](https://raw.githubusercontent.com/doramirdor/NadirClaw/main/docs/devto-v060.md) | ⭐ 236 | `development` |
| [Permission Presets](development/devops/permission-presets_562dc9de/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/project-health/references/permission-presets.md) | ⭐ 543 | `development` |

### Development/Testing (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/testing/002-name-skill_0af17ce7/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/tdd-orchestrator/SKILL.md) | 🔥 15.4k | `development` |
| [Oauth Plugin Architecture](development/testing/084-oauth_plugin_architecture_8d4366d2/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/OAUTH_PLUGIN_ARCHITECTURE.md) | ⭐ 177 | `development` |
| [Auth Providers](development/testing/085-auth-providers_346cccf7/) | [CaddyGlow/ccproxy-api](https://raw.githubusercontent.com/CaddyGlow/ccproxy-api/main/docs/development/auth-providers.md) | ⭐ 177 | `development` |

### Development/Tools (19 skills)

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
| [Skill](development/tools/name-skill_d043a7e0/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/skills/gemini-peer-review/SKILL.md) | ⭐ 543 | `development` |
| [Project Health](development/tools/project-health_3eadc1b7/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/commands/project-health.md) | ⭐ 543 | `development` |
| [Release](development/tools/release_e34bdc1d/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/dev-tools/commands/release.md) | ⭐ 543 | `development` |

### Investment (8 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](investment/021-name-skill_3e0f379b/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/data-scientist/SKILL.md) | 🔥 15.4k | `investment` |
| [Skill](investment/021-name-skill_0db988ba/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/SKILL.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/049-decision-frameworks_95ad10a0/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/carrier-relationship-management/references/decision-frameworks.md) | 🔥 15.4k | `investment` |
| [Decision Frameworks](investment/049-decision-frameworks_30ceb8b9/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/energy-procurement/references/decision-frameworks.md) | 🔥 15.4k | `investment` |
| [Seven Of Cups](investment/052-seven-of-cups_852b35c0/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/cups/seven-of-cups.md) | ⭐ 2.9k | `investment` |
| [Nine Of Pentacles](investment/053-nine-of-pentacles_272b240f/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/pentacles/nine-of-pentacles.md) | ⭐ 2.9k | `investment` |
| [Skill](investment/021-name-skill_ea8c4627/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/connectors-available/SKILL.md) | ⭐ 11 | `investment` |
| [Server Functions](investment/server-functions_535f27be/) | [jezweb/claude-skills](https://raw.githubusercontent.com/jezweb/claude-skills/main/plugins/cloudflare/skills/tanstack-start/references/server-functions.md) | ⭐ 543 | `investment` |

### Other (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [04 The Emperor](other/036-04-the-emperor_1b48cba1/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/04-the-emperor.md) | ⭐ 2.9k | `other` |
| [Eight Of Swords](other/037-eight-of-swords_24e30c87/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/swords/eight-of-swords.md) | ⭐ 2.9k | `other` |

### Productivity (6 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Four Of Cups](productivity/174-four-of-cups_e7c80a02/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/cups/four-of-cups.md) | ⭐ 2.9k | `productivity` |
| [07 The Chariot](productivity/175-07-the-chariot_fd9f65cd/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/major/07-the-chariot.md) | ⭐ 2.9k | `productivity` |
| [Two Of Pentacles](productivity/176-two-of-pentacles_95edc988/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/pentacles/two-of-pentacles.md) | ⭐ 2.9k | `productivity` |
| [Knight Of Wands](productivity/177-knight-of-wands_4ae17136/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/wands/knight-of-wands.md) | ⭐ 2.9k | `productivity` |
| [Page Of Wands](productivity/178-page-of-wands_6da5cc15/) | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/let-fate-decide/skills/let-fate-decide/cards/wands/page-of-wands.md) | ⭐ 2.9k | `productivity` |
| [Lp Rebalancer Guide](productivity/174-lp_rebalancer_guide_63b97f02/) | [hummingbot/skills](https://raw.githubusercontent.com/hummingbot/skills/main/skills/lp-agent/references/lp_rebalancer_guide.md) | ⭐ 11 | `productivity` |

### Research (15 skills)

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
| [Completeness Report Plugin Creator](research/completeness-report-plugin-creator_bf91f9a2/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/audits/completeness-report-plugin-creator.md) | ⭐ 20 | `research` |

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

*Last updated: 2026-02-27 10:14:11 UTC*
*Automatically maintained by SkillFlow*
