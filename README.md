# X-Skills

A curated collection of **1115 AI-powered skills** organized into 15 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (7 skills)
- **Automation/Workflow** (62 skills)
- **Commercial** (60 skills)
- **Communication** (30 skills)
- **Content Creation** (88 skills)
- **Daily Assistant** (49 skills)
- **Data Analysis** (96 skills)
- **Development** (355 skills)
- **Development/Devops** (172 skills)
- **Development/Testing** (25 skills)
- **Development/Tools** (117 skills)
- **Investment** (9 skills)
- **Other** (1 skill)
- **Productivity** (10 skills)
- **Research** (34 skills)

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


### Automation/Scripting (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agent Calibration](automation/scripting/085-agent-calibration_a87c080f/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/agent-calibration.md) | ⭐ 61 | `automation` |
| [Skill](automation/scripting/003-name-skill_570982e2/) | [ALBEDO-TABAI/lets-go-rss](https://raw.githubusercontent.com/ALBEDO-TABAI/lets-go-rss/main/SKILL.md) | ⭐ 31 | `automation` |
| [Skill](automation/scripting/003-name-skill_28ff4f04/) | [artwist-polyakov/polyakov-claude-skills](https://raw.githubusercontent.com/artwist-polyakov/polyakov-claude-skills/main/plugins/yandex-search-api/skills/yandex-search-api/SKILL.md) | ⭐ 38 | `automation` |
| [01 Configuration](automation/scripting/080-01-configuration_398bf016/) | [volcengine/OpenViking](https://raw.githubusercontent.com/volcengine/OpenViking/main/docs/en/guides/01-configuration.md) | ⭐ 3.2k | `automation` |
| [Ecomode](automation/scripting/086-ecomode_f42e5d09/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/ecomode.md) | ⭐ 10 | `automation` |
| [Fleet](automation/scripting/086-fleet_1a9c6299/) | [TechNickAI/openclaw-config](https://raw.githubusercontent.com/TechNickAI/openclaw-config/main/.claude/commands/fleet.md) | ⭐ 13 | `automation` |
| [Vulnerability Databases](automation/scripting/086-vulnerability_databases_01d8449e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/time-aware-dependency-cve-scanner/references/vulnerability_databases.md) | ⭐ 10 | `automation` |

### Automation/Workflow (62 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Vfunc Sig](automation/workflow/136-vfunc_sig_9dfecfa2/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.serena/memories/vfunc_sig.md) | ⭐ 19 | `automation` |
| [Skill](automation/workflow/002-name-skill_2854fe0b/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-deploy/SKILL.md) | ⭐ 80 | `automation` |
| [Cli Guide](automation/workflow/137-cli-guide_69513c0f/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-ai-agentscript/references/cli-guide.md) | ⭐ 80 | `automation` |
| [Export Import Architecture](automation/workflow/148-export-import-architecture_1fc68b51/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/export-import-architecture.md) | ⭐ 3.3k | `automation` |
| [Export Import Tutorial](automation/workflow/149-export-import-tutorial_372778e0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import-tutorial.md) | ⭐ 3.3k | `automation` |
| [Export Import](automation/workflow/150-export-import_af809418/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import.md) | ⭐ 3.3k | `automation` |
| [Mcp Cli](automation/workflow/151-mcp-cli_dbea27bc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/mcp-cli.md) | ⭐ 3.3k | `automation` |
| [Libreoffice Server](automation/workflow/152-libreoffice-server_2c0ac43b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/libreoffice-server.md) | ⭐ 3.3k | `automation` |
| [Skill](automation/workflow/002-name-skill_15b7286f/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/gas-fee-optimizer/skills/optimizing-gas-fees/SKILL.md) | ⭐ 1.4k | `automation` |
| [Readme.Ja](automation/workflow/040-readmeja_c8bfa185/) | [japan1988/multi-agent-mediation](https://raw.githubusercontent.com/japan1988/multi-agent-mediation/main/README.ja.md) | ⭐ 28 | `automation` |
| [Custom Workflows](automation/workflow/133-custom-workflows_4c3bac3c/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/advanced/custom-workflows.md) | 🔥 9.7k | `automation` |
| [04 Next Steps](automation/workflow/134-04-next-steps_12451d6e/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/04-next-steps.md) | 🔥 9.7k | `automation` |
| [04 Packaging](automation/workflow/135-04-packaging_49de9086/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/04-packaging.md) | 🔥 9.7k | `automation` |
| [Detecting Llm Hallucinations In Ci](automation/workflow/133-detecting-llm-hallucinations-in-ci_da061733/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/guides/detecting-llm-hallucinations-in-ci.md) | ⭐ 43 | `automation` |
| [Github Workflows](automation/workflow/135-github_workflows_50507981/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/GITHUB_WORKFLOWS.md) | ⭐ 151 | `automation` |
| [Release](automation/workflow/054-release_3485bb27/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/.claude/commands/github/release.md) | ⭐ 151 | `automation` |
| [Skill](automation/workflow/002-name-skill_80557270/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/complete-milestone/SKILL.md) | ⭐ 18 | `automation` |
| [Skill](automation/workflow/002-name-skill_927abc66/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/create-milestone/SKILL.md) | ⭐ 18 | `automation` |
| [Skill](automation/workflow/002-name-skill_a46079fe/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/start-milestone/SKILL.md) | ⭐ 18 | `automation` |
| [Skill](automation/workflow/002-name-skill_e93ca9fa/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/development/complete-implementation/SKILL.md) | ⭐ 18 | `automation` |
| [Skill](automation/workflow/002-name-skill_79eacb2d/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-bootstrap-fewshot/SKILL.md) | ⭐ 38 | `automation` |
| [Skill](automation/workflow/002-name-skill_e7c99d7c/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-gepa-reflective/SKILL.md) | ⭐ 38 | `automation` |
| [Skill](automation/workflow/002-name-skill_cf8e98c4/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-miprov2-optimizer/SKILL.md) | ⭐ 38 | `automation` |
| [Skill](automation/workflow/002-name-skill_36fd08e8/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-rag-pipeline/SKILL.md) | ⭐ 38 | `automation` |
| [Skill](automation/workflow/002-name-skill_91a7f888/) | [tradermonty/claude-trading-skills](https://raw.githubusercontent.com/tradermonty/claude-trading-skills/main/skills/edge-candidate-agent/SKILL.md) | ⭐ 41 | `automation` |
| [Skill](automation/workflow/002-name-skill_5b34f272/) | [tradermonty/claude-trading-skills](https://raw.githubusercontent.com/tradermonty/claude-trading-skills/main/skills/vcp-screener/SKILL.md) | ⭐ 41 | `automation` |
| [Skill](automation/workflow/002-name-skill_7f260aa5/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-memory/skills/past-conversations/SKILL.md) | ⭐ 34 | `automation` |
| [Skill](automation/workflow/002-name-skill_d59eab83/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-utilities/skills/youtube-research/SKILL.md) | ⭐ 34 | `automation` |
| [Frontmatter Options](automation/workflow/133-frontmatter-options_58425037/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-skills/skills/skill-repair/references/frontmatter-options.md) | ⭐ 34 | `automation` |
| [Cli Reference](automation/workflow/134-cli-reference_ef10a01a/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-utilities/skills/youtube-research/references/cli-reference.md) | ⭐ 34 | `automation` |
| [Readme Flat Skills Created](automation/workflow/138-readme_flat_skills_created_e05fd193/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_CREATED.md) | 🔥 24.5k | `automation` |
| [Skill](automation/workflow/002-name-skill_1209e1d0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/bug-history-summarizer/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_f0b22b17/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cd-pipeline-generator/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_790d366d/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/ci-pipeline-synthesizer/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_af396e21/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/critical-interval-security-checker/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_9d2db4ac/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/design-smell-detector/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_7a87c712/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interval-difference-analyzer/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_c200f5ba/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/markdown-document-structurer/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_70f25ca2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-carrying-code-generator/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_25b1e455/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-failure-explainer/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_0c5cc123/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-skeleton-generator/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_ea28e9c9/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/regression-consistency-checker/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_6f6be860/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-coverage-checker/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_529d2618/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rtl-property-inference/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_9de1f003/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-szz-analyzer/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_e18f1be4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-reasoning-verifier/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_920f39fb/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verified-spec-code-mapper/SKILL.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_cc708543/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/vulnerability-pattern-matcher/SKILL.md) | ⭐ 10 | `automation` |
| [Interval Analysis](automation/workflow/133-interval_analysis_38a233b2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interval-difference-analyzer/references/interval_analysis.md) | ⭐ 10 | `automation` |
| [Isabelle Lemmas](automation/workflow/134-isabelle_lemmas_b858c64e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/lemma-discovery-assistant/references/isabelle_lemmas.md) | ⭐ 10 | `automation` |
| [Formal Tools](automation/workflow/135-formal_tools_dea54acb/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rtl-equivalence-checker/references/formal_tools.md) | ⭐ 10 | `automation` |
| [Skill](automation/workflow/002-name-skill_d13fa0ae/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/design-orchestration/SKILL.md) | 🔥 13.4k | `automation` |
| [Readme Flat Skills Updated](automation/workflow/135-readme_flat_skills_updated_7e94df01/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_UPDATED.md) | 🔥 24.4k | `automation` |
| [Readme Flat Skills Az](automation/workflow/135-readme_flat_skills_az_51532ca9/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_AZ.md) | 🔥 24.6k | `automation` |
| [Kernel Sh](automation/workflow/153-kernel-sh_b42fc39c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/agent-infrastructure/kernel-sh.md) | ⭐ 18 | `automation` |
| [Skill](automation/workflow/002-name-skill_8fb34807/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/development/implement-feature/SKILL.md) | ⭐ 18 | `automation` |
| [Automation](automation/workflow/032-automation_3dfb4b4c/) | [doobidoo/MCP-Context-Provider](https://raw.githubusercontent.com/doobidoo/MCP-Context-Provider/main/AUTOMATION.md) | ⭐ 24 | `automation` |
| [Build](automation/workflow/136-build_34d67a2d/) | [doobidoo/MCP-Context-Provider](https://raw.githubusercontent.com/doobidoo/MCP-Context-Provider/main/BUILD.md) | ⭐ 24 | `automation` |
| [Learning Guide](automation/workflow/137-learning_guide_9fc21bb7/) | [doobidoo/MCP-Context-Provider](https://raw.githubusercontent.com/doobidoo/MCP-Context-Provider/main/docs/phase3/LEARNING_GUIDE.md) | ⭐ 24 | `automation` |
| [Skill](automation/workflow/002-name-skill_d8e49039/) | [jim60105/copilot-prompt](https://raw.githubusercontent.com/jim60105/copilot-prompt/master/skills/nanobanana-restore/SKILL.md) | ⭐ 17 | `automation` |
| [Claude](automation/workflow/claude_f54cc682/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/CLAUDE.md) | ⭐ 10 | `automation` |
| [Tinyfish](automation/workflow/tinyfish_03b55251/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/agent-infrastructure/tinyfish.md) | ⭐ 20 | `automation` |

### Commercial (60 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Express Adapter](commercial/387-express-adapter_6db5888c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/express-adapter.md) | ⭐ 102 | `commercial` |
| [Passkey](commercial/388-passkey_b16def64/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/passkey.md) | ⭐ 102 | `commercial` |
| [Routing Patterns](commercial/389-routing-patterns_959a3406/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/routing-patterns.md) | ⭐ 102 | `commercial` |
| [Caching Strategies](commercial/390-caching-strategies_3decca32/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/caching-strategies.md) | ⭐ 102 | `products` |
| [Server Components](commercial/179-server-components_249b33ee/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/server-components.md) | ⭐ 102 | `commercial` |
| [Advanced](commercial/002-advanced_87d29bd8/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/advanced.md) | ⭐ 102 | `commercial` |
| [Migration 0.7.0](commercial/391-migration-070_0d56b710/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/MIGRATION-0.7.0.md) | ⭐ 3.3k | `commercial` |
| [Kubernetes](commercial/344-kubernetes_0493cab9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/kubernetes.md) | ⭐ 3.3k | `commercial` |
| [Developer Workstation](commercial/060-developer-workstation_2f62c709/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/developer-workstation.md) | ⭐ 3.3k | `commercial` |
| [Oauth Troubleshooting](commercial/392-oauth-troubleshooting_2786d167/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/oauth-troubleshooting.md) | ⭐ 3.3k | `commercial` |
| [Supported Databases](commercial/393-supported-databases_d8a56a4d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/supported-databases.md) | ⭐ 3.3k | `commercial` |
| [Index](commercial/102-index_833de47b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/index.md) | ⭐ 3.3k | `commercial` |
| [Grpc Services](commercial/394-grpc-services_a969a94f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/grpc-services.md) | ⭐ 3.3k | `commercial` |
| [Query Param Auth](commercial/395-query-param-auth_e000b191/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/query-param-auth.md) | ⭐ 3.3k | `commercial` |
| [Cards](commercial/378-cards_6c04640b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/cards.md) | ⭐ 117 | `commercial` |
| [Skill](commercial/210-name-skill_1b74b4fe/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/stored-procedure-generator/skills/generating-stored-procedures/SKILL.md) | ⭐ 1.4k | `commercial` |
| [Tutorials](commercial/368-tutorials_9acb070e/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TUTORIALS.md) | ⭐ 43 | `commercial` |
| [Ui Store Model Db Setting](commercial/367-ui_store_model_db_setting_e86381b4/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/ui_store_model_db_setting.md) | 🔥 36.1k | `commercial` |
| [Skill](commercial/210-name-skill_24f220a9/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-haystack-integration/SKILL.md) | ⭐ 38 | `commercial` |
| [Prompt Extraction](commercial/369-prompt-extraction_571db0db/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-haystack-integration/references/prompt-extraction.md) | ⭐ 38 | `commercial` |
| [00 Overview](commercial/368-00-overview_bf57f153/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/00-overview.md) | ⭐ 10 | `setup` `test` |
| [Cco](commercial/369-cco_e8197f0c/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/cco.md) | ⭐ 10 | `commercial` |
| [kc](commercial/370-kc_75e40a38/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/kc.md) | ⭐ 10 | `commercial` |
| [me](commercial/371-me_15e356ae/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/me.md) | ⭐ 10 | `commercial` |
| [ta](commercial/372-ta_4704bd09/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/ta.md) | ⭐ 10 | `commercial` |
| [Continue](commercial/373-continue_f379530a/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/continue.md) | ⭐ 10 | `commercial` |
| [Knowledge Copilot](commercial/374-knowledge-copilot_33ced865/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/knowledge-copilot.md) | ⭐ 10 | `commercial` |
| [Memory](commercial/375-memory_f09c9b79/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/memory.md) | ⭐ 10 | `commercial` |
| [Pause](commercial/376-pause_b44ab317/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/pause.md) | ⭐ 10 | `commercial` |
| [Setup Copilot](commercial/377-setup-copilot_be32540f/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/setup-copilot.md) | ⭐ 10 | `commercial` |
| [Setup Knowledge Sync](commercial/378-setup-knowledge-sync_9688f5ba/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/setup-knowledge-sync.md) | ⭐ 10 | `commercial` |
| [Setup Project](commercial/379-setup-project_aebf9fc8/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/setup-project.md) | ⭐ 10 | `commercial` |
| [Setup](commercial/380-setup_d9762a9a/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/setup.md) | ⭐ 10 | `commercial` |
| [Skills Approve](commercial/381-skills-approve_58a41431/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/skills-approve.md) | ⭐ 10 | `commercial` |
| [Update Copilot](commercial/382-update-copilot_89b98d8f/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/update-copilot.md) | ⭐ 10 | `commercial` |
| [Protocol Injection](commercial/383-protocol-injection_65b58405/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/hooks/protocol-injection.md) | ⭐ 10 | `commercial` |
| [00 Overview](commercial/368-00-overview_6550b0a9/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/10-architecture/00-overview.md) | ⭐ 10 | `commercial` |
| [04 Token Efficiency Playbook](commercial/384-04-token-efficiency-playbook_f25854fb/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/30-operations/04-token-efficiency-playbook.md) | ⭐ 10 | `commercial` |
| [00 Enhancement Features](commercial/385-00-enhancement-features_2f6b098e/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/00-enhancement-features.md) | ⭐ 10 | `agent-improvement` `typescript` `compilation` |
| [03 Auto Checkpoint Hooks](commercial/386-03-auto-checkpoint-hooks_9b1442c1/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/03-auto-checkpoint-hooks.md) | ⭐ 10 | `commercial` |
| [Web Security](commercial/387-web-security_71e93443/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/security/web-security.md) | ⭐ 10 | `security` `owasp` `web` |
| [Style Reference](commercial/388-style-reference_1406cd42/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-content/skills/image-generation/references/style-reference.md) | ⭐ 34 | `commercial` |
| [Exceptions](commercial/249-exceptions_3dbd3b44/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/api/exceptions.md) | ⭐ 106 | `commercial` |
| [Cooldown Enforcement](commercial/379-cooldown-enforcement_6206b24f/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/docs/development/cooldown-enforcement.md) | 🔥 24.5k | `commercial` |
| [Skill](commercial/210-name-skill_1689a2dc/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/ambiguity-detector/SKILL.md) | ⭐ 10 | `commercial` |
| [Skill](commercial/210-name-skill_4de8311f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/assertion-synthesizer/SKILL.md) | ⭐ 10 | `commercial` |
| [Skill](commercial/210-name-skill_d32d6cde/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-comment-generator/SKILL.md) | ⭐ 10 | `commercial` |
| [Skill](commercial/210-name-skill_7a71ebfa/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rollback-strategy-advisor/SKILL.md) | ⭐ 10 | `commercial` |
| [Skill](commercial/210-name-skill_cf7a2136/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/scenario-generator/SKILL.md) | ⭐ 10 | `commercial` |
| [Abstract Domains](commercial/368-abstract_domains_5d58dce8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-domain-explorer/references/abstract_domains.md) | ⭐ 10 | `commercial` |
| [Loop Invariants](commercial/369-loop_invariants_62b05ab3/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-invariant-generator/references/loop_invariants.md) | ⭐ 10 | `commercial` |
| [Smell Patterns](commercial/370-smell-patterns_372f30a8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-smell-detector/references/smell-patterns.md) | ⭐ 10 | `commercial` |
| [Conflict Patterns](commercial/371-conflict_patterns_b5ce72e0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/conflict-analyzer/references/conflict_patterns.md) | ⭐ 10 | `commercial` |
| [Selection Guide](commercial/372-selection_guide_aae206f0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/design-pattern-suggestor/references/selection_guide.md) | ⭐ 10 | `commercial` |
| [Smell Catalog](commercial/373-smell_catalog_2cd2eff1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/design-smell-detector/references/smell_catalog.md) | ⭐ 10 | `commercial` |
| [Invariant Patterns](commercial/374-invariant-patterns_e28ad843/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/invariant-inference/references/invariant-patterns.md) | ⭐ 10 | `commercial` |
| [Architecture Patterns](commercial/375-architecture_patterns_abdb78e2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/legacy-code-summarizer/references/architecture_patterns.md) | ⭐ 10 | `commercial` |
| [Isabelle Syntax](commercial/376-isabelle_syntax_ba669e84/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/program-to-model-extractor/references/isabelle_syntax.md) | ⭐ 10 | `commercial` |
| [Platform Guides](commercial/377-platform_guides_5a3e2d4f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rollback-strategy-advisor/references/platform_guides.md) | ⭐ 10 | `commercial` |
| [Distributed Patterns](commercial/378-distributed_patterns_e5145e30/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tlaplus-spec-generator/references/distributed_patterns.md) | ⭐ 10 | `commercial` |

### Communication (30 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Schema](communication/267-schema_737de31c/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/SCHEMA.md) | ⭐ 102 | `communication` |
| [Nestjs Setup](communication/268-nestjs-setup_476bb950/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/nestjs-setup.md) | ⭐ 102 | `communication` |
| [Nextjs Setup](communication/269-nextjs-setup_c674b3e1/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/better-auth/references/nextjs-setup.md) | ⭐ 102 | `communication` |
| [Server Actions](communication/270-server-actions_4b85acf1/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/server-actions.md) | ⭐ 102 | `communication` |
| [Database Adapter](communication/271-database-adapter_0b083935/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-authentication/references/database-adapter.md) | ⭐ 102 | `communication` |
| [Rfc9728 Compliance](communication/272-rfc9728-compliance_0ee721f6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/rfc9728-compliance.md) | ⭐ 3.3k | `communication` |
| [Password Management](communication/273-password-management_bc5d8451/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/password-management.md) | ⭐ 3.3k | `communication` |
| [Inputs](communication/269-inputs_0e8fceb4/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/inputs.md) | ⭐ 117 | `communication` |
| [Toasts](communication/270-toasts_88d5e31e/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/toasts.md) | ⭐ 117 | `communication` |
| [Skill](communication/127-name-skill_5efa5b4b/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/SKILL.md) | ⭐ 61 | `team:support` |
| [Import Export](communication/251-import-export_9220633d/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/import-export.md) | ⭐ 61 | `communication` |
| [Sdk Usage](communication/252-sdk-usage_484ce649/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/sdk-usage.md) | ⭐ 61 | `tenant:acme` |
| [Cost Tracking](communication/251-cost_tracking_1bfec4b3/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/COST_TRACKING.md) | ⭐ 43 | `communication` |
| [Openclaw Integration](communication/253-openclaw_integration_9a5f9478/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/tutorials/openclaw_integration.md) | 🔥 36.1k | `communication` |
| [Tool Reference](communication/253-tool-reference_51413bad/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-memory/skills/past-conversations/references/tool-reference.md) | ⭐ 34 | `communication` |
| [Skill](communication/127-name-skill_933b4179/) | [TechNickAI/openclaw-config](https://raw.githubusercontent.com/TechNickAI/openclaw-config/main/skills/workflow-builder/SKILL.md) | ⭐ 13 | `communication` |
| [Agent](communication/253-agent_dcf1a7b4/) | [TechNickAI/openclaw-config](https://raw.githubusercontent.com/TechNickAI/openclaw-config/main/workflows/email-steward/AGENT.md) | ⭐ 13 | `communication` |
| [Skill](communication/127-name-skill_6ad0ed8f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/api-design-assistant/SKILL.md) | ⭐ 10 | `communication` |
| [Ambiguity Patterns](communication/253-ambiguity_patterns_6fbc78bd/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/ambiguity-detector/references/ambiguity_patterns.md) | ⭐ 10 | `communication` |
| [Best Practices](communication/254-best-practices_b5adde16/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/api-design-assistant/references/best-practices.md) | ⭐ 10 | `communication` |
| [Mutation Operators](communication/255-mutation_operators_eee9b805/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/behavioral-mutation-analyzer/references/mutation_operators.md) | ⭐ 10 | `communication` |
| [Time Intervals](communication/256-time_intervals_ff219d1f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/critical-interval-security-checker/references/time_intervals.md) | ⭐ 10 | `communication` |
| [Refactoring Strategies](communication/257-refactoring_strategies_3786b869/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/design-smell-detector/references/refactoring_strategies.md) | ⭐ 10 | `communication` |
| [Openapi Patterns](communication/258-openapi-patterns_7c7d0e06/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interface-specification-generator/references/openapi-patterns.md) | ⭐ 10 | `communication` |
| [Python Patterns](communication/259-python-patterns_2f31991f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/module-component-generator/references/python-patterns.md) | ⭐ 10 | `communication` |
| [Constraint Patterns](communication/260-constraint_patterns_3f4286c8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/nl-to-constraints/references/constraint_patterns.md) | ⭐ 10 | `communication` |
| [Best Practices](communication/261-best_practices_3273e09c/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/security-sensitive-path-instrumenter/references/best_practices.md) | ⭐ 10 | `communication` |
| [Security Events](communication/262-security_events_e42caca5/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/security-sensitive-path-instrumenter/references/security_events.md) | ⭐ 10 | `communication` |
| [Todowrite Usage Guide](communication/253-todowrite-usage-guide_04def9e0/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/docs/todowrite-usage-guide.md) | ⭐ 600 | `communication` |
| [Snapshot Refs](communication/274-snapshot-refs_edfb3fd8/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.agents/skills/agent-browser/references/snapshot-refs.md) | ⭐ 18 | `communication` |

### Content Creation (88 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](content-creation/049-name-skill_df94dda4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/drizzle-orm-patterns/SKILL.md) | ⭐ 102 | `drizzle` `orm` `database` |
| [Skill](content-creation/049-name-skill_feaebcdb/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/SKILL.md) | ⭐ 102 | `nextjs` `next.js` `app-router` |
| [Skill](content-creation/049-name-skill_2cfa63f4/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-data-fetching/SKILL.md) | ⭐ 102 | `posts` |
| [App Router Fundamentals](content-creation/405-app-router-fundamentals_5db2d98d/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/app-router-fundamentals.md) | ⭐ 102 | `content creation` |
| [Metadata Api](content-creation/406-metadata-api_464de7c9/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/metadata-api.md) | ⭐ 102 | `content creation` |
| [Nextjs16 Migration](content-creation/407-nextjs16-migration_807e73c6/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-app-router/references/nextjs16-migration.md) | ⭐ 102 | `content creation` |
| [React Query](content-creation/408-react-query_66bf81c6/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-data-fetching/references/REACT-QUERY.md) | ⭐ 102 | `content creation` |
| [Core Web Vitals](content-creation/409-core-web-vitals_c726a13b/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/core-web-vitals.md) | ⭐ 102 | `content creation` |
| [Image Optimization](content-creation/410-image-optimization_8a4da78f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/image-optimization.md) | ⭐ 102 | `content creation` |
| [Metadata Seo](content-creation/411-metadata-seo_eea42aa7/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/metadata-seo.md) | ⭐ 102 | `content creation` |
| [Nextjs 16 Patterns](content-creation/412-nextjs-16-patterns_1faa29b3/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-performance/references/nextjs-16-patterns.md) | ⭐ 102 | `content creation` |
| [Nextjs Config](content-creation/413-nextjs-config_d4dbc259/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/nextjs-config.md) | ⭐ 102 | `content creation` |
| [Troubleshooting](content-creation/110-troubleshooting_7c93e51b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/troubleshooting.md) | ⭐ 3.3k | `content creation` |
| [Passthrough](content-creation/414-passthrough_6dcf2b83/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/passthrough.md) | ⭐ 3.3k | `content creation` |
| [Tool Annotations](content-creation/415-tool-annotations_9f3db504/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/tool-annotations.md) | ⭐ 3.3k | `content creation` |
| [038 Experimental Rust Transport Backend](content-creation/416-038-experimental-rust-transport-backend_0eab4f20/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/038-experimental-rust-transport-backend.md) | ⭐ 3.3k | `content creation` |
| [Index](content-creation/019-index_c4edf73d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/kit/index.md) | ⭐ 3.3k | `content creation` |
| [Index](content-creation/019-index_1e66c9e7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/press/index.md) | ⭐ 3.3k | `content creation` |
| [Chunker Server](content-creation/417-chunker-server_334bf908/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/chunker-server.md) | ⭐ 3.3k | `content creation` |
| [Eval Server](content-creation/418-eval-server_35d409dc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/eval-server.md) | ⭐ 3.3k | `content creation` |
| [Url To Markdown Server](content-creation/419-url-to-markdown-server_55cb28c7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/url-to-markdown-server.md) | ⭐ 3.3k | `content creation` |
| [Skill](content-creation/049-name-skill_6b88bec6/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/SKILL.md) | ⭐ 117 | `package` |
| [Skill](content-creation/049-name-skill_e7bf1e39/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/SKILL.md) | ⭐ 117 | `content creation` |
| [Content Guidelines](content-creation/360-content-guidelines_b4f2960a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/references/content-guidelines.md) | ⭐ 117 | `content creation` |
| [Shiny Formatting](content-creation/361-shiny-formatting_30dd05e6/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/references/shiny-formatting.md) | ⭐ 117 | `content creation` |
| [Conversion Blogdown](content-creation/362-conversion-blogdown_038f094b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-blogdown.md) | ⭐ 117 | `data` |
| [Conversion Distill](content-creation/363-conversion-distill_69039bbc/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-distill.md) | ⭐ 117 | `content creation` |
| [Features](content-creation/109-features_40be1130/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/docs/FEATURES.md) | ⭐ 829 | `content creation` |
| [Skill](content-creation/049-name-skill_aa320b97/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/scientific-writing/SKILL.md) | ⭐ 829 | `content creation` |
| [V1.81.12](content-creation/351-v18112_b76af1d8/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/release_notes/v1.81.12.md) | 🔥 36.4k | `content creation` |
| [Skill](content-creation/049-name-skill_232ae8b1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/work-backlog-item/SKILL.md) | ⭐ 18 | `content creation` |
| [Whats New](content-creation/354-whats-new_6a65e564/) | [NTCoding/claude-skillz](https://raw.githubusercontent.com/NTCoding/claude-skillz/main/claude-code-updates/commands/whats-new.md) | ⭐ 239 | `content creation` |
| [Skill](content-creation/049-name-skill_a0bce514/) | [jim60105/copilot-prompt](https://raw.githubusercontent.com/jim60105/copilot-prompt/master/skills/create-blog-post/SKILL.md) | ⭐ 17 | `content creation` |
| [Writing Guidelines](content-creation/355-writing-guidelines_49868b1b/) | [jim60105/copilot-prompt](https://raw.githubusercontent.com/jim60105/copilot-prompt/master/skills/create-blog-post/references/writing-guidelines.md) | ⭐ 17 | `content creation` |
| [P5 Guide](content-creation/353-p5-guide_f110bc04/) | [dmccreary/claude-skills](https://raw.githubusercontent.com/dmccreary/claude-skills/main/skills/microsim-generator/references/p5-guide.md) | ⭐ 49 | `content creation` |
| [02 Scraping](content-creation/354-02-scraping_eb8fe18d/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/02-scraping.md) | 🔥 9.7k | `content creation` |
| [Usage](content-creation/042-usage_cd9a1109/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/archive/legacy/USAGE.md) | 🔥 9.7k | `content creation` |
| [Index](content-creation/019-index_a0e2ee39/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/index.md) | ⭐ 36 | `content creation` |
| [Toolsets](content-creation/360-toolsets_5e962330/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/api/toolsets.md) | ⭐ 36 | `content creation` |
| [Console Toolset](content-creation/420-console-toolset_a2a7b2fa/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/console-toolset.md) | ⭐ 36 | `content creation` |
| [Agents](content-creation/185-agents_7984f1c9/) | [dirnbauer/webconsulting-skills](https://raw.githubusercontent.com/dirnbauer/webconsulting-skills/main/AGENTS.md) | ⭐ 12 | `content creation` |
| [Skill](content-creation/049-name-skill_9adf125e/) | [dirnbauer/webconsulting-skills](https://raw.githubusercontent.com/dirnbauer/webconsulting-skills/main/skills/firecrawl/SKILL.md) | ⭐ 12 | `content creation` |
| [Skill Solrfal](content-creation/356-skill-solrfal_37d65441/) | [dirnbauer/webconsulting-skills](https://raw.githubusercontent.com/dirnbauer/webconsulting-skills/main/skills/typo3-solr/SKILL-SOLRFAL.md) | ⭐ 12 | `content creation` |
| [Migration V2](content-creation/360-migration_v2_e4def2aa/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/MIGRATION_V2.md) | ⭐ 151 | `content creation` |
| [Config Settings](content-creation/361-config_settings_6990ffd4/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/config_settings.md) | 🔥 36.1k | `cache_hit` `cache_key` `proxy_base_url` |
| [Skill](content-creation/049-name-skill_047ceb54/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/session-historian/SKILL.md) | ⭐ 18 | `content creation` |
| [Skill](content-creation/049-name-skill_0237a6a0/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/work-backlog-item/SKILL.md) | ⭐ 18 | `content creation` |
| [Development Guidelines](content-creation/364-development-guidelines_8c46d36d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/fastmcp-creator/skills/fastmcp-creator/references/development-guidelines.md) | ⭐ 18 | `content creation` |
| [cw](content-creation/355-cw_50bb228b/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/cw.md) | ⭐ 10 | `content creation` |
| [Sec](content-creation/357-sec_c729bce6/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/sec.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/049-name-skill_aff47338/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-content/skills/video-social/SKILL.md) | ⭐ 34 | `content creation` |
| [Skill](content-creation/049-name-skill_880c972f/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-skills/skills/skill-creator/SKILL.md) | ⭐ 34 | `content creation` |
| [Skill](content-creation/049-name-skill_e0c8b935/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-advanced-module-composition/SKILL.md) | ⭐ 38 | `content creation` |
| [Frontmatter Options](content-creation/362-frontmatter-options_2e27656a/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-skills/skills/skill-creator/references/frontmatter-options.md) | ⭐ 34 | `content creation` |
| [Platforms](content-creation/212-platforms_1737d684/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-utilities/skills/youtube-research/references/platforms.md) | ⭐ 34 | `content creation` |
| [Skill](content-creation/049-name-skill_15483b26/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/cfb-data/SKILL.md) | ⭐ 24 | `content creation` |
| [Skill](content-creation/049-name-skill_4213150f/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/sports-news/SKILL.md) | ⭐ 24 | `content creation` |
| [Skill](content-creation/049-name-skill_e26dd851/) | [TechNickAI/openclaw-config](https://raw.githubusercontent.com/TechNickAI/openclaw-config/main/skills/create-great-prompts/SKILL.md) | ⭐ 13 | `content creation` |
| [Skill](content-creation/049-name-skill_db245f9b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/git-master/SKILL.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/049-name-skill_35f1170f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/lemma-discovery-assistant/SKILL.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/049-name-skill_24937c5a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-trace-summarizer/SKILL.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/049-name-skill_8609e66e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/reference-searcher/SKILL.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/049-name-skill_730e35d5/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/strategic-planner/SKILL.md) | ⭐ 10 | `content creation` |
| [Resolution Strategies](content-creation/353-resolution_strategies_56d0b692/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/conflict-analyzer/references/resolution_strategies.md) | ⭐ 10 | `content creation` |
| [Python Deprecations](content-creation/354-python_deprecations_ce2aa423/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/deprecated-api-updater/references/python_deprecations.md) | ⭐ 10 | `content creation` |
| [Coq Lemmas](content-creation/355-coq_lemmas_7f1b2e55/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/lemma-discovery-assistant/references/coq_lemmas.md) | ⭐ 10 | `content creation` |
| [Proof Patterns](content-creation/356-proof_patterns_f7f352c2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/lemma-discovery-assistant/references/proof_patterns.md) | ⭐ 10 | `content creation` |
| [Failure Patterns](content-creation/357-failure_patterns_47486de5/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-failure-explainer/references/failure_patterns.md) | ⭐ 10 | `content creation` |
| [Coq Tactics](content-creation/358-coq_tactics_1b27841a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-skeleton-generator/references/coq_tactics.md) | ⭐ 10 | `content creation` |
| [Tactic Interpretation](content-creation/359-tactic_interpretation_6843b67e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-trace-summarizer/references/tactic_interpretation.md) | ⭐ 10 | `content creation` |
| [Common Violations](content-creation/360-common_violations_08c1ab71/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rtl-specification-consistency-checker/references/common_violations.md) | ⭐ 10 | `content creation` |
| [Ambiguity Resolution](content-creation/361-ambiguity_resolution_054c3383/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/specification-to-temporal-logic-generator/references/ambiguity_resolution.md) | ⭐ 10 | `content creation` |
| [Ctl Patterns](content-creation/362-ctl_patterns_c1f0cf77/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/specification-to-temporal-logic-generator/references/ctl_patterns.md) | ⭐ 10 | `content creation` |
| [Cwe Patterns](content-creation/363-cwe_patterns_ffe73f52/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-vulnerability-detector/references/cwe_patterns.md) | ⭐ 10 | `content creation` |
| [Coq Tactics](content-creation/358-coq_tactics_d29b3b60/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tactic-suggestion-assistant/references/coq_tactics.md) | ⭐ 10 | `content creation` |
| [Proof Patterns](content-creation/356-proof_patterns_877a5e99/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tactic-suggestion-assistant/references/proof_patterns.md) | ⭐ 10 | `content creation` |
| [Skill](content-creation/049-name-skill_4759d125/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/skills/LLM/SKILL.md) | ⭐ 600 | `content creation` |
| [Skill](content-creation/049-name-skill_13b3ef93/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/skills/TTS/SKILL.md) | ⭐ 600 | `content creation` |
| [Skill](content-creation/049-name-skill_7e56cadf/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/skills/video-generation/SKILL.md) | ⭐ 600 | `content creation` |
| [Skill](content-creation/049-name-skill_3b4f2543/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/skills/web-reader/SKILL.md) | ⭐ 600 | `content creation` |
| [Topic Specialist](content-creation/421-topic-specialist_df4d546d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/agents/topic-specialist.md) | ⭐ 18 | `content creation` |
| [Skill](content-creation/049-name-skill_4c8c6a8a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/work-backlog-item/SKILL.md) | ⭐ 18 | `content creation` |
| [Relief On Station](content-creation/422-relief-on-station_cde95eb3/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/damage-control/relief-on-station.md) | ⭐ 143 | `content creation` |
| [Skill](content-creation/049-name-skill_8d3e47e8/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/design-media/image-forge/SKILL.md) | ⭐ 10 | `content creation` |
| [Magick Reference](content-creation/353-magick-reference_7febcf3f/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/design-media/image-forge/references/magick-reference.md) | ⭐ 10 | `content creation` |
| [Recipes](content-creation/354-recipes_24c9be07/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/design-media/image-forge/references/recipes.md) | ⭐ 10 | `content creation` |
| [Claude](content-creation/007-claude_fae8645d/) | [calderbuild/agentcut](https://raw.githubusercontent.com/calderbuild/agentcut/main/CLAUDE.md) | ⭐ 12 | `content creation` |
| [Paperdraw](content-creation/paperdraw_8e90ac11/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/developer-tools/paperdraw.md) | ⭐ 20 | `content creation` |

### Daily Assistant (49 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Package Configs](daily-assistant/294-package-configs_66bf674a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/package-configs.md) | ⭐ 102 | `daily assistant` |
| [Crewai](daily-assistant/295-crewai_c93ce999/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/crewai.md) | ⭐ 3.3k | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_165a91fb/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/create-release-checklist/SKILL.md) | ⭐ 117 | `daily assistant` |
| [Claude](daily-assistant/037-claude_60cb35f0/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/CLAUDE.md) | ⭐ 1.4k | `daily assistant` |
| [Codex Models](daily-assistant/270-codex-models_4cf972bb/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-orchestrator/references/codex-models.md) | ⭐ 10 | `daily assistant` |
| [Claude](daily-assistant/037-claude_53ac8506/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/CLAUDE.md) | ⭐ 151 | `daily assistant` |
| [Claude Md Guardian](daily-assistant/265-claude-md-guardian_3c6331d8/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/agent/claude-md-guardian.md) | ⭐ 151 | `daily assistant` |
| [Project Workflow.Draft](daily-assistant/271-project_workflowdraft_208d2544/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/project_workflow.draft.md) | ⭐ 18 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_ef7d33fa/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/development/add-new-feature/SKILL.md) | ⭐ 18 | `daily assistant` |
| [Config](daily-assistant/262-config_1bbda31c/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/config.md) | ⭐ 10 | `ecomode_config` |
| [Reflect](daily-assistant/263-reflect_d329d123/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/reflect.md) | ⭐ 10 | `daily assistant` |
| [Update Project](daily-assistant/264-update-project_cb2cb2b8/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/update-project.md) | ⭐ 10 | `daily assistant` |
| [03 Time Free Language Guide](daily-assistant/265-03-time-free-language-guide_9ea7470f/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/60-qa/03-time-free-language-guide.md) | ⭐ 10 | `daily assistant` |
| [Kubernetes](daily-assistant/266-kubernetes_9cfcbaab/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/devops/kubernetes.md) | ⭐ 10 | `kubernetes` `k8s` `containers` |
| [Skill](daily-assistant/032-name-skill_d5f7ff6c/) | [tradermonty/claude-trading-skills](https://raw.githubusercontent.com/tradermonty/claude-trading-skills/main/skills/kanchi-dividend-sop/SKILL.md) | ⭐ 41 | `daily assistant` |
| [Security Analysis V070](daily-assistant/267-security_analysis_v070_14e0aad1/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/SECURITY_ANALYSIS_v070.md) | ⭐ 34 | `daily assistant` |
| [Quality Checklist](daily-assistant/268-quality-checklist_aed25b68/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-skills/skills/skill-repair/references/quality-checklist.md) | ⭐ 34 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_3480447c/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/fastf1/SKILL.md) | ⭐ 24 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_35f731fa/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/football-data/SKILL.md) | ⭐ 24 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_a23716e1/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/mlb-data/SKILL.md) | ⭐ 24 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_df9bd73a/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/nba-data/SKILL.md) | ⭐ 24 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_c4d6be83/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/nfl-data/SKILL.md) | ⭐ 24 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_cc67ac0f/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/nhl-data/SKILL.md) | ⭐ 24 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_9307b4b7/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/tennis-data/SKILL.md) | ⭐ 24 | `daily assistant` |
| [Commands](daily-assistant/271-commands_29378f0a/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/fastf1/references/commands.md) | ⭐ 24 | `daily assistant` |
| [Commands](daily-assistant/271-commands_5dfbe33f/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/football-data/references/commands.md) | ⭐ 24 | `daily assistant` |
| [Schemas](daily-assistant/272-schemas_34ab5339/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/football-data/references/schemas.md) | ⭐ 24 | `daily assistant` |
| [Update Model](daily-assistant/262-update-model_2796dd59/) | [TechNickAI/openclaw-config](https://raw.githubusercontent.com/TechNickAI/openclaw-config/main/.claude/commands/update-model.md) | ⭐ 13 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_1882d3aa/) | [TechNickAI/openclaw-config](https://raw.githubusercontent.com/TechNickAI/openclaw-config/main/skills/librarian/SKILL.md) | ⭐ 13 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_76f22eb9/) | [TechNickAI/openclaw-config](https://raw.githubusercontent.com/TechNickAI/openclaw-config/main/skills/smart-delegation/SKILL.md) | ⭐ 13 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_c007a5d0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interface-specification-generator/SKILL.md) | ⭐ 10 | `daily assistant` |
| [Skill](daily-assistant/032-name-skill_90ef8c73/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/program-to-tlaplus-spec-generator/SKILL.md) | ⭐ 10 | `daily assistant` |
| [Verification Languages](daily-assistant/262-verification_languages_11c5b1da/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-invariant-generator/references/verification_languages.md) | ⭐ 10 | `daily assistant` |
| [Build Systems](daily-assistant/263-build-systems_4f54dc2b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/build-ci-migration-assistant/references/build-systems.md) | ⭐ 10 | `daily assistant` |
| [Kubernetes Patterns](daily-assistant/264-kubernetes_patterns_551806a2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/containerization-assistant/references/kubernetes_patterns.md) | ⭐ 10 | `daily assistant` |
| [Java Edge Cases](daily-assistant/265-java_edge_cases_f163cab6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/edge-case-generator/references/java_edge_cases.md) | ⭐ 10 | `daily assistant` |
| [Clarification Patterns](daily-assistant/266-clarification_patterns_8eea096a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-to-tlaplus-property-generator/references/clarification_patterns.md) | ⭐ 10 | `daily assistant` |
| [Requirement Patterns](daily-assistant/267-requirement_patterns_5931ff4c/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-to-tlaplus-property-generator/references/requirement_patterns.md) | ⭐ 10 | `daily assistant` |
| [Step1 Customer Problems](daily-assistant/263-step1-customer-problems_9c1de879/) | [RafaelGorski/Problem-Based-SRS](https://raw.githubusercontent.com/RafaelGorski/Problem-Based-SRS/main/skills/problem-based-srs/references/step1-customer-problems.md) | ⭐ 10 | `daily assistant` |
| [Tasks 5 Sam Error Recovery](daily-assistant/296-tasks-5-sam-error-recovery_f6d8494e/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/tasks-5-sam-error-recovery.md) | ⭐ 18 | `daily assistant` |
| [Profiling](daily-assistant/297-profiling_896fbd9e/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.agents/skills/agent-browser/references/profiling.md) | ⭐ 18 | `daily assistant` |
| [Reddit Reply Graduated Discipline](daily-assistant/298-reddit-reply-graduated-discipline_dc3c1c09/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/docs/reddit-reply-graduated-discipline.md) | ⭐ 143 | `daily assistant` |
| [Commendations](daily-assistant/299-commendations_e16b6ff4/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/commendations.md) | ⭐ 143 | `daily assistant` |
| [Squadron Composition](daily-assistant/300-squadron-composition_d5a66e2c/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/squadron-composition.md) | ⭐ 143 | `daily assistant` |
| [Escalation](daily-assistant/301-escalation_806cb18a/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/damage-control/escalation.md) | ⭐ 143 | `daily assistant` |
| [Hull Integrity](daily-assistant/302-hull-integrity_530e2263/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/damage-control/hull-integrity.md) | ⭐ 143 | `daily assistant` |
| [Man Overboard](daily-assistant/303-man-overboard_6d08a8b7/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/damage-control/man-overboard.md) | ⭐ 143 | `daily assistant` |
| [Scuttle And Reform](daily-assistant/304-scuttle-and-reform_496947da/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/damage-control/scuttle-and-reform.md) | ⭐ 143 | `daily assistant` |
| [Session Resumption](daily-assistant/305-session-resumption_9df13ad4/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/damage-control/session-resumption.md) | ⭐ 143 | `daily assistant` |

### Data Analysis (96 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](data-analysis/226-name-skill_2414c4db/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-core/skills/drawio-logical-diagrams/SKILL.md) | ⭐ 102 | `data analysis` |
| [Raw Typescript Lambda](data-analysis/498-raw-typescript-lambda_c2fffb80/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/raw-typescript-lambda.md) | ⭐ 102 | `data analysis` |
| [Performance Architecture](data-analysis/499-performance-architecture_f977d60e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/performance-architecture.md) | ⭐ 3.3k | `data analysis` |
| [Dcr](data-analysis/500-dcr_14b0e9e4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/dcr.md) | ⭐ 3.3k | `data analysis` |
| [Postgres Upgrade Process](data-analysis/501-postgres-upgrade-process_60a48da6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/postgres-upgrade-process.md) | ⭐ 3.3k | `data analysis` |
| [007 Pluggable Cache Backend](data-analysis/502-007-pluggable-cache-backend_6ee2c5f9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/007-pluggable-cache-backend.md) | ⭐ 3.3k | `data analysis` |
| [Index](data-analysis/113-index_974fad2a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/index.md) | ⭐ 3.3k | `data analysis` |
| [Internal Observability](data-analysis/503-internal-observability_328ec0d8/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/internal-observability.md) | ⭐ 3.3k | `data analysis` |
| [Plugins](data-analysis/396-plugins_5f80238a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/plugins.md) | ⭐ 3.3k | `data analysis` |
| [Csv Pandas Chat Server](data-analysis/504-csv-pandas-chat-server_454d219f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/csv-pandas-chat-server.md) | ⭐ 3.3k | `data analysis` |
| [Data Analysis Server](data-analysis/505-data-analysis-server_cef0d00a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/data-analysis-server.md) | ⭐ 3.3k | `data analysis` |
| [Docx Server](data-analysis/506-docx-server_a79fb0e0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/docx-server.md) | ⭐ 3.3k | `data analysis` |
| [Graphviz Server](data-analysis/507-graphviz-server_b48cf167/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/graphviz-server.md) | ⭐ 3.3k | `data analysis` |
| [Plotly Server](data-analysis/508-plotly-server_3921f69a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/plotly-server.md) | ⭐ 3.3k | `data analysis` |
| [Pptx Server](data-analysis/509-pptx-server_9c6d336a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/pptx-server.md) | ⭐ 3.3k | `data analysis` |
| [Python Sandbox Server](data-analysis/510-python-sandbox-server_563cf18c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/python-sandbox-server.md) | ⭐ 3.3k | `data analysis` |
| [Brand Yml In R](data-analysis/479-brand-yml-in-r_780fb3bb/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/brand-yml-in-r.md) | ⭐ 117 | `data analysis` |
| [Quarto](data-analysis/480-quarto_81427cd4/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/quarto.md) | ⭐ 117 | `data analysis` |
| [Shiny R](data-analysis/481-shiny-r_0109581d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/shiny-r.md) | ⭐ 117 | `data analysis` |
| [Skill](data-analysis/226-name-skill_533d0b3b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/SKILL.md) | ⭐ 117 | `data analysis` |
| [Callouts](data-analysis/482-callouts_e4d4e5e1/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/callouts.md) | ⭐ 117 | `data analysis` |
| [Conversion Bookdown](data-analysis/483-conversion-bookdown_d994073b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-bookdown.md) | ⭐ 117 | `data analysis` |
| [Conversion Rmarkdown](data-analysis/484-conversion-rmarkdown_aee6105e/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-rmarkdown.md) | ⭐ 117 | `data analysis` |
| [Cross References](data-analysis/485-cross-references_202cb6fe/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/cross-references.md) | ⭐ 117 | `data analysis` |
| [Diagrams](data-analysis/442-diagrams_c3d9252d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/diagrams.md) | ⭐ 117 | `data analysis` |
| [Shortcodes](data-analysis/486-shortcodes_3dd7acd7/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/shortcodes.md) | ⭐ 117 | `data analysis` |
| [Tables](data-analysis/487-tables_15e26f8d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/tables.md) | ⭐ 117 | `data analysis` |
| [Accordions](data-analysis/488-accordions_c86bfde7/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/accordions.md) | ⭐ 117 | `data analysis` |
| [Navigation](data-analysis/489-navigation_d4984c12/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/navigation.md) | ⭐ 117 | `data analysis` |
| [Page Layouts](data-analysis/490-page-layouts_7aa85757/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/page-layouts.md) | ⭐ 117 | `data analysis` |
| [Theming](data-analysis/491-theming_77a63bbc/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/theming.md) | ⭐ 117 | `data analysis` |
| [Value Boxes](data-analysis/492-value-boxes_3b530b15/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/value-boxes.md) | ⭐ 117 | `data analysis` |
| [Claude](data-analysis/036-claude_8dccc45b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/CLAUDE.md) | ⭐ 18 | `data analysis` |
| [Orchestrator Discipline Patterns](data-analysis/486-orchestrator-discipline-patterns_95fa0f32/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/codebase/orchestrator-discipline-patterns.md) | ⭐ 18 | `data analysis` |
| [Research Scripts Refs](data-analysis/487-research-scripts-refs_9880003d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/plan/the-rewrite-room/research-scripts-refs.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_3854a340/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/research-curator/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_0637d133/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/swarm-primitives/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_a73903d5/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agentskill-kaizen/skills/kaizen-improvement/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_4667e05f/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/dasel/skills/dasel-reference/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_0f557dc7/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/audit-skill-lifecycle/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_5d22706c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/plugin-creator/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_7bfb55e1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/summarizer/skills/file-summarization/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_7326fbc9/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/summarizer/skills/image-summarization/SKILL.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_58a16c4c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/development-harness/skills/workflows/execution/SKILL.md) | ⭐ 18 | `data analysis` |
| [Prompt Optimization](data-analysis/488-prompt-optimization_f32f3366/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/the-rewrite-room/skills/the-rewrite-room/workflows/prompt-optimization.md) | ⭐ 18 | `data analysis` |
| [Research Utilities](data-analysis/489-research-utilities_a32adef1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/the-rewrite-room/skills/the-rewrite-room/workflows/research-utilities.md) | ⭐ 18 | `data analysis` |
| [Contextual Ai Documentation Optimizer Adapter](data-analysis/490-contextual-ai-documentation-optimizer-adapter_394973bd/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/the-rewrite-room/skills/the-rewrite-room/workflows/adapters/contextual-ai-documentation-optimizer-adapter.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/226-name-skill_bf9acbc8/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/mempool-analyzer/skills/analyzing-mempool/SKILL.md) | ⭐ 1.4k | `data analysis` |
| [Skill](data-analysis/226-name-skill_fecbac2d/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/nft-rarity-analyzer/skills/analyzing-nft-rarity/SKILL.md) | ⭐ 1.4k | `data analysis` |
| [Skill](data-analysis/226-name-skill_3ce7b1ff/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/core-development/claude-usage/SKILL.md) | ⭐ 10 | `data analysis` |
| [014 Security Headers Cors Middleware](data-analysis/481-014-security-headers-cors-middleware_f2c11a7a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/014-security-headers-cors-middleware.md) | ⭐ 3.3k | `data analysis` |
| [Skill](data-analysis/226-name-skill_95cf7107/) | [bowenliang123/md_exporter](https://raw.githubusercontent.com/bowenliang123/md_exporter/main/SKILL.md) | ⭐ 182 | `data analysis` |
| [Docker](data-analysis/511-docker_d9e87c50/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/docker.md) | ⭐ 36 | `data analysis` |
| [Skill](data-analysis/226-name-skill_ef75fc67/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/seed_skills/code-to-diagram/SKILL.md) | ⭐ 867 | `data analysis` |
| [Mermaid Patterns](data-analysis/478-mermaid-patterns_17eac38b/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/seed_skills/code-to-diagram/references/mermaid-patterns.md) | ⭐ 867 | `data analysis` |
| [Skill](data-analysis/226-name-skill_510e04d3/) | [bowenliang123/md_exporter](https://raw.githubusercontent.com/bowenliang123/md_exporter/main/SKILL.md) | ⭐ 182 | `data analysis` |
| [Skill Frontend](data-analysis/479-skill-frontend_03f46ca5/) | [dirnbauer/webconsulting-skills](https://raw.githubusercontent.com/dirnbauer/webconsulting-skills/main/skills/typo3-solr/SKILL-FRONTEND.md) | ⭐ 12 | `data analysis` |
| [Getting Started](data-analysis/480-getting_started_c7386cf9/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/GETTING_STARTED.md) | ⭐ 43 | `data analysis` |
| [Skill](data-analysis/226-name-skill_cf8ea07d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/group-items-to-milestone/SKILL.md) | ⭐ 18 | `data analysis` |
| [Readme Cn](data-analysis/436-readme_cn_29975317/) | [volcengine/OpenViking](https://raw.githubusercontent.com/volcengine/OpenViking/main/README_CN.md) | ⭐ 3.2k | `data analysis` |
| [Skill](data-analysis/226-name-skill_7a78b921/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-optimize-anything/SKILL.md) | ⭐ 38 | `data analysis` |
| [Orchestration Guide](data-analysis/478-orchestration_guide_2a0a46be/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/orchestrator/ORCHESTRATION_GUIDE.md) | ⭐ 10 | `data analysis` |
| [02 Interaction Design](data-analysis/479-02-interaction-design_55738d76/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/tasks/interactive-tui/02-interaction-design.md) | ⭐ 10 | `data analysis` |
| [Python Idioms](data-analysis/480-python-idioms_4902bf23/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/code/python-idioms.md) | ⭐ 10 | `python` `idioms` `patterns` |
| [React Patterns](data-analysis/481-react-patterns_fd65e047/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/code/react-patterns.md) | ⭐ 10 | `react` `hooks` `components` |
| [Skill](data-analysis/226-name-skill_3bb374d5/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-coding/skills/clean-branches/SKILL.md) | ⭐ 34 | `data analysis` |
| [Script Patterns](data-analysis/482-script-patterns_842143ab/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-skills/skills/skill-creator/references/script-patterns.md) | ⭐ 34 | `data analysis` |
| [Skill](data-analysis/226-name-skill_3ac3579a/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/cbb-data/SKILL.md) | ⭐ 24 | `data analysis` |
| [Skill](data-analysis/226-name-skill_154869da/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/golf-data/SKILL.md) | ⭐ 24 | `data analysis` |
| [Advanced](data-analysis/334-advanced_8b673c79/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/advanced.md) | ⭐ 106 | `data analysis` |
| [Quickstart](data-analysis/253-quickstart_25901d75/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/docs/docs/quickstart.md) | ⭐ 858 | `data analysis` |
| [Skill](data-analysis/226-name-skill_15ac6ffa/) | [bowenliang123/md_exporter](https://raw.githubusercontent.com/bowenliang123/markdown-exporter/main/SKILL.md) | ⭐ 183 | `data analysis` |
| [Skill](data-analysis/226-name-skill_e9b319b6/) | [wwwzhouhui/skills_collection](https://raw.githubusercontent.com/wwwzhouhui/skills_collection/main/github-readme-generator/SKILL.md) | ⭐ 107 | `data analysis` |
| [Skill](data-analysis/226-name-skill_cd720ae1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/build-ci-migration-assistant/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_cdcc854e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/component-boundary-identifier/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_b0a4f9f8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/conflict-analyzer/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_836250f1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/invariant-inference/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_70c3d5a4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/metamorphic-property-extractor/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_e195bd17/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/nl-to-constraints/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_bba613c7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/pseudocode-to-java-code/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_714a343b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/python-to-dafny-translator/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_6edd9188/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/reproduction-trace-instrumenter/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_7bb11b2b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-summarizer/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_e98d823e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rtl-specification-consistency-checker/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_801d35ad/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/specification-to-temporal-logic-generator/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_30b268f4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/spring-mvc-to-boot-migrator/SKILL.md) | ⭐ 10 | `data analysis` |
| [Skill](data-analysis/226-name-skill_23dd6bfa/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tlaplus-model-reduction/SKILL.md) | ⭐ 10 | `data analysis` |
| [Function Contracts](data-analysis/478-function_contracts_b53b7484/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-invariant-generator/references/function_contracts.md) | ⭐ 10 | `data analysis` |
| [Mutation Analysis Report](data-analysis/479-mutation_analysis_report_a17c88c6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/behavioral-mutation-analyzer/assets/mutation_analysis_report.md) | ⭐ 10 | `data analysis` |
| [Vulnerability Patterns](data-analysis/480-vulnerability_patterns_cbd37ad3/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/critical-interval-security-checker/references/vulnerability_patterns.md) | ⭐ 10 | `data analysis` |
| [Bug Patterns](data-analysis/481-bug_patterns_8f60e2c0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-bug-detector/references/bug_patterns.md) | ⭐ 10 | `data analysis` |
| [Semantic Analysis](data-analysis/482-semantic_analysis_f2c05cc4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-szz-analyzer/references/semantic_analysis.md) | ⭐ 10 | `data analysis` |
| [Detection Patterns](data-analysis/483-detection_patterns_888fe31a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-bug-detector/references/detection_patterns.md) | ⭐ 10 | `data analysis` |
| [Extraction Patterns](data-analysis/484-extraction_patterns_84b5b7b6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/traceability-matrix-generator/references/extraction_patterns.md) | ⭐ 10 | `data analysis` |
| [Claude](data-analysis/036-claude_671c19ce/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/CLAUDE.md) | ⭐ 18 | `data analysis` |
| [Skill](data-analysis/294-description-skill_db9c540a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/daily-releases/SKILL.md) | ⭐ 18 | `data analysis` |

### Development (355 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Multitenancy](development/1203-multitenancy_d26020b1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/multitenancy.md) | ⭐ 3.3k | `development` |
| [Plugins](development/698-plugins_b81593c7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins.md) | ⭐ 3.3k | `openai` `moderation` `content-safety` |
| [Selecting An Mcp Gateway](development/955-selecting-an-mcp-gateway_cec67c55/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/selecting-an-mcp-gateway.md) | ⭐ 3.3k | `development` |
| [Postgresql Schema Configuration](development/2977-postgresql-schema-configuration_9f6cd11d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/postgresql-schema-configuration.md) | ⭐ 3.3k | `development` |
| [Mcp Developer Guide Json Rpc](development/2978-mcp-developer-guide-json-rpc_2bfea8c1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/mcp-developer-guide-json-rpc.md) | ⭐ 3.3k | `development` |
| [Profiling](development/2979-profiling_d8b52edd/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/profiling.md) | ⭐ 3.3k | `development` |
| [Sso Entra Role Mapping](development/2980-sso-entra-role-mapping_c6294cc1/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-entra-role-mapping.md) | ⭐ 3.3k | `development` |
| [Sso Github Tutorial](development/2981-sso-github-tutorial_68c131fe/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-github-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Ibm Tutorial](development/2982-sso-ibm-tutorial_61b0c795/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-ibm-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Microsoft Entra Id Tutorial](development/2983-sso-microsoft-entra-id-tutorial_e579028c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-microsoft-entra-id-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Okta Tutorial](development/2984-sso-okta-tutorial_76fd8f99/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-okta-tutorial.md) | ⭐ 3.3k | `development` |
| [Teams](development/2985-teams_c57a4406/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/teams.md) | ⭐ 3.3k | `development` |
| [003 Expose Multi Transport Endpoints](development/2986-003-expose-multi-transport-endpoints_9637d88e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/003-expose-multi-transport-endpoints.md) | ⭐ 3.3k | `development` |
| [005 Vscode Devcontainer Support](development/2987-005-vscode-devcontainer-support_6e77ea8b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/005-vscode-devcontainer-support.md) | ⭐ 3.3k | `development` |
| [008 Federation Discovery](development/2988-008-federation-discovery_753ace5b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/008-federation-discovery.md) | ⭐ 3.3k | `development` |
| [010 Observability Prometheus](development/2989-010-observability-prometheus_fcec712d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/010-observability-prometheus.md) | ⭐ 3.3k | `development` |
| [011 Tool Federation](development/2990-011-tool-federation_24908a26/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/011-tool-federation.md) | ⭐ 3.3k | `development` |
| [016 Plugin Framework Ai Middleware](development/2991-016-plugin-framework-ai-middleware_f40a7f2c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/016-plugin-framework-ai-middleware.md) | ⭐ 3.3k | `development` |
| [017 Adopt Orjson Json Serialization](development/2992-017-adopt-orjson-json-serialization_7a66024e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/017-adopt-orjson-json-serialization.md) | ⭐ 3.3k | `development` |
| [022 Elicitation Passthrough Implementation](development/2993-022-elicitation-passthrough-implementation_288c5524/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/022-elicitation-passthrough-implementation.md) | ⭐ 3.3k | `development` |
| [024 Uvicorn Standard Extras](development/2994-024-uvicorn-standard-extras_6d24b0a8/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/024-uvicorn-standard-extras.md) | ⭐ 3.3k | `development` |
| [040 Flexible Admin Ui Sections](development/2995-040-flexible-admin-ui-sections_361375d2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/040-flexible-admin-ui-sections.md) | ⭐ 3.3k | `development` |
| [Gateway Hooks](development/2996-gateway-hooks_9aa359bb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins/gateway-hooks.md) | ⭐ 3.3k | `development` |
| [Security Hooks](development/2997-security-hooks_aa01ca2a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/plugins/security-hooks.md) | ⭐ 3.3k | `development` |
| [A2A](development/2998-a2a_9a2ca7a0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/a2a.md) | ⭐ 3.3k | `development` |
| [Cline](development/2999-cline_8713195d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/cline.md) | ⭐ 3.3k | `development` |
| [Llm Chat](development/3000-llm-chat_b4fec4c4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/llm-chat.md) | ⭐ 3.3k | `development` |
| [Index](development/468-index_79b6cf6d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/index.md) | ⭐ 3.3k | `custom` `filter` |
| [Lifecycle](development/2173-lifecycle_17a38751/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/lifecycle.md) | ⭐ 3.3k | `plugin` |
| [Rust Plugins](development/3001-rust-plugins_b6ecb04f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/rust-plugins.md) | ⭐ 3.3k | `development` |
| [Index](development/468-index_817e60d4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/index.md) | ⭐ 3.3k | `development` |
| [Calculator Server](development/3002-calculator-server_84a343ea/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/calculator-server.md) | ⭐ 3.3k | `development` |
| [Box](development/3003-box_1e5cdfee/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/box/box.md) | ⭐ 3.3k | `development` |
| [Monday Mcp](development/3004-monday-mcp_91796fb7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/monday/monday-mcp.md) | ⭐ 3.3k | `development` |
| [Skill](development/1178-name-skill_17e6beab/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/SKILL.md) | ⭐ 117 | `development` |
| [Brand Yml Spec](development/2903-brand-yml-spec_4c1c8800/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/brand-yml-spec.md) | ⭐ 117 | `development` |
| [Shiny Python](development/2904-shiny-python_6d034f3b/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/brand-yml/references/shiny-python.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_d8e79aed/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/posit-dev/critical-code-reviewer/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_e5da82b8/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/posit-dev/describe-design/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_06ff523f/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/posit-dev/pr-create/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_b817f91f/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_faade469/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cran-extrachecks/SKILL.md) | ⭐ 117 | `development` |
| [Skill](development/1178-name-skill_1a9af13a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/lifecycle/SKILL.md) | ⭐ 117 | `development` |
| [Tidyverse Formatting](development/2905-tidyverse-formatting_6f9d8aae/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/open-source/release-post/references/tidyverse-formatting.md) | ⭐ 117 | `package-name` `category` |
| [Code Cells](development/2906-code-cells_4883683d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/code-cells.md) | ⭐ 117 | `development` |
| [Conditional Content](development/2907-conditional-content_3c229e8a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conditional-content.md) | ⭐ 117 | `development` |
| [Conversion Xaringan](development/2908-conversion-xaringan_39d57dc1/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/conversion-xaringan.md) | ⭐ 117 | `development` |
| [Extensions](development/2909-extensions_7560b3bd/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/extensions.md) | ⭐ 117 | `development` |
| [Figures](development/2910-figures_c0e7f57e/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/figures.md) | ⭐ 117 | `development` |
| [Markdown Linting](development/2911-markdown-linting_0dff392d/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/markdown-linting.md) | ⭐ 117 | `development` |
| [Ansi Operations](development/2912-ansi-operations_3b86a218/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/ansi-operations.md) | ⭐ 117 | `development` |
| [Inline Markup](development/2913-inline-markup_1ed5cb03/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/inline-markup.md) | ⭐ 117 | `development` |
| [Progress](development/714-progress_d7998082/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/progress.md) | ⭐ 117 | `development` |
| [Lifecycle Stages](development/2914-lifecycle-stages_4f024f57/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/lifecycle/references/lifecycle-stages.md) | ⭐ 117 | `development` |
| [Migration](development/2204-migration_55087718/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/migration.md) | ⭐ 117 | `development` |
| [Feature Context Validate Orchestrator Discipline](development/2977-feature-context-validate-orchestrator-discipline_1964559c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/feature-context-validate-orchestrator-discipline.md) | ⭐ 18 | `development` |
| [Task File Format](development/2797-task_file_format_c1ef2af8/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/docs/TASK_FILE_FORMAT.md) | ⭐ 18 | `development` |
| [Claude](development/140-claude_76000a30/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/CLAUDE.md) | ⭐ 18 | `development` |
| [Research Agents](development/2978-research-agents_1d3d3c9d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/plan/the-rewrite-room/research-agents.md) | ⭐ 18 | `development` |
| [Research Architecture Pattern](development/2979-research-architecture-pattern_edc4bb85/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/plan/the-rewrite-room/research-architecture-pattern.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_3bc11637/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/agent-creator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_771f59ba/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/create-merge-request-changelog/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_4a8cf475/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/external-pattern-integrator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_45901b1b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/orchestrating-swarms/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_013658d1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/scientific-thinking/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_d1022c7b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/swarm-operations/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_d132f8a1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/swarm-patterns/SKILL.md) | ⭐ 18 | `development` |
| [Claude](development/140-claude_f69e49c5/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/scripts/CLAUDE.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_798aebc4/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agent-orchestration/skills/how-to-delegate/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_fde04101/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/bash-development/skills/bash-53-features/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_bbf386d9/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/bash-development/skills/bash-lint/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_8e1c6987/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/bash-development/skills/bash-logging/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_e6f15e2d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/commitlint/skills/commitlint/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_67bd8e87/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/holistic-linting/skills/holistic-linting-orchestrator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_fe04c1e3/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/holistic-linting/skills/holistic-linting/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_6e18624e/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/litellm/skills/litellm/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_4fd6044a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/llamafile/skills/llamafile/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_ff9fe8fb/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/agent-creator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_c856efe0/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/audit-skill-completeness/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_1640461c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/lint/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_5f0e587a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/memory-and-rules/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_becad5b3/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/token-launch-tracker/skills/tracking-token-launches/SKILL.md) | ⭐ 1.4k | `development` |
| [Response Api](development/1220-response_api_7468ad4b/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/response_api.md) | 🔥 36.4k | `development` |
| [Auto Routing](development/2892-auto_routing_e7ca64fd/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/auto_routing.md) | 🔥 36.4k | `development` |
| [Grooming 2026 02 21](development/2869-grooming-2026-02-21_f58716b3/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/grooming-reports/grooming-2026-02-21.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_7063d1fd/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/gh/SKILL.md) | ⭐ 18 | `development` |
| [Cto](development/2917-cto_087484d0/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-cto/agents/cto.md) | ⭐ 10 | `development` |
| [Subagent Patterns](development/2918-subagent-patterns_4b0a81e2/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-orchestrator/references/subagent-patterns.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_4a33313e/) | [artwist-polyakov/polyakov-claude-skills](https://raw.githubusercontent.com/artwist-polyakov/polyakov-claude-skills/main/plugins/codex-review/skills/codex-review/SKILL.md) | ⭐ 38 | `development` |
| [Documentation Updates Summary](development/2874-documentation_updates_summary_665ee62f/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/DOCUMENTATION_UPDATES_SUMMARY.md) | 🔥 9.7k | `development` |
| [Multi Source](development/2875-multi-source_9f256d80/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/advanced/multi-source.md) | 🔥 9.7k | `development` |
| [Unified Scraping](development/1108-unified_scraping_1d019192/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/features/UNIFIED_SCRAPING.md) | 🔥 9.7k | `development` |
| [01 Installation](development/2876-01-installation_39803a31/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/01-installation.md) | 🔥 9.7k | `development` |
| [03 Your First Skill](development/2877-03-your-first-skill_c9fe1c7e/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/03-your-first-skill.md) | 🔥 9.7k | `development` |
| [Cli Reference](development/1593-cli_reference_ba8e75e8/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/CLI_REFERENCE.md) | 🔥 9.7k | `development` |
| [Config Format](development/2878-config_format_e39c2327/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/CONFIG_FORMAT.md) | 🔥 9.7k | `development` |
| [06 Troubleshooting](development/2879-06-troubleshooting_78ddfebd/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/06-troubleshooting.md) | 🔥 9.7k | `development` |
| [Cli Reference](development/1593-cli_reference_2e61f4f8/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/zh-CN/reference/CLI_REFERENCE.md) | 🔥 9.7k | `development` |
| [Config Format](development/2878-config_format_af066045/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/zh-CN/reference/CONFIG_FORMAT.md) | 🔥 9.7k | `development` |
| [Claude](development/140-claude_8e8685ee/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/CLAUDE.md) | ⭐ 36 | `development` |
| [Abbreviations](development/3076-abbreviations_cdb8d3d2/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/includes/abbreviations.md) | ⭐ 36 | `development` |
| [Skill](development/1178-name-skill_d55098c1/) | [NTCoding/claude-skillz](https://raw.githubusercontent.com/NTCoding/claude-skillz/main/separation-of-concerns/SKILL.md) | ⭐ 240 | `development` |
| [Prompting Guide](development/730-prompting_guide_e84adb81/) | [promptdriven/pdd](https://raw.githubusercontent.com/promptdriven/pdd/main/docs/prompting_guide.md) | ⭐ 478 | `development` |
| [Prompting Guide](development/730-prompting_guide_7e3a5594/) | [promptdriven/pdd](https://raw.githubusercontent.com/promptdriven/pdd/main/pdd/docs/prompting_guide.md) | ⭐ 478 | `development` |
| [Setup With Gemini](development/970-setup_with_gemini_14e7ae1e/) | [promptdriven/pdd](https://raw.githubusercontent.com/promptdriven/pdd/main/SETUP_WITH_GEMINI.md) | ⭐ 478 | `development` |
| [Onboarding](development/635-onboarding_e59a6011/) | [promptdriven/pdd](https://raw.githubusercontent.com/promptdriven/pdd/main/docs/ONBOARDING.md) | ⭐ 478 | `development` |
| [Skill](development/1178-name-skill_58a38edd/) | [dirnbauer/webconsulting-skills](https://raw.githubusercontent.com/dirnbauer/webconsulting-skills/main/skills/refactor/SKILL.md) | ⭐ 12 | `development` |
| [Skill](development/1178-name-skill_541f68b0/) | [dirnbauer/webconsulting-skills](https://raw.githubusercontent.com/dirnbauer/webconsulting-skills/main/skills/typo3-solr/SKILL.md) | ⭐ 12 | `development` |
| [Adapters](development/2870-adapters_98b83b97/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/ADAPTERS.md) | ⭐ 43 | `development` |
| [Backend Requirements](development/2871-backend_requirements_4ae9bbe4/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/BACKEND_REQUIREMENTS.md) | ⭐ 43 | `development` |
| [Database Setup](development/2872-database_setup_b57196ca/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DATABASE_SETUP.md) | ⭐ 43 | `development` |
| [Faq](development/360-faq_97bfa995/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FAQ.md) | ⭐ 43 | `development` |
| [Framework Support](development/2873-framework_support_3e670e71/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FRAMEWORK_SUPPORT.md) | ⭐ 43 | `development` |
| [Trace Spec](development/2874-trace_spec_5795a259/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TRACE_SPEC.md) | ⭐ 43 | `development` |
| [Enhance Claude Md](development/2888-enhance-claude-md_3a2fa8ac/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/command/enhance-claude-md.md) | ⭐ 151 | `development` |
| [Branching Strategy](development/2889-branching_strategy_a98c62b9/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/BRANCHING_STRATEGY.md) | ⭐ 151 | `development` |
| [Ci Cd Fix Validation](development/2890-ci_cd_fix_validation_1fcdc8df/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/CI_CD_FIX_VALIDATION.md) | ⭐ 151 | `development` |
| [Installation](development/474-installation_65d262bf/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/INSTALLATION.md) | ⭐ 151 | `development` |
| [Quick Start](development/758-quick_start_beef973f/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/QUICK_START.md) | ⭐ 151 | `development` |
| [Troubleshooting](development/1097-troubleshooting_532ef3e4/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/TROUBLESHOOTING.md) | ⭐ 151 | `development` |
| [How To Use](development/451-how_to_use_abac3394/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/skill/HOW_TO_USE.md) | ⭐ 151 | `development` |
| [Skill](development/1178-name-skill_901ff1ab/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/skill/SKILL.md) | ⭐ 151 | `development` |
| [Commit Smart](development/2449-commit-smart_792ffa0e/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/.claude/commands/github/commit-smart.md) | ⭐ 151 | `development` |
| [Localai](development/2923-localai_42b73c5b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/llm-infrastructure/localai.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_92a46630/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/create-backlog-item/SKILL.md) | ⭐ 18 | `development` |
| [Issue Stories](development/2924-issue-stories_38ea1d07/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/gh/references/issue-stories.md) | ⭐ 18 | `development` |
| [Github Integration](development/2925-github-integration_38beca7d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/work-backlog-item/references/github-integration.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_c08c967a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agent-orchestration/skills/how-to-delegate/SKILL.md) | ⭐ 18 | `development` |
| [Community Practices](development/2926-community-practices_97872c86/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/fastmcp-creator/skills/fastmcp-creator/references/community-practices.md) | ⭐ 18 | `development` |
| [Evaluation Guide](development/2927-evaluation-guide_8b06b6e9/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/fastmcp-creator/skills/fastmcp-creator/references/evaluation-guide.md) | ⭐ 18 | `development` |
| [Mcp Best Practices](development/2928-mcp-best-practices_8ff83b58/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/fastmcp-creator/skills/fastmcp-creator/references/mcp-best-practices.md) | ⭐ 18 | `development` |
| [01 Overview](development/2703-01-overview_f7adeb76/) | [volcengine/OpenViking](https://raw.githubusercontent.com/volcengine/OpenViking/main/docs/zh/api/01-overview.md) | ⭐ 3.2k | `development` |
| [Faq](development/360-faq_1dfa4d66/) | [volcengine/OpenViking](https://raw.githubusercontent.com/volcengine/OpenViking/main/docs/zh/faq/faq.md) | ⭐ 3.2k | `development` |
| [Agents](development/028-agents_bf83e48f/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/AGENTS.md) | ⭐ 38 | `development` |
| [Dspy Framework](development/2880-dspy-framework_71928957/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/docs/dspy-framework.md) | ⭐ 38 | `development` |
| [Skill](development/1178-name-skill_e9cd0dfe/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-output-refinement-constraints/SKILL.md) | ⭐ 38 | `development` |
| [Skill](development/1178-name-skill_9babc8cf/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-react-agent-builder/SKILL.md) | ⭐ 38 | `development` |
| [Skill](development/1178-name-skill_65d28443/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/skills/dspy-signature-designer/SKILL.md) | ⭐ 38 | `development` |
| [Claude Reference](development/2885-claude_reference_e90f0b78/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/CLAUDE_REFERENCE.md) | ⭐ 10 | `development` |
| [Code Review Overview](development/2886-code_review_overview_7b5d10f4/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/CODE_REVIEW_OVERVIEW.md) | ⭐ 10 | `development` |
| [Project Context](development/2887-project-context_3c779b49/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/project-context.md) | ⭐ 10 | `development` |
| [Service Designer](development/2888-service-designer_3269e6ee/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/service-designer.md) | ⭐ 10 | `development` |
| [Ux Designer](development/1303-ux-designer_19082b95/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/ux-designer.md) | ⭐ 10 | `development` |
| [Development Rules](development/2889-development-rules_01de8421/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/reference/development-rules.md) | ⭐ 10 | `development` |
| [01 Agents](development/2890-01-agents_692ac406/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/10-architecture/01-agents.md) | ⭐ 10 | `development` |
| [02 Philosophy](development/2891-02-philosophy_ae785666/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/10-architecture/02-philosophy.md) | ⭐ 10 | `development` |
| [03 Decision Guide](development/2892-03-decision-guide_556a6b6d/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/10-architecture/03-decision-guide.md) | ⭐ 10 | `development` |
| [02 Customization](development/2893-02-customization_1fd255ef/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/20-configuration/02-customization.md) | ⭐ 10 | `development` |
| [02 Documentation Guide](development/2894-02-documentation-guide_28037300/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/30-operations/02-documentation-guide.md) | ⭐ 10 | `development` |
| [02 Api Reference](development/2895-02-api-reference_abbd76cb/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/02-api-reference.md) | ⭐ 10 | `development` |
| [03 Knowledge Sync](development/2896-03-knowledge-sync_0850e4bf/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/03-knowledge-sync.md) | ⭐ 10 | `development` |
| [04 Goal Driven Agents](development/2897-04-goal-driven-agents_75941d63/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/04-goal-driven-agents.md) | ⭐ 10 | `development` |
| [06 Opus 46 Capabilities](development/2898-06-opus-46-capabilities_9717c785/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/06-opus-46-capabilities.md) | ⭐ 10 | `development` |
| [Magic Keywords](development/2899-magic-keywords_518a7c2b/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/magic-keywords.md) | ⭐ 10 | `development` |
| [Progress Hud](development/2900-progress-hud_7ed6d284/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/progress-hud.md) | ⭐ 10 | `development` |
| [01 Framework Validation Strategy](development/2901-01-framework-validation-strategy_183331a6/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/60-qa/01-framework-validation-strategy.md) | ⭐ 10 | `development` |
| [00 Quick Reference](development/2902-00-quick-reference_4ff3f348/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/70-reference/00-quick-reference.md) | ⭐ 10 | `development` |
| [01 Usage Guide](development/2903-01-usage-guide_f94a5e9c/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/70-reference/01-usage-guide.md) | ⭐ 10 | `development` |
| [Triggers](development/2904-triggers_eefca6aa/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/mcp-servers/skills-copilot/TRIGGERS.md) | ⭐ 10 | `development` |
| [Deliverables](development/2905-deliverables_45b207a3/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/scripts/knowledge-sync/DELIVERABLES.md) | ⭐ 10 | `development` |
| [Quick Start](development/756-quick-start_5839d418/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/scripts/knowledge-sync/QUICK-START.md) | ⭐ 10 | `development` |
| [00 Plan Overview](development/2906-00-plan-overview_b35a7edc/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/tasks/interactive-tui/00-plan-overview.md) | ⭐ 10 | `development` |
| [01 Experience Design](development/2907-01-experience-design_717b1214/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/tasks/interactive-tui/01-experience-design.md) | ⭐ 10 | `development` |
| [04 Phase Tasks](development/2908-04-phase-tasks_60addeaf/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/tasks/interactive-tui/04-phase-tasks.md) | ⭐ 10 | `development` |
| [Javascript Patterns](development/2909-javascript-patterns_8aa06c4c/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/code/javascript-patterns.md) | ⭐ 10 | `javascript` `typescript` `nodejs` |
| [Api Docs](development/2759-api-docs_cd3abd1c/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/documentation/api-docs.md) | ⭐ 10 | `api` `documentation` `openapi` |
| [Crypto Patterns](development/2910-crypto-patterns_f252f989/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/security/crypto-patterns.md) | ⭐ 10 | `security` `cryptography` `encryption` |
| [Claude](development/140-claude_873a74f2/) | [tradermonty/claude-trading-skills](https://raw.githubusercontent.com/tradermonty/claude-trading-skills/main/CLAUDE.md) | ⭐ 41 | `development` |
| [Pipeline If V1](development/2911-pipeline_if_v1_190b7f52/) | [tradermonty/claude-trading-skills](https://raw.githubusercontent.com/tradermonty/claude-trading-skills/main/skills/edge-candidate-agent/references/pipeline_if_v1.md) | ⭐ 41 | `development` |
| [Claude](development/140-claude_e9228d7d/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/CLAUDE.md) | ⭐ 34 | `development` |
| [Skill](development/1178-name-skill_379c3eee/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-coding/skills/commit/SKILL.md) | ⭐ 34 | `development` |
| [Skill](development/1178-name-skill_9cba2968/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-memory/skills/extract-learnings/SKILL.md) | ⭐ 34 | `development` |
| [Badges](development/096-badges_7a989932/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-coding/skills/readme-maker/references/badges.md) | ⭐ 34 | `development` |
| [Sqlcipher Install](development/1539-sqlcipher_install_a9f7b6cf/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/SQLCIPHER_INSTALL.md) | ⭐ 4.0k | `development` |
| [Agents](development/028-agents_fa9578ce/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/AGENTS.md) | ⭐ 106 | `development` |
| [Warp](development/2871-warp_d095774e/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/seed_skills/humanizer/WARP.md) | ⭐ 858 | `development` |
| [Claude](development/140-claude_560359e5/) | [frmoretto/hardstop](https://raw.githubusercontent.com/frmoretto/hardstop/main/CLAUDE.md) | ⭐ 16 | `development` |
| [Readme Zh](development/2870-readme-zh_1d655b75/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/README-zh.md) | ⭐ 10 | `development` |
| [Deployment](development/272-deployment_c894faca/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skill-manager/DEPLOYMENT.md) | ⭐ 10 | `development` |
| [Github Pages Setup](development/2871-github-pages-setup_3ea86591/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skill-manager/GITHUB-PAGES-SETUP.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_4cb1e073/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-domain-explorer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_41240558/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-invariant-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_8aaa9d19/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-state-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_d131277b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-trace-summarizer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_cf753f0e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/api-documentation-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_552fdaca/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/behavior-preservation-checker/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_039c7923/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/behavioral-mutation-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_3e035e64/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/bisect-aware-instrumentation/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_719e8f6a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/bug-localization/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_9f69aa96/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/bug-to-patch-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_3e1e46b9/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/c-cpp-to-lean4-translator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_2ec9c676/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-change-summarizer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_1088e1e6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-completion-semantic-constraints/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_9f522f28/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-pattern-extractor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_e03f470f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-refactoring-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_34483783/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-repair-generation-combo/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_f7f5b1b1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-review-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_0bec5f8f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-search-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_9232fa36/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-summarizer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_5ec17124/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-translation/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_9df9a74d/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cve-reachability-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_09a577ca/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/dead-code-removal/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_b1cc8add/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/dependency-resolver/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_e3f28d67/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/design-pattern-suggestor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_97859a59/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/edge-case-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_74dcf989/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/exploitability-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_2ea0b7b3/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/failure-oriented-instrumentation/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_4273f85a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/formal-spec-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_6e888f33/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/frontend-ui-ux/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_ab3b61e3/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/function-class-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_c38f99c7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/fuzzing-input-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_0c83b772/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/git-bisect-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_6d4ef9c3/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/github-triage/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_abb965fe/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/imperative-to-coq-model-extractor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_4c5c35eb/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/incremental-java-programmer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_75bac2fd/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interface-contract-verifier/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_17cb7ac9/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/issue-report-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_b5f041fc/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/lsp-refactoring/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_e30ce380/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/model-guided-code-repair/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_dc35d689/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/modular-code-enforcement/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_7257094b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/module-level-code-translator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_84894089/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/playwright-automation/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_6fac8170/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/program-to-model-extractor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_85d7bf18/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-refactoring-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_081380af/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/pseudocode-extractor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_e246bb36/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/python-repo-quickstart/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_36c47dfe/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/python-to-lean4-translator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_97af31c7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/readme-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_b885c309/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/release-change-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_6687e21f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/release-notes-writer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_c74442e6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/replay-oriented-instrumentation/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_382004ca/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-comparison-reporter/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_45ae72f0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-enhancer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_328e3976/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-summary/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_b1cae492/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rtl-equivalence-checker/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_0747239a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/runtime-error-explainer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_647ce833/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/security-patch-advisor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_725edb80/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/security-sensitive-path-instrumenter/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_8977aa19/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-bug-detector/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_a99fc336/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/session-handoff/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_7d7e61cb/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/smart-mutation-operator-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_d6ab2ea2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/smv-model-extractor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_7c185aa2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/specification-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_7fdc7c5a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-bug-detector/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_b2bb0fdd/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-vulnerability-detector/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_53e00a5f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/symbolic-execution-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_11102597/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/szz-bug-identifier/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_af3b8b54/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/technical-debt-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_635bd5c1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/time-aware-dependency-cve-scanner/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_eed33985/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tlaplus-guided-code-repair/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_b8705431/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tlaplus-spec-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_e113998e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/traceability-matrix-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_28a1c731/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verification-boundary-reporter/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_00f2dfd7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verified-pseudocode-extractor/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/1178-name-skill_f78acac1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/vulnerability-root-cause-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Api Reference](development/051-api_reference_47f35914/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-domain-explorer/references/api_reference.md) | ⭐ 10 | `development` |
| [Abstract Interpretation](development/2872-abstract_interpretation_dc1d8434/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/abstract-trace-summarizer/references/abstract_interpretation.md) | ⭐ 10 | `development` |
| [Frama C Integration](development/2873-frama_c_integration_3d5496e9/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/acsl-annotation-assistant/references/frama_c_integration.md) | ⭐ 10 | `development` |
| [Comparison Techniques](development/2874-comparison_techniques_00532933/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/behavior-preservation-checker/references/comparison_techniques.md) | ⭐ 10 | `development` |
| [Difference Patterns](development/2875-difference_patterns_003cae23/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/behavior-preservation-checker/references/difference_patterns.md) | ⭐ 10 | `development` |
| [Exit Codes](development/2876-exit_codes_0f6c6de3/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/bisect-aware-instrumentation/references/exit_codes.md) | ⭐ 10 | `development` |
| [Git Bisect Guide](development/2877-git_bisect_guide_24b6e3a4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/bisect-aware-instrumentation/references/git_bisect_guide.md) | ⭐ 10 | `development` |
| [Style Guides](development/2878-style_guides_958e560a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-comment-generator/references/style_guides.md) | ⭐ 10 | `development` |
| [Refactoring Patterns](development/2879-refactoring-patterns_3ffdddfd/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-smell-detector/references/refactoring-patterns.md) | ⭐ 10 | `development` |
| [Cve Analysis](development/2880-cve_analysis_66dbec59/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cve-reachability-analyzer/references/cve_analysis.md) | ⭐ 10 | `development` |
| [Language Guide](development/2881-language_guide_46ebc092/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cve-reachability-analyzer/references/language_guide.md) | ⭐ 10 | `development` |
| [Reachability Patterns](development/2882-reachability_patterns_aaeeebc9/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cve-reachability-analyzer/references/reachability_patterns.md) | ⭐ 10 | `development` |
| [Action Guidelines](development/2883-action_guidelines_db08b980/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cve-watchlist-action-recommendation-generator/references/action_guidelines.md) | ⭐ 10 | `development` |
| [Risk Scoring](development/2884-risk_scoring_b916321e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cve-watchlist-action-recommendation-generator/references/risk_scoring.md) | ⭐ 10 | `development` |
| [Java Deprecations](development/2885-java_deprecations_14ed814d/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/deprecated-api-updater/references/java_deprecations.md) | ⭐ 10 | `development` |
| [Javascript Deprecations](development/2886-javascript_deprecations_72e1bfae/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/deprecated-api-updater/references/javascript_deprecations.md) | ⭐ 10 | `development` |
| [C Cpp Edge Cases](development/2887-c_cpp_edge_cases_ca826cd2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/edge-case-generator/references/c_cpp_edge_cases.md) | ⭐ 10 | `development` |
| [Javascript Edge Cases](development/2888-javascript_edge_cases_9dc16b42/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/edge-case-generator/references/javascript_edge_cases.md) | ⭐ 10 | `development` |
| [Platform Instructions](development/2889-platform_instructions_817b3b55/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/environment-setup-assistant/references/platform_instructions.md) | ⭐ 10 | `development` |
| [Assessment Criteria](development/2890-assessment_criteria_53e1884e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/exploitability-analyzer/references/assessment_criteria.md) | ⭐ 10 | `development` |
| [Framework Comparison](development/2891-framework_comparison_cf5d2ba4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/framework-migration-assistant/references/framework_comparison.md) | ⭐ 10 | `development` |
| [Fuzzing Patterns](development/2892-fuzzing-patterns_882c12c7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/fuzzing-input-generator/references/fuzzing-patterns.md) | ⭐ 10 | `development` |
| [Extraction Patterns](development/2893-extraction_patterns_2f3afd0a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/imperative-to-coq-model-extractor/references/extraction_patterns.md) | ⭐ 10 | `development` |
| [Implementation Patterns](development/2076-implementation-patterns_290fdc50/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/incremental-python-programmer/references/implementation-patterns.md) | ⭐ 10 | `development` |
| [Api Reference](development/051-api_reference_ea6d127a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interface-specification-generator/references/api_reference.md) | ⭐ 10 | `development` |
| [Abstract Interpretation](development/2872-abstract_interpretation_fb585db4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interval-difference-analyzer/references/abstract_interpretation.md) | ⭐ 10 | `development` |
| [Optimization Patterns](development/2894-optimization-patterns_9f396a8e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/interval-profiling-performance-analyzer/references/optimization-patterns.md) | ⭐ 10 | `development` |
| [Report Patterns](development/2895-report_patterns_89b708c1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/issue-report-generator/references/report_patterns.md) | ⭐ 10 | `development` |
| [Dependency Analysis](development/2896-dependency_analysis_290ff993/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/legacy-code-summarizer/references/dependency_analysis.md) | ⭐ 10 | `development` |
| [Temporal Logic Patterns](development/2897-temporal_logic_patterns_7c96bb0c/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/model-guided-code-repair/references/temporal_logic_patterns.md) | ⭐ 10 | `development` |
| [Documentation Patterns](development/2898-documentation-patterns_2d211ed6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/module-component-generator/references/documentation-patterns.md) | ⭐ 10 | `development` |
| [Java Patterns](development/2899-java-patterns_605d10e6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/module-component-generator/references/java-patterns.md) | ⭐ 10 | `development` |
| [Language Patterns](development/2900-language_patterns_b94de45b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/program-to-tlaplus-spec-generator/references/language_patterns.md) | ⭐ 10 | `development` |
| [Coq Pcc](development/2901-coq_pcc_a8152856/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-carrying-code-generator/references/coq_pcc.md) | ⭐ 10 | `development` |
| [Safety Properties](development/2902-safety_properties_5b9928f0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-carrying-code-generator/references/safety_properties.md) | ⭐ 10 | `development` |
| [Refactoring Patterns](development/2903-refactoring_patterns_745b5a7e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/proof-refactoring-assistant/references/refactoring_patterns.md) | ⭐ 10 | `development` |
| [Pseudocode Patterns](development/2904-pseudocode-patterns_d74fefe1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/pseudocode-to-python-code/references/pseudocode-patterns.md) | ⭐ 10 | `development` |
| [Coq Refinement](development/2905-coq_refinement_f2931f85/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/refinement-step-generator/references/coq_refinement.md) | ⭐ 10 | `development` |
| [Isabelle Refinement](development/2906-isabelle_refinement_3f8663a1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/refinement-step-generator/references/isabelle_refinement.md) | ⭐ 10 | `development` |
| [Refinement Patterns](development/2907-refinement_patterns_ce1ff05f/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/refinement-step-generator/references/refinement_patterns.md) | ⭐ 10 | `development` |
| [Detection Strategies](development/2908-detection_strategies_a741df25/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/regression-consistency-checker/references/detection_strategies.md) | ⭐ 10 | `development` |
| [Instrumentation Techniques](development/2909-instrumentation_techniques_9337f605/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/reproduction-trace-instrumenter/references/instrumentation_techniques.md) | ⭐ 10 | `development` |
| [Analysis Patterns](development/2910-analysis-patterns_4825f0f8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-comparison-reporter/references/analysis-patterns.md) | ⭐ 10 | `development` |
| [Analysis Framework](development/2911-analysis_framework_521328e8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/requirement-enhancer/references/analysis_framework.md) | ⭐ 10 | `development` |
| [Debugging Guide](development/1383-debugging_guide_823bbbef/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/runtime-error-explainer/references/debugging_guide.md) | ⭐ 10 | `development` |
| [Remediation Strategies](development/2912-remediation_strategies_6c76acc2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/security-patch-advisor/references/remediation_strategies.md) | ⭐ 10 | `development` |
| [Language Patterns](development/2900-language_patterns_d5190632/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/security-sensitive-path-instrumenter/references/language_patterns.md) | ⭐ 10 | `development` |
| [Formal Verification](development/2913-formal_verification_752d3796/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-equivalence-verifier/references/formal_verification.md) | ⭐ 10 | `development` |
| [Symbolic Execution](development/2914-symbolic_execution_64987aba/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-equivalence-verifier/references/symbolic_execution.md) | ⭐ 10 | `development` |
| [Language Support](development/2915-language_support_c5de8ba0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-szz-analyzer/references/language_support.md) | ⭐ 10 | `development` |
| [Szz Algorithm](development/2916-szz_algorithm_7389292c/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/semantic-szz-analyzer/references/szz_algorithm.md) | ⭐ 10 | `development` |
| [Extraction Patterns](development/2893-extraction_patterns_3fa08400/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/smv-model-extractor/references/extraction_patterns.md) | ⭐ 10 | `development` |
| [Instrumentation Guide](development/2917-instrumentation_guide_7d5fc1cb/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/state-snapshot-instrumenter/references/instrumentation_guide.md) | ⭐ 10 | `development` |
| [Snapshot Format](development/2918-snapshot_format_adcc0b5a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/state-snapshot-instrumenter/references/snapshot_format.md) | ⭐ 10 | `development` |
| [Use Cases](development/1779-use_cases_a8936535/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/state-snapshot-instrumenter/references/use_cases.md) | ⭐ 10 | `development` |
| [Constraint Solving](development/2919-constraint_solving_c7029da5/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/symbolic-execution-assistant/references/constraint_solving.md) | ⭐ 10 | `development` |
| [Path Exploration](development/2920-path_exploration_34b7ff48/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/symbolic-execution-assistant/references/path_exploration.md) | ⭐ 10 | `development` |
| [Tool Integration](development/2921-tool_integration_72b88b4e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/symbolic-execution-assistant/references/tool_integration.md) | ⭐ 10 | `development` |
| [Szz Algorithm](development/2916-szz_algorithm_106b9607/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/szz-bug-identifier/references/szz_algorithm.md) | ⭐ 10 | `development` |
| [Api Reference](development/051-api_reference_9d07b088/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/taint-instrumentation-assistant/references/api_reference.md) | ⭐ 10 | `development` |
| [Repair Patterns](development/2922-repair_patterns_14367a84/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tlaplus-guided-code-repair/references/repair_patterns.md) | ⭐ 10 | `development` |
| [Code Scanning](development/2923-code_scanning_a8773842/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/traceability-matrix-generator/references/code_scanning.md) | ⭐ 10 | `development` |
| [Boundary Patterns](development/2924-boundary_patterns_0cf177d7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verification-boundary-reporter/references/boundary_patterns.md) | ⭐ 10 | `development` |
| [Coq Analysis](development/2925-coq_analysis_abfbdbfd/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verification-boundary-reporter/references/coq_analysis.md) | ⭐ 10 | `development` |
| [Dafny Analysis](development/2926-dafny_analysis_cdbf3d93/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verification-boundary-reporter/references/dafny_analysis.md) | ⭐ 10 | `development` |
| [Isabelle Analysis](development/2927-isabelle_analysis_b7eb79b4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verification-boundary-reporter/references/isabelle_analysis.md) | ⭐ 10 | `development` |
| [Mapping Patterns](development/2928-mapping_patterns_ad93136e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/verified-spec-code-mapper/references/mapping_patterns.md) | ⭐ 10 | `development` |
| [Api Reference](development/051-api_reference_988443f4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/vulnerability-pattern-matcher/references/api_reference.md) | ⭐ 10 | `development` |
| [Vulnerability Patterns](development/2929-vulnerability_patterns_2d771fb2/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/vulnerability-pattern-matcher/references/vulnerability_patterns.md) | ⭐ 10 | `development` |
| [Analysis Strategies](development/2930-analysis_strategies_b9283726/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/vulnerability-root-cause-analyzer/references/analysis_strategies.md) | ⭐ 10 | `development` |
| [Root Causes](development/2931-root_causes_720b1819/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/vulnerability-root-cause-analyzer/references/root_causes.md) | ⭐ 10 | `development` |
| [Vulnerability Patterns](development/2929-vulnerability_patterns_fb8a1472/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/vulnerability-root-cause-analyzer/references/vulnerability_patterns.md) | ⭐ 10 | `development` |
| [Catalog](development/126-catalog_15bc0fce/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/CATALOG.md) | 🔥 13.4k | `development` |
| [Skill](development/1178-name-skill_96c5cffb/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/code-reviewer/SKILL.md) | 🔥 13.4k | `development` |
| [Readme Flat Skills Releases](development/799-readme_flat_skills_releases_9c4e255d/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_SKILLS_RELEASES.md) | 🔥 24.4k | `development` |
| [Streaming Sse Guide](development/2873-streaming-sse-guide_56799ef2/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/docs/streaming-sse-guide.md) | ⭐ 600 | `development` |
| [Optimization Summary](development/2874-optimization_summary_829204c9/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/skills/frontend-design/OPTIMIZATION_SUMMARY.md) | ⭐ 600 | `development` |
| [Skill](development/1178-name-skill_b2af3b48/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/skills/web-search/SKILL.md) | ⭐ 600 | `development` |
| [Skill](development/1178-name-skill_53f1cebe/) | [RafaelGorski/Problem-Based-SRS](https://raw.githubusercontent.com/RafaelGorski/Problem-Based-SRS/main/skills/problem-based-srs/SKILL.md) | ⭐ 10 | `development` |
| [Cooldown](development/2920-cooldown_4267ed35/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/docs/COOLDOWN.md) | 🔥 24.6k | `development` |
| [Ai Agents Skills For Embedded Engineers Cursor](development/3078-ai-agents-skills-for-embedded-engineers-cursor_1284f72d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/workshops/ai-agents-skills-for-embedded-engineers-cursor.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_692faca7/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.agents/skills/agent-browser/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/1178-name-skill_993fd8f5/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/create-merge-request-changelog/SKILL.md) | ⭐ 18 | `development` |
| [Proxy Support](development/3079-proxy-support_4c44ba60/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.agents/skills/agent-browser/references/proxy-support.md) | ⭐ 18 | `development` |
| [Skill](development/1530-description-skill_4c9f2d09/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/workshops/.cursor/skills/c-embedded-standards/SKILL.md) | ⭐ 18 | `development` |
| [Claude](development/140-claude_867c1488/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/CLAUDE.md) | ⭐ 143 | `development` |
| [Action Stations](development/3080-action-stations_d64772d5/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/action-stations.md) | ⭐ 143 | `development` |
| [Crew Roles](development/3081-crew-roles_ff61add7/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/crew-roles.md) | ⭐ 143 | `development` |
| [Royal Marines](development/3082-royal-marines_92007cea/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/references/royal-marines.md) | ⭐ 143 | `development` |
| [Context Guide](development/2884-context_guide_789638a4/) | [doobidoo/MCP-Context-Provider](https://raw.githubusercontent.com/doobidoo/MCP-Context-Provider/main/docs/guides/CONTEXT_GUIDE.md) | ⭐ 24 | `development` |
| [Workflow Troubleshooting](development/2875-workflow_troubleshooting_d48aaa23/) | [doobidoo/MCP-Context-Provider](https://raw.githubusercontent.com/doobidoo/MCP-Context-Provider/main/docs/WORKFLOW_TROUBLESHOOTING.md) | ⭐ 24 | `development` |
| [Configuration](development/191-configuration_41aedf4b/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/CONFIGURATION.md) | ⭐ 4.0k | `development` |
| [Session Summary](development/2877-session_summary_a4a22d16/) | [doobidoo/MCP-Context-Provider](https://raw.githubusercontent.com/doobidoo/MCP-Context-Provider/main/docs/phase3/SESSION_SUMMARY.md) | ⭐ 24 | `development` |
| [Skill](development/1178-name-skill_8f7ac62b/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/seed_skills/skill-finder/SKILL.md) | ⭐ 915 | `development` |
| [Skill](development/1178-name-skill_a93c5015/) | [agentic-community/mcp-gateway-registry](https://raw.githubusercontent.com/agentic-community/mcp-gateway-registry/main/.claude/skills/new-feature-design/SKILL.md) | ⭐ 449 | `development` |
| [Skill](development/1178-name-skill_1c3e32a8/) | [agentic-community/mcp-gateway-registry](https://raw.githubusercontent.com/agentic-community/mcp-gateway-registry/main/.claude/skills/pr-review/SKILL.md) | ⭐ 449 | `development` |
| [Copilotkit](development/copilotkit_9ebf6b0d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/agent-frameworks/copilotkit.md) | ⭐ 20 | `development` |
| [Cline](development/cline_883e86e1/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/coding-agents/cline.md) | ⭐ 20 | `development` |
| [Google Ai Studio](development/google-ai-studio_293ed116/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/developer-tools/google-ai-studio.md) | ⭐ 20 | `development` |
| [Google Ai Studio](development/google-ai-studio_fc156345/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/prompt-engineering/google-ai-studio.md) | ⭐ 20 | `development` |

### Development/Devops (172 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agent Deployment Guide](development/devops/229-agent-deployment-guide_fc15347e/) | [Jaganpro/sf-skills](https://raw.githubusercontent.com/Jaganpro/sf-skills/main/skills/sf-deploy/references/agent-deployment-guide.md) | ⭐ 80 | `development` |
| [Guide Skills Frontend](development/devops/376-guide-skills-frontend_2f63a088/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/docs/guide-skills-frontend.md) | ⭐ 102 | `development` |
| [Guide Skills Monorepo](development/devops/377-guide-skills-monorepo_83c4b53f/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/docs/guide-skills-monorepo.md) | ⭐ 102 | `development` |
| [Skill](development/devops/014-name-skill_1ec7fc18/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `java` |
| [Skill](development/devops/014-name-skill_9bc6e00e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `php` |
| [Skill](development/devops/014-name-skill_3c0620fc/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `python` |
| [Skill](development/devops/014-name-skill_2658dda5/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/SKILL.md) | ⭐ 102 | `aws` `lambda` `typescript` |
| [Skill](development/devops/014-name-skill_62d193fa/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/SKILL.md) | ⭐ 102 | `nextjs` `next.js` `deployment` |
| [Skill](development/devops/014-name-skill_d0e8129e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/SKILL.md) | ⭐ 102 | `nx` `monorepo` `typescript` |
| [Skill](development/devops/014-name-skill_4133613a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/SKILL.md) | ⭐ 102 | `development` |
| [Raw Java Lambda](development/devops/378-raw-java-lambda_bc77ac25/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/raw-java-lambda.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_7a94e770/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-java/skills/aws-lambda-java-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_0963c52b/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-php/skills/aws-lambda-php-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Chalice Lambda](development/devops/380-chalice-lambda_73723400/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/chalice-lambda.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_64158df7/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-python/skills/aws-lambda-python-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Serverless Deployment](development/devops/379-serverless-deployment_5ae8bdfe/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/aws-lambda-typescript-integration/references/serverless-deployment.md) | ⭐ 102 | `development` |
| [Docker Patterns](development/devops/381-docker-patterns_6cc95fae/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nextjs-deployment/references/docker-patterns.md) | ⭐ 102 | `development` |
| [Basics](development/devops/382-basics_94f6201e/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/basics.md) | ⭐ 102 | `development` |
| [Mcp Server Python](development/devops/383-mcp-server-python_5f654215/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/mcp-server-python.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_ccaab12c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/index.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_3c572c12/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/index.md) | ⭐ 3.3k | `development` |
| [Roadmap](development/devops/097-roadmap_c9aa8183/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/roadmap.md) | ⭐ 3.3k | `development` |
| [Security Features](development/devops/384-security-features_012c8460/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/security-features.md) | ⭐ 3.3k | `development` |
| [Input Validation](development/devops/228-input-validation_01d87216/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/input-validation.md) | ⭐ 3.3k | `development` |
| [Mcp Architecture Patterns](development/devops/385-mcp-architecture-patterns_230d8859/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/best-practices/mcp-architecture-patterns.md) | ⭐ 3.3k | `development` |
| [Argocd](development/devops/386-argocd_31a49405/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/argocd.md) | ⭐ 3.3k | `development` |
| [Cforge Gateway](development/devops/387-cforge-gateway_ac5a8cdb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/cforge-gateway.md) | ⭐ 3.3k | `security` `filter` |
| [Compose](development/devops/388-compose_366cfc28/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/compose.md) | ⭐ 3.3k | `development` |
| [Fly Io](development/devops/389-fly-io_e5ccc7e2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/fly-io.md) | ⭐ 3.3k | `development` |
| [Google Cloud Run](development/devops/390-google-cloud-run_8be8613f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/google-cloud-run.md) | ⭐ 3.3k | `development` |
| [Helm](development/devops/391-helm_2945bf9d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/helm.md) | ⭐ 3.3k | `development` |
| [Ibm Code Engine](development/devops/392-ibm-code-engine_fe3cb3dc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/ibm-code-engine.md) | ⭐ 3.3k | `development` |
| [Local](development/devops/393-local_93c0d987/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/local.md) | ⭐ 3.3k | `development` |
| [Minikube](development/devops/394-minikube_0c8a4b03/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/minikube.md) | ⭐ 3.3k | `development` |
| [Proxy Auth](development/devops/395-proxy-auth_7d43d3de/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/proxy-auth.md) | ⭐ 3.3k | `development` |
| [Tls Configuration](development/devops/201-tls-configuration_54ff4560/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/deployment/tls-configuration.md) | ⭐ 3.3k | `development` |
| [Developer Onboarding](development/devops/396-developer-onboarding_61358b2f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/developer-onboarding.md) | ⭐ 3.3k | `development` |
| [Github](development/devops/397-github_487e80b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/github.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_a20a1f7e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/index.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_0010150a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/faq/index.md) | ⭐ 3.3k | `development` |
| [Admin Ui Customization](development/devops/398-admin-ui-customization_b4697c47/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/admin-ui-customization.md) | ⭐ 3.3k | `development` |
| [Api Usage](development/devops/399-api-usage_930cc30b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/api-usage.md) | ⭐ 3.3k | `development` |
| [Configuration](development/devops/009-configuration_ed70f974/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/configuration.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_7c9efd7f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/index.md) | ⭐ 3.3k | `development` |
| [Logging](development/devops/400-logging_7601ca0a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/logging.md) | ⭐ 3.3k | `development` |
| [Proxy](development/devops/401-proxy_5a83e05a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/proxy.md) | ⭐ 3.3k | `development` |
| [Rbac](development/devops/010-rbac_46494f27/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/rbac.md) | ⭐ 3.3k | `development` |
| [Scale](development/devops/402-scale_5375b413/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/scale.md) | ⭐ 3.3k | `development` |
| [Securing](development/devops/403-securing_264f96d4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/securing.md) | ⭐ 3.3k | `development` |
| [Self Signed Certificates](development/devops/404-self-signed-certificates_c29006e4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/self-signed-certificates.md) | ⭐ 3.3k | `development` |
| [Sso Generic Oidc Tutorial](development/devops/405-sso-generic-oidc-tutorial_1cf7119b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-generic-oidc-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Google Tutorial](development/devops/406-sso-google-tutorial_ec0b87b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-google-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso Keycloak Tutorial](development/devops/407-sso-keycloak-tutorial_bc429daa/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso-keycloak-tutorial.md) | ⭐ 3.3k | `development` |
| [Sso](development/devops/408-sso_3ff6196e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso.md) | ⭐ 3.3k | `development` |
| [Tuning](development/devops/409-tuning_b13077bc/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/tuning.md) | ⭐ 3.3k | `development` |
| [Upgrade](development/devops/410-upgrade_b3150af6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/upgrade.md) | ⭐ 3.3k | `development` |
| [Well Known Uris](development/devops/411-well-known-uris_b7a6cdf4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/well-known-uris.md) | ⭐ 3.3k | `development` |
| [Config Validation](development/devops/412-config-validation_dc2663b4/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/operations/config-validation.md) | ⭐ 3.3k | `development` |
| [Cpu Spin Loop Mitigation](development/devops/413-cpu-spin-loop-mitigation_fdb3b575/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/operations/cpu-spin-loop-mitigation.md) | ⭐ 3.3k | `development` |
| [Features](development/devops/361-features_4c3fbb11/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/features.md) | ⭐ 3.3k | `development` |
| [Admin UI](development/devops/414-ui_223c1cd2/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/ui.md) | ⭐ 3.3k | `development` |
| [Argocd Helm Deployment Ibm Cloud Iks](development/devops/415-argocd-helm-deployment-ibm-cloud-iks_8a45b06f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/argocd-helm-deployment-ibm-cloud-iks.md) | ⭐ 3.3k | `development` |
| [Openwebui Tutorial](development/devops/416-openwebui-tutorial_ace5026d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/tutorials/openwebui-tutorial.md) | ⭐ 3.3k | `development` |
| [Mcpgateway Translate](development/devops/417-mcpgateway-translate_199775f5/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/mcpgateway-translate.md) | ⭐ 3.3k | `development` |
| [Multi Auth Headers](development/devops/418-multi-auth-headers_d4462e80/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/multi-auth-headers.md) | ⭐ 3.3k | `development` |
| [Reverse Proxy](development/devops/419-reverse-proxy_a49bea67/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/reverse-proxy.md) | ⭐ 3.3k | `development` |
| [014 Security Headers Cors Middleware](development/devops/420-014-security-headers-cors-middleware_009499c3/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/014-security-headers-cors-middleware.md) | ⭐ 3.3k | `development` |
| [015 Well Known Uri Handler](development/devops/421-015-well-known-uri-handler_198ce9fb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/015-well-known-uri-handler.md) | ⭐ 3.3k | `development` |
| [018 Built In Response Compression](development/devops/422-018-built-in-response-compression_161e1ccb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/018-built-in-response-compression.md) | ⭐ 3.3k | `development` |
| [019 Modular Architecture Split](development/devops/423-019-modular-architecture-split_606af208/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/019-modular-architecture-split.md) | ⭐ 3.3k | `development` |
| [025 Granian Http Server](development/devops/424-025-granian-http-server_7849918e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/025-granian-http-server.md) | ⭐ 3.3k | `development` |
| [032 Mcp Session Pool](development/devops/425-032-mcp-session-pool_fd39527c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/032-mcp-session-pool.md) | ⭐ 3.3k | `development` |
| [036 Bootstrap Custom Roles](development/devops/426-036-bootstrap-custom-roles_98fdf334/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/036-bootstrap-custom-roles.md) | ⭐ 3.3k | `development` |
| [038 Multi Worker Session Affinity](development/devops/193-038-multi-worker-session-affinity_de1547c3/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/038-multi-worker-session-affinity.md) | ⭐ 3.3k | `development` |
| [Observability](development/devops/427-observability_358deba7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/observability.md) | ⭐ 3.3k | `development` |
| [Phoenix](development/devops/428-phoenix_a6b52052/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/observability/phoenix.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_49786343/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/media/social/index.md) | ⭐ 3.3k | `development` |
| [Bee](development/devops/429-bee_bd521c8a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/bee.md) | ⭐ 3.3k | `development` |
| [Index](development/devops/050-index_6da5357e/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/index.md) | ⭐ 3.3k | `development` |
| [Claude Desktop](development/devops/430-claude-desktop_33f0496c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/claude-desktop.md) | ⭐ 3.3k | `development` |
| [Continue](development/devops/431-continue_4f371801/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/clients/continue.md) | ⭐ 3.3k | `development` |
| [Grpc Transport](development/devops/226-grpc-transport_c50d7487/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/grpc-transport.md) | ⭐ 3.3k | `development` |
| [Mtls](development/devops/432-mtls_54eac01a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/mtls.md) | ⭐ 3.3k | `development` |
| [Terraform](development/devops/433-terraform_86d5b2a9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/hashicorp/terraform.md) | ⭐ 3.3k | `development` |
| [Langflow Server](development/devops/434-langflow-server_4739c90c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/third-party/langflow-server.md) | ⭐ 3.3k | `development` |
| [Github](development/devops/397-github_cda2891c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/microsoft/github.md) | ⭐ 3.3k | `development` |
| [Rollback](development/devops/364-rollback_4f482aad/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/troubleshoot/rollback.md) | 🔥 36.4k | `development` |
| [Cli Commands](development/devops/319-cli-commands_4a42a248/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/cli-commands.md) | ⭐ 61 | `development` |
| [Multi Tenancy](development/devops/361-multi-tenancy_677a145e/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/multi-tenancy.md) | ⭐ 61 | ``tenant:${tenantId}`` `role:support` |
| [Skill](development/devops/014-name-skill_7c48cfe9/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-orchestrator/SKILL.md) | ⭐ 10 | `development` |
| [Readme Es](development/devops/435-readme_es_1474cf17/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_es.md) | ⭐ 820 | `development` |
| [Readme Ja](development/devops/436-readme_ja_a026d9c8/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_ja.md) | ⭐ 820 | `development` |
| [Readme Pt Br](development/devops/437-readme_pt-br_698660d2/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_pt-BR.md) | ⭐ 820 | `development` |
| [Readme Zh Cn](development/devops/438-readme_zh-cn_b5fce0bf/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_zh-CN.md) | ⭐ 820 | `development` |
| [Mcp](development/devops/025-mcp_b3a0b417/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/docs/docs/concepts/mcp.md) | ⭐ 820 | `development` |
| [Securing](development/devops/361-securing_0ab63085/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/securing.md) | ⭐ 3.3k | `development` |
| [Agents](development/devops/053-agents_05228376/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/AGENTS.md) | 🔥 9.7k | `development` |
| [Environment Variables](development/devops/362-environment_variables_588fc08a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/ENVIRONMENT_VARIABLES.md) | 🔥 9.7k | `development` |
| [05 Workflows](development/devops/363-05-workflows_75f7af4a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/05-workflows.md) | 🔥 9.7k | `development` |
| [Quick Reference](development/devops/364-quick_reference_07ac890f/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/archive/legacy/QUICK_REFERENCE.md) | 🔥 9.7k | `development` |
| [05 Workflows](development/devops/363-05-workflows_f0c92b2a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/zh-CN/user-guide/05-workflows.md) | 🔥 9.7k | `development` |
| [Readme Es](development/devops/361-readme_es_853948d5/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_es.md) | ⭐ 867 | `development` |
| [Readme Ja](development/devops/362-readme_ja_c5cebf0e/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_ja.md) | ⭐ 867 | `development` |
| [Readme Pt Br](development/devops/363-readme_pt-br_bb8b45ba/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_pt-BR.md) | ⭐ 867 | `development` |
| [Readme Zh Cn](development/devops/364-readme_zh-cn_a711b2d7/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_zh-CN.md) | ⭐ 867 | `development` |
| [Golden Traces](development/devops/365-golden_traces_f3a5dd38/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/GOLDEN_TRACES.md) | ⭐ 43 | `development` |
| [Mcp Contracts](development/devops/366-mcp_contracts_a2cd7e77/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/MCP_CONTRACTS.md) | ⭐ 43 | `development` |
| [Quickstart Langgraph](development/devops/367-quickstart_langgraph_89171e4f/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/QUICKSTART_LANGGRAPH.md) | ⭐ 43 | `development` |
| [V1.81.14](development/devops/370-v18114_726abb57/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/release_notes/v1.81.14.md) | 🔥 36.1k | `development` |
| [Index](development/devops/050-index_89352769/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/blog/server_root_path/index.md) | 🔥 36.1k | `incident-report` `ui` `stability` |
| [Skill](development/devops/014-name-skill_48cceb05/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agent-orchestration/skills/agent-orchestration/SKILL.md) | ⭐ 18 | `development` |
| [Gitlab Ci Local Guide](development/devops/123-gitlab-ci-local-guide_cc9fc740/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/gitlab-skill/skills/gitlab-skill/references/gitlab-ci-local-guide.md) | ⭐ 18 | `development` |
| [Keyword Updates V2.8](development/devops/362-keyword_updates_v28_ec0aa0de/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/KEYWORD_UPDATES_V2.8.md) | ⭐ 10 | `development` |
| [Stream E Implementation Summary](development/devops/363-stream-e-implementation-summary_a69ff379/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/STREAM-E-IMPLEMENTATION-SUMMARY.md) | ⭐ 10 | `development` |
| [02 Learning Roadmap](development/devops/364-02-learning-roadmap_5da15d75/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/01-getting-started/02-learning-roadmap.md) | ⭐ 10 | `development` |
| [01 Configuration](development/devops/314-01-configuration_5d02e38c/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/20-configuration/01-configuration.md) | ⭐ 10 | `development` |
| [01 Working Protocol](development/devops/365-01-working-protocol_c4d03ee5/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/30-operations/01-working-protocol.md) | ⭐ 10 | `development` |
| [Lifecycle Hooks](development/devops/366-lifecycle-hooks_b2dfa81a/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/lifecycle-hooks.md) | ⭐ 10 | `development` |
| [Skill Evaluation](development/devops/367-skill-evaluation_75680f21/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/skill-evaluation.md) | ⭐ 10 | `frontend` `ui` |
| [Ci Cd Patterns](development/devops/368-ci-cd-patterns_80a16cb1/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/devops/ci-cd-patterns.md) | ⭐ 10 | `ci-cd` `github-actions` `pipeline` |
| [Git Workflows](development/devops/369-git-workflows_870cd392/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/skills/devops/git-workflows.md) | ⭐ 10 | `git` `version-control` `branching` |
| [Developing](development/devops/371-developing_551cf265/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/developing.md) | ⭐ 4.0k | `development` |
| [Readme Es](development/devops/362-readme_es_ae0fbdc4/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_es.md) | ⭐ 858 | `development` |
| [Readme Ja](development/devops/363-readme_ja_d27afd34/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_ja.md) | ⭐ 858 | `development` |
| [Readme Zh Cn](development/devops/364-readme_zh-cn_7b293e68/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/README_zh-CN.md) | ⭐ 858 | `development` |
| [Security Features](development/devops/365-security-features_62041f2c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/security-features.md) | ⭐ 3.3k | `development` |
| [Configuration](development/devops/009-configuration_befa1307/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/configuration.md) | ⭐ 3.3k | `development` |
| [Proxy](development/devops/366-proxy_ca57609f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/proxy.md) | ⭐ 3.3k | `development` |
| [Securing](development/devops/367-securing_37c36d58/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/securing.md) | ⭐ 3.3k | `development` |
| [Sso](development/devops/368-sso_078f40c7/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/sso.md) | ⭐ 3.3k | `development` |
| [Reverse Proxy](development/devops/369-reverse-proxy_047a8541/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/reverse-proxy.md) | ⭐ 3.3k | `development` |
| [Readme Extra](development/devops/156-readme_extra_7bd9e4c4/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_EXTRA.md) | 🔥 24.5k | `development` |
| [Readme Flat All Created](development/devops/158-readme_flat_all_created_dd2a8c4f/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_CREATED.md) | 🔥 24.5k | `development` |
| [Readme Flat All Updated](development/devops/160-readme_flat_all_updated_4ea6dcff/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_UPDATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Claude Md Created](development/devops/378-readme_flat_claude-md_created_6adcd3ea/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_CREATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Claude Md Updated](development/devops/379-readme_flat_claude-md_updated_8b4b36b6/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_UPDATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Commands Az](development/devops/380-readme_flat_commands_az_4e6b10d3/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_AZ.md) | 🔥 24.5k | `development` |
| [Readme Flat Commands Updated](development/devops/381-readme_flat_commands_updated_051d5cf1/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_UPDATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Tooling Created](development/devops/162-readme_flat_tooling_created_66311c52/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_CREATED.md) | 🔥 24.5k | `development` |
| [Skill](development/devops/014-name-skill_dc6361a6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/change-log-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/devops/014-name-skill_9a58dcd1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/configuration-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/devops/014-name-skill_eadf1bfd/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/containerization-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/devops/014-name-skill_3d634c6a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/cve-watchlist-action-recommendation-generator/SKILL.md) | ⭐ 10 | `development` |
| [Ci Platforms](development/devops/362-ci-platforms_999207f0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/build-ci-migration-assistant/references/ci-platforms.md) | ⭐ 10 | `linux` |
| [Project Migration](development/devops/363-project_migration_37a98145/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-translation/references/project_migration.md) | ⭐ 10 | `development` |
| [Infra Configs](development/devops/364-infra_configs_015d228e/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/configuration-generator/references/infra_configs.md) | ⭐ 10 | `development` |
| [Dockerfile Best Practices](development/devops/365-dockerfile_best_practices_fbac54ba/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/containerization-assistant/references/dockerfile_best_practices.md) | ⭐ 10 | `development` |
| [Tool Guides](development/devops/366-tool_guides_4609e245/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/environment-setup-assistant/references/tool_guides.md) | ⭐ 10 | `development` |
| [Database Rollback Patterns](development/devops/367-database_rollback_patterns_1204166d/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/rollback-strategy-advisor/references/database_rollback_patterns.md) | ⭐ 10 | `development` |
| [Framework Comparison](development/devops/368-framework_comparison_5ee1caae/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/spring-mvc-to-boot-migrator/references/framework_comparison.md) | ⭐ 10 | `development` |
| [Migration Guide](development/devops/011-migration_guide_6ad9289a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/spring-mvc-to-boot-migrator/references/migration_guide.md) | ⭐ 10 | `development` |
| [Graphviz Patterns](development/devops/369-graphviz-patterns_ddba74c0/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/system-diagram-generator/references/graphviz-patterns.md) | ⭐ 10 | `development` |
| [Skill](development/devops/014-name-skill_a1cc67a6/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/architect-review/SKILL.md) | 🔥 13.4k | `development` |
| [Skill](development/devops/014-name-skill_f9416fde/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/performance-engineer/SKILL.md) | 🔥 13.4k | `development` |
| [Readme Flat Tooling Releases](development/devops/163-readme_flat_tooling_releases_7e2219ec/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_RELEASES.md) | 🔥 24.4k | `development` |
| [Readme Awesome](development/devops/154-readme_awesome_8f88676e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_AWESOME.md) | 🔥 24.6k | `development` |
| [Readme Classic](development/devops/155-readme_classic_58de1a77/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_CLASSIC.md) | 🔥 24.6k | `development` |
| [Readme Flat All Az](development/devops/157-readme_flat_all_az_46e4168e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_AZ.md) | 🔥 24.6k | `development` |
| [Readme Flat All Releases](development/devops/159-readme_flat_all_releases_1b29b908/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_ALL_RELEASES.md) | 🔥 24.6k | `development` |
| [Readme Flat Claude Md Az](development/devops/370-readme_flat_claude-md_az_12d1ad30/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_AZ.md) | 🔥 24.6k | `development` |
| [Readme Flat Claude Md Releases](development/devops/371-readme_flat_claude-md_releases_b8e41ed9/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLAUDE-MD_RELEASES.md) | 🔥 24.6k | `development` |
| [Readme Flat Tooling Az](development/devops/161-readme_flat_tooling_az_38ad2f96/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_TOOLING_AZ.md) | 🔥 24.6k | `development` |
| [Skill](development/devops/085-description-skill_dc4dba10/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/workshops/.cursor/skills/ha-zigbee2mqtt-docker/SKILL.md) | ⭐ 18 | `development` |
| [Developing](development/devops/365-developing_f9dbbab5/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/developing.md) | ⭐ 4.0k | `development` |
| [Smr Sdk Mcp Control Proposal 2026 02 19](development/devops/372-smr-sdk-mcp-control-proposal-2026-02-19_d3b17850/) | [synth-laboratories/synth-ai](https://raw.githubusercontent.com/synth-laboratories/synth-ai/main/contracts/smr-sdk-mcp-control-proposal-2026-02-19.md) | ⭐ 75 | `development` |
| [Rbac](development/devops/010-rbac_941810a9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/rbac.md) | ⭐ 3.3k | `development` |
| [Securing](development/devops/366-securing_f01095ce/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/securing.md) | ⭐ 3.3k | `development` |
| [Reverse Proxy](development/devops/378-reverse-proxy_c536397f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/reverse-proxy.md) | ⭐ 3.3k | `development` |
| [Entra Id Setup](development/devops/240-entra-id-setup_33e50d9e/) | [agentic-community/mcp-gateway-registry](https://raw.githubusercontent.com/agentic-community/mcp-gateway-registry/main/docs/entra-id-setup.md) | ⭐ 449 | `development` |
| [Security Scanner](development/devops/004-security-scanner_23f0cff7/) | [agentic-community/mcp-gateway-registry](https://raw.githubusercontent.com/agentic-community/mcp-gateway-registry/main/docs/security-scanner.md) | ⭐ 449 | `development` |
| [Fly Io](development/devops/fly-io_8538885d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/agent-infrastructure/fly-io.md) | ⭐ 20 | `development` |
| [Motia](development/devops/motia_5fa96f76/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/research/api-frameworks/motia.md) | ⭐ 20 | `development` |

### Development/Testing (25 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](development/testing/002-name-skill_187b30c5/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nestjs-drizzle-crud-generator/SKILL.md) | ⭐ 102 | `development` |
| [Nestjs Config](development/testing/089-nestjs-config_ab544f6a/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/turborepo-monorepo/references/nestjs-config.md) | ⭐ 102 | `development` |
| [Db Performance](development/testing/090-db-performance_9aa0f078/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/development/db-performance.md) | ⭐ 3.3k | `development` |
| [Backends](development/testing/091-backends_3471597a/) | [vstorm-co/pydantic-ai-backend](https://raw.githubusercontent.com/vstorm-co/pydantic-ai-backend/main/docs/concepts/backends.md) | ⭐ 36 | `development` |
| [Architecture](development/testing/085-architecture_40b6849a/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/docs/ARCHITECTURE.md) | ⭐ 151 | `development` |
| [Create Pr](development/testing/086-create-pr_a6287a1c/) | [alirezarezvani/ClaudeForge](https://raw.githubusercontent.com/alirezarezvani/ClaudeForge/dev/.claude/commands/github/create-pr.md) | ⭐ 151 | `development` |
| [Patterns](development/testing/059-patterns_aed86128/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/patterns.md) | ⭐ 106 | `development` |
| [Programmatic Skills](development/testing/085-programmatic-skills_3ca2713a/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/programmatic-skills.md) | ⭐ 106 | `development` |
| [Skill](development/testing/002-name-skill_0eb6112a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/agent-browser/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/testing/002-name-skill_21ba0401/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/coverage-enhancer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/testing/002-name-skill_ada59ba6/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/deprecated-api-updater/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/testing/002-name-skill_fa978da1/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/error-explanation-generator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/testing/002-name-skill_94e280c7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/framework-migration-assistant/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/testing/002-name-skill_0f0dd6fd/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/incremental-python-programmer/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/testing/002-name-skill_76eee166/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/init-deep/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/testing/002-name-skill_52e52268/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/regression-root-cause-analyzer/SKILL.md) | ⭐ 10 | `development` |
| [Pattern Catalog](development/testing/082-pattern_catalog_d1a14b18/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/design-pattern-suggestor/references/pattern_catalog.md) | ⭐ 10 | `development` |
| [Python Edge Cases](development/testing/083-python_edge_cases_c20c67e3/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/edge-case-generator/references/python_edge_cases.md) | ⭐ 10 | `development` |
| [Debugging Strategies](development/testing/084-debugging_strategies_d90ef8c8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/error-explanation-generator/references/debugging_strategies.md) | ⭐ 10 | `development` |
| [Error Patterns](development/testing/085-error_patterns_f4a78eb7/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/error-explanation-generator/references/error_patterns.md) | ⭐ 10 | `development` |
| [Migration Guide](development/testing/039-migration_guide_395a5136/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/framework-migration-assistant/references/migration_guide.md) | ⭐ 10 | `development` |
| [Failure Patterns](development/testing/086-failure-patterns_d331f78a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/regression-root-cause-analyzer/references/failure-patterns.md) | ⭐ 10 | `development` |
| [Python Replay](development/testing/087-python-replay_5cd4a143/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/replay-oriented-instrumentation/references/python-replay.md) | ⭐ 10 | `development` |
| [Dependency Formats](development/testing/088-dependency_formats_fb90d6f4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/time-aware-dependency-cve-scanner/references/dependency_formats.md) | ⭐ 10 | `development` |
| [Commands](development/testing/092-commands_cc6cc178/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.agents/skills/agent-browser/references/commands.md) | ⭐ 18 | `development` |

### Development/Tools (117 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Typescript](development/tools/349-typescript_d85d7b93/) | [giuseppe-trisciuoglio/developer-kit](https://raw.githubusercontent.com/giuseppe-trisciuoglio/developer-kit/main/plugins/developer-kit-typescript/skills/nx-monorepo/references/typescript.md) | ⭐ 102 | `development` |
| [Skill](development/tools/002-name-skill_034a479c/) | [netresearch/jira-skill](https://raw.githubusercontent.com/netresearch/jira-skill/main/skills/jira-communication/SKILL.md) | ⭐ 24 | `development` |
| [Developing](development/tools/350-developing_4eef14b0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/DEVELOPING.md) | ⭐ 3.3k | `development` |
| [Enable Payload Logging](development/tools/351-enable_payload_logging_c75da92d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/enable_payload_logging.md) | ⭐ 3.3k | `development` |
| [Mcpgateway](development/tools/010-mcpgateway_f02247c0/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/mcpgateway.md) | ⭐ 3.3k | `development` |
| [Plugins Llms](development/tools/352-plugins-llms_8fca01fe/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/llms/plugins-llms.md) | ⭐ 3.3k | `development` |
| [Oauth Authorization Code Ui Design](development/tools/011-oauth-authorization-code-ui-design_ec055b8d/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/oauth-authorization-code-ui-design.md) | ⭐ 3.3k | `development` |
| [Oauth Design](development/tools/012-oauth-design_dd0d2a41/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/oauth-design.md) | ⭐ 3.3k | `development` |
| [Bulk Import](development/tools/353-bulk-import_8f7144ba/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/bulk-import.md) | ⭐ 3.3k | `development` |
| [Export Import Reference](development/tools/354-export-import-reference_cd3b304c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/export-import-reference.md) | ⭐ 3.3k | `development` |
| [Metadata Tracking](development/tools/355-metadata-tracking_93c516b6/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/metadata-tracking.md) | ⭐ 3.3k | `development` |
| [Oauth](development/tools/356-oauth_2ff3a5f9/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/oauth.md) | ⭐ 3.3k | `development` |
| [Tags](development/tools/357-tags_0898f285/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/overview/tags.md) | ⭐ 3.3k | `development` |
| [006 Gateway Tool Rate Limiting](development/tools/358-006-gateway-tool-rate-limiting_bb786b5a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/006-gateway-tool-rate-limiting.md) | ⭐ 3.3k | `development` |
| [021 Built In Proxy Vs Service Mesh](development/tools/359-021-built-in-proxy-vs-service-mesh_87238005/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/021-built-in-proxy-vs-service-mesh.md) | ⭐ 3.3k | `development` |
| [023 One Time Authentication Servers](development/tools/360-023-one-time-authentication-servers_6f4caffb/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/023-one-time-authentication-servers.md) | ⭐ 3.3k | `development` |
| [027 Migrate Psycopg3](development/tools/361-027-migrate-psycopg3_ae99362f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/027-migrate-psycopg3.md) | ⭐ 3.3k | `development` |
| [035 Query Parameter Authentication](development/tools/362-035-query-parameter-authentication_7da5eb7c/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/035-query-parameter-authentication.md) | ⭐ 3.3k | `development` |
| [Llamaindex](development/tools/363-llamaindex_332803db/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/agents/llamaindex.md) | ⭐ 3.3k | `development` |
| [Http Auth Hooks](development/tools/364-http-auth-hooks_18dede4f/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/plugins/http-auth-hooks.md) | ⭐ 3.3k | `development` |
| [Fast Time Server](development/tools/365-fast-time-server_c35f9972/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/fast-time-server.md) | ⭐ 3.3k | `development` |
| [Instana](development/tools/366-instana_c81a89af/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/ibm/instana.md) | ⭐ 3.3k | `development` |
| [Mermaid Server](development/tools/367-mermaid-server_f0be2d55/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/mermaid-server.md) | ⭐ 3.3k | `development` |
| [Index](development/tools/062-index_ec82d631/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/external/open/index.md) | ⭐ 3.3k | `development` |
| [Skill](development/tools/002-name-skill_fe5d0c93/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib-theming/SKILL.md) | ⭐ 117 | `development` |
| [Conditions](development/tools/329-conditions_fbf1013a/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/r-lib/cli/references/conditions.md) | ⭐ 117 | `development` |
| [Sass And Css Variables](development/tools/330-sass-and-css-variables_db885c80/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib-theming/references/sass-and-css-variables.md) | ⭐ 117 | `development` |
| [Backlog](development/tools/336-backlog_5f579c4c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/BACKLOG.md) | ⭐ 18 | `development` |
| [Copilot Instructions](development/tools/337-copilot-instructions_727ce275/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.github/copilot-instructions.md) | ⭐ 18 | `development` |
| [Architect Validate Orchestrator Discipline](development/tools/338-architect-validate-orchestrator-discipline_2adb21c5/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/architect-validate-orchestrator-discipline.md) | ⭐ 18 | `development` |
| [Tasks 4 Validate Orchestrator Discipline](development/tools/339-tasks-4-validate-orchestrator-discipline_55f8f536/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plan/tasks-4-validate-orchestrator-discipline.md) | ⭐ 18 | `development` |
| [Orchestrator Discipline Grooming 2026 02 20](development/tools/340-orchestrator-discipline-grooming-2026-02-20_115fd850/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/grooming-reports/orchestrator-discipline-grooming-2026-02-20.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_f8e7b13a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/commit-staged/SKILL.md) | ⭐ 18 | `development` |
| [2026 02 20 Duckdb Lock Scope Flag](development/tools/341-2026-02-20-duckdb-lock-scope-flag_18db1838/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agentskill-kaizen/docs/plans/2026-02-20-duckdb-lock-scope-flag.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_0fffb33b/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/clang-format/skills/clang-format/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_98d88f7e/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/conventional-commits/skills/conventional-commits/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_97365666/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/development-harness/skills/implementation-manager/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_b65d9a7d/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/fastmcp-creator/skills/fastmcp-creator/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_900d15ae/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/agentskills/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_da240963/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/plugin-creator/skills/claude-plugins-reference-2026/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_ab8de90a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/hatchling/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_7c7e3416/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/implementation-manager/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_2be1a33c/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/python3-development/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_d5053e84/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/python3-development/skills/toml-python/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_6f3596ac/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/workshops/.claude/skills/embedded-debug-tools/SKILL.md) | ⭐ 18 | `development` |
| [Fact Checker](development/tools/329-fact-checker_7699dca6/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/agents/fact-checker.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_f380706f/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/fact-check/SKILL.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_07dab008/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/skills/integration-automation/codex-cto/SKILL.md) | ⭐ 10 | `development` |
| [Agents](development/tools/015-agents_d9ec1552/) | [bowenliang123/md_exporter](https://raw.githubusercontent.com/bowenliang123/md_exporter/main/AGENTS.md) | ⭐ 182 | `development` |
| [Claude](development/tools/017-claude_69bd09db/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/CLAUDE.md) | 🔥 9.7k | `development` |
| [Architecture](development/tools/051-architecture_97bc174a/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/ARCHITECTURE.md) | 🔥 9.7k | `development` |
| [Mcp Server](development/tools/330-mcp-server_11035d67/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/advanced/mcp-server.md) | 🔥 9.7k | `development` |
| [Mcp Reference](development/tools/331-mcp_reference_9b2d06aa/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/reference/MCP_REFERENCE.md) | 🔥 9.7k | `development` |
| [Debugging](development/tools/331-debugging_7b3de054/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DEBUGGING.md) | ⭐ 43 | `development` |
| [Evaluation Metrics](development/tools/332-evaluation_metrics_131ff394/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/EVALUATION_METRICS.md) | ⭐ 43 | `development` |
| [Tool Categories](development/tools/333-tool_categories_7ce43e0f/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TOOL_CATEGORIES.md) | ⭐ 43 | `development` |
| [Troubleshooting](development/tools/205-troubleshooting_50903b7f/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TROUBLESHOOTING.md) | ⭐ 43 | `development` |
| [Backlog](development/tools/331-backlog_385378a6/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/BACKLOG.md) | ⭐ 18 | `development` |
| [Labels](development/tools/332-labels_8c835f7a/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/gh/references/labels.md) | ⭐ 18 | `development` |
| [Milestones](development/tools/333-milestones_b673eaf0/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/gh/references/milestones.md) | ⭐ 18 | `development` |
| [Projects V2](development/tools/334-projects-v2_67064eaf/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/skills/gh/references/projects-v2.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_39c40974/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/clang-format/skills/clang-format/SKILL.md) | ⭐ 18 | `development` |
| [Typescript Mcp Server](development/tools/096-typescript-mcp-server_292a1388/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/fastmcp-creator/skills/fastmcp-creator/references/typescript-mcp-server.md) | ⭐ 18 | `development` |
| [Skill](development/tools/002-name-skill_1b580663/) | [seojoonkim/prompt-guard](https://raw.githubusercontent.com/seojoonkim/prompt-guard/main/SKILL.md) | ⭐ 94 | `development` |
| [01 Overview](development/tools/334-01-overview_dfb3df0b/) | [volcengine/OpenViking](https://raw.githubusercontent.com/volcengine/OpenViking/main/docs/en/api/01-overview.md) | ⭐ 3.2k | `development` |
| [Skill](development/tools/002-name-skill_27b948e8/) | [OmidZamani/dspy-skills](https://raw.githubusercontent.com/OmidZamani/dspy-skills/master/.claude/skills-skill-perfection/SKILL.md) | ⭐ 38 | `development` |
| [Claude](development/tools/017-claude_94ff948e/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/CLAUDE.md) | ⭐ 10 | `development` |
| [Engineer](development/tools/331-engineer_be08fbfc/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/engineer.md) | ⭐ 10 | `development` |
| [Ui Designer](development/tools/332-ui-designer_4f74bc0e/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/agents/ui-designer.md) | ⭐ 10 | `development` |
| [01 User Journey](development/tools/333-01-user-journey_31174cba/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/01-getting-started/01-user-journey.md) | ⭐ 10 | `development` |
| [03 Agent Guide](development/tools/334-03-agent-guide_096afe72/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/30-operations/03-agent-guide.md) | ⭐ 10 | `development` |
| [02 Orchestration Workflow](development/tools/335-02-orchestration-workflow_564b5875/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/02-orchestration-workflow.md) | ⭐ 10 | `development` |
| [04 Orchestration Troubleshooting](development/tools/336-04-orchestration-troubleshooting_fe428adf/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/04-orchestration-troubleshooting.md) | ⭐ 10 | `development` |
| [05 Worktree Isolation](development/tools/337-05-worktree-isolation_1c408e7e/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/05-worktree-isolation.md) | ⭐ 10 | `development` |
| [Zero Config Installation](development/tools/338-zero-config-installation_417d1ca7/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/zero-config-installation.md) | ⭐ 10 | `development` |
| [02 Upgrade Guide](development/tools/339-02-upgrade-guide_df24dcca/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/70-reference/02-upgrade-guide.md) | ⭐ 10 | `development` |
| [Changelog Slim](development/tools/340-changelog-slim_aca6472f/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/mcp-servers/copilot-memory/CHANGELOG-SLIM.md) | ⭐ 10 | `development` |
| [Implementation](development/tools/341-implementation_5d198633/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/scripts/knowledge-sync/IMPLEMENTATION.md) | ⭐ 10 | `development` |
| [03 Technical Architecture](development/tools/342-03-technical-architecture_eb886e76/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/tasks/interactive-tui/03-technical-architecture.md) | ⭐ 10 | `development` |
| [Skill](development/tools/002-name-skill_5cc9ba34/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-coding/skills/push-pr/SKILL.md) | ⭐ 34 | `development` |
| [Lenses](development/tools/343-lenses_d6bf568b/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-memory/skills/past-conversations/references/lenses.md) | ⭐ 34 | `development` |
| [Skill Anatomy](development/tools/344-skill-anatomy_b3b8b3b0/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-skills/skills/skill-repair/references/skill-anatomy.md) | ⭐ 34 | `development` |
| [Types](development/tools/331-types_bae5a72e/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/api/types.md) | ⭐ 106 | `development` |
| [Skill](development/tools/002-name-skill_41f381c2/) | [frmoretto/hardstop](https://raw.githubusercontent.com/frmoretto/hardstop/main/.claude/skills/hs/SKILL.md) | ⭐ 16 | `development` |
| [Readme Flat Clients Created](development/tools/157-readme_flat_clients_created_b24ab63b/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLIENTS_CREATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Docs Releases](development/tools/163-readme_flat_docs_releases_82fd8d60/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_DOCS_RELEASES.md) | 🔥 24.5k | `development` |
| [Readme Flat Hooks Created](development/tools/166-readme_flat_hooks_created_04ada004/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_HOOKS_CREATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Statusline Created](development/tools/171-readme_flat_statusline_created_26ae1ec4/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STATUSLINE_CREATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Styles Created](development/tools/175-readme_flat_styles_created_41e22aa8/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STYLES_CREATED.md) | 🔥 24.5k | `development` |
| [Readme Flat Workflows Created](development/tools/179-readme_flat_workflows_created_7a42507e/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_WORKFLOWS_CREATED.md) | 🔥 24.5k | `development` |
| [Skill](development/tools/002-name-skill_c3528e8b/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-smell-detector/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/tools/002-name-skill_248c9a31/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/dead-code-eliminator/SKILL.md) | ⭐ 10 | `development` |
| [Skill](development/tools/002-name-skill_13bb49ce/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/legacy-code-summarizer/SKILL.md) | ⭐ 10 | `development` |
| [Dead Code Patterns](development/tools/330-dead-code-patterns_6ae799f4/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/dead-code-eliminator/references/dead-code-patterns.md) | ⭐ 10 | `development` |
| [Migration Patterns](development/tools/331-migration_patterns_d4bc0348/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/framework-migration-assistant/references/migration_patterns.md) | ⭐ 10 | `development` |
| [Code Quality Checklist](development/tools/332-code_quality_checklist_8e709a3a/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/legacy-code-summarizer/references/code_quality_checklist.md) | ⭐ 10 | `development` |
| [Python Contracts](development/tools/333-python_contracts_a8a56814/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/static-reasoning-verifier/references/python_contracts.md) | ⭐ 10 | `development` |
| [Readme Flat Hooks Releases](development/tools/167-readme_flat_hooks_releases_75b26790/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_HOOKS_RELEASES.md) | 🔥 24.4k | `development` |
| [Readme Flat Docs Az](development/tools/161-readme_flat_docs_az_3c66a52c/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_DOCS_AZ.md) | 🔥 24.4k | `development` |
| [Readme Flat Styles Releases](development/tools/176-readme_flat_styles_releases_3f7c42c6/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STYLES_RELEASES.md) | 🔥 24.4k | `development` |
| [Readme Flat Workflows Releases](development/tools/180-readme_flat_workflows_releases_676983da/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_WORKFLOWS_RELEASES.md) | 🔥 24.4k | `development` |
| [Readme Flat Statusline Updated](development/tools/173-readme_flat_statusline_updated_202be4cb/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STATUSLINE_UPDATED.md) | 🔥 24.4k | `development` |
| [Readme Flat Workflows Updated](development/tools/181-readme_flat_workflows_updated_f60e4311/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_WORKFLOWS_UPDATED.md) | 🔥 24.4k | `development` |
| [Readme Flat Clients Az](development/tools/156-readme_flat_clients_az_34ec63b2/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLIENTS_AZ.md) | 🔥 24.6k | `development` |
| [Readme Flat Clients Updated](development/tools/159-readme_flat_clients_updated_038bed5f/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_CLIENTS_UPDATED.md) | 🔥 24.6k | `development` |
| [Readme Flat Commands Releases](development/tools/160-readme_flat_commands_releases_34ade65c/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_COMMANDS_RELEASES.md) | 🔥 24.6k | `development` |
| [Readme Flat Docs Created](development/tools/162-readme_flat_docs_created_fc64977c/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_DOCS_CREATED.md) | 🔥 24.6k | `development` |
| [Readme Flat Hooks Az](development/tools/165-readme_flat_hooks_az_687deee8/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_HOOKS_AZ.md) | 🔥 24.6k | `development` |
| [Readme Flat Statusline Az](development/tools/170-readme_flat_statusline_az_fb5537d0/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STATUSLINE_AZ.md) | 🔥 24.6k | `development` |
| [Readme Flat Styles Az](development/tools/174-readme_flat_styles_az_ef4a07bb/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_STYLES_AZ.md) | 🔥 24.6k | `development` |
| [Readme Flat Workflows Az](development/tools/178-readme_flat_workflows_az_ec5961ff/) | [hesreallyhim/awesome-claude-code](https://raw.githubusercontent.com/hesreallyhim/awesome-claude-code/main/README_ALTERNATIVES/README_FLAT_WORKFLOWS_AZ.md) | 🔥 24.6k | `development` |
| [Backlog](development/tools/368-backlog_23d86cd2/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/.claude/BACKLOG.md) | ⭐ 18 | `development` |
| [Skill](development/tools/086-description-skill_aebb20d3/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/workshops/.cursor/skills/embedded-debug-tools/SKILL.md) | ⭐ 18 | `development` |
| [Claude](development/tools/017-claude_2b3b5dfe/) | [doobidoo/MCP-Context-Provider](https://raw.githubusercontent.com/doobidoo/MCP-Context-Provider/main/CLAUDE.md) | ⭐ 24 | `development` |
| [Api Surface](development/tools/345-api-surface_956ec1cb/) | [synth-laboratories/synth-ai](https://raw.githubusercontent.com/synth-laboratories/synth-ai/main/skills/synth-smr-control/references/api-surface.md) | ⭐ 75 | `development` |
| [Iam Settings Ui](development/tools/335-iam-settings-ui_2f16319a/) | [agentic-community/mcp-gateway-registry](https://raw.githubusercontent.com/agentic-community/mcp-gateway-registry/main/docs/iam-settings-ui.md) | ⭐ 449 | `development` |
| [Tmux Claude Workflow](development/tools/tmux-claude-workflow_ec924085/) | [tdimino/claude-code-minoan](https://raw.githubusercontent.com/tdimino/claude-code-minoan/main/docs/tmux-claude-workflow.md) | ⭐ 10 | `development` |

### Investment (9 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](investment/021-name-skill_65b23c64/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/whale-alert-monitor/skills/monitoring-whale-activity/SKILL.md) | ⭐ 1.4k | `investment` |
| [Input Schema](investment/048-input-schema_cad1183d/) | [tradermonty/claude-trading-skills](https://raw.githubusercontent.com/tradermonty/claude-trading-skills/main/skills/kanchi-dividend-review-monitor/references/input-schema.md) | ⭐ 41 | `investment` |
| [Theme Detection Methodology](investment/049-theme_detection_methodology_bfb5c432/) | [tradermonty/claude-trading-skills](https://raw.githubusercontent.com/tradermonty/claude-trading-skills/main/skills/theme-detector/references/theme_detection_methodology.md) | ⭐ 41 | `investment` |
| [Skill](investment/050-description-skill_8b0ec79b/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-thinking/skills/thinking-partner/SKILL.md) | ⭐ 34 | `investment` |
| [Skill](investment/021-name-skill_76ec2b17/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/kalshi/SKILL.md) | ⭐ 24 | `investment` |
| [Skill](investment/021-name-skill_b71219b9/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/polymarket/SKILL.md) | ⭐ 24 | `investment` |
| [Skill](investment/021-name-skill_9fd4c364/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/wnba-data/SKILL.md) | ⭐ 24 | `investment` |
| [Commands](investment/050-commands_b3537439/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/kalshi/references/commands.md) | ⭐ 24 | `investment` |
| [Series Tickers](investment/051-series-tickers_fb5fea49/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/kalshi/references/series-tickers.md) | ⭐ 24 | `investment` |

### Other (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Team Ids](other/036-team-ids_bb32d02b/) | [machina-sports/sports-skills](https://raw.githubusercontent.com/machina-sports/sports-skills/main/skills/mlb-data/references/team-ids.md) | ⭐ 24 | `other` |

### Productivity (10 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Parallel Session Cleanup](productivity/179-parallel-session-cleanup_d83f616a/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/manage/parallel-session-cleanup.md) | ⭐ 3.3k | `productivity` |
| [Dark Mode](productivity/176-dark-mode_77760605/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib-theming/references/dark-mode.md) | ⭐ 117 | `productivity` |
| [Tooltips Popovers](productivity/177-tooltips-popovers_f123d4bd/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/shiny/shiny-bslib/references/tooltips-popovers.md) | ⭐ 117 | `productivity` |
| [Protocol](productivity/173-protocol_d5300386/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/.claude/commands/protocol.md) | ⭐ 10 | `productivity` |
| [Correction Detection](productivity/174-correction-detection_42a7a5e8/) | [Everyone-Needs-A-Copilot/claude-copilot](https://raw.githubusercontent.com/Everyone-Needs-A-Copilot/claude-copilot/main/docs/50-features/correction-detection.md) | ⭐ 10 | `correction` `user-feedback` |
| [Skill](productivity/093-name-skill_a23fac79/) | [gupsammy/Claudest](https://raw.githubusercontent.com/gupsammy/Claudest/main/plugins/claude-coding/skills/updateclaudemd/SKILL.md) | ⭐ 34 | `productivity` |
| [Skill](productivity/093-name-skill_34f178e8/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/code-optimizer/SKILL.md) | ⭐ 10 | `productivity` |
| [Skill](productivity/093-name-skill_1f0c7074/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/refinement-step-generator/SKILL.md) | ⭐ 10 | `productivity` |
| [Skill](productivity/093-name-skill_7d504321/) | [sickn33/antigravity-awesome-skills](https://raw.githubusercontent.com/sickn33/antigravity-awesome-skills/main/skills/multi-agent-brainstorming/SKILL.md) | 🔥 13.4k | `productivity` |
| [Skill](productivity/093-name-skill_3e8d7700/) | [harrymunro/nelson](https://raw.githubusercontent.com/harrymunro/nelson/main/skills/nelson/SKILL.md) | ⭐ 143 | `productivity` |

### Research (34 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Pandoc Server](research/265-pandoc-server_a65d4d68/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/go/pandoc-server.md) | ⭐ 3.3k | `research` |
| [Latex Server](research/266-latex-server_e6f7093b/) | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/using/servers/python/latex-server.md) | ⭐ 3.3k | `research` |
| [Citations](research/261-citations_00e0d6bf/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/citations.md) | ⭐ 117 | `research` |
| [Divs And Spans](research/262-divs-and-spans_a11251da/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/divs-and-spans.md) | ⭐ 117 | `research` |
| [Layout](research/263-layout_f4e230b1/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/layout.md) | ⭐ 117 | `research` |
| [Yaml Front Matter](research/264-yaml-front-matter_8c124133/) | [posit-dev/skills](https://raw.githubusercontent.com/posit-dev/skills/main/quarto/quarto-authoring/references/yaml-front-matter.md) | ⭐ 117 | `research` |
| [Claude](research/015-claude_b0ed45bf/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/CLAUDE.md) | ⭐ 829 | `research` |
| [Writer](research/144-writer_1664b873/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/.claude/WRITER.md) | ⭐ 829 | `research` |
| [Api](research/146-api_a79b5e1a/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/docs/API.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_83163c0c/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/citation-management/SKILL.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_995bd408/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/SKILL.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_ed983e74/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/research-lookup/SKILL.md) | ⭐ 829 | `research` |
| [Citation Validation](research/014-citation_validation_5901b5ff/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/citation-management/references/citation_validation.md) | ⭐ 829 | `research` |
| [Api Reference](research/007-api_reference_bd3cf235/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/api_reference.md) | ⭐ 829 | `research` |
| [Deep Research Guide](research/258-deep_research_guide_c3d92b15/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/deep_research_guide.md) | ⭐ 829 | `research` |
| [Extraction Patterns](research/259-extraction_patterns_49ce48cf/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/extraction_patterns.md) | ⭐ 829 | `research` |
| [Search Best Practices](research/260-search_best_practices_8ba3e7fe/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/search_best_practices.md) | ⭐ 829 | `research` |
| [Workflow Recipes](research/261-workflow_recipes_0676ea23/) | [K-Dense-AI/claude-scientific-writer](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-writer/main/skills/parallel-web/references/workflow_recipes.md) | ⭐ 829 | `research` |
| [Skill](research/139-name-skill_41583fad/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/cross-chain-bridge-monitor/skills/monitoring-cross-chain-bridges/SKILL.md) | ⭐ 1.4k | `research` |
| [Fleet Config](research/257-fleet-config_ee28fcec/) | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/fleet-management/reference/fleet-config.md) | ⭐ 61 | `]                            # Optional: key:value strings for filtering
    memory_blocks: [` |
| [Introduction](research/267-introduction_9d11f1af/) | [MooseGoose0701/skill-compose](https://raw.githubusercontent.com/MooseGoose0701/skill-compose/main/docs/docs/introduction.md) | ⭐ 820 | `research` |
| [02 Quick Start](research/258-02-quick-start_aa356bca/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/getting-started/02-quick-start.md) | 🔥 9.7k | `research` |
| [01 Core Concepts](research/259-01-core-concepts_8eaa4c6e/) | [yusufkaraaslan/Skill_Seekers](https://raw.githubusercontent.com/yusufkaraaslan/Skill_Seekers/development/docs/user-guide/01-core-concepts.md) | 🔥 9.7k | `research` |
| [Clear Framework](research/168-clear-framework_08305a60/) | [Jamie-BitFlight/claude_skills](https://raw.githubusercontent.com/Jamie-BitFlight/claude_skills/main/plugins/agent-orchestration/skills/agent-orchestration/clear-framework.md) | ⭐ 18 | `research` |
| [Concepts](research/265-concepts_6e75c0f3/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/concepts.md) | ⭐ 106 | `research` |
| [Creating Skills](research/266-creating-skills_2fd29c76/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/creating-skills.md) | ⭐ 106 | `csv` `data` `analysis` |
| [Index](research/267-index_e98d9570/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/index.md) | ⭐ 106 | `research` |
| [Quick Start](research/268-quick-start_543c444a/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/quick-start.md) | ⭐ 106 | `research` |
| [Toolset](research/269-toolset_53992f79/) | [DougTrajano/pydantic-ai-skills](https://raw.githubusercontent.com/DougTrajano/pydantic-ai-skills/main/docs/api/toolset.md) | ⭐ 106 | `research` |
| [Skill](research/139-name-skill_b96d1e1c/) | [ThepExcel/agent-skills](https://raw.githubusercontent.com/ThepExcel/agent-skills/main/deep-research/SKILL.md) | ⭐ 15 | `research` |
| [Tlaplus Syntax](research/258-tlaplus_syntax_cf31ead9/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/program-to-tlaplus-spec-generator/references/tlaplus_syntax.md) | ⭐ 10 | `research` |
| [Tlaplus Syntax](research/258-tlaplus_syntax_54ea0eea/) | [ArabelaTso/Skills-4-SE](https://raw.githubusercontent.com/ArabelaTso/Skills-4-SE/main/skills/tlaplus-spec-generator/references/tlaplus_syntax.md) | ⭐ 10 | `research` |
| [Docx Js](research/258-docx-js_0d24521f/) | [jjyaoao/HelloAgents](https://raw.githubusercontent.com/jjyaoao/HelloAgents/main/skills/docx/docx-js.md) | ⭐ 600 | `research` |
| [Faq](research/153-faq_3a909ea1/) | [LearningCircuit/local-deep-research](https://raw.githubusercontent.com/LearningCircuit/local-deep-research/main/docs/faq.md) | ⭐ 4.0k | `research` |

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

*Last updated: 2026-02-23 08:47:10 UTC*
*Automatically maintained by SkillFlow*
