# X-Skills

A curated collection of **84 AI-powered skills** organized into 10 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (3 skills)
- **Automation/Workflow** (4 skills)
- **Communication** (9 skills)
- **Content Creation** (10 skills)
- **Daily Assistant** (1 skill)
- **Data Analysis** (23 skills)
- **Development** (24 skills)
- **Development/Testing** (2 skills)
- **Development/Tools** (7 skills)
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
| [Automation Agent](automation-agent/) | 40 | automation | `xskills patch install automation-agent` |
| [Communication Agent](communication-agent/) | 40 | communication | `xskills patch install communication-agent` |
| [Content Creator](content-creator/) | 60 | content-creation | `xskills patch install content-creator` |
| [Data Analyst](data-analyst/) | 60 | data-analysis | `xskills patch install data-analyst` |
| [DevOps Engineer](devops-engineer/) | 60 | development | `xskills patch install devops-engineer` |
| [Productivity Assistant](productivity-assistant/) | 50 | productivity, daily-assistant | `xskills patch install productivity-assistant` |
| [Python Developer](python-dev/) | 50 | development | `xskills patch install python-dev` |
| [Research Agent](research-agent/) | 30 | research | `xskills patch install research-agent` |
| [Web Development Agent](web-dev-agent/) | 80 | development | `xskills patch install web-dev-agent` |

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


### Automation/Scripting (3 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Workflow Nodes Memories](automation/scripting/workflow_nodes_memories_dbd21b95/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Workflow_Nodes_Memories.md) | ⭐ 803 | `automation` |
| [Chatsummarysummarizer](automation/scripting/chatsummarysummarizer_ca7a6b00/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/ChatSummarySummarizer.md) | ⭐ 803 | `automation` |
| [Workflows](automation/scripting/workflows_5e8f8d55/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Workflows.md) | ⭐ 803 | `automation` |

### Automation/Workflow (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Timestamps And Time Tracking](automation/workflow/timestamps_and_time_tracking_882ee04e/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/TimeStamps_And_Time_Tracking.md) | ⭐ 803 | `automation` |
| [Customworkflow](automation/workflow/customworkflow_0379bd6e/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/CustomWorkflow.md) | ⭐ 803 | `automation` |
| [Getcustomfile](automation/workflow/getcustomfile_87d8b837/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/GetCustomFile.md) | ⭐ 803 | `automation` |
| [Nested Workflows](automation/workflow/nested_workflows_397f680a/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Nested_Workflows.md) | ⭐ 803 | `automation` |

### Communication (9 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Adaptable Front End Api](communication/adaptable_front_end_api_924b1f8e/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Adaptable_Front_End_Api.md) | ⭐ 803 | `communication` |
| [Apitype](communication/apitype_0d6a86e5/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Configuration_Files/ApiType.md) | ⭐ 803 | `communication` |
| [Endpoint](communication/endpoint_0de3946c/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Configuration_Files/Endpoint.md) | ⭐ 803 | `communication` |
| [Preset](communication/preset_469fe489/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Configuration_Files/Preset.md) | ⭐ 803 | `communication` |
| [User](communication/user_d2c5e1bf/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Configuration_Files/User.md) | ⭐ 803 | `communication` |
| [Workflowlock](communication/workflowlock_c725e3d9/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/WorkflowLock.md) | ⭐ 803 | `communication` |
| [Fullchatsummary](communication/fullchatsummary_28dcdaf5/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/FullChatSummary.md) | ⭐ 803 | `communication` |
| [Recentmemorysummarizertool](communication/recentmemorysummarizertool_1589406d/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/RecentMemorySummarizerTool.md) | ⭐ 803 | `communication` |
| [Vectormemorysearch](communication/vectormemorysearch_09ff1270/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/VectorMemorySearch.md) | ⭐ 803 | `communication` |

### Content Creation (10 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Prompt Routing](content-creation/prompt_routing_571b41bd/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Prompt_Routing.md) | ⭐ 803 | `content creation` |
| [Workflow Prompting Methodologies Socg](content-creation/workflow_prompting_methodologies_socg_27eedd83/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/LLM_Assisted_Workflow_Generation/Workflow_Prompting_Methodologies_Socg.md) | ⭐ 803 | `content creation` |
| [Open Webui](content-creation/open-webui_7a12e1af/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Open-WebUI.md) | ⭐ 803 | `content creation` |
| [Getting Start Wilmer Api](content-creation/_getting-start_wilmer-api_770a3f14/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/_Getting-Start_Wilmer-Api.md) | ⭐ 803 | `content creation` |
| [Savecustomfile](content-creation/savecustomfile_75540104/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/SaveCustomFile.md) | ⭐ 803 | `content creation` |
| [Standard Conversational](content-creation/standard_conversational_cca1bd67/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Standard_Conversational.md) | ⭐ 803 | `content creation` |
| [Vector Qualitymemory](content-creation/vector_qualitymemory_572fe7a4/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/Vector_QualityMemory.md) | ⭐ 803 | `content creation` |
| [Custom Python Scripts](content-creation/custom_python_scripts_8a0ada32/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Custom_Python_Scripts.md) | ⭐ 803 | `content creation` |
| [Memories](content-creation/memories_7d437633/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Memories.md) | ⭐ 803 | `content creation` |
| [Example Guide Wikipedia Search](content-creation/example_guide_wikipedia_search_799fa025/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/LLM_Assisted_Workflow_Generation/Example_Guide_Wikipedia_Search.md) | ⭐ 803 | `content creation` |

### Daily Assistant (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Default Endpoints And Presets](daily-assistant/default_endpoints_and_presets_abc32ef8/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/LLM_Assisted_Workflow_Generation/Default_Endpoints_And_Presets.md) | ⭐ 803 | `daily assistant` |

### Data Analysis (23 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Geepers Business Plan](data-analysis/423-geepers_business_plan_c0dd00dd/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/agents/geepers_business_plan.md) | ⭐ 1.3k | `data analysis` |
| [Geepers Corpus Ux](data-analysis/424-geepers_corpus_ux_e56c7603/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/agents/geepers_corpus_ux.md) | ⭐ 1.3k | `data analysis` |
| [Geepers Design](data-analysis/425-geepers_design_d56e9e29/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/agents/geepers_design.md) | ⭐ 1.3k | `data analysis` |
| [Geepers Orchestrator Fullstack](data-analysis/426-geepers_orchestrator_fullstack_5e055b9f/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/agents/geepers_orchestrator_fullstack.md) | ⭐ 1.3k | `data analysis` |
| [Geepers Corpus Ux](data-analysis/424-geepers_corpus_ux_b0a49b87/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/corpus/geepers_corpus_ux.md) | ⭐ 1.3k | `data analysis` |
| [Geepers Design](data-analysis/425-geepers_design_c93ca317/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/fullstack/geepers_design.md) | ⭐ 1.3k | `data analysis` |
| [Geepers Orchestrator Fullstack](data-analysis/426-geepers_orchestrator_fullstack_d00b9e50/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/fullstack/geepers_orchestrator_fullstack.md) | ⭐ 1.3k | `data analysis` |
| [Geepers Business Plan](data-analysis/423-geepers_business_plan_1cf8d17c/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/geepers-agents/product/geepers_business_plan.md) | ⭐ 1.3k | `data analysis` |
| [004 Gemini Supervised Fine Tuning Predictive Maintenance](data-analysis/427-004-gemini-supervised-fine-tuning-predictive-maintenance_d16113fc/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/community/jeremy-firebase/000-usermanuals/004-gemini-supervised-fine-tuning-predictive-maintenance.md) | ⭐ 1.3k | `data analysis` |
| [Analyze Sentiment](data-analysis/428-analyze-sentiment_6d0a94b5/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/crypto/market-sentiment-analyzer/commands/analyze-sentiment.md) | ⭐ 1.3k | `data analysis` |
| [Db Diff](data-analysis/429-db-diff_5e7c5b4f/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/database-diff-tool/commands/db-diff.md) | ⭐ 1.3k | `data analysis` |
| [Db Docs](data-analysis/430-db-docs_a382d1fe/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/database-documentation-gen/commands/db-docs.md) | ⭐ 1.3k | `data analysis` |
| [Sharding](data-analysis/431-sharding_66505e24/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/database-sharding-manager/commands/sharding.md) | ⭐ 1.3k | `data analysis` |
| [Transactions](data-analysis/432-transactions_5321f762/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/database/database-transaction-monitor/commands/transactions.md) | ⭐ 1.3k | `data analysis` |
| [Openbb Portfolio](data-analysis/433-openbb-portfolio_7ad8614f/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/finance/openbb-terminal/commands/openbb-portfolio.md) | ⭐ 1.3k | `data analysis` |
| [Currency](data-analysis/434-currency_358fbd92/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/productivity/travel-assistant/commands/currency.md) | ⭐ 1.3k | `data analysis` |
| [Check Deps](data-analysis/435-check-deps_6e67b65b/) | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/security/dependency-checker/commands/check-deps.md) | ⭐ 1.3k | `data analysis` |
| [Wilmer Prompt Flow Beginning To End](data-analysis/wilmer_prompt_flow_beginning_to_end_777dcc2d/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Wilmer_Prompt_Flow_Beginning_To_End.md) | ⭐ 803 | `data analysis` |
| [Prompttemplates](data-analysis/prompttemplates_e34645fb/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Configuration_Files/PromptTemplates.md) | ⭐ 803 | `data analysis` |
| [Arithmetic Processor](data-analysis/arithmetic_processor_9d59a18e/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Arithmetic_Processor.md) | ⭐ 803 | `data analysis` |
| [Jsonextractor](data-analysis/jsonextractor_d1957e17/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/JsonExtractor.md) | ⭐ 803 | `data analysis` |
| [Offline Wikipedia](data-analysis/offline_wikipedia_62ac0241/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Offline_Wikipedia.md) | ⭐ 803 | `data analysis` |
| [Pythonmodule](data-analysis/pythonmodule_95a9134b/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/PythonModule.md) | ⭐ 803 | `data analysis` |

### Development (24 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Readme Mcp Tools](development/readme_mcp_tools_031f87ef/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Public/modules/README_MCP_TOOLS.md) | ⭐ 803 | `development` |
| [Custom Workflows](development/custom_workflows_042ab23e/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Custom_Workflows.md) | ⭐ 803 | `development` |
| [Llm Apis](development/llm_apis_58ed0d2e/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/LLM_Apis.md) | ⭐ 803 | `development` |
| [Llm Apis Remove Unwanted Text From Llm Stream](development/llm_apis_remove_unwanted_text_from_llm_stream_c1a60822/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/LLM_Apis_Remove_Unwanted_Text_From_Llm_Stream.md) | ⭐ 803 | `development` |
| [Memories](development/memories_29aca277/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Memories.md) | ⭐ 803 | `development` |
| [Workflow Variables](development/workflow_variables_c92aa314/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Workflow_Variables.md) | ⭐ 803 | `development` |
| [Adaptable Llm Backend Apis](development/adaptable_llm_backend_apis_651d6bc3/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Adaptable_LLM_Backend_Apis.md) | ⭐ 803 | `development` |
| [Offline Wikipedia Support](development/offline_wikipedia_support_bc9a006d/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Offline_Wikipedia_Support.md) | ⭐ 803 | `development` |
| [Prompt Routing In Workflow](development/prompt_routing_in_workflow_a60b1458/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Prompt_Routing_In_Workflow.md) | ⭐ 803 | `development` |
| [Workflow Jinja2 Support](development/workflow_jinja2_support_c421c6a0/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Core_Features/Workflow_Jinja2_Support.md) | ⭐ 803 | `development` |
| [Recursive Workflows](development/recursive_workflows_0659acaa/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/LLM_Assisted_Workflow_Generation/Recursive_Workflows.md) | ⭐ 803 | `development` |
| [Sillytavern](development/sillytavern_252c7120/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/SillyTavern.md) | ⭐ 803 | `development` |
| [Routing](development/routing_b62af615/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Configuration_Files/Routing.md) | ⭐ 803 | `development` |
| [Workflow Features](development/workflow_features_82cae8e7/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Workflow_Features.md) | ⭐ 803 | `development` |
| [Workflow Nodes](development/workflow_nodes_81a3fb81/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Workflow_Nodes.md) | ⭐ 803 | `development` |
| [Workflow Variables](development/workflow_variables_ecdd2861/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Workflow_Variables.md) | ⭐ 803 | `development` |
| [Conditional](development/conditional_005838a8/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Conditional.md) | ⭐ 803 | `development` |
| [Conditionalcustomworkflow](development/conditionalcustomworkflow_27739597/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/ConditionalCustomWorkflow.md) | ⭐ 803 | `development` |
| [Image Processor](development/image_processor_8ffed74d/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Image_Processor.md) | ⭐ 803 | `development` |
| [Staticresponse](development/staticresponse_72c49d77/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/StaticResponse.md) | ⭐ 803 | `development` |
| [Tagtextextractor](development/tagtextextractor_32898c82/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/TagTextExtractor.md) | ⭐ 803 | `development` |
| [Getcurrentsummaryfromfile](development/getcurrentsummaryfromfile_e4c10994/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/GetCurrentSummaryFromFile.md) | ⭐ 803 | `development` |
| [Recentmemory](development/recentmemory_7f32f8f0/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/RecentMemory.md) | ⭐ 803 | `development` |
| [Conversationmemory](development/conversationmemory_ed41173d/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/Legacy_And_Less_Used_Nodes/ConversationMemory.md) | ⭐ 803 | `development` |

### Development/Testing (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Llm Assisted Setup](development/testing/llm_assisted_setup_4f7a18d7/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/LLM_Assisted_Setup.md) | ⭐ 803 | `development` |
| [Unit Tests Llm Apis](development/testing/unit_tests_llm_apis_91d947a5/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Unit_Tests_LLM_Apis.md) | ⭐ 803 | `development` |

### Development/Tools (7 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Api](development/tools/api_1513c7a4/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Api.md) | ⭐ 803 | `development` |
| [Unit Tests](development/tools/unit_tests_74064428/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Unit_Tests.md) | ⭐ 803 | `development` |
| [Workflows](development/tools/workflows_46d35683/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/Developer_Docs/Features_And_Packages/Workflows.md) | ⭐ 803 | `development` |
| [Workflows](development/tools/workflows_479d3243/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Workflows.md) | ⭐ 803 | `development` |
| [String Concatenator](development/tools/string_concatenator_a058f360/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/String_Concatenator.md) | ⭐ 803 | `development` |
| [Release V3.1.0](development/tools/release-v310_fecf7cf6/) | [seojoonkim/prompt-guard](https://raw.githubusercontent.com/seojoonkim/prompt-guard/main/RELEASE-v3.1.0.md) | ⭐ 64 | `development` |
| [Skill](development/tools/name-skill_9d61a891/) | [seojoonkim/prompt-guard](https://raw.githubusercontent.com/seojoonkim/prompt-guard/main/SKILL.md) | ⭐ 64 | `development` |

### Research (1 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Chatsummarymemorygatheringtool](research/chatsummarymemorygatheringtool_f7c70128/) | [SomeOddCodeGuy/WilmerAI](https://raw.githubusercontent.com/SomeOddCodeGuy/WilmerAI/master/Docs/User_Documentation/Setup/Workflow_Details/Nodes/Memory_Nodes/Legacy_And_Less_Used_Nodes/ChatSummaryMemoryGatheringTool.md) | ⭐ 803 | `research` |

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

*Last updated: 2026-02-09 04:31:59 UTC*
*Automatically maintained by SkillFlow*
