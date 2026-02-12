# X-Skills

A curated collection of **326 AI-powered skills** organized into 14 categories.

## Overview

This repository contains automatically aggregated skills from various open-source projects. Each skill is designed to work with AI assistants like Claude Code to automate specific tasks.

## Categories

- **Automation/Scripting** (43 skills)
- **Automation/Workflow** (2 skills)
- **Commercial** (2 skills)
- **Communication** (20 skills)
- **Content Creation** (18 skills)
- **Daily Assistant** (26 skills)
- **Data Analysis** (21 skills)
- **Development** (82 skills)
- **Development/Devops** (8 skills)
- **Development/Testing** (2 skills)
- **Development/Tools** (17 skills)
- **Investment** (76 skills)
- **Other** (4 skills)
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


### Automation/Scripting (43 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Skill](automation/scripting/name-skill_70005990/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_EmitSoundFilter/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_5df4f174/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_EmitSoundParams/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_c2eca301/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_IsAlive-AND-CBaseEntity_GetEyePosition-AND-CBasePlayerPawn_GetEyePosition/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_e144da4e/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_SetStateChanged/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_6d3ce05f/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_TakeDamageOld/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_df204dab/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerController_SetPawn/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_61ec8170/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerPawn_CommitSuicide/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_39f4fb87/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerPawn_DropActivePlayerWeapon-AND-CCSPlayer_ItemServices_DropActivePlayerWeapon-AND-CCSPlayer_WeaponServices_DropWeapon/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_25b8d3ce/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerPawn_GetEyeAngles/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_7d2de804/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayerController_ChangeTeam/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_65ef3acb/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayerController_InventoryUpdateThink/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_b7dc2f45/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayerController_Respawn/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_9061b023/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayerPawnBase_PostThink/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_f80fd8dc/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayer_ItemServices_CanAcquire/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_90165ae9/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayer_ItemServices_RemoveWeapons/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_8cf1a6bd/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayer_WeaponServices_SelectItem/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_68711a03/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CEntitySystem_AddEntityIOEvent/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_ea59362f/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CGameResourceService_BuildResourceManifest/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_95d1d913/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CPointTeleport_Teleport-AND-CBaseEntity_Teleport/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_393f49c0/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CSource2GameEntities_CheckTransmit-AND-CCheckTransmitInfo/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_fc79156d/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-NetworkStateChanged/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_cc90869c/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_SetGroundEntity/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_c992fdee/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseModelEntity_SetModel-AND-CBaseEntity_SetGravityScale/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_13d4116b/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerController_HandleCommand_JoinTeam/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_d345296d/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerController_SetPlayerName/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_1bb94861/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSGameRules_FindPickerEntity/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_65a8d466/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayerController_SwitchTeam/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_232be25b/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayer_WeaponServices_CanUse/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_5b1a9239/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CEntityIOOutput_FireOutputInternal/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_72c8d4e4/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CEntityIdentity_AcceptInput/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_cbed7c1b/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CGameEntitySystem_FindEntityByName/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_dc671043/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CGameRules_ClientSettingsChanged/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_73dfc733/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CServerSideClient_IsHearingClient/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_034ba7f8/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CTakeDamageInfo_GetWeaponName/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_2e135196/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CTriggerGravity_GravityTouch/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_46cf78ac/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-ClientPrint/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_ba2bb923/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-GetCSWeaponDataFromKey/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_47d42b1d/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-Host_Say-AND-UTIL_SayTextFilter-AND-UTIL_SayTextFilter2/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_ac2fb943/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-UTIL_ClientPrintAll/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_9241fc69/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-UTIL_CreateEntityByName/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_0814ab81/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-UTIL_PlayerSlotToPlayerController/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_c949942b/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-UTIL_Remove/SKILL.md) | ⭐ 16 | `automation` |
| [Skill](automation/scripting/name-skill_cccebfba/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-WeaponBuy/SKILL.md) | ⭐ 16 | `automation` |

### Automation/Workflow (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Executive Pitch](automation/workflow/067-executive-pitch_61fed1da/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/docs/presenter/executive-pitch.md) | ⭐ 60 | `automation` |
| [Detecting Llm Hallucinations In Ci](automation/workflow/138-detecting-llm-hallucinations-in-ci_471c775c/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/guides/detecting-llm-hallucinations-in-ci.md) | ⭐ 43 | `automation` |

### Commercial (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Tasks](commercial/tasks_b3f6b25e/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-custom-ai-workflow/tasks.md) | ⭐ 19 | `commercial` |
| [Spec](commercial/spec_b955f4f4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-news-performance/specs/news-performance/spec.md) | ⭐ 19 | `commercial` |

### Communication (20 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Help](communication/251-help_dcae5407/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/help.md) | 🔥 40.5k | `communication` |
| [Model Warnings](communication/252-model-warnings_ba625717/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/model-warnings.md) | 🔥 40.5k | `communication` |
| [Multi Line](communication/253-multi-line_5264f142/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_includes/multi-line.md) | 🔥 40.5k | `communication` |
| [Git](communication/254-git_1469e351/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/git.md) | 🔥 40.5k | `communication` |
| [Docker](communication/255-docker_24d9f59d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/install/docker.md) | 🔥 40.5k | `communication` |
| [Analytics](communication/256-analytics_073a41ba/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/analytics.md) | 🔥 40.5k | `communication` |
| [Edit Errors](communication/257-edit-errors_29cca0f5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/edit-errors.md) | 🔥 40.5k | `communication` |
| [Caching](communication/258-caching_50e5e286/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/caching.md) | 🔥 40.5k | `communication` |
| [Notifications](communication/144-notifications_38a3472f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/notifications.md) | 🔥 40.5k | `communication` |
| [Cost Tracking](communication/259-cost_tracking_b8c7b063/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/COST_TRACKING.md) | ⭐ 43 | `communication` |
| [Debugging](communication/260-debugging_6f2c48f9/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DEBUGGING.md) | ⭐ 43 | `communication` |
| [Trace Spec](communication/261-trace_spec_24289c84/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TRACE_SPEC.md) | ⭐ 43 | `communication` |
| [Design](communication/design_6d130448/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-system-logs-ui/design.md) | ⭐ 19 | `communication` |
| [Design](communication/design_c8e52461/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-scoped-features/design.md) | ⭐ 19 | `communication` |
| [Tasks](communication/tasks_600ce05f/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-scoped-features/tasks.md) | ⭐ 19 | `communication` |
| [Spec](communication/spec_22f7467d/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-scoped-features/specs/user-chat-history/spec.md) | ⭐ 19 | `communication` |
| [Proposal](communication/proposal_62e59c2b/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-auth/proposal.md) | ⭐ 19 | `communication` |
| [Tasks](communication/tasks_d0b1d3c4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-auth/tasks.md) | ⭐ 19 | `communication` |
| [Proposal](communication/proposal_82b4bcfe/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-scoped-features/proposal.md) | ⭐ 19 | `communication` |
| [Spec](communication/spec_d157ed37/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-auth/specs/user-auth/spec.md) | ⭐ 19 | `communication` |

### Content Creation (18 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [No Heredoc.Instructions](content-creation/357-no-heredocinstructions_a1bdea6d/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/no-heredoc.instructions.md) | ⭐ 60 | `content creation` |
| [Skill](content-creation/049-name-skill_83d53e39/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/skills/docs-writer/SKILL.md) | ⭐ 60 | `content creation` |
| [2023 10 22 Repomap](content-creation/358-2023-10-22-repomap_89a9b7a9/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2023-10-22-repomap.md) | 🔥 40.5k | `content creation` |
| [2024 05 02 Browser](content-creation/359-2024-05-02-browser_bd617fca/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-02-browser.md) | 🔥 40.5k | `content creation` |
| [2024 05 24 Self Assembly](content-creation/360-2024-05-24-self-assembly_37b50370/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-05-24-self-assembly.md) | 🔥 40.5k | `content creation` |
| [2024 07 01 Sonnet Not Lazy](content-creation/361-2024-07-01-sonnet-not-lazy_630e8d66/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-07-01-sonnet-not-lazy.md) | 🔥 40.5k | `content creation` |
| [2024 11 21 Quantization](content-creation/362-2024-11-21-quantization_7b1bdc60/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-11-21-quantization.md) | 🔥 40.5k | `content creation` |
| [2024 12 03 Qwq](content-creation/363-2024-12-03-qwq_223f2745/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-12-03-qwq.md) | 🔥 40.5k | `content creation` |
| [2025 01 28 Deepseek Down](content-creation/364-2025-01-28-deepseek-down_b46d65b3/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-01-28-deepseek-down.md) | 🔥 40.5k | `content creation` |
| [2025 05 08 Qwen3](content-creation/365-2025-05-08-qwen3_b8c3e75c/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-05-08-qwen3.md) | 🔥 40.5k | `content creation` |
| [Benchmarks 0125](content-creation/366-benchmarks-0125_c20592e5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/benchmarks-0125.md) | 🔥 40.5k | `content creation` |
| [Ctags](content-creation/367-ctags_9dd05206/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/ctags.md) | 🔥 40.5k | `content creation` |
| [Faq](content-creation/368-faq_af3fbd08/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/faq.md) | 🔥 40.5k | `content creation` |
| [Unified Diffs](content-creation/369-unified-diffs_16b627c7/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/unified-diffs.md) | 🔥 40.5k | `content creation` |
| [Index](content-creation/019-index_2e12c12a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/index.md) | 🔥 40.5k | `content creation` |
| [Browser](content-creation/370-browser_27e56b0f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/browser.md) | 🔥 40.5k | `content creation` |
| [Modes](content-creation/371-modes_77c3d350/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/modes.md) | 🔥 40.5k | `content creation` |
| [Tasks](content-creation/tasks_7ac00ad6/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-plugin-data-explorer/tasks.md) | ⭐ 19 | `content creation` |

### Daily Assistant (26 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Development Guide](daily-assistant/development_guide_08eb06aa/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/DEVELOPMENT_GUIDE.md) | ⭐ 19 | `daily assistant` |
| [Agents](daily-assistant/agents_7e46a546/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/AGENTS.md) | ⭐ 19 | `daily assistant` |
| [Proposal](daily-assistant/proposal_907f2a23/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-etf-data-plugins/proposal.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_5cdc3e3a/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-financial-plugins/tasks.md) | ⭐ 19 | `daily assistant` |
| [Proposal](daily-assistant/proposal_48e1f616/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-predefined-plugin-groups/proposal.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_8f87cfac/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/enhance-etf-index-ui/tasks.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_3515ad40/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/fix-portfolio-user-isolation/tasks.md) | ⭐ 19 | `daily assistant` |
| [Design](daily-assistant/design_9331dc5b/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/design.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_0628ed7c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/tasks.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_9f728d7c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/refactor-redis-queue-sync/tasks.md) | ⭐ 19 | `daily assistant` |
| [Design](daily-assistant/design_574ae00a/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-data-scheduler/design.md) | ⭐ 19 | `daily assistant` |
| [Spec](daily-assistant/spec_fd1b54a2/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/specs/data-management/spec.md) | ⭐ 19 | `daily assistant` |
| [Design](daily-assistant/design_7f4154fc/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-02-06-enhance-data-management/design.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_bb0c59e9/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-02-06-enhance-data-management/tasks.md) | ⭐ 19 | `daily assistant` |
| [Spec](daily-assistant/spec_da992562/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-multi-agent-strategy-arena/specs/multi-agent-arena/spec.md) | ⭐ 19 | `daily assistant` |
| [Spec](daily-assistant/spec_d112eb02/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/specs/trade-calendar-service/spec.md) | ⭐ 19 | `daily assistant` |
| [Spec](daily-assistant/spec_f6686104/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/refactor-redis-queue-sync/specs/redis-sync-queue/spec.md) | ⭐ 19 | `daily assistant` |
| [Spec](daily-assistant/spec_0f9eb0f9/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-data-scheduler/specs/unified-data-scheduler/spec.md) | ⭐ 19 | `daily assistant` |
| [Spec](daily-assistant/spec_792c5a0e/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-02-06-enhance-data-management/specs/data-management/spec.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_9efbf991/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-chat-orchestrator-agent/tasks.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_8500eebb/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-stock-screener/tasks.md) | ⭐ 19 | `daily assistant` |
| [Design](daily-assistant/design_04ede0bf/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-news-performance/design.md) | ⭐ 19 | `daily assistant` |
| [Proposal](daily-assistant/proposal_2650e42c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/proposal.md) | ⭐ 19 | `daily assistant` |
| [Proposal](daily-assistant/proposal_c4ee0eed/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-data-scheduler/proposal.md) | ⭐ 19 | `daily assistant` |
| [Tasks](daily-assistant/tasks_57da6f25/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-news-center/tasks.md) | ⭐ 19 | `daily assistant` |
| [Proposal](daily-assistant/proposal_7d888f42/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-02-06-enhance-data-management/proposal.md) | ⭐ 19 | `daily assistant` |

### Data Analysis (21 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Tasks](data-analysis/tasks_ae85bca4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-etf-data-plugins/tasks.md) | ⭐ 19 | `data analysis` |
| [Design](data-analysis/design_8e31eeb0/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-analysis/design.md) | ⭐ 19 | `data analysis` |
| [Design](data-analysis/design_8efe06fb/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-predefined-plugin-groups/design.md) | ⭐ 19 | `data analysis` |
| [Spec](data-analysis/spec_97812b5e/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/specs/financial-report-analysis/spec.md) | ⭐ 19 | `data analysis` |
| [Baseplugin Quick Reference](data-analysis/baseplugin_quick_reference_eae05a91/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/BASEPLUGIN_QUICK_REFERENCE.md) | ⭐ 19 | `data analysis` |
| [Spec](data-analysis/spec_e1c1d262/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-analysis/specs/market-analysis/spec.md) | ⭐ 19 | `data analysis` |
| [Spec](data-analysis/spec_1baaadfa/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-predefined-plugin-groups/specs/predefined-plugin-groups/spec.md) | ⭐ 19 | `data analysis` |
| [Intelligent Strategy System Archive](data-analysis/intelligent-strategy-system-archive_4ac63d50/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/archives/intelligent-strategy-system-archive.md) | ⭐ 19 | `data analysis` |
| [Tasks](data-analysis/tasks_14f1b9d6/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-plugins/tasks.md) | ⭐ 19 | `data analysis` |
| [Proposal](data-analysis/proposal_c3b9fd03/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-screener/proposal.md) | ⭐ 19 | `data analysis` |
| [Tasks](data-analysis/tasks_211feb97/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-screener/tasks.md) | ⭐ 19 | `data analysis` |
| [Proposal](data-analysis/proposal_373d9599/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-strategy-system/proposal.md) | ⭐ 19 | `data analysis` |
| [Proposal](data-analysis/proposal_bdb892c2/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-analysis/proposal.md) | ⭐ 19 | `data analysis` |
| [Tasks](data-analysis/tasks_b092e131/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-analysis/tasks.md) | ⭐ 19 | `data analysis` |
| [Tasks](data-analysis/tasks_0d436ddf/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-system-logs-ui/tasks.md) | ⭐ 19 | `data analysis` |
| [Design](data-analysis/design_310dce7d/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-financial-report-agent/design.md) | ⭐ 19 | `data analysis` |
| [Tasks](data-analysis/tasks_af099e54/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-financial-report-agent/tasks.md) | ⭐ 19 | `data analysis` |
| [Design](data-analysis/design_dd5cfbc5/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-portfolio-management/design.md) | ⭐ 19 | `data analysis` |
| [Spec](data-analysis/spec_367bdaff/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-screener/specs/index-screener/spec.md) | ⭐ 19 | `data analysis` |
| [Spec](data-analysis/spec_95e0a164/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/specs/plugin-dependency-management/spec.md) | ⭐ 19 | `data analysis` |
| [Spec](data-analysis/spec_9c48bfd5/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-financial-report-agent/specs/financial-report-analysis/spec.md) | ⭐ 19 | `data analysis` |

### Development (82 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Scripting](development/2875-scripting_163d09dc/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/scripting.md) | 🔥 40.5k | `development` |
| [Usage](development/2876-usage_8a5d1b8d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage.md) | 🔥 40.5k | `development` |
| [Index](development/468-index_3f9d5a0f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/share/index.md) | 🔥 40.5k | `development` |
| [Adv Model Settings](development/2877-adv-model-settings_45ec2acd/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/adv-model-settings.md) | 🔥 40.5k | `development` |
| [Reasoning](development/2878-reasoning_a05be33d/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/reasoning.md) | 🔥 40.5k | `development` |
| [Contrib](development/2879-contrib_40379351/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/contrib.md) | 🔥 40.5k | `development` |
| [Edit](development/2880-edit_92b57853/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/edit.md) | 🔥 40.5k | `development` |
| [Notes](development/2881-notes_aebac3fc/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/notes.md) | 🔥 40.5k | `development` |
| [Refactor](development/827-refactor_1a6d950f/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/leaderboards/refactor.md) | 🔥 40.5k | `development` |
| [Privacy](development/2882-privacy_a4277601/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/legal/privacy.md) | 🔥 40.5k | `development` |
| [Anthropic](development/2185-anthropic_4e802ab8/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/anthropic.md) | 🔥 40.5k | `development` |
| [Gemini](development/1371-gemini_48668ea1/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/gemini.md) | 🔥 40.5k | `development` |
| [Github](development/2317-github_bc027382/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/github.md) | 🔥 40.5k | `development` |
| [Ollama](development/1338-ollama_22770bc8/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/llms/ollama.md) | 🔥 40.5k | `development` |
| [Edit Formats](development/2883-edit-formats_0d3c5081/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/edit-formats.md) | 🔥 40.5k | `development` |
| [Auto Accept Architect](development/2884-auto-accept-architect_4d97fe45/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/auto-accept-architect.md) | 🔥 40.5k | `development` |
| [Dont Drop Original Read Files](development/2885-dont-drop-original-read-files_e4e2aebf/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/dont-drop-original-read-files.md) | 🔥 40.5k | `development` |
| [Model Accepts Settings](development/2886-model-accepts-settings_4bf52d6e/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/model-accepts-settings.md) | 🔥 40.5k | `development` |
| [Tree Sitter Language Pack](development/2887-tree-sitter-language-pack_d7fc4c53/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/recordings/tree-sitter-language-pack.md) | 🔥 40.5k | `development` |
| [Support](development/2888-support_2a9ff792/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/support.md) | 🔥 40.5k | `development` |
| [Token Limits](development/2889-token-limits_f51ac904/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/troubleshooting/token-limits.md) | 🔥 40.5k | `development` |
| [Copypaste](development/2890-copypaste_845f243a/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/copypaste.md) | 🔥 40.5k | `development` |
| [Tips](development/2891-tips_b26427c6/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/tips.md) | 🔥 40.5k | `development` |
| [Tutorials](development/1100-tutorials_64a073b5/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/tutorials.md) | 🔥 40.5k | `development` |
| [Voice](development/2892-voice_06d5fba4/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/voice.md) | 🔥 40.5k | `development` |
| [Watch](development/2893-watch_395175e2/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/watch.md) | 🔥 40.5k | `development` |
| [Adapters](development/2894-adapters_beef3be8/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/ADAPTERS.md) | ⭐ 43 | `development` |
| [Backend Requirements](development/2895-backend_requirements_51bd5367/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/BACKEND_REQUIREMENTS.md) | ⭐ 43 | `development` |
| [Database Setup](development/2896-database_setup_061e7377/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/DATABASE_SETUP.md) | ⭐ 43 | `development` |
| [Framework Support](development/2897-framework_support_cd166e68/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FRAMEWORK_SUPPORT.md) | ⭐ 43 | `development` |
| [Generic Guardrail Api](development/2871-generic_guardrail_api_60da4bb2/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/adding_provider/generic_guardrail_api.md) | 🔥 35.8k | `development` |
| [Project](development/project_852f50da/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/project.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_c67b1b6e/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-chat-orchestrator-agent/proposal.md) | ⭐ 19 | `development` |
| [Design](development/design_497a09d9/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-etf-data-plugins/design.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_8df4efd1/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-financial-statement-plugins/proposal.md) | ⭐ 19 | `development` |
| [Design](development/design_6d97a46c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-strategy-system/design.md) | ⭐ 19 | `development` |
| [Tasks](development/tasks_da032476/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-strategy-system/tasks.md) | ⭐ 19 | `development` |
| [Design](development/design_02377a64/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-plugin-data-explorer/design.md) | ⭐ 19 | `development` |
| [Idx Factor Pro](development/idx_factor_pro_fbfdf148/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_idx_factor_pro/idx_factor_pro.md) | ⭐ 19 | `development` |
| [Tushare Ths Index](development/tushare_ths_index_7cf1bc79/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_ths_index/tushare_ths_index.md) | ⭐ 19 | `development` |
| [Skill](development/name-skill_d4daeece/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_StartTouch-AND-CBaseEntity_Touch-CBaseEntity_EndTouch/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_937d2e20/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBasePlayerPawn_RemovePlayerItem/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_5d9c767a/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseTrigger_StartTouch-AND-CBaseTrigger_EndTouch/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_390ec490/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayer_ItemServices_GiveNamedItem/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_3f900628/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CTriggerPush_Touch/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_897bdac8/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-FireBullets-AND-TraceAttack-AND-CTakeDamageInfo/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_07c27589/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-IGameSystem_LoopPostInitAllSystems_pEventDispatcher-AND-IGameSystem_LoopDestroyAllSystems_s_GameSystems/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_df3e6ebd/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/generate-signature-for-globalvar/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_1b2e0cfd/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/generate-signature-for-vfuncoffset/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_82e1e8ce/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/write-func-as-yaml/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_b1d9bd95/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/write-vfunc-as-yaml/SKILL.md) | ⭐ 16 | `development` |
| [Proposal](development/proposal_c8958610/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-custom-ai-workflow/proposal.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_83fcace3/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-extended-index-plugins/proposal.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_2d399f90/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-financial-plugins/proposal.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_72176259/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-plugins/proposal.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_2a2f99f3/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-research-data-ui/proposal.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_a492a65f/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/refactor-redis-queue-sync/proposal.md) | ⭐ 19 | `development` |
| [Proposal](development/proposal_ce8fec95/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-financial-report-agent/proposal.md) | ⭐ 19 | `development` |
| [Daily](development/daily_0f172f6a/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_daily/daily.md) | ⭐ 19 | `development` |
| [Daily Basic](development/daily_basic_aaf4d69d/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_daily_basic/daily_basic.md) | ⭐ 19 | `development` |
| [Etf Basic](development/etf_basic_7906c298/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_etf_basic/etf_basic.md) | ⭐ 19 | `development` |
| [Fund Adj](development/fund_adj_51d40d45/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_etf_fund_adj/fund_adj.md) | ⭐ 19 | `development` |
| [Etf Index](development/etf_index_df6df2d4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_etf_index/etf_index.md) | ⭐ 19 | `development` |
| [Stk Mins](development/stk_mins_23b772e2/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_etf_stk_mins/stk_mins.md) | ⭐ 19 | `development` |
| [Finace Indicator](development/finace_indicator_af395750/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_finace_indicator/finace_indicator.md) | ⭐ 19 | `development` |
| [Index Weight](development/index_weight_99125714/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_index_weight/index_weight.md) | ⭐ 19 | `development` |
| [Spec](development/spec_df6acc73/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-docker-compose-deployment/specs/redis-cache/spec.md) | ⭐ 19 | `development` |
| [Spec](development/spec_fa3f4473/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-etf-data-plugins/specs/etf-data-plugins/spec.md) | ⭐ 19 | `development` |
| [Spec](development/spec_46f83991/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-extended-index-plugins/specs/extended-index-data-sources/spec.md) | ⭐ 19 | `development` |
| [Spec](development/spec_a4e77e5c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-financial-statement-plugins/specs/financial-statement-plugins/spec.md) | ⭐ 19 | `development` |
| [Spec](development/spec_de0d2d7c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-financial-plugins/specs/hk-financial-plugins/spec.md) | ⭐ 19 | `development` |
| [Spec](development/spec_5a24b642/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-top-list-tracking/specs/top-list-tracking/spec.md) | ⭐ 19 | `development` |
| [Spec](development/spec_6baf46a2/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/enhance-etf-index-ui/specs/etf-filter-enhancement/spec.md) | ⭐ 19 | `development` |
| [Spec](development/spec_0b511433/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/enhance-etf-index-ui/specs/index-filter-enhancement/spec.md) | ⭐ 19 | `development` |
| [Skill](development/name-skill_8415e994/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CBaseEntity_SetMoveType/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_f4f0b363/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSGameRules_TerminateRound-AND-CEntityInstance_AcceptInput/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_305b5ac0/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CCSPlayer_WeaponServices_Weapon_GetSlot/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_6c652dda/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-CGameEntitySystem_FindEntityByClassName/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_ed04afdc/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-GiveNamedItem-AND-UTIL_CreateEntityByName-AND-CBaseEntity_DispatchSpawn/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_9301d740/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-IGameSystem_InitAllSystems-AND-IGameSystem_InitAllSystems_pFirst/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_abbe367f/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/find-LegacyGameEventListener/SKILL.md) | ⭐ 16 | `development` |
| [Skill](development/name-skill_345b149d/) | [hzqst/CS2_VibeSignatures](https://raw.githubusercontent.com/hzqst/CS2_VibeSignatures/main/.claude/skills/get-func-from-yaml/SKILL.md) | ⭐ 16 | `development` |

### Development/Devops (8 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Configs](development/devops/034-configs_40652549/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/configs.md) | 🔥 35.8k | `development` |
| [Prod](development/devops/041-prod_27e95b4b/) | [BerriAI/litellm](https://raw.githubusercontent.com/BerriAI/litellm/main/docs/my-website/docs/proxy/prod.md) | 🔥 35.8k | `development` |
| [Arena Api](development/devops/arena_api_3dbfee82/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/ARENA_API.md) | ⭐ 19 | `development` |
| [Design](development/devops/design_d5e07314/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-docker-compose-deployment/design.md) | ⭐ 19 | `development` |
| [Tasks](development/devops/tasks_ffb87591/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-startup-workers/tasks.md) | ⭐ 19 | `development` |
| [Proposal](development/devops/proposal_bb5b8093/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-docker-compose-deployment/proposal.md) | ⭐ 19 | `development` |
| [Tasks](development/devops/tasks_8edf371f/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-docker-compose-deployment/tasks.md) | ⭐ 19 | `development` |
| [Spec](development/devops/spec_22b873ac/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-docker-compose-deployment/specs/docker-deployment/spec.md) | ⭐ 19 | `development` |

### Development/Testing (2 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [2024 04 09 Gpt 4 Turbo](development/testing/082-2024-04-09-gpt-4-turbo_9e2cf0f0/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2024-04-09-gpt-4-turbo.md) | 🔥 40.5k | `development` |
| [Design](development/testing/design_1f1bbfe1/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-stock-screener/design.md) | ⭐ 19 | `development` |

### Development/Tools (17 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Docs.Instructions](development/tools/229-docsinstructions_c7327425/) | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/.github/instructions/docs.instructions.md) | ⭐ 60 | `development` |
| [2025 01 15 Uv](development/tools/319-2025-01-15-uv_954a0942/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/_posts/2025-01-15-uv.md) | 🔥 40.5k | `development` |
| [Aider Conf](development/tools/320-aider_conf_49df9413/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/aider_conf.md) | 🔥 40.5k | `development` |
| [Dotenv](development/tools/321-dotenv_0088c147/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/dotenv.md) | 🔥 40.5k | `development` |
| [Options](development/tools/322-options_6bd20c61/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/config/options.md) | 🔥 40.5k | `development` |
| [Commands](development/tools/252-commands_57549608/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/usage/commands.md) | 🔥 40.5k | `development` |
| [Agents](development/tools/015-agents_1bbd0b4d/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/AGENTS.md) | ⭐ 43 | `development` |
| [Evaluation Metrics](development/tools/323-evaluation_metrics_0d6d74bb/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/EVALUATION_METRICS.md) | ⭐ 43 | `development` |
| [Faq](development/tools/324-faq_1839ca4b/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/FAQ.md) | ⭐ 43 | `development` |
| [Mcp Contracts](development/tools/325-mcp_contracts_53288f40/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/MCP_CONTRACTS.md) | ⭐ 43 | `development` |
| [Tool Categories](development/tools/326-tool_categories_6c5054ba/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TOOL_CATEGORIES.md) | ⭐ 43 | `development` |
| [Troubleshooting](development/tools/205-troubleshooting_271bcac1/) | [hidai25/eval-view](https://raw.githubusercontent.com/hidai25/eval-view/main/docs/TROUBLESHOOTING.md) | ⭐ 43 | `development` |
| [Skill](development/tools/name-skill_75938eda/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/skills/tushare-plugin-builder/SKILL.md) | ⭐ 19 | `development` |
| [Design](development/tools/design_8a7bb14c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-custom-ai-workflow/design.md) | ⭐ 19 | `development` |
| [Design](development/tools/design_504151b9/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-financial-statement-plugins/design.md) | ⭐ 19 | `development` |
| [Proposal](development/tools/proposal_a6b56f16/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-stock-screener/proposal.md) | ⭐ 19 | `development` |
| [Tasks](development/tools/tasks_686ee296/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-extended-index-plugins/tasks.md) | ⭐ 19 | `development` |

### Investment (76 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Plugin Quick Start](investment/plugin_quick_start_db0b7f1c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/PLUGIN_QUICK_START.md) | ⭐ 19 | `investment` |
| [Ai Stock Platform Design](investment/ai_stock_platform_design_77dcc7a7/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/AI_STOCK_PLATFORM_DESIGN.md) | ⭐ 19 | `investment` |
| [Cli Guide](investment/cli_guide_07b09d25/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/CLI_GUIDE.md) | ⭐ 19 | `investment` |
| [Database Setup](investment/database_setup_419ba41a/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/DATABASE_SETUP.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_426a57b4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/spec.md) | ⭐ 19 | `investment` |
| [Ods Stock Basic](investment/ods_stock_basic_7ecb27d3/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/tables/ods_stock_basic.md) | ⭐ 19 | `investment` |
| [Design](investment/design_8b69661b/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-chat-orchestrator-agent/design.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_b02292e4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-financial-statement-plugins/tasks.md) | ⭐ 19 | `investment` |
| [Design](investment/design_8197a4dd/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-financial-plugins/design.md) | ⭐ 19 | `investment` |
| [Design](investment/design_0c3c9097/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-stock-support/design.md) | ⭐ 19 | `investment` |
| [Design](investment/design_fafbfdfa/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-plugins/design.md) | ⭐ 19 | `investment` |
| [Design](investment/design_eb8a3615/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-screener/design.md) | ⭐ 19 | `investment` |
| [Design](investment/design_b591d76c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-overview-ui/design.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_1bb615e3/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-overview-ui/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_66480e91/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-overview-ui/tasks.md) | ⭐ 19 | `investment` |
| [Design](investment/design_8759c38a/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-multi-agent-strategy-arena/design.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_a2eaae60/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-multi-agent-strategy-arena/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_4c951f08/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-multi-agent-strategy-arena/tasks.md) | ⭐ 19 | `investment` |
| [Design](investment/design_7087e0f6/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-analyst-agent/design.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_3e835521/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-analyst-agent/tasks.md) | ⭐ 19 | `investment` |
| [Design](investment/design_30ae1ca4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-frontend-ui/design.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_e3873265/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-plugin-data-explorer/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_a7f63c1e/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-predefined-plugin-groups/tasks.md) | ⭐ 19 | `investment` |
| [Design](investment/design_e172ac32/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-research-data-ui/design.md) | ⭐ 19 | `investment` |
| [Design](investment/design_4c153dd7/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-top-list-tracking/design.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_eda05418/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-top-list-tracking/proposal.md) | ⭐ 19 | `investment` |
| [Design](investment/design_a0120fcf/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-auth/design.md) | ⭐ 19 | `investment` |
| [Design](investment/design_c6d3697c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/enhance-etf-index-ui/design.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_151801e7/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/enhance-etf-index-ui/proposal.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_750f6e67/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-news-performance/proposal.md) | ⭐ 19 | `investment` |
| [Design](investment/design_67bca546/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/refactor-redis-queue-sync/design.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_b1b18f60/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-data-scheduler/tasks.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_505c9af2/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-index-plugins/specs/index-data-sources/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_7b190a0f/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-strategy-system/specs/intelligent-strategy-system/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_cf8674d0/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-overview-ui/specs/market-overview/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_7de31a7d/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-multi-agent-strategy-arena/specs/agent-discussion/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_f3c770ab/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-analyst-agent/specs/news-analysis/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_50cf2cf9/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/enhance-etf-index-ui/specs/backend-api/spec.md) | ⭐ 19 | `investment` |
| [Quick Load Guide](investment/quick_load_guide_e87670ae/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/QUICK_LOAD_GUIDE.md) | ⭐ 19 | `investment` |
| [Toplist Features](investment/toplist_features_19892a6c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/TOPLIST_FEATURES.md) | ⭐ 19 | `investment` |
| [Test Report](investment/test_report_3323f76f/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/TEST_REPORT.md) | ⭐ 19 | `investment` |
| [Portfolio Partition Strategy](investment/portfolio_partition_strategy_ce3c21e4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/tables/portfolio_partition_strategy.md) | ⭐ 19 | `investment` |
| [Design](investment/design_ac1db0fd/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-extended-index-plugins/design.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_3b6c3237/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-stock-support/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_57c1aeb2/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-stock-support/tasks.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_c91d0465/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-analyst-agent/proposal.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_d34d7cea/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-frontend-ui/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_e8e9ed15/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-frontend-ui/tasks.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_eb849698/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-research-data-ui/tasks.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_711fd4e0/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-system-logs-ui/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_949f724a/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-top-list-tracking/tasks.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_a9969f64/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/fix-portfolio-user-isolation/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_40b94f39/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-news-performance/tasks.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_15d862f2/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-startup-workers/proposal.md) | ⭐ 19 | `investment` |
| [Plugin Conventions](investment/plugin_conventions_fce59ddc/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/skills/tushare-plugin-builder/references/plugin_conventions.md) | ⭐ 19 | `investment` |
| [Proposal](investment/proposal_42f169f4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-portfolio-management/proposal.md) | ⭐ 19 | `investment` |
| [Tasks](investment/tasks_3a828af4/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-portfolio-management/tasks.md) | ⭐ 19 | `investment` |
| [Tushare Hk Adjfactor](investment/tushare_hk_adjfactor_fbdb8643/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_hk_adjfactor/tushare_hk_adjfactor.md) | ⭐ 19 | `investment` |
| [Tushare Hk Basic](investment/tushare_hk_basic_1ac82266/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_hk_basic/tushare_hk_basic.md) | ⭐ 19 | `investment` |
| [Tushare Hk Daily](investment/tushare_hk_daily_e18a0531/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_hk_daily/tushare_hk_daily.md) | ⭐ 19 | `investment` |
| [Tushare Hk Daily Adj](investment/tushare_hk_daily_adj_1f2df216/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_hk_daily_adj/tushare_hk_daily_adj.md) | ⭐ 19 | `investment` |
| [Tushare Hk Tradecal](investment/tushare_hk_tradecal_6ed41bd5/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_hk_tradecal/tushare_hk_tradecal.md) | ⭐ 19 | `investment` |
| [Stock Basic](investment/stock_basic_292fd149/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/src/stock_datasource/plugins/tushare_stock_basic/stock_basic.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_0b4ff484/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-chat-orchestrator-agent/specs/chat-orchestration/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_2e26449c/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-custom-ai-workflow/specs/custom-ai-workflow/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_80d4a6d6/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-stock-support/specs/stock-detail/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_7f90cfac/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-intelligent-stock-screener/specs/stock-screener/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_c71776f5/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-market-overview-ui/specs/etf-display/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_535e1583/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-multi-agent-strategy-arena/specs/strategy-competition/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_fb297c9a/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-research-data-ui/specs/research-data-ui/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_fd9843b3/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-top-list-tracking/specs/frontend-integration/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_833c72fc/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-user-scoped-features/specs/user-data-isolation/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_d4735bee/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/specs/schedule-management/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_7f92c7e0/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-portfolio-management/specs/agent-system/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_734c3985/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-portfolio-management/specs/api-endpoints/spec.md) | ⭐ 19 | `investment` |
| [Spec](investment/spec_ee96166e/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/archive/2026-01-11-add-portfolio-management/specs/portfolio-management/spec.md) | ⭐ 19 | `investment` |

### Other (4 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Agents](other/agents_98937111/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/AGENTS.md) | ⭐ 19 | `other` |
| [Spec](other/spec_df36458b/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-news-frontend-ui/specs/news-frontend/spec.md) | ⭐ 19 | `other` |
| [Spec](other/spec_f9ae8b65/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/enhance-etf-index-ui/specs/frontend-navigation/spec.md) | ⭐ 19 | `other` |
| [Spec](other/spec_64e6f326/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-startup-workers/specs/data-management/spec.md) | ⭐ 19 | `other` |

### Research (5 skills)

| Skill | Source | Popularity | Tags |
|-------|--------|------------|------|
| [Infinite Output](research/259-infinite-output_a05bef43/) | [Aider-AI/aider](https://raw.githubusercontent.com/Aider-AI/aider/main/aider/website/docs/more/infinite-output.md) | 🔥 40.5k | `research` |
| [Design](research/design_596f4d17/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-news-center/design.md) | ⭐ 19 | `research` |
| [Proposal](research/proposal_e6644861/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-news-center/proposal.md) | ⭐ 19 | `research` |
| [Spec](research/spec_715c123e/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-system-logs-ui/specs/system-logs-viewer/spec.md) | ⭐ 19 | `research` |
| [Spec](research/spec_5c986edb/) | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-news-center/specs/news-center/spec.md) | ⭐ 19 | `research` |

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

*Last updated: 2026-02-12 18:22:28 UTC*
*Automatically maintained by SkillFlow*
