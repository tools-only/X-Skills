# SF-Skills Orchestration Architecture

> How skill recommendations and workflow orchestration currently work

---

## Overview

The orchestration system has **two entry points** that suggest skills at different moments:

```
┌─────────────────────────────────────────────────────────────────┐
│                                                                 │
│   USER TYPES PROMPT                                             │
│         │                                                       │
│         ▼                                                       │
│   ┌─────────────────┐                                          │
│   │ UserPromptSubmit│  ◄── Hook #1: BEFORE Claude responds     │
│   │ Hook            │      "What skill should handle this?"    │
│   └────────┬────────┘                                          │
│            │                                                    │
│            ▼                                                    │
│   ┌─────────────────┐                                          │
│   │ CLAUDE WORKS    │                                          │
│   │ (writes files)  │                                          │
│   └────────┬────────┘                                          │
│            │                                                    │
│            ▼                                                    │
│   ┌─────────────────┐                                          │
│   │ PostToolUse     │  ◄── Hook #2: AFTER file is written      │
│   │ Hook            │      "What skill should come next?"      │
│   └─────────────────┘                                          │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Data Flow

### Single Source of Truth

All orchestration data lives in one file:

```
shared/hooks/skills-registry.json
```

This file contains:
- **16 skill definitions** with keywords, patterns, triggers
- **5 workflow chains** (full_feature, integration, agentforce, triangle, troubleshooting)
- **Confidence levels** (REQUIRED=3, RECOMMENDED=2, OPTIONAL=1)

---

## Hook #1: UserPromptSubmit

**File:** `shared/hooks/skill-activation-prompt.py`

**When:** Before Claude sees the user's message

**Purpose:** Suggest which skill(s) to use based on the prompt

### Scoring Algorithm

```python
score = 0

# 1. Keyword matches (+2 each, max 3 counted)
for keyword in skill["keywords"]:
    if keyword in prompt:
        score += 2

# 2. Intent pattern match (+3)
for pattern in skill["intentPatterns"]:
    if regex_match(pattern, prompt):
        score += 3
        break

# 3. File pattern match (+2)
for pattern in skill["filePatterns"]:
    if any_active_file_matches(pattern):
        score += 2
        break
```

### Confidence Assignment

```
score >= 7  →  REQUIRED    (★★★)
score >= 4  →  RECOMMENDED (★★)
score >= 2  →  OPTIONAL    (★)
score <  2  →  No suggestion
```

### Chain Detection

Checks if prompt contains trigger phrases:

```python
chains = {
    "full_feature": ["build feature", "complete feature", "end to end"],
    "agentforce": ["build agent", "create agentforce", "copilot"],
    "integration": ["integrate api", "external integration"],
    "triangle": ["apex lwc flow", "triangle"],
    "troubleshooting": ["debug", "troubleshoot", "fix failing"]
}

for chain_name, triggers in chains.items():
    for trigger in triggers:
        if trigger in prompt.lower():
            return chain_name  # First match wins
```

### Example Output

```
═══════════════════════════════════════════════════════
💡 SKILL SUGGESTIONS (based on your request)
═══════════════════════════════════════════════════════

📋 DETECTED WORKFLOW: agentforce
   Order: sf-metadata → sf-apex → sf-flow → sf-deploy → sf-ai-agentscript
   ⭐ START WITH: /sf-metadata

⭐⭐⭐ /sf-apex - REQUIRED
   └─ Apex code development with validation
⭐⭐ /sf-flow - RECOMMENDED
   └─ Flow Builder automation

───────────────────────────────────────────────────────
```

---

## Hook #2: PostToolUse

**File:** `shared/hooks/suggest-related-skills.py`

**When:** After Claude writes or edits a file

**Purpose:** Suggest next steps based on what was just created

### Detection Logic

```python
# 1. Which skill owns this file?
def detect_skill_from_file(file_path):
    for skill, config in registry["skills"].items():
        for pattern in config["filePatterns"]:
            if regex_match(pattern, file_path):
                return skill
    return None

# 2. What patterns are in the file content?
def detect_content_triggers(content, skill_config):
    suggestions = []
    for pattern, suggestion in skill_config["contentTriggers"].items():
        if regex_search(pattern, content):
            suggestions.append(suggestion)
    return suggestions
```

### Orchestration Relationships

Each skill defines three relationship types:

```json
{
  "sf-apex": {
    "orchestration": {
      "prerequisites": [
        {
          "skill": "sf-metadata",
          "condition": "custom_field_reference",
          "message": "Ensure custom fields exist before Apex references them",
          "confidence": 2
        }
      ],
      "next_steps": [
        {
          "skill": "sf-testing",
          "condition": "always",
          "message": "Run tests to validate Apex code",
          "confidence": 3
        }
      ],
      "commonly_with": [
        {
          "skill": "sf-flow",
          "trigger": "@InvocableMethod",
          "message": "Create Flow to call this @InvocableMethod",
          "confidence": 2
        }
      ]
    }
  }
}
```

### Context Persistence

Saves state to track workflow progress:

```
/tmp/sf-skills-context.json
```

```json
{
  "timestamp": "2026-01-14T10:30:00",
  "last_skill": "sf-apex",
  "last_file": "/path/to/MyClass.cls",
  "detected_triggers": ["@InvocableMethod"],
  "suggested_next": ["sf-flow", "sf-testing"]
}
```

### Example Output

```
═════════════════════════════════════════════════════════
🔗 SKILL SUGGESTIONS (working with sf-apex)
═════════════════════════════════════════════════════════

📋 DETECTED WORKFLOW: triangle
   Step 2 of 4: sf-apex
   Next: sf-lwc → sf-flow

⚠️ PREREQUISITE: /sf-metadata *** REQUIRED
   └─ Ensure custom fields exist before Apex references them
➡️ NEXT STEP: /sf-testing *** REQUIRED
   └─ Run tests to validate Apex code
🔄 RELATED: /sf-flow ** RECOMMENDED
   └─ Create Flow to call this @InvocableMethod

─────────────────────────────────────────────────────────
```

---

## Content Triggers

Pattern detection in file content that suggests related skills:

| Source Skill | Pattern | Suggests | Why |
|--------------|---------|----------|-----|
| sf-apex | `@InvocableMethod` | sf-flow | Flow can call this method |
| sf-apex | `@AuraEnabled` | sf-lwc | LWC needs this controller |
| sf-apex | `@IsTest` | sf-testing | Test class needs execution |
| sf-flow | `actionCalls.*apex` | sf-apex | Flow calls Apex action |
| sf-flow | `ComponentInstance` | sf-lwc | Screen flow uses LWC |
| sf-lwc | `import.*@salesforce/apex` | sf-apex | LWC imports Apex |
| sf-ai-agentscript | `flow://` | sf-flow | Agent uses Flow target |

---

## Workflow Chains

Pre-defined sequences for common development patterns:

### full_feature
```
sf-metadata → sf-apex → sf-flow → sf-lwc → sf-deploy → sf-testing → sf-data
```
Triggered by: "build feature", "complete feature", "full implementation"

### agentforce
```
sf-metadata → sf-apex → sf-flow → sf-deploy → sf-ai-agentscript → sf-ai-agentforce-testing
```
Triggered by: "build agent", "create agentforce", "copilot"

### integration
```
sf-connected-apps → sf-integration → sf-flow → sf-deploy
```
Triggered by: "integrate api", "external integration"

### triangle
```
sf-apex → sf-lwc → sf-flow → sf-deploy
```
Triggered by: "apex lwc flow", "triangle"

### troubleshooting
```
sf-testing → sf-debug → sf-apex → sf-deploy → sf-testing
```
Triggered by: "debug", "troubleshoot", "fix failing"

---

## File Structure

```
sf-skills/
├── .claude/
│   └── hooks.json                    # Global hook config (UserPromptSubmit)
│
├── shared/hooks/
│   ├── skills-registry.json          # ★ SINGLE SOURCE OF TRUTH ★
│   ├── skill-activation-prompt.py    # Hook #1 implementation
│   ├── suggest-related-skills.py     # Hook #2 implementation
│   └── README.md                     # Hook documentation
│
├── sf-apex/
│   ├── SKILL.md                      # Skill instructions
│   ├── hooks/
│   │   └── hooks.json                # PostToolUse hook config
│   └── ...
│
└── [15 more sf-* skills with same structure]
```

---

## Key Observations

### What Works Well
1. **Centralized registry** - One file to understand all relationships
2. **Two-stage suggestions** - Before (intent) and after (content)
3. **Content triggers** - Code patterns drive next-skill suggestions
4. **Chain awareness** - Shows progress through workflows

### Current Limitations

1. **Exact matching** - "create apex" works, "make a class" doesn't
2. **Fixed thresholds** - score≥7 always means REQUIRED, regardless of context
3. **Volatile state** - `/tmp/` context lost on reboot
4. **Manual sync** - Adding a skill requires editing skills-registry.json separately from SKILL.md

---

## How to Test

### Test Hook #1 (UserPromptSubmit)
```bash
echo '{"prompt":"I need to create an apex trigger for Account"}' | \
  python3 shared/hooks/skill-activation-prompt.py
```

### Test Hook #2 (PostToolUse)
```bash
echo '{"tool":"Write","toolInput":{"file_path":"test.cls"}}' | \
  python3 shared/hooks/suggest-related-skills.py sf-apex
```

### View Current Context
```bash
cat /tmp/sf-skills-context.json | jq .
```

---

## Registry Schema Reference

```json
{
  "version": "3.0.0",
  "skills": {
    "<skill-name>": {
      "keywords": ["keyword1", "keyword2"],
      "intentPatterns": ["regex1", "regex2"],
      "filePatterns": ["regex1", "regex2"],
      "contentTriggers": {
        "pattern": { "suggests": "skill", "message": "reason" }
      },
      "priority": "high|medium|low",
      "description": "...",
      "orchestration": {
        "prerequisites": [{ "skill": "...", "condition": "...", "message": "...", "confidence": 1-3 }],
        "next_steps": [{ "skill": "...", "condition": "...", "message": "...", "confidence": 1-3 }],
        "commonly_with": [{ "skill": "...", "trigger": "...", "message": "...", "confidence": 1-3 }]
      }
    }
  },
  "chains": {
    "<chain-name>": {
      "description": "...",
      "trigger_phrases": ["phrase1", "phrase2"],
      "order": ["skill1", "skill2", "skill3"]
    }
  },
  "confidence_levels": {
    "3": { "label": "REQUIRED", "icon": "***" },
    "2": { "label": "RECOMMENDED", "icon": "**" },
    "1": { "label": "OPTIONAL", "icon": "*" }
  }
}
```
