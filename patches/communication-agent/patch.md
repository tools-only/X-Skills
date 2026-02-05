# Communication Agent

**Version:** 1.0.0
**Generated:** 2026-02-05T03:01:56.250992
**Skills:** 40

## Description

Email and messaging skills

## Use Case

Agents handling email, messaging, and communication tasks

## Contents

This patch includes **40 skills** organized into **1 categories**.


### Communication (40 skills)

| Skill | Subcategory | Tags | Status |
|-------|-------------|------|--------|
| [Chat Environment](../communication/023-claude_810b05a9) | - | `data analysis` | **Required** |
| [Setup](../communication/1297-setup_55e51e1e) | - | `communication` | **Required** |
| [Maintenance](../communication/071-postgresql_791c00df) | - | `communication` | **Required** |
| [Workspaces](../communication/124-workspaces_b3e4d3a9) | - | `communication` | **Required** |
| [Skill](../communication/003-name-skill_7a1bf3b5) | - | `communication` | **Required** |
| [Skill](../communication/003-name-skill_881b1d20) | - | `communication` | **Required** |
| [Active Record](../communication/002-ruby_b64643ce) | - | `communication` | **Required** |
| [Skill](../communication/003-name-skill_07c008b4) | - | `communication` | **Required** |
| [Shadcn Accessibility](../communication/106-aria_a22cd4cc) | - | `communication` | **Required** |
| [Riverpod State](../communication/099-dart_ba711b61) | - | `communication` | **Required** |
| [Security Checklists](../communication/102-strong_3be0bd62) | - | `communication` | **Required** |
| [Email Password Auth](../communication/051-email_a96376f5) | - | `communication` | **Required** |
| [Alerting Rules](../communication/005-yaml_0b8d94c8) | - | `communication` | **Required** |
| [Optimization](../communication/078-sql_dd197a02) | - | `communication` | **Required** |
| [Skill](../communication/003-name-skill_a40ef0cd) | - | `communication` | **Required** |
| [Userguide](../communication/1041-userguide_a2542e04) | - | `communication` | **Required** |
| [Standup Notes](../communication/109-ticket_afc86e1b) | - | `communication` | **Required** |
| [Marketplace Schema](../communication/072-marketplace_5a9b0859) | - | `communication` | **Required** |
| [Llms Txt](../communication/070-purpose_0bb1b665) | - | `communication` | **Required** |
| [Query Patterns](../communication/090-sql_19c3b2e8) | - | `communication` | **Required** |
| [Employee Onboarding Offboarding Process](../communication/1058-employee_onboarding_offboarding_process_b697ed21) | - | `communication` | **Required** |
| [CLI Guide](../communication/026-tier_f807443d) | - | `communication` | **Required** |
| [Dtos Validation](../communication/049-typescript_e4282304) | - | `communication` | **Required** |
| [Checklist](../communication/1137-checklist_7165eb4b) | - | `communication` | **Required** |
| [Reference Set 02 Idor Prevention](../communication/445-reference-set-02-idor-prevention_802750c5) | - | `communication` | **Required** |
| [Notifications](../communication/612-notifications_d933f28e) | - | `communication` | **Required** |
| [Orchestration Conditional](../communication/079-conditional_bbec71f3) | - | `communication` | **Required** |
| [Conversation Framework](../communication/038-three_7bde1cef) | - | `daily assistant` | **Required** |
| [Brand Guidelines](../communication/021-framework_6d6e9dc4) | - | `content creation` | **Required** |
| [Skill](../communication/107-ai-sdk-5_bb521bed) | - | `content creation` | **Required** |
| [Endpoints Routing](../communication/052-python_352d08a0) | - | `communication` | **Required** |
| [Migration Checklist](../communication/257-migration-checklist_2f4b8456) | - | `communication` | **Required** |
| [Skill](../communication/003-name-skill_6ec2bf56) | - | `communication` | **Required** |
| [Reference Set 02 Idor Prevention](../communication/096-threat_8d00c0dc) | - | `communication` | **Required** |
| [Hidden Conversations Sidebar Click Fix](../communication/1054-hidden_conversations_sidebar_click_fix_788ea95e) | - | `communication` | **Required** |
| [Authentication](../communication/016-python_1c8cb677) | - | `communication` | **Required** |
| [Faqs](../communication/1038-faqs_51096eab) | - | `communication` | **Required** |
| [Skill](../communication/003-name-skill_9f6e9d65) | - | `communication` | **Required** |
| [Components](../communication/031-vue_91388a4b) | - | `communication` | **Required** |
| [Component Patterns](../communication/1279-component-patterns_81c1b135) | - | `communication` | **Required** |


## Installation

To install this patch to your Claude Code skills directory:

```bash
python -m src.patch_installer install communication-agent
```

Or using the skillflow CLI (if available):

```bash
skillflow patch install communication-agent
```

## Dependencies

None

## Manifest

```json
{
  "spec": {
    "id": "communication-agent",
    "name": "Communication Agent",
    "description": "Email and messaging skills",
    "use_case": "Agents handling email, messaging, and communication tasks",
    "categories": [
      "communication"
    ],
    "subcategories": [],
    "tags": [],
    "exclude": [
      "*test*",
      "*example*",
      "*demo*",
      "*template*"
    ],
    "min_stars": null,
    "max_skills": 40,
    "required_skills": [],
    "optional_skills": [],
    "dependencies": [],
    "version": "1.0.0"
  },
  "skills": [
    {
      "local_path": "communication/023-claude_810b05a9",
      "display_name": "Chat Environment",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "data analysis"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/1297-setup_55e51e1e",
      "display_name": "Setup",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/071-postgresql_791c00df",
      "display_name": "Maintenance",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/124-workspaces_b3e4d3a9",
      "display_name": "Workspaces",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/003-name-skill_7a1bf3b5",
      "display_name": "Skill",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/003-name-skill_881b1d20",
      "display_name": "Skill",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/002-ruby_b64643ce",
      "display_name": "Active Record",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/003-name-skill_07c008b4",
      "display_name": "Skill",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/106-aria_a22cd4cc",
      "display_name": "Shadcn Accessibility",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/099-dart_ba711b61",
      "display_name": "Riverpod State",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/102-strong_3be0bd62",
      "display_name": "Security Checklists",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/051-email_a96376f5",
      "display_name": "Email Password Auth",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/005-yaml_0b8d94c8",
      "display_name": "Alerting Rules",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/078-sql_dd197a02",
      "display_name": "Optimization",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/003-name-skill_a40ef0cd",
      "display_name": "Skill",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/1041-userguide_a2542e04",
      "display_name": "Userguide",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/109-ticket_afc86e1b",
      "display_name": "Standup Notes",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/072-marketplace_5a9b0859",
      "display_name": "Marketplace Schema",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/070-purpose_0bb1b665",
      "display_name": "Llms Txt",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/090-sql_19c3b2e8",
      "display_name": "Query Patterns",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/1058-employee_onboarding_offboarding_process_b697ed21",
      "display_name": "Employee Onboarding Offboarding Process",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/026-tier_f807443d",
      "display_name": "CLI Guide",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/049-typescript_e4282304",
      "display_name": "Dtos Validation",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/1137-checklist_7165eb4b",
      "display_name": "Checklist",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/445-reference-set-02-idor-prevention_802750c5",
      "display_name": "Reference Set 02 Idor Prevention",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/612-notifications_d933f28e",
      "display_name": "Notifications",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/079-conditional_bbec71f3",
      "display_name": "Orchestration Conditional",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/038-three_7bde1cef",
      "display_name": "Conversation Framework",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "daily assistant"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/021-framework_6d6e9dc4",
      "display_name": "Brand Guidelines",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "content creation"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/107-ai-sdk-5_bb521bed",
      "display_name": "Skill",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "content creation"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/052-python_352d08a0",
      "display_name": "Endpoints Routing",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/257-migration-checklist_2f4b8456",
      "display_name": "Migration Checklist",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/003-name-skill_6ec2bf56",
      "display_name": "Skill",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/096-threat_8d00c0dc",
      "display_name": "Reference Set 02 Idor Prevention",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/1054-hidden_conversations_sidebar_click_fix_788ea95e",
      "display_name": "Hidden Conversations Sidebar Click Fix",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/016-python_1c8cb677",
      "display_name": "Authentication",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/1038-faqs_51096eab",
      "display_name": "Faqs",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/003-name-skill_9f6e9d65",
      "display_name": "Skill",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/031-vue_91388a4b",
      "display_name": "Components",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    },
    {
      "local_path": "communication/1279-component-patterns_81c1b135",
      "display_name": "Component Patterns",
      "category": "communication",
      "subcategory": "",
      "tags": [
        "communication"
      ],
      "required": true,
      "reason": "Matches filters"
    }
  ],
  "total_count": 40,
  "generated_at": "2026-02-05T03:01:56.250992",
  "checksum": "d720a103240731d5c56b8bb8ea5d955ceea8efbebfd5bda422056ea11bc47888"
}
```

---

*Generated by [SkillFlow](https://github.com/tools-only/SkillFlow)*
