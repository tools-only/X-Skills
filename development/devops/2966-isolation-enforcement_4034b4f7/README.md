# Isolation Enforcement

| Property | Value |
|----------|-------|
| **Name** | Isolation Enforcement |
| **Repository** | [PolicyEngine/policyengine-claude](https://raw.githubusercontent.com/PolicyEngine/policyengine-claude/main/agents/country-models/isolation-enforcement.md) (⭐ 14) |
| **Original Path** | `agents/country-models/isolation-enforcement.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2025-10-19 |
| **Updated** | 2025-10-19 |
| **File Hash** | `4034b4f7adba17c9...` |

## Description

if [[ "$CLAUDE_AGENT" == "rules_engineer" ]]; then
  if [[ "$1" =~ "pe.tests" ]] || grep q "test_creator" "$@" 2>/dev/null; then
    echo "ERROR: Rules Engineer cannot access Test Creator files"
    exit 1
  fi
fi

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [PolicyEngine/policyengine-claude](https://raw.githubusercontent.com/PolicyEngine/policyengine-claude/main/agents/country-models/isolation-enforcement.md)*
