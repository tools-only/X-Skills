# Architecture

| Property | Value |
|----------|-------|
| **Name** | Architecture |
| **Repository** | [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/scenarios/S04-service-validation/scenario/architecture.md) (⭐ 57) |
| **Original Path** | `scenarios/S04-service-validation/scenario/architecture.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2026-01-19 |
| **Updated** | 2026-01-19 |
| **File Hash** | `723c39951ffd7884...` |

## Description

subgraph "Data Tier"
                SQL_SRV[SQL Server<br/>sqlsaifswc01xxx<br/>Entra ID Admin Only]
                SQL_DB[SQL Database<br/>sqldbsaifswc01<br/>Basic Tier]
            end

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jonathan-vella/azure-agentic-infraops](https://raw.githubusercontent.com/jonathan-vella/azure-agentic-infraops/main/scenarios/S04-service-validation/scenario/architecture.md)*
