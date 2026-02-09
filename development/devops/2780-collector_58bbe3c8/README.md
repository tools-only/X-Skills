# Collector

| Property | Value |
|----------|-------|
| **Name** | Collector |
| **Repository** | [julianobarbosa/claude-code-skills](https://raw.githubusercontent.com/julianobarbosa/claude-code-skills/main/skills/opentelemetry-skill/references/COLLECTOR.md) (⭐ 14) |
| **Original Path** | `skills/opentelemetry-skill/references/COLLECTOR.md` |
| **Category** | development |
| **Subcategory** | devops |
| **Tags** | development |
| **Created** | 2025-12-16 |
| **Updated** | 2025-12-16 |
| **File Hash** | `58bbe3c81033d75c...` |

## Description

yaml
receivers:
  otlp:
    protocols:
      grpc:
        endpoint: ${env:MY_POD_IP}:4317
        max_recv_msg_size_mib: 4
      http:
        endpoint: ${env:MY_POD_IP}:4318
        cors:
          allowed_origins: [""]

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [julianobarbosa/claude-code-skills](https://raw.githubusercontent.com/julianobarbosa/claude-code-skills/main/skills/opentelemetry-skill/references/COLLECTOR.md)*
