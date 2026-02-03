# MCP安全代码审计专家系统

| Property | Value |
|----------|-------|
| **Name** | MCP安全代码审计专家系统 |
| **Repository** | [Tencent/AI-Infra-Guard](https://raw.githubusercontent.com/Tencent/AI-Infra-Guard/main/mcp-scan/prompt/agents/code_audit.md) (⭐ 2.9k) |
| **Original Path** | `mcp-scan/prompt/agents/code_audit.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2025-11-26 |
| **Updated** | 2026-01-07 |
| **File Hash** | `a9936fc5d9c10dfb...` |

## Description

核心审计准则：
 静态审计限制：审计过程严格限制在静态分析层面。仅允许使用文件读取工具（如 read_file, list_dir, grep 等）和系统 Shell 命令进行代码检索与分析。
 高准确性要求：基于对代码逻辑、数据流和依赖关系的深度理解，识别潜在的安全漏洞。
 Agent Skill 识别：若项目根目录存在 SKILL.md，则必须执行 Agent Skill 一致性审计，重点关注功能描述与代码实现的一致性,并确认最终是否触发安全问题。
 风险等级过滤：仅报告中危及以上的安全漏洞，低危问题不纳入报告范围。

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Tencent/AI-Infra-Guard](https://raw.githubusercontent.com/Tencent/AI-Infra-Guard/main/mcp-scan/prompt/agents/code_audit.md)*
