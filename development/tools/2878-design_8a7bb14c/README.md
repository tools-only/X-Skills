# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-custom-ai-workflow/design.md) (⭐ 19) |
| **Original Path** | `openspec/changes/add-custom-ai-workflow/design.md` |
| **Category** | development |
| **Subcategory** | tools |
| **Tags** | development |
| **Created** | 2026-01-19 |
| **Updated** | 2026-01-19 |
| **File Hash** | `8a7bb14cfd3eef5f...` |

## Description

当前系统架构：
 MCP Server (services/mcp_server.py): 动态发现并注册所有插件的工具，暴露HTTP接口
 Orchestrator (agents/orchestrator.py): 智能路由用户请求到对应Agent，支持MCP工具调用
 插件系统 (plugins/): 丰富的数据插件（日线、估值、财务、ETF、指数等）
 现有策略工作台 (views/StrategyWorkbench.vue): 策略管理和回测界面

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-custom-ai-workflow/design.md)*
