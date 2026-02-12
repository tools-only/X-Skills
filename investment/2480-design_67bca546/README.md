# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/refactor-redis-queue-sync/design.md) (⭐ 19) |
| **Original Path** | `openspec/changes/refactor-redis-queue-sync/design.md` |
| **Category** | investment |
| **Subcategory** | trading |
| **Tags** | investment |
| **Created** | 2026-01-26 |
| **Updated** | 2026-01-26 |
| **File Hash** | `67bca546308c5003...` |

## Description

此外，项目对 TuShare 数据拉取存在 HTTP 代理需求，代理通过 stock_datasource/core/proxy.py 的 proxy_context() 在“数据拉取上下文”内设置进程环境变量，从而影响 requests/tushare。

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/refactor-redis-queue-sync/design.md)*
