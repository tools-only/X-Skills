# Spec

| Property | Value |
|----------|-------|
| **Name** | Spec |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/spec.md) (⭐ 19) |
| **Original Path** | `docs/spec.md` |
| **Category** | investment |
| **Subcategory** | trading |
| **Tags** | investment |
| **Created** | 2025-10-24 |
| **Updated** | 2025-10-24 |
| **File Hash** | `426a57b42dc649e1...` |

## Description

> 版本：v1.1（含 SchemaonAPI 动态表构建 / ODS 自动演进）
> 决策要点：A股使用原子接口；历史回填 20200101→今日；日更调度 18:00 Asia/Shanghai；不复权为主表；并发仅做速率限制（不超过积分档限制）；失败任务入表记录；与交易日/公告对齐行数并做涨跌停/停复牌一致性 DQ；港股沿用 TuShare ts_code 规范且当前仅占位；权限控制暂缓（后续对外 API 再启用）。

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/docs/spec.md)*
