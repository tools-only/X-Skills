# Tasks

| Property | Value |
|----------|-------|
| **Name** | Tasks |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-data-scheduler/tasks.md) (⭐ 19) |
| **Original Path** | `openspec/changes/unify-data-scheduler/tasks.md` |
| **Category** | investment |
| **Subcategory** | analysis |
| **Tags** | investment |
| **Created** | 2026-02-06 |
| **Updated** | 2026-02-06 |
| **File Hash** | `b1b18f60a875ff6c...` |

## Description

[x] 3.1 在 UnifiedScheduler 中注册每日 16:00 缺失检测任务（missing_check_job）
 [x] 3.2 实现 missing_check_job() 方法
   [x] 检查是否为交易日，非交易日跳过
   [x] 调用 detect_missing_data() 获取缺失信息
   [x] 将检测结果记录到日志（INFO / WARNING 级别）
 [x] 3.3 缺失检测时间可通过 missing_check_time 配置项动态调整

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/unify-data-scheduler/tasks.md)*
