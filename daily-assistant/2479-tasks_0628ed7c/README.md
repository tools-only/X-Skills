# Tasks

| Property | Value |
|----------|-------|
| **Name** | Tasks |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/tasks.md) (⭐ 19) |
| **Original Path** | `openspec/changes/optimize-plugin-dependencies/tasks.md` |
| **Category** | daily-assistant |
| **Subcategory** | scheduling |
| **Tags** | daily assistant |
| **Created** | 2026-01-13 |
| **Updated** | 2026-01-23 |
| **File Hash** | `0628ed7c399e7ebc...` |

## Description

1. ✅ TradeCalendarService 可全局访问，提供统一的交易日期查询
2. ✅ 所有 daily 类插件正确声明对应的 basic 依赖
3. ✅ 执行插件前自动检查依赖是否满足
4. ✅ 依赖未满足时返回清晰的错误信息，包含缺失的依赖列表
5. ✅ 现有功能不受影响（向后兼容）
6. ✅ 支持按类别筛选插件（A股/港股/指数/ETF）
7. ✅ 支持批量触发同步任务
8. ✅ 可选依赖（如复权因子）可关联同步
9. ✅ 定时调度可配置并自动按依赖顺序执行
10. ✅ 界面显示操作说明和帮助提示

**Tags:** `daily assistant`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/optimize-plugin-dependencies/tasks.md)*
