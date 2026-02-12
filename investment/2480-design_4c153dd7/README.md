# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-top-list-tracking/design.md) (⭐ 19) |
| **Original Path** | `openspec/changes/add-top-list-tracking/design.md` |
| **Category** | investment |
| **Subcategory** | trading |
| **Tags** | investment |
| **Created** | 2026-01-13 |
| **Updated** | 2026-01-13 |
| **File Hash** | `4c153dd7ea604451...` |

## Description

2. 持仓管理页面集成：在 PortfolioView.vue 中显示持仓股票的龙虎榜状态
vue
<! 在持仓列表中添加龙虎榜状态列 >
<template topliststatus="{ row }">
  <TopListStatusBadge :tscode="row.ts_code" />
</template>

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-top-list-tracking/design.md)*
