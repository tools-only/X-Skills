# Design

| Property | Value |
|----------|-------|
| **Name** | Design |
| **Repository** | [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-stock-support/design.md) (⭐ 19) |
| **Original Path** | `openspec/changes/add-hk-stock-support/design.md` |
| **Category** | investment |
| **Subcategory** | trading |
| **Tags** | investment |
| **Created** | 2026-02-12 |
| **Updated** | 2026-02-12 |
| **File Hash** | `0c3c909703b5b97d...` |

## Description

def detect_market_type(ts_code: str) > MarketType:
    if ts_code.endswith('.HK'):
        return MarketType.HK_STOCK
    elif ts_code.endswith(('.SZ', '.SH')):
        return MarketType.A_SHARE
    return MarketType.UNKNOWN

**Tags:** `investment`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Yourdaylight/stock_datasource](https://raw.githubusercontent.com/Yourdaylight/stock_datasource/main/openspec/changes/add-hk-stock-support/design.md)*
