# Control Center Metrics Caching

| Property | Value |
|----------|-------|
| **Name** | Control Center Metrics Caching |
| **Repository** | [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/features/v0.235.001/CONTROL_CENTER_METRICS_CACHING.md) (⭐ 110) |
| **Original Path** | `docs/explanation/features/v0.235.001/CONTROL_CENTER_METRICS_CACHING.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-13 |
| **Updated** | 2026-01-13 |
| **File Hash** | `699d023f61ad1291...` |

## Description

async function loadRefreshStatus() {
    const response = await fetch('/api/admin/controlcenter/refreshstatus');
    const result = await response.json();
    
    document.getElementById('lastRefreshTime').textContent = 
        result.last_refresh_formatted || 'Never';
}

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [microsoft/simplechat](https://raw.githubusercontent.com/microsoft/simplechat/main/docs/explanation/features/v0.235.001/CONTROL_CENTER_METRICS_CACHING.md)*
