# Graph Viewer Sidebar Updates

**Date:** 2025-12-21

## Summary

Updated the `/docs/sims/graph-viewer/` MicroSim to match the new template in `skills/installer/references/learning-graph-viewer.md`. Added color swatches to the legend and moved the search form into the sidebar.

## Changes Made

### 1. Added glightbox Plugin to MkDocs

- Installed `mkdocs-glightbox` (v0.5.2)
- Added `glightbox` to plugins in `mkdocs.yml`

### 2. Legend Redesign (local.css)

Replaced table-based legend with flexbox layout:

**Removed:**
- `.legend-table` styles
- `.color-indicator` styles (second column with color name)

**Added:**
```css
.legend-item {
  display: flex;
  align-items: center;
  padding: 6px 0;
  gap: 8px;
}

.color-box {
  width: 20px;
  height: 20px;
  border-radius: 4px;
  border: 1px solid rgba(0,0,0,0.1);
  flex-shrink: 0;
}

.legend-item label {
  font-size: 13px;
  color: #555;
  cursor: pointer;
  flex: 1;
}

#legend {
  max-height: 400px;
  overflow-y: auto;
}
```

### 3. Legend Generation (script.js)

Updated `generateLegend()` function:

**Before:** Created table rows with two columns (checkbox+label, color indicator with color name)

**After:** Creates div elements with:
- Checkbox
- Color swatch (`<span class="color-box">`)
- Label text

Removed the second column that displayed the color name.

### 4. HTML Structure (main.html)

**Legend container:** Changed from `<table class="legend-table" id="legend-table">` to `<div id="legend">`

**Search form:** Moved from `#main` area into `#sidebar`:
- Added "Search" heading
- Renamed "Legend & Controls" to "Categories"
- Updated placeholder to "Search concepts..."

### New Sidebar Layout

```
[≡] Toggle button
━━━━━━━━━━━━━━━━━━━━━━
Search
[Search concepts...    ]

Categories
[Check All] [Uncheck All]
[✓] [■] Foundation Concepts
[✓] [■] Bloom's Taxonomy
[✓] [■] Visualization Types
...

Graph Statistics
Nodes: 200
Edges: 450
Orphans: 0
```

## Files Modified

| File | Changes |
|------|---------|
| `mkdocs.yml` | Added glightbox plugin |
| `docs/sims/graph-viewer/main.html` | Moved search to sidebar, changed legend to div |
| `docs/sims/graph-viewer/local.css` | Added color-box styles, removed table styles |
| `docs/sims/graph-viewer/script.js` | Updated generateLegend() for new structure |

## Testing

```bash
mkdocs serve
```

Navigate to: `http://127.0.0.1:8000/claude-skills/sims/graph-viewer/main.html`
