# Updates to Learning Graph Viewer Installer Skill

**Date:** 2025-12-19
**File:** `/skills/installer/references/learning-graph-viewer.md`
**Source:** `/Users/dan/Documents/ws/automating-instructional-design/docs/sims/graph-viewer/`

## Summary

Updated the learning-graph-viewer installer skill to reflect the design from the automating-instructional-design project's graph-viewer. Added comprehensive documentation for all CSS classes and vis-network options.

## Changes Made

### 1. Complete File Templates

Added full source code for all four viewer files:
- `main.html` - Complete HTML structure with sidebar, search, legend, stats, and graph container
- `local.css` - Full CSS with section comments explaining each block
- `script.js` - Complete JavaScript with all functions documented
- `index.md` - Documentation page template

### 2. CSS Documentation

Added detailed CSS reference section with tables documenting:

| Section | Classes Documented |
|---------|-------------------|
| Layout | `.container`, `.sidebar`, `.sidebar.collapsed`, `.sidebar-content`, `.graph-container`, `#network` |
| Sidebar Header | `.sidebar-header`, `.sidebar-header h4`, `.toggle-btn` |
| Search | `.search-container`, `.search-results`, `.search-result-item`, `.result-label`, `.result-category` |
| Legend | `.legend-container`, `.legend-controls`, `.legend-btn`, `#legend`, `.legend-item`, `.color-box` |
| Statistics | `.stats-container`, `#stats`, `#visible-nodes`, `#visible-edges`, `#foundational-nodes` |

Added key colors/variables table:

| Element | Color | Usage |
|---------|-------|-------|
| Sidebar header | `#2196F3` | Blue Material Design primary |
| Focus ring | `rgba(33, 150, 243, 0.1)` | Light blue focus glow |
| Stats numbers | `#2196F3` | Blue accent for statistics |
| Graph background | `aliceblue` | Light blue canvas background |
| Borders | `#ddd` | Light gray for borders |
| Text primary | `#333` | Dark gray for headings |
| Text secondary | `#555` / `#666` | Medium gray for labels |

### 3. vis-network Options Reference

Added comprehensive documentation for all vis-network configuration:

#### Groups Configuration
```javascript
visGroups[groupId] = {
    color: {
        background: groupInfo.color,
        border: groupInfo.color,
        highlight: { background: groupInfo.color, border: '#333' },
        hover: { background: groupInfo.color, border: '#666' }
    },
    font: { color: groupInfo.font?.color || 'black' }
};
```

#### Layout Options
```javascript
layout: {
    randomSeed: 42,        // Consistent positions across reloads
    improvedLayout: true   // Better initial spread
}
```

#### Physics Options (forceAtlas2Based)
| Parameter | Value | Purpose |
|-----------|-------|---------|
| `gravitationalConstant` | -50 | Repulsion force (negative = push apart) |
| `centralGravity` | 0.01 | Pull toward center (low = spread out) |
| `springLength` | 100 | Ideal edge length in pixels |
| `springConstant` | 0.08 | Edge spring stiffness |
| `damping` | 0.4 | Velocity damping |
| `avoidOverlap` | 0.5 | Node overlap prevention (0-1) |

#### Node Options
```javascript
nodes: {
    shape: 'box',
    margin: 4,
    font: { size: 14, face: 'Arial' },
    borderWidth: 2,
    shadow: true
}
```

#### Edge Options
```javascript
edges: {
    smooth: {
        type: 'cubicBezier',
        forceDirection: 'horizontal',
        roundness: 0.4
    },
    width: 1.5
}
```

#### Interaction Options
```javascript
interaction: {
    hover: true,
    tooltipDelay: 200,
    zoomView: true,
    dragView: true
}
```

### 4. CSS Section Comments

Added descriptive section headers throughout the CSS:
- `/* RESET AND BASE STYLES */`
- `/* MAIN LAYOUT CONTAINER */`
- `/* SIDEBAR STYLES */`
- `/* SIDEBAR HEADER */`
- `/* TOGGLE BUTTON */`
- `/* SIDEBAR CONTENT */`
- `/* SEARCH CONTAINER STYLES */`
- `/* SEARCH RESULTS DROPDOWN */`
- `/* LEGEND CONTAINER STYLES */`
- `/* STATISTICS CONTAINER STYLES */`
- `/* GRAPH CONTAINER STYLES */`
- `/* RESPONSIVE ADJUSTMENTS */`

### 5. New Troubleshooting Issue

Added Issue #4: Graph takes too long to stabilize

**Fix suggestions:**
- Increase `damping` (0.4 → 0.6) for faster settling
- Decrease `stabilization.iterations` for quicker initial render
- Increase `springConstant` for stiffer edges

## Files Referenced

Source files from automating-instructional-design project:
- `docs/sims/graph-viewer/main.html` (52 lines)
- `docs/sims/graph-viewer/script.js` (355 lines)
- `docs/sims/graph-viewer/local.css` (266 lines)
- `docs/sims/graph-viewer/index.md` (36 lines)

## Result

The skill file grew from ~270 lines to ~1189 lines with complete documentation for reproducing the graph-viewer design in any intelligent textbook project.
