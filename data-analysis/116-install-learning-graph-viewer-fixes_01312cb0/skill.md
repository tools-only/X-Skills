# Learning Graph Viewer Installer Fixes

**Date:** 2025-12-17
**Project:** claude-skills/skills/installer

## Summary

Fixed three critical issues in the learning graph viewer that were causing incorrect behavior when installed via the installer skill.

## Issues Encountered and Fixed

### Issue 1: Hierarchical Layout Not Working

**Symptom:** The graph layout was incorrect and nodes were not positioned properly.

**Root Cause:** The script.js was using vis-network's hierarchical layout, which doesn't work well with DAGs that have complex dependency structures.

**Fix:** Changed from hierarchical layout to physics-based layout using forceAtlas2Based solver:

```javascript
// BEFORE (broken)
layout: {
    hierarchical: {
        enabled: true,
        direction: 'LR',
        sortMethod: 'directed'
    }
},
physics: {
    enabled: false
}

// AFTER (fixed)
layout: {
    randomSeed: 42,
    improvedLayout: true
},
physics: {
    enabled: true,
    solver: 'forceAtlas2Based',
    forceAtlas2Based: {
        gravitationalConstant: -50,
        centralGravity: 0.01,
        springLength: 100,
        springConstant: 0.08,
        damping: 0.4,
        avoidOverlap: 0.5
    }
}
```

### Issue 2: Legend Labels Showing Group IDs

**Symptom:** The legend displayed "BLOOM", "VISUA", "AUDIE" instead of human-readable names like "Bloom's Taxonomy", "Visualization Types", "Audience Adaptation".

**Root Cause:** The `classifierName` values in learning-graph.json were set to the group ID instead of human-readable names. This was caused by the csv-to-json.py script not having mappings for custom taxonomy IDs.

**Fix:** Updated the groups section in learning-graph.json with proper classifierName values:

```json
// BEFORE (broken)
"BLOOM": {
    "classifierName": "BLOOM",
    "color": "PeachPuff"
}

// AFTER (fixed)
"BLOOM": {
    "classifierName": "Bloom's Taxonomy",
    "color": "PeachPuff"
}
```

### Issue 3: Legend Colors Not Matching Node Colors

**Symptom:** When filtering to show only "Foundation Concepts" (red in legend), blue nodes appeared instead.

**Root Cause:** The script was setting node colors directly, but vis-network was overriding them with its default group color scheme because the `groups` option wasn't passed to the network configuration.

**Fix:** Build a visGroups configuration object from the JSON groups and pass it to vis-network:

```javascript
// Build vis-network groups configuration from JSON groups
const visGroups = {};
Object.entries(groups).forEach(([groupId, groupInfo]) => {
    visGroups[groupId] = {
        color: {
            background: groupInfo.color || 'lightgray',
            border: groupInfo.color || 'lightgray',
            highlight: {
                background: groupInfo.color || 'lightgray',
                border: '#333'
            },
            hover: {
                background: groupInfo.color || 'lightgray',
                border: '#666'
            }
        },
        font: {
            color: groupInfo.font?.color || 'black'
        }
    };
});

const options = {
    groups: visGroups,  // THIS IS REQUIRED!
    // ... other options
};
```

Also simplified the node creation since vis-network now handles colors via groups:

```javascript
// BEFORE (trying to set colors manually)
const nodes = new vis.DataSet(allNodes.map(node => ({
    ...node,
    color: groups[node.group]?.color || 'lightgray',
    font: { color: groups[node.group]?.font?.color || 'black' }
})));

// AFTER (let vis-network handle it via groups option)
const nodes = new vis.DataSet(allNodes);
```

## Files Modified

### In automating-instructional-design project:
- `docs/sims/graph-viewer/script.js` - Fixed layout and groups configuration
- `docs/learning-graph/learning-graph.json` - Fixed classifierName values

### In claude-skills installer:
- `skills/installer/references/learning-graph-viewer.md` - Updated documentation with:
  - New Step 3.5: Verify classifierName Values
  - New "Common Issues and Fixes" section
  - Updated prerequisites
- `skills/installer/references/assets/main.html` - Created template file
- `skills/installer/references/assets/script.js` - Created with all fixes applied
- `skills/installer/references/assets/local.css` - Created stylesheet
- `skills/installer/references/assets/index.md` - Created documentation page

## Key Lessons Learned

1. **vis-network groups option is required** - When nodes have a `group` property, vis-network will use its default color scheme unless you explicitly pass a `groups` configuration in the options.

2. **Don't use hierarchical layout for learning graphs** - The hierarchical layout in vis-network doesn't work well with DAGs that have complex cross-dependencies. Use physics-based layout (forceAtlas2Based) instead.

3. **classifierName must be human-readable** - The csv-to-json.py script may not have mappings for custom taxonomy IDs, resulting in classifierName being set to the group ID. Always verify and update these to human-readable names.

## Taxonomy Mapping Reference

| Group ID | classifierName |
|----------|----------------|
| FOUND | Foundation Concepts |
| BLOOM | Bloom's Taxonomy |
| VISUA | Visualization Types |
| LIBRA | Libraries & Tools |
| SPECI | Specification |
| COGNI | Cognitive Science |
| AUDIE | Audience Adaptation |
| EVALU | Evaluation & Testing |
| ITERA | Iteration & Workflow |
| ACCES | Accessibility |
| DEPLO | Deployment |
| CAPST | Capstone |
