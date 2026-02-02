# Mermaid Hover Tooltip Fix

**Date:** 2024-12-14
**Issue:** Hover tooltips not working in iframes
**Status:** Resolved

## Problem Description

When testing the Mermaid workflow diagram at `http://localhost:8000/microsims/sims/basic/#workflow-example-microsim-creation`, the hover tooltips were not appearing. However, the same diagram worked correctly when viewed directly at `http://localhost:8000/microsims/sims/microsim-creation-workflow/`.

Both pages used iframes to embed the same `main.html` file, so the issue was not with the embedding method itself.

## Root Cause Analysis

The original tooltip initialization used a fixed 500ms timeout:

```javascript
// Original problematic code
setTimeout(setupTooltips, 500);
```

This approach is unreliable because:

1. **Variable rendering time**: Mermaid.js rendering time varies depending on diagram complexity, page load order, and browser performance
2. **Iframe context**: When loaded in an iframe on a longer page (like `basic.md` with multiple iframes), the Mermaid diagram may take longer to render
3. **Race condition**: If Mermaid hasn't finished rendering by 500ms, the `.node` elements don't exist yet, and `setupTooltips()` attaches event listeners to nothing

## Solution Implemented

Replaced the fixed timeout with a polling-based approach that waits for Mermaid to actually finish rendering:

```javascript
// Robust polling approach - waits for Mermaid to finish rendering
function waitForMermaid() {
    const mermaidDiv = document.querySelector('.mermaid');
    const svg = mermaidDiv.querySelector('svg');
    if (svg && document.querySelectorAll('.node').length > 0) {
        setupTooltips();
    } else {
        // Check again after a short delay
        setTimeout(waitForMermaid, 100);
    }
}

// Start checking after initial load
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', () => setTimeout(waitForMermaid, 100));
} else {
    setTimeout(waitForMermaid, 100);
}
```

This approach:
- Polls every 100ms until both conditions are met
- Checks for SVG element presence (Mermaid creates this)
- Checks for `.node` elements (these are the interactive diagram nodes)
- Works regardless of how long Mermaid takes to render

## Files Modified

### 1. microsims repository

**File:** `/Users/dan/Documents/ws/microsims/docs/sims/microsim-creation-workflow/main.html`

**Changes:**
- Replaced `setTimeout(setupTooltips, 500)` with polling-based `waitForMermaid()` function
- Added `font-family: Arial, Helvetica, sans-serif` to tooltip CSS

### 2. claude-skills repository

**File:** `/Users/dan/Documents/ws/claude-skills/skills/microsim-generator/references/mermaid-guide.md`

**Changes:**

1. **Added new section: "Interactive Tooltips (Required Feature)"**
   - Tooltip HTML structure (`<div id="tooltip"></div>`)
   - Complete CSS styling with Arial font family
   - Tooltip data structure pattern
   - Robust polling initialization code
   - Full main.html template with tooltips
   - Tooltip content writing guidelines

2. **Updated Step 4.1 (Create main.html)**
   - Added `{{TOOLTIP_DATA}}` placeholder to list
   - Noted that tooltips are a required default feature

3. **Updated Step 6 (Validate and Test)**
   - Added tooltip verification checklist items:
     - Confirm every node has a tooltip entry
     - Verify Arial font family
     - Test iframe compatibility

4. **Updated Troubleshooting section**
   - Added "Tooltips not appearing" with solutions
   - Added "Tooltips not working in iframes" with solutions
   - Added "Tooltip positioning incorrect" with solutions

5. **Updated Step 7 (Inform the User)**
   - Added "Interactive hover tooltips on all nodes" to features list

6. **Updated Example 1**
   - Added tooltip data example alongside Mermaid code

## Tooltip CSS Specification

The tooltip must use this CSS (note the required font-family):

```css
#tooltip {
    position: absolute;
    background-color: #333;
    color: #fff;
    padding: 8px 12px;
    border-radius: 6px;
    font-family: Arial, Helvetica, sans-serif;
    font-size: 14px;
    max-width: 280px;
    pointer-events: none;
    opacity: 0;
    transition: opacity 0.2s;
    z-index: 1000;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}
#tooltip.visible {
    opacity: 1;
}
.node {
    cursor: pointer;
}
```

## Testing Verification

After the fix:
- Tooltips work on the dedicated page: `http://localhost:8000/microsims/sims/microsim-creation-workflow/`
- Tooltips work when embedded via iframe: `http://localhost:8000/microsims/sims/basic/#workflow-example-microsim-creation`
- Tooltips appear with Arial font
- Tooltip positioning respects viewport boundaries

## Key Takeaways

1. **Never use fixed timeouts for DOM readiness** - Always poll or use observers when waiting for dynamically generated content
2. **Test in both standalone and iframe contexts** - Behavior can differ due to load timing
3. **Document the pattern** - Added comprehensive documentation to the skill so future Mermaid diagrams include robust tooltips by default

## Related Links

- [Mermaid.js Documentation](https://mermaid.js.org/)
- Microsims basic examples page: `/docs/sims/basic.md`
- Mermaid generator skill: `/skills/microsim-generator/references/mermaid-guide.md`
