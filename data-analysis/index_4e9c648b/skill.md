---
title: Skill Development Priority Matrix Visualization
description: Interactive Chart.js visualization showing skill development priority matrix visualization
image: /sims/skill-impact-chart/skill-impact-chart.png
og:image: /sims/skill-impact-chart/skill-impact-chart.png
quality_score: 70
---


# Skill Development Priority Matrix Visualization


An interactive Chart.js bubble chart that visualizes the priority matrix for skill development in creating educational visualizations. This chart helps identify which visualization types should be prioritized based on their impact (frequency of use) versus effort (similarity to existing MicroSims).

## Interactive Chart



```html


```

[View Chart Fullscreen](./main.html){ .md-button .md-button--primary }

## Overview

This visualization uses a **bubble chart** (scatter plot with variable bubble sizes) to display the relationship between:

- **X-axis (Effort)**: Development effort required (0-10 scale), based on dissimilarity to existing MicroSims
- **Y-axis (Impact)**: Frequency of use in the textbook (0-10 scale)
- **Bubble Size**: Number of instances of each visualization type
- **Color**: Status indicator (Red = Build needed, Green = Already exists)

### Quadrant Interpretation

The chart is divided into four quadrants (using 5 as the midpoint):

1. **High Impact, Low Effort** (Green shading) - **Priority 1**: Build these first for maximum ROI
2. **High Impact, High Effort** (Yellow shading) - **Priority 2**: Important but resource-intensive
3. **Low Impact, Low Effort** (Blue shading) - **Priority 3**: Quick wins but limited impact
4. **Low Impact, High Effort** (Red shading) - **Priority 4**: Avoid unless strategic necessity

## Features

### Interactive Elements

- **Hover tooltips**: Display detailed information including type, count, impact, effort, priority score, and status
- **Labeled bubbles**: Each bubble shows the visualization type and instance count
- **Quadrant backgrounds**: Subtle color shading indicates priority zones
- **Legend**: Distinguishes between "Build" and "Exists" items
- **Quadrant lists**: Below the chart, items are categorized by quadrant

### Visual Design

- **Bubble sizes**: Scaled proportionally from smallest (1 instance) to largest (15 instances)
- **Status colors**: Red for items that need to be built, green for existing items
- **Grid lines**: Light grid for easy value reading
- **Quadrant dividers**: Dashed lines at x=5 and y=5

## Customization Guide

### Adjusting Chart Margins and Layout

To prevent bubbles from being clipped at the edges, you can modify several parameters in [main.html](main.html):

#### 1. Chart Scale Ranges

In the `options.scales` section (around line 100-130), adjust the `min` and `max` values:

```javascript
scales: {
    x: {
        min: -0.5,  // Add padding before 0 (default: 0)
        max: 10.5,  // Add padding after 10 (default: 10)
        // ... other options
    },
    y: {
        min: -0.5,  // Add padding before 0 (default: 0)
        max: 10.5,  // Add padding after 10 (default: 10)
        // ... other options
    }
}
```

**Effect**: Extends the visible range beyond 0-10, providing space for edge bubbles.

#### 2. Chart Padding

Add layout padding in the `options` section:

```javascript
options: {
    layout: {
        padding: {
            top: 20,
            right: 20,
            bottom: 20,
            left: 20
        }
    },
    // ... other options
}
```

**Effect**: Adds space around the entire chart area.

#### 3. Aspect Ratio

Modify the aspect ratio to change chart dimensions (around line 80):

```javascript
options: {
    maintainAspectRatio: true,
    aspectRatio: 1.5,  // Width:Height ratio (default: 1.5)
    // Try 1.8 for wider, 1.2 for taller
}
```

**Effect**: Changes the width-to-height ratio of the chart.

#### 4. Container Height

In [style.css](style.css), adjust the chart container height (around line 18):

```css
.chart-container {
    height: 600px;  /* Default: 600px */
    /* Try 700px or 800px for more vertical space */
}
```

**Effect**: Increases the vertical space available for the chart.

### Adjusting Bubble Sizes

In [main.html](main.html), find the `getBubbleSize` function (around line 65):

```javascript
function getBubbleSize(count) {
    const minSize = 8;   // Minimum bubble radius (default: 8)
    const maxSize = 30;  // Maximum bubble radius (default: 30)
    return minSize + ((count - minCount) / (maxCount - minCount)) * (maxSize - minSize);
}
```

**Parameters to adjust**:
- `minSize`: Smaller values create tinier bubbles for low-count items
- `maxSize`: Larger values create bigger bubbles for high-count items

**Recommended ranges**:
- For more subtle differences: `minSize: 10, maxSize: 25`
- For more dramatic differences: `minSize: 5, maxSize: 40`

### Changing Colors

#### Status Colors

Modify the color scheme (around line 50):

```javascript
const colors = {
    'Build': 'rgba(220, 53, 69, 0.8)',  // Red with 80% opacity
    'Exists': 'rgba(40, 167, 69, 0.8)'  // Green with 80% opacity
};

const borderColors = {
    'Build': 'rgb(220, 53, 69)',   // Solid red border
    'Exists': 'rgb(40, 167, 69)'   // Solid green border
};
```

**Format**: `rgba(red, green, blue, alpha)` where values are 0-255, alpha is 0-1

#### Quadrant Background Colors

In the `afterDraw` plugin (around line 180):

```javascript
ctx.fillStyle = 'green';      // High Impact, Low Effort
ctx.fillStyle = 'yellow';     // High Impact, High Effort
ctx.fillStyle = 'lightblue';  // Low Impact, Low Effort
ctx.fillStyle = 'red';        // Low Impact, High Effort
```

The opacity is controlled by `ctx.globalAlpha = 0.05` (5% opacity).

### Modifying Data

To update the visualization data, edit the `data` array (around line 36):

```javascript
const data = [
    {
        type: 'graph-model',        // Visualization type name
        count: 10,                  // Number of instances
        impact: 6.7,                // Impact score (0-10)
        effort: 1,                  // Effort score (0-10)
        priority: 6.67,             // Calculated priority score
        status: 'Build'             // 'Build' or 'Exists'
    },
    // ... more items
];
```

**Note**: Priority score is typically calculated as `impact / effort`.

### Adjusting Label Positioning

To move the labels above bubbles (around line 210):

```javascript
ctx.fillText(
    `${dataPoint.label} (n=${dataPoint.count})`,
    element.x,
    element.y - element.options.radius - 8  // Offset from bubble (default: -8)
);
```

**Adjust the offset**:
- Larger negative values (e.g., `-12`) move labels further from bubbles
- Smaller values (e.g., `-5`) move labels closer

### Customizing Tooltips

Modify tooltip content in the `options.plugins.tooltip.callbacks` section (around line 105):

```javascript
tooltip: {
    callbacks: {
        label: function(context) {
            const item = context.raw;
            return [
                `Type: ${item.label}`,
                `Count: ${item.count}`,
                `Impact: ${item.y}`,
                `Effort: ${item.x}`,
                `Priority Score: ${item.priority}`,
                `Status: ${item.status}`
                // Add more lines here
            ];
        }
    }
}
```

## Technical Details

### Dependencies

- **Chart.js v4.4.0**: Loaded from CDN (`https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js`)
- **No other external dependencies required**

### Browser Compatibility

- Modern browsers with Canvas support (Chrome, Firefox, Safari, Edge)
- Responsive design works on mobile and tablet devices

### File Structure

```
skill-impact-chart/
├── index.md          # This documentation page
├── main.html         # Interactive Chart.js visualization
└── style.css         # Styling for chart and quadrant lists
```

### Data Flow

1. **Data source**: Hardcoded array in [main.html](main.html) (line 36-47)
2. **Processing**: Bubbles are sized, colored, and positioned based on data
3. **Rendering**: Chart.js renders to HTML5 Canvas element
4. **Interaction**: Hover events trigger tooltips, quadrant lists are auto-populated

## Use Cases

This visualization pattern is useful for:

- **Prioritization matrices**: Any 2D priority framework (Eisenhower Matrix, Risk/Impact, etc.)
- **Portfolio analysis**: Product portfolios, project portfolios
- **Resource allocation**: Identifying high-ROI opportunities
- **Technology selection**: Comparing options on multiple dimensions
- **Strategic planning**: Visualizing trade-offs between competing factors

## Future Enhancements

Potential improvements for this visualization:

1. **Dynamic data loading**: Load data from JSON file instead of hardcoded array
2. **Interactive filtering**: Checkboxes to show/hide specific statuses or quadrants
3. **Click-to-explore**: Click bubbles to see detailed breakdown or examples
4. **Animation**: Animate bubbles on load for visual impact
5. **Export functionality**: Download chart as PNG or PDF
6. **Comparison mode**: Compare two different time periods or scenarios

## References

- [Chart.js Documentation](https://www.chartjs.org/docs/latest/)
- [Bubble Chart Guide](https://www.chartjs.org/docs/latest/charts/bubble.html)
- [Priority Matrix Theory](https://en.wikipedia.org/wiki/Priority_Matrix)

## License

This visualization is part of the IT Management Graph educational project and is licensed under Creative Commons ShareAlike Attribution Noncommercial (CC BY-NC-SA 4.0).