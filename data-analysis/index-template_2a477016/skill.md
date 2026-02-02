---
title: {{TITLE}}
description: {{DESCRIPTION}}
quality_score: 85
image: /sims/{{MICROSIM_NAME}}/{{MICROSIM_NAME}}.png
og:image: /sims/{{MICROSIM_NAME}}/{{MICROSIM_NAME}}.png
twitter:image: /sims/{{MICROSIM_NAME}}/{{MICROSIM_NAME}}.png
social:
   cards: false
---
# {{TITLE}}

<iframe src="main.html" height="500px" width="100%" scrolling="no"></iframe>

[View {{TITLE}} Fullscreen](./main.html){ .md-button .md-button--primary }

## Embed This Visualization

Place the following line in your website to include this visualization:

```html
<iframe src="https://dmccreary.github.io/{{REPO_NAME}}/sims/{{MICROSIM_NAME}}/main.html" height="500px" width="100%" scrolling="no"></iframe>
```

## Overview

{{OVERVIEW}}

## How to Use

This interactive plot allows you to explore {{FUNCTION_DESCRIPTION}}:

1. **Drag the slider** to move the point along the curve
2. **Hover over the curve** to see exact x,y values at any point
3. **Use the toolbar** to zoom, pan, or download the plot as an image

## Mathematical Background

{{MATHEMATICAL_BACKGROUND}}

## Key Concepts

{{KEY_CONCEPTS}}

## Lesson Plan

### Learning Objectives

After using this visualization, students will be able to:

{{LEARNING_OBJECTIVES}}

### Discussion Questions

{{DISCUSSION_QUESTIONS}}

### Activities

1. **Explore the Function**: Move the slider from minimum to maximum. What patterns do you observe?
2. **Find Key Points**: Identify where the function crosses zero, reaches its maximum, and reaches its minimum.
3. **Make Predictions**: Before moving the slider, predict what the y-value will be for a given x-value.

## Technical Notes

This MicroSim uses:

- **Plotly.js 2.27.0**: A high-level JavaScript charting library built on D3.js and stack.gl
- Plotly provides built-in interactivity including zoom, pan, hover tooltips, and image export

### Plotly Features Used

- `Plotly.newPlot()`: Creates the initial visualization
- `Plotly.restyle()`: Efficiently updates trace data without full re-render
- `Plotly.Plots.resize()`: Handles responsive resizing
- Mode bar with zoom, pan, and download tools

## References

{{REFERENCES}}
