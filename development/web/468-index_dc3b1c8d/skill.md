# Learning Graph Color Test

<iframe src="main.html" width="100%" height="880px" scrolling="no"></iframe>

## Overview

This interactive visualization tests the effectiveness of 17 pastel web-safe colors for differentiating nodes in a learning graph. The simulation displays 50 nodes distributed across 17 color groups, each with 1-3 random connections to other nodes.

**[Open Full Screen Version](main.html){ target="_blank" }**

## Purpose

The goal is to evaluate whether the selected pastel color palette provides sufficient visual distinction between different taxonomy groups in a learning graph visualization. This helps ensure that users can easily identify and distinguish between different concept categories.

## Interactive Features

- **50 Nodes**: Numbered 1-50, distributed evenly across 17 color groups
- **Random Connections**: Each node has 1-3 edges to other nodes (black color)
- **Ellipse Shape**: All nodes use the ellipse shape for consistency
- **Interactive**: Hover over nodes to highlight, drag to rearrange
- **Physics Simulation**: Nodes arrange themselves automatically using force-directed layout

## Color Palette

The visualization uses 17 web-safe pastel colors:

1. **MistyRose** - Soft pink/rose
2. **PeachPuff** - Light peach
3. **LightYellow** - Pale yellow
4. **Honeydew** - Light green
5. **PaleTurquoise** - Soft cyan
6. **AliceBlue** - Very light blue
7. **Lavender** - Light purple
8. **LavenderBlush** - Pink lavender
9. **Thistle** - Soft mauve/purple
10. **MintCream** - Very light mint
11. **LightCoral** - Soft coral
12. **Plum** - Medium purple
13. **Gainsboro** - Light gray
14. **PowderBlue** - Soft blue
15. **PaleGreen** - Soft green
16. **Aquamarine** - Blue-green
17. **LightPink** - Soft pink

## How to Use

1. **Observe Color Differentiation**: Can you easily distinguish between different color groups?
2. **Identify Patterns**: Notice how nodes of the same color group are numbered (1, 18, 35... are all Group 1)
3. **Interact**: Drag nodes around to see connections more clearly
4. **Check the Legend**: The color legend at the bottom shows all 17 colors with their names

## Technical Details

- **Library**: vis-network.js
- **Node Distribution**: Nodes are distributed round-robin across groups (Node 1 → Group 1, Node 2 → Group 2, ..., Node 18 → Group 1, etc.)
- **Edge Generation**: Random connections ensure realistic graph structure
- **Font Color**: Black text on light backgrounds, white text on darker colors (Plum, LightCoral)

## Design Considerations

This color palette was chosen for:

- **Accessibility**: All colors are light pastels with good contrast
- **Web-Safe**: Named CSS colors that work across all browsers
- **Distinctiveness**: Colors span the spectrum (reds, oranges, yellows, greens, blues, purples, grays)
- **No Hex Codes**: Uses only named web-safe colors for simplicity

## View Source

- [main.html](main.html) - Standalone visualization file
- [Learning Graph Generator Skill](../../skills/learning-graph-generator/) - Uses these colors

## Related Resources

- [vis-network Documentation](https://visjs.github.io/vis-network/docs/network/)
- [Web-Safe Color Names](https://www.w3schools.com/colors/colors_names.asp)
- [Learning Graph Schema](https://github.com/dmccreary/learning-graphs)
