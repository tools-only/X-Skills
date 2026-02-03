# MicroSim Add Icons

The microsim-add-icons skill adds clickable Creative Commons license and fullscreen
navigation icons to the control region of an existing p5.js MicroSim. The icons
appear in the lower right corner using distance-based click detection.

## Key Capabilities

This skill enhances MicroSims with:

- **Creative Commons Icon**: Links to license information
- **Fullscreen Icon**: Enables fullscreen viewing mode
- **Distance-based Click Detection**: Responsive icon interaction
- **Consistent Positioning**: Lower right corner placement

## When to Use

Use this skill when:

- Adding license information access to a MicroSim
- Enabling fullscreen mode via an icon
- Enhancing a MicroSim with clickable UI elements
- Following the icons demo pattern from the MicroSims repository
- Standardizing icon behavior across multiple MicroSims

## What Gets Added

The skill adds approximately 40 lines of code including:

1. **Icon Variables**: Position and size declarations
2. **Draw Functions**: Icon rendering in the draw loop
3. **Click Handler**: `mousePressed()` function for icon interaction
4. **Position Updates**: Icon repositioning in `windowResized()`

## Prerequisites

The target file should be a p5.js MicroSim following the standard pattern with:

- Global variables section
- `setup()` function
- `draw()` function
- `windowResized()` function

## Workflow

1. Identify the target MicroSim JavaScript file
2. Read the existing file to understand its structure
3. Add icon variables to the global section
4. Add icon drawing code to the `draw()` function
5. Add or update the `mousePressed()` function
6. Update icon positions in `windowResized()`

## Icon Appearance

The icons use simple geometric shapes:

- **CC Icon**: Circle with "CC" text
- **Fullscreen Icon**: Four corner brackets

Both icons change appearance on hover to indicate interactivity.

## Integration

This skill is typically used after a MicroSim is functionally complete
to add polish and standard navigation elements before publication.
