# MicroSim Screen Capture

The microsim-screen-capture skill automates the capture of high-quality screenshots
for MicroSim visualizations using Chrome headless mode. It properly handles dynamic
JavaScript content and ensures visualizations fully render before capture.

## Key Capabilities

This skill provides:

- **Headless Chrome Capture**: Automated browser-based screenshots
- **JavaScript Rendering**: Waits for dynamic content to load
- **CDN Library Support**: Handles external libraries (p5.js, vis-network, etc.)
- **Configurable Delay**: Adjustable wait time for complex animations
- **Consistent Output**: Standard image dimensions and format

## When to Use

Use this skill when:

- Creating preview images for MicroSim social media metadata
- Generating screenshots for documentation
- Capturing visualizations for quality assessment
- Working with microsim-standardization for 100/100 quality scores
- Creating thumbnails for MicroSim galleries

## Typical User Requests

- "Create a screenshot of this MicroSim"
- "I need a preview image for the org-chart MicroSim"
- "Generate a social media preview for this visualization"
- "Capture a screenshot of the main.html file"

## Workflow

1. **Validate MicroSim**: Verify directory exists with main.html
2. **Run Capture Script**: Execute shell script with path argument
3. **Wait for Render**: Allow time for JavaScript execution
4. **Save Screenshot**: Output PNG to MicroSim directory

## Prerequisites

- Chrome or Chromium browser installed
- MicroSim directory with main.html file
- MicroSim name following kebab-case convention

## Script Usage

```bash
bash scripts/capture-screenshot.sh <microsim-directory-path>
```

Example:
```bash
bash scripts/capture-screenshot.sh docs/sims/org-chart
```

## Output

The script generates:

- **screenshot.png**: Full-resolution capture in MicroSim directory
- Dimensions: Typically 1200x630 (optimized for social sharing)

## Integration

This skill works with the microsim-standardization skill to ensure
all MicroSims have preview images for social media sharing and documentation.
The screenshot is referenced in the MicroSim's index.md frontmatter.
