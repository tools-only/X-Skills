---
description: Describe images and screenshots by viewing content with multimodal Read tool and documenting only visible elements. Activates on what's in this image, describe this screenshot, summarize this diagram, image summary, what does this screenshot show, explain this diagram, break down this chart. Handles UI screenshots, architecture diagrams, charts, photos, code screenshots, and terminal output. Never infers from filenames.
---

# Image Summarization

Methodology for summarizing visual content including screenshots, diagrams, charts, photos, code screenshots, and terminal output.

## Reading Images

The model MUST use the Read tool to view images. Claude Code is multimodal and can directly read image files.

**Supported formats**:

- .png, .jpg, .jpeg, .gif, .webp, .svg

**Special handling for SVG**:

- Read the image with Read tool to view it visually
- ALSO read the SVG file as text to extract structured data (labels, IDs, embedded text)

**Example**:

```text
Read("/path/to/diagram.svg")        # Visual view
Read("/path/to/diagram.svg")        # Text view for extracting labels
```

## Image Type Strategies

The summarization approach varies by image type. The model MUST identify the image type and apply the corresponding strategy.

### Screenshots (UI)

The model MUST describe:

- Layout structure (header, sidebar, main content, footer)
- Visible text (labels, headings, button text, menu items)
- Interactive elements (buttons, dropdowns, checkboxes, input fields)
- Application state indicators (error messages, success banners, loading spinners)
- Visual hierarchy (what draws attention first)

**Example output structure**:

```text
UI screenshot showing [application name if visible] with three-column layout.
Left sidebar contains navigation menu with items: [list items].
Main content area displays [describe primary content].
Top banner shows [error/success/info message if present].
```

### Architecture Diagrams

The model MUST describe:

- Components (boxes, circles, shapes with labels)
- Connections (arrows, lines) with directionality
- Data flow direction (source → destination)
- Labels on components and connections
- Groupings or boundaries (dotted lines, colored regions)

**Example output structure**:

```text
Architecture diagram with 5 components:
- [Component A] connects to [Component B] via [labeled arrow]
- Data flows from [source] through [intermediate] to [destination]
- [Component X] is grouped with [Component Y] within [boundary label]
```

### Charts and Graphs

The model MUST describe:

- Chart type (bar chart, line graph, pie chart, scatter plot, histogram)
- Axes labels and scale (X-axis, Y-axis, units)
- Data trends (increasing, decreasing, cyclical, outliers)
- Notable data points (peaks, valleys, intersections)
- Legend entries if present

**Example output structure**:

```text
Line graph titled [title] showing [Y-axis label] over [X-axis label].
Trend: [describe pattern - e.g., "steady increase from 10 to 50 over time period"].
Notable points: peak at [X value] with [Y value], lowest point at [X value].
```

### Photos

The model MUST describe:

- Subject (person, object, scene)
- Setting (location type, environment)
- Relevant contextual details (time of day if visible, weather, activity)
- Focus: describe what appears to be the intended subject

**Do NOT**:

- Speculate about intent or meaning beyond what is visible
- Identify people by name unless text labels are present
- Infer relationships or emotions without visible evidence

### Code Screenshots

The model MUST:

- Extract visible code text verbatim
- Describe programming language if identifiable (from syntax, filename, or IDE indicators)
- Describe apparent purpose from visible context (function names, comments, file structure)
- Note if code is truncated (top/bottom cut off)

**Example output structure**:

```text
Code screenshot showing [language] file [filename if visible].
Visible code defines [function/class name] that [describe from visible code].
Lines [X-Y] are visible; code appears truncated at [top/bottom/both].
```

### Terminal Output Screenshots

The model MUST:

- Extract visible command(s) from prompt
- Extract visible output text
- Describe command purpose from visible context
- Note if output is truncated

**Example output structure**:

```text
Terminal screenshot showing command: [extract exact command]
Output: [extract visible output text]
[If truncated] Output appears truncated; [X] lines visible.
```

## Fidelity Rules for Images

The model MUST describe only what IS visible in the image.

**Prohibited behaviors**:

- Inferring functionality from UI elements ("this button saves the file" vs "button labeled Save")
- Guessing obscured text ("the label probably says..." vs "label partially obscured, visible text: [readable portion]")
- Assuming diagram connections represent specific protocols without labels
- Identifying chart data points beyond what is readable

**Required behaviors**:

- State when text is partially obscured: "partially visible text: [what's readable]"
- State when labels are unclear: "label appears to read [best interpretation] (uncertain)"
- Describe visible state indicators: "red error icon next to field X" not "field X has an error"
- Use conditional language for uncertain interpretations

## Output Format

The model MUST produce a structured summary following the format defined in [Structured Summary](../summarizer/templates/structured.md).

**Image-specific frontmatter**:

```yaml
---
source_type: image
source_path: "/absolute/path/to/image.png"
summarized_at: "2026-02-06T10:30:00Z"
method: abstractive
word_count_source: null
word_count_summary: <integer>
compression_ratio: null
confidence: high | medium | low
confidence_notes: "Image clearly visible with high resolution" | "Some labels obscured" | "Low resolution, text partially unreadable"
---
```

**Required sections**:

1. **Summary** - BLUF description of what the image shows
2. **What Was Found** - itemized list of visible elements
3. **What Was NOT Found** - expected elements that are absent or not visible
4. **Uncertain** - ambiguous or partially obscured elements
5. **Sources** - image file path with read timestamp

## Confidence Assessment for Images

**High confidence**:

- Image is high resolution with clear, readable text
- Labels and UI elements are fully visible
- No ambiguity in interpretation

**Medium confidence**:

- Some labels or text are small but readable
- Parts of the image are cropped but main content is clear
- Interpretation required for unlabeled elements

**Low confidence**:

- Low resolution or blurry image
- Text is partially unreadable
- Significant portions cropped or obscured
- Diagram lacks labels, requiring guesswork

## Output Rendering

1. **Read template** - Load the template file at `../summarizer/templates/{format_id}.md` (default: `structured`). The template defines the schema, required sections, and fidelity constraints for the selected format.
2. **Render** - Produce output following the template's Schema section. Use the template's Example as a reference for structure and style.
3. **Verify fidelity** - Confirm the output satisfies the template's Fidelity Constraints and all applicable [Fidelity Rules](../summarizer/references/fidelity-rules.md).

## Integration with Fidelity Rules

All image summaries MUST comply with the shared fidelity rules defined in [Fidelity Rules](../summarizer/references/fidelity-rules.md).

**Key rules for images**:

- **Rule 1: Read Before Summarizing** - Use Read tool to view the image, never describe from filename alone
- **Rule 2: Extract Before Abstracting** - Extract visible text (labels, code, terminal output, error messages, headings) verbatim before describing the image. For SVGs, also extract text from the markup
- **Rule 3: Preserve Counts and Specifics** - Count visible elements exactly ("5 buttons" not "several buttons")
- **Rule 4: Distinguish Absence from Nonexistence** - "No visible label on component X" not "component X is unlabeled"
- **Rule 6: State Confidence Explicitly** - Assess based on image clarity, resolution, and completeness

## Examples

### Screenshot Summary

```yaml
---
source_type: image
source_path: "/home/user/screenshot.png"
summarized_at: "2026-02-06T10:45:00Z"
method: abstractive
word_count_source: null
word_count_summary: 85
compression_ratio: null
confidence: high
confidence_notes: "High resolution screenshot with all UI elements clearly visible"
---
```

**Summary**: Login page for application "ExampleApp" showing username and password fields, "Remember me" checkbox, and "Sign In" button. Blue banner at top displays "Welcome back!" message.

**What Was Found**:

- Application name "ExampleApp" in header logo
- Two input fields labeled "Username" and "Password"
- Checkbox labeled "Remember me"
- Blue primary button labeled "Sign In"
- Link below button: "Forgot password?"
- Welcome message banner with blue background

**What Was NOT Found**:

- No visible error messages or validation warnings
- No "Sign Up" or registration option visible

**Uncertain**: N/A

**Sources**:

- Screenshot image /home/user/screenshot.png (read 2026-02-06 at 10:45 UTC)

### Architecture Diagram Summary

```yaml
---
source_type: image
source_path: "/home/user/architecture.svg"
summarized_at: "2026-02-06T11:00:00Z"
method: abstractive
word_count_source: null
word_count_summary: 120
compression_ratio: null
confidence: medium
confidence_notes: "Most components labeled clearly, one connection label partially obscured"
---
```

**Summary**: Three-tier architecture diagram showing web frontend, API gateway, and database layer with message queue for async processing.

**What Was Found**:

- Frontend: "React Web App" component
- Middle tier: "API Gateway" labeled component with bidirectional arrow to frontend
- Backend: "PostgreSQL Database" with arrow from API Gateway
- "Redis Cache" component connected to API Gateway
- "Message Queue" component with arrow to "Background Worker"
- Data flow: Frontend → API Gateway → Database and Cache

**What Was NOT Found**:

- No authentication or security layer shown
- No load balancer component visible

**Uncertain**:

- Connection between API Gateway and Message Queue has partially obscured label (visible text: "async...jobs")
- Background Worker box shows partial text, appears truncated on right edge

**Sources**:

- Architecture diagram /home/user/architecture.svg (read 2026-02-06 at 11:00 UTC)
- SVG text extraction from same file
