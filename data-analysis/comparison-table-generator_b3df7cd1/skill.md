# Comparison Table Generator

## Overview

The **Comparison Table Generator** skill creates interactive comparison table MicroSims for educational content. These tables allow students to compare multiple items across several criteria using color-coded star ratings, difficulty badges, logos, and hover tooltips. The skill generates complete MicroSim packages following all microsim-standardization rules, including proper documentation and mkdocs.yml navigation updates.

This skill uses pure CSS (no JavaScript library required) for lightweight, fast-loading interactive tables.

## When to Use This Skill

Use the Comparison Table Generator when you need to create:

- **Educational Comparisons** - Linux distributions, programming languages, frameworks, tools
- **Decision-Making Aids** - Help students choose between multiple options with clear criteria
- **Multi-Criteria Analysis** - Compare items across 2-4 rating dimensions
- **Visual Summaries** - Quick-reference tables with ratings and categories
- **Interactive Learning** - Tables with hover tooltips for deeper explanations

## Key Features

### Star Ratings (1-5 Scale)

Color-coded ratings using Unicode stars:

- **5 stars** (Green #22c55e) - Excellent
- **4 stars** (Yellow-green #84cc16) - Very Good
- **3 stars** (Orange #f59e0b) - Good/Average
- **2 stars** (Red-orange #f97316) - Below Average
- **1 star** (Red #ef4444) - Poor

### Difficulty/Category Badges

Pill-shaped badges with semantic colors:

- **Easy** - Green background, dark green text
- **Medium** - Yellow background, dark orange text
- **Hard** - Red background, dark red text

Custom categories can be added following the same pattern.

### Hover Tooltips

Pure CSS tooltips that appear on row hover:

- Smooth fade transitions
- First row tooltip appears BELOW to avoid header overlap
- Maximum 400px width for readability
- Prefixed with item name for context

### Responsive Design

- Mobile breakpoint at 700px
- Horizontal scrolling for narrow screens
- Adjusted font sizes and logo dimensions

### Logo Support

- SVG format recommended (32x32px)
- Consistent row heights with flexbox alignment
- Stored in `logos/` subdirectory

## Workflow

The skill follows a 9-step process:

1. **Gather Requirements** - Table title, items, ratings, logos, badges
2. **Create Directory Structure** - `docs/sims/[microsim-name]/`
3. **Generate main.html** - From template with star ratings and tooltips
4. **Generate style.css** - Copy template CSS with all patterns
5. **Create index.md** - Documentation with YAML frontmatter
6. **Create metadata.json** - Dublin Core metadata
7. **Add Logo Files** - SVG files in logos/ subdirectory
8. **Update mkdocs.yml Navigation** - Add entry to nav section
9. **Validate and Report** - Verify files and suggest preview

## File Structure

Each generated comparison table creates:

```
docs/sims/[table-name]/
├── index.md          # Documentation with iframe embed
├── main.html         # Interactive comparison table
├── style.css         # Styling (star colors, badges, tooltips)
├── metadata.json     # Dublin Core metadata
└── logos/            # SVG logo files
    ├── item1.svg
    ├── item2.svg
    └── ...
```

## Example Use Cases

### 1. Linux Distribution Comparison

**Scenario**: Help students choose the right Linux distribution.

**Columns**:
- Distribution (logo + name)
- User Friendly (1-5 stars)
- Stability (1-5 stars)
- Fresh Software (1-5 stars)
- Learning Curve (Easy/Medium/Hard badge)
- Best For (text description)

**Items**: Debian, Ubuntu, Fedora, Arch Linux, Linux Mint

### 2. Programming Language Comparison

**Scenario**: Compare languages for a beginner's course.

**Columns**:
- Language (logo + name)
- Ease of Learning (1-5 stars)
- Job Market (1-5 stars)
- Performance (1-5 stars)
- Difficulty (Easy/Medium/Hard)
- Best For (text)

**Items**: Python, JavaScript, Java, C++, Go

### 3. Database Comparison

**Scenario**: Help students choose a database technology.

**Columns**:
- Database (logo + name)
- Query Flexibility (1-5 stars)
- Scalability (1-5 stars)
- Learning Curve (Easy/Medium/Hard)
- Best For (text)

**Items**: PostgreSQL, MySQL, MongoDB, Redis, SQLite

### 4. Cloud Provider Comparison

**Scenario**: Compare major cloud platforms.

**Columns**:
- Provider (logo + name)
- Ease of Use (1-5 stars)
- Feature Breadth (1-5 stars)
- Documentation (1-5 stars)
- Cost (Low/Medium/High)
- Best For (text)

**Items**: AWS, Google Cloud, Azure, DigitalOcean

## Technical Details

### HTML Structure

```html
<tr data-tooltip="Item Name: Description for tooltip">
    <td class="distro-cell">
        <img src="logos/item.svg" alt="Item" class="distro-logo">
        <span class="distro-name">Item Name</span>
    </td>
    <td class="rating">
        <span class="stars stars-4">★★★★</span>
        <span class="stars-empty">★</span>
    </td>
    <td class="difficulty medium">Medium</td>
    <td class="best-for">Description text</td>
</tr>
```

### Star Rating Pattern

```html
<!-- 4 out of 5 stars -->
<td class="rating">
    <span class="stars stars-4">★★★★</span>
    <span class="stars-empty">★</span>
</td>
```

### Difficulty Badge Pattern

```html
<td class="difficulty easy">Easy</td>
<td class="difficulty medium">Medium</td>
<td class="difficulty hard">Hard</td>
```

### First Row Tooltip Fix

The CSS automatically positions the first row's tooltip below (instead of above) to avoid being hidden by the table header:

```css
.comparison-table tbody tr:first-child[data-tooltip]::after {
    bottom: auto;
    top: calc(100% + 10px);
}
```

## Customization

### Custom Badge Categories

To add categories beyond Easy/Medium/Hard:

```css
.difficulty.beginner {
    background-color: #dbeafe;  /* Light blue */
    color: #1e40af;             /* Dark blue */
}

.difficulty.advanced {
    background-color: #fae8ff;  /* Light purple */
    color: #86198f;             /* Dark purple */
}
```

### Custom Star Colors

Modify in style.css:

```css
.stars-5 { color: #22c55e; }  /* Green */
.stars-4 { color: #84cc16; }  /* Yellow-green */
.stars-3 { color: #f59e0b; }  /* Orange */
.stars-2 { color: #f97316; }  /* Red-orange */
.stars-1 { color: #ef4444; }  /* Red */
```

### Iframe Height Calculation

Approximately 60px per row + 150px for header/legend:

- 3 rows: ~330px
- 5 rows: ~450px (default 470px)
- 8 rows: ~630px
- 10 rows: ~750px

## Educational Framework Integration

### Bloom's Taxonomy Alignment

Comparison tables support multiple cognitive levels:

- **Remember** - Identify items and their characteristics
- **Understand** - Interpret rating scales and categories
- **Apply** - Use comparisons to make decisions
- **Analyze** - Compare trade-offs between options
- **Evaluate** - Assess which option best fits specific needs

### Learning Objectives Template

```markdown
After reviewing this comparison, students should be able to:

1. **Identify** the major [items] and their characteristics
2. **Compare** [items] based on different criteria
3. **Evaluate** which [item] is best for specific use cases
4. **Analyze** the trade-offs between [criteria]
```

## Best Practices

### Content Guidelines

1. **Limit items** - 3-8 items per table for readability
2. **Limit rating columns** - 2-4 rating criteria maximum
3. **Consistent tooltips** - Always prefix with item name
4. **Clear badges** - Use semantic colors (green=easy, red=hard)
5. **Meaningful descriptions** - "Best For" column should be actionable

### Logo Guidelines

1. **SVG format** - Scalable, small file size
2. **32x32px** - Consistent dimensions
3. **Simple designs** - Clear at small sizes
4. **Official sources** - Check brand guidelines or vectorlogo.zone

### Accessibility

1. **Semantic colors** - Green=positive, Red=challenging
2. **Sufficient contrast** - WCAG AA compliant
3. **Clear labels** - Descriptive column headers
4. **Text alternatives** - Logo alt text, tooltip descriptions

## Comparison with Other Skills

### vs. ChartJS Generator

- **Comparison Table**: Categorical comparisons with ratings
- **ChartJS**: Numerical data visualization (bar, line, pie charts)

**Use Comparison Table for**: Side-by-side feature comparisons
**Use ChartJS for**: Quantitative data trends and distributions

### vs. Radar Chart (ChartJS)

- **Comparison Table**: Multiple items, multiple criteria, badges
- **Radar Chart**: Multi-dimensional profile visualization

**Use Comparison Table for**: Discrete ratings with categories
**Use Radar for**: Continuous metrics with overlapping profiles

### vs. Vis-Network

- **Comparison Table**: Static comparisons in table format
- **Vis-Network**: Interactive network graph diagrams

**Use Comparison Table for**: Linear comparisons
**Use Vis-Network for**: Relationship networks and dependencies

## Troubleshooting

### Tooltip Hidden Behind Header

**Cause**: First row tooltip positioned above table header
**Solution**: CSS automatically handles this - verify `tr:first-child` rules are present

### Row Heights Misaligned

**Cause**: Logos with different aspect ratios
**Solution**: Set explicit height on `.distro-cell` (default 60px)

### Stars Not Colored

**Cause**: Missing star color classes
**Solution**: Ensure `stars-N` class matches the rating value (1-5)

### Badges Not Styled

**Cause**: Class name doesn't match CSS
**Solution**: Use lowercase class names: `easy`, `medium`, `hard`

## Related Skills

- **ChartJS Generator** - Data visualization charts
- **Microsim-P5** - Custom interactive simulations
- **MicroSim Standardization** - Quality validation for MicroSims
- **MicroSim Matcher** - Select the right generator for specifications

## References

- [CSS Flexbox Guide](https://css-tricks.com/snippets/css/a-guide-to-flexbox/) - Layout for logo cells
- [Pure CSS Tooltips](https://css-tricks.com/a-complete-guide-to-tooltip-patterns/) - Tooltip implementation
- [Vector Logo Zone](https://www.vectorlogo.zone/) - Source for SVG logos
- [Dublin Core Metadata](https://www.dublincore.org/specifications/dublin-core/dcmi-terms/) - Metadata standards

## Quick Start Example

To generate a comparison table:

```
Input Requirements:
- Title: Linux Distribution Comparison
- MicroSim Name: linux-distro-comparison
- Items: Debian, Ubuntu, Fedora, Arch, Mint (with logos and tooltips)
- Rating Columns: User Friendly, Stability, Fresh Software
- Badge Column: Learning Curve (Easy/Medium/Hard)
- Description Column: Best For
```

**Output**: Complete MicroSim package with:
- Interactive table with hover tooltips
- Color-coded star ratings
- Difficulty badges with legend
- Responsive design
- Documentation with learning objectives
- Dublin Core metadata
- mkdocs.yml navigation entry

---

*Comparison tables integrate seamlessly into MkDocs Material textbooks as iframe-embedded MicroSims with professional styling and educational metadata.*
