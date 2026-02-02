# Causal Loop MicroSim Generator - Skill Documentation Session

**Date:** 2025-11-27
**Task:** Create detailed skill description documentation for the causal-loop-microsim-generator skill

## Summary

Created comprehensive documentation for the `causal-loop-microsim-generator` skill, which generates interactive Causal Loop Diagram (CLD) MicroSims using the vis-network JavaScript library for systems thinking education.

## Files Created

### 1. Detailed Skill Description
**Path:** `docs/skill-descriptions/microsims/causal-loop-microsim-generator.md`

A 360+ line comprehensive documentation file covering:

- **Overview**: Purpose and capabilities of the skill
- **Key Features**: Interactive CLDs, polarity indicators (+/-), loop markers (R/B), systems archetypes
- **When to Use**: Trigger phrases and use cases
- **MicroSim Architecture**: Folder structure and output files
- **CLD Concepts**: Causal relationships, feedback loop types, loop determination rules
- **Systems Archetypes**: Table of common patterns (limits-to-growth, fixes-that-fail, etc.)
- **JSON Data Schema**: Complete schema for nodes, edges, loops, and educational content
- **Node Positioning Guidelines**: Layout patterns for 3, 4, and 5-node loops
- **Visual Configuration**: Polarity colors, loop indicators, edge curves
- **Interactive Features**: Click events, URL parameters
- **Educational Content Features**: Bloom's Taxonomy alignment, discussion questions
- **Best Practices**: Diagram design, educational value, accessibility
- **Output Files**: index.md, main.html, JavaScript, data.json, style.css
- **Integration**: Related skills (microsim-p5, vis-network, learning-graph-generator)
- **Troubleshooting**: Common issues and solutions
- **References**: vis-network docs, systems thinking literature

## Files Updated

### 2. MicroSims Index Page
**Path:** `docs/skill-descriptions/microsims/index.md`

Changes made:

1. **Added to numbered generator list** (now item 5):
   ```
   5. **Causal Loop Diagram Generator** - creates interactive causal loop diagrams for systems thinking education
   ```

2. **Added detailed section** after Vis-Network section with:
   - Skill name: `causal-loop-microsim-generator`
   - Width responsive: Yes
   - Key features bullet list
   - Description of use cases
   - Link to full documentation

3. **Renumbered items** 6-10 (Math Function Plotter through Bubble Chart Generator)

## Source Files Analyzed

The following skill files were read to understand the skill's capabilities:

| File | Purpose |
|------|---------|
| `skills/causal-loop-microsim-generator/SKILL.md` | Main skill definition and workflow |
| `skills/causal-loop-microsim-generator/assets/rules.md` | CLD generation rules and JSON schema |
| `skills/causal-loop-microsim-generator/assets/templates/microsim.js` | JavaScript template for vis-network |
| `skills/causal-loop-microsim-generator/assets/templates/main.html` | HTML template structure |
| `skills/causal-loop-microsim-generator/assets/templates/data.json` | Example JSON data schema |
| `docs/skill-descriptions/microsims/vis-network.md` | Reference for documentation style |

## Skill Capabilities Documented

The causal-loop-microsim-generator skill:

1. **Gathers requirements** from user (name, title, nodes, edges, loops)
2. **Generates 5 files** per MicroSim:
   - `index.md` - Documentation with iframe embed
   - `main.html` - HTML container with vis-network CDN
   - `[name].js` - JavaScript for CLD rendering
   - `data.json` - Node, edge, loop definitions
   - `style.css` - Layout and legend styling
3. **Supports systems archetypes**: limits-to-growth, fixes-that-fail, shifting-the-burden, success-to-the-successful, tragedy-of-the-commons, escalation, drifting-goals
4. **Provides educational content**: Learning objectives, discussion questions, key insights, common misconceptions
5. **Updates mkdocs.yml** navigation in alphabetical order

## Technical Details

- **Library:** vis-network.js (loaded from CDN)
- **Canvas:** 600x600 pixels standard
- **Polarity Colors:** Green (#28a745) for positive, Red (#dc3545) for negative
- **Loop Indicators:** Red ellipse for Reinforcing (R), Green ellipse for Balancing (B)
- **Iframe Height:** 500px default (width responsive)

## Next Steps

- Take screenshots of example CLDs for documentation
- Update `mkdocs.yml` if new documentation pages need navigation entries
- Test iframe embedding in MkDocs Material theme
