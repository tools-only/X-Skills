# Session Log: Create comparison-table-generator Skill

**Date:** 2025-12-07
**Task:** Create a new skill called `comparison-table-generator` based on the linux-distro-comparison MicroSim template

## Summary

Created a new Claude Code skill that generates interactive comparison table MicroSims for educational content. The skill creates side-by-side comparisons with star ratings, difficulty badges, logos, hover tooltips, and updates mkdocs.yml navigation.

## Input Sources

### Template MicroSim
- **Location:** `/Users/dan/Documents/ws/learning-linux/docs/sims/linux-distro-comparison/`
- **Files analyzed:**
  - `index.md` - Documentation structure
  - `main.html` - HTML template with comprehensive comments
  - `style.css` - CSS with star ratings, badges, tooltips, responsive design
  - `logos/` - SVG logo files (debian.svg, ubuntu.svg, fedora.svg, arch.svg, mint.svg)

### Standardization Rules
- **Location:** `/Users/dan/Documents/ws/claude-skills/skills/microsim-standardization/`
- **Key files referenced:**
  - `SKILL.md` - Quality checklist and standardization workflow
  - `assets/index-template.md` - Documentation template
  - `assets/metadata-schema.json` - Dublin Core JSON schema
  - `assets/metadata-template.json` - Sample metadata

## Steps Performed

1. **Read template files** - Analyzed the linux-distro-comparison MicroSim structure
2. **Read microsim-standardization rules** - Understood quality requirements
3. **Initialized skill** - Used `init_skill.py` script from skill-creator
4. **Removed unneeded directories** - Deleted `scripts/` and `references/` (not needed for this skill)
5. **Created SKILL.md** - 9-step workflow with comprehensive instructions
6. **Created asset templates:**
   - `main-template.html` - HTML template with placeholders and documentation
   - `style-template.css` - Complete CSS from original template
   - `index-template.md` - Documentation template with placeholders
7. **Copied metadata files** - Schema and template from microsim-standardization
8. **Validated and packaged** - Used `package_skill.py` to verify skill structure

## Created Files

### Skill Directory Structure
```
/Users/dan/Documents/ws/claude-skills/skills/comparison-table-generator/
├── SKILL.md                      # 8,769 bytes - Main skill definition
└── assets/
    ├── main-template.html        # 7,398 bytes - HTML template
    ├── style-template.css        # 12,268 bytes - CSS template
    ├── index-template.md         # 1,300 bytes - Documentation template
    ├── metadata-schema.json      # 4,166 bytes - JSON schema
    ├── metadata-template.json    # 967 bytes - Sample metadata
    └── logos/                    # Empty directory for logo files
```

### SKILL.md Description
```
This skill generates interactive comparison table MicroSims for educational
content. Use this skill when users need to create side-by-side comparisons
of items with star ratings (1-5 scale), difficulty badges (Easy/Medium/Hard),
logos, hover tooltips, and description columns. The skill creates a complete
MicroSim package with HTML, CSS, logos directory, index.md documentation,
and metadata.json, then updates mkdocs.yml navigation.
```

### Workflow Steps in SKILL.md
1. **Gather Requirements** - Table title, items, ratings, logos, badges
2. **Create Directory Structure** - `docs/sims/[microsim-name]/`
3. **Generate main.html** - From template with star ratings and tooltips
4. **Generate style.css** - Copy template CSS
5. **Create index.md** - Documentation with YAML frontmatter
6. **Create metadata.json** - Dublin Core metadata
7. **Add Logo Files** - SVG files in logos/ subdirectory
8. **Update mkdocs.yml Navigation** - Add entry to nav section
9. **Validate and Report** - Verify files and suggest preview

## Key Features of Generated MicroSims

### Star Ratings
- 5-color scale: green (5) → yellow-green (4) → orange (3) → red-orange (2) → red (1)
- Unicode star character (★) for compatibility
- Separate spans for filled and empty stars

### Difficulty Badges
- Easy: Green background (#dcfce7), dark green text (#166534)
- Medium: Yellow background (#fef3c7), dark orange text (#92400e)
- Hard: Red background (#fee2e2), dark red text (#991b1b)

### Hover Tooltips
- Pure CSS implementation using ::after pseudo-elements
- First row tooltip appears BELOW to avoid header overlap
- Smooth fade transitions with opacity/visibility

### Responsive Design
- Mobile breakpoint at 700px
- Horizontal scrolling for narrow screens
- Adjusted font sizes and logo dimensions

## Validation Result

```
📦 Packaging skill: comparison-table-generator

🔍 Validating skill...
✅ Skill is valid!

  Added: comparison-table-generator/SKILL.md
  Added: comparison-table-generator/assets/index-template.md
  Added: comparison-table-generator/assets/metadata-template.json
  Added: comparison-table-generator/assets/main-template.html
  Added: comparison-table-generator/assets/metadata-schema.json
  Added: comparison-table-generator/assets/style-template.css

✅ Successfully packaged skill to: comparison-table-generator.zip
```

## Usage

To use this skill, invoke it when you need to create a comparison table MicroSim:

```
/skill comparison-table-generator
```

The skill will guide you through:
1. Defining what items to compare
2. Setting up rating criteria
3. Providing logos
4. Generating all required files
5. Updating mkdocs.yml navigation

## Notes

- The skill conforms to all microsim-standardization rules
- Templates include comprehensive HTML/CSS comments explaining each section
- Iframe height should be calculated: ~60px per row + 150px for header/legend
- First row tooltip fix is built into the CSS template
- Logo files should be 32x32px SVG format for consistency
