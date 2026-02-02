# Install MkDocs Template

The install-mkdocs-template skill creates a complete MkDocs Material project
structure optimized for intelligent textbooks. It sets up a Conda virtual
environment, installs dependencies, generates all configuration files, and
deploys to GitHub Pages.

## Key Capabilities

This skill provides end-to-end project setup:

1. **Conda Environment** - Creates `mkdocs` virtual environment with Python 3.11
2. **Dependencies** - Installs MkDocs, Material theme, and required packages
3. **Complete mkdocs.yml** - All Material theme options pre-configured
4. **Custom CSS** - Brand color customization with CSS variables
5. **Social Override Plugin** - Per-page custom social media card images
6. **Build & Deploy** - Builds site and deploys to GitHub Pages

## What Gets Created

```
project-root/
├── docs/
│   ├── css/extra.css           # Brand colors and custom styles
│   ├── img/                    # Logo and favicon placeholders
│   ├── chapters/index.md       # Chapter section starter
│   ├── learning-graph/index.md # Learning graph section starter
│   ├── sims/index.md           # MicroSims section starter
│   └── index.md                # Home page template
├── plugins/
│   ├── __init__.py
│   └── social_override.py      # Custom social media plugin
├── mkdocs.yml                  # Full configuration
└── setup.py                    # Plugin installation
```

## MkDocs Material Features Included

### Navigation Features
- Expandable sections, breadcrumbs, section indexes
- Back to top button, footer navigation
- Table of contents that follows scroll

### Content Features
- Code copy button with syntax highlighting
- Edit on GitHub button
- Mermaid diagram support

### Markdown Extensions
- Admonitions (note, warning, tip, etc.)
- Collapsible details blocks
- Content tabs
- Math support (LaTeX via MathJax)
- Task lists, emoji, and more

### Plugins
- Full-text search with suggestions
- Social media card generation
- Custom social card override per page

## Required Information

When running this skill, provide:

1. **site_name** - Textbook title
2. **site_description** - Brief description for SEO
3. **site_author** - Author name(s)
4. **site_url** - Deployment URL (e.g., https://username.github.io/repo/)
5. **repo_url** - GitHub repository URL
6. **primary_color_rgb** - Brand color (default: 218, 120, 87)
7. **google_analytics_id** - Optional analytics property

## Social Override Plugin

The included plugin allows custom social media images per page:

```markdown
---
title: My Page Title
image: img/my-custom-social-card.png
---
```

This overrides the auto-generated social card for that specific page.

## Prerequisites

- Conda (Miniconda or Anaconda) installed
- Git repository initialized with remote origin configured
- GitHub repository created

## Workflow Summary

The skill executes these steps:

1. Create Conda environment: `conda create -n mkdocs python=3.11 -y`
2. Activate and install: `conda activate mkdocs && pip install mkdocs mkdocs-material ...`
3. Create all project files (mkdocs.yml, CSS, plugins, starter content)
4. Install social_override plugin: `pip install -e .`
5. Build site: `mkdocs build`
6. Deploy to GitHub Pages: `mkdocs gh-deploy`
7. Provide the live site URL

## After Deployment

The skill provides the GitHub Pages URL for testing:

```
https://<username>.github.io/<repo-name>/
```

User should:

1. Add logo to `docs/img/logo.png` (50x50px recommended)
2. Add favicon to `docs/img/favicon.ico`
3. Verify the live site loads correctly

## Reactivating the Environment

When returning to work on the textbook:

```bash
conda activate mkdocs
```

## Integration

This skill is typically the **first step** in the intelligent textbook creation
workflow. After setting up the MkDocs structure, proceed with:

1. Course Description Analyzer - Create/validate course description
2. Learning Graph Generator - Generate concept dependencies
3. Book Chapter Generator - Design chapter structure
4. Continue with content generation skills
