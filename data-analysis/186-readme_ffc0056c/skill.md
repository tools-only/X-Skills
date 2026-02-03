# README Generator

The readme-generator skill creates or updates comprehensive README.md files for
GitHub repositories following best practices. It generates all essential sections
including badges, project overview, site metrics, and getting started instructions.

## Key Capabilities

This skill generates README.md files with:

- **Technology Badges**: Python, MkDocs, Material theme, etc.
- **Project Overview**: Description, purpose, key features
- **Site Metrics**: From book-metrics.md if available
- **Getting Started**: Installation and setup instructions
- **Project Structure**: Directory layout explanation
- **Contributing Guidelines**: How to contribute
- **License and Contact**: Standard footer sections

## When to Use

Use this skill when:

- Starting a new GitHub repository that needs a README
- Updating an existing README to follow best practices
- After significant project changes that need documentation
- Before publishing or sharing a repository
- Migrating from another documentation system
- After adding new technologies or dependencies

## Workflow

The skill follows these steps:

1. **Analyze Repository**: Check for existing README, identify technologies
2. **Gather Metadata**: Read mkdocs.yml for site info, check for docs/
3. **Generate Badges**: Create shields.io badges for all technologies
4. **Extract Metrics**: Pull statistics from book-metrics.md if available
5. **Write README**: Generate complete markdown with all sections

## Prerequisites

For best results, the repository should have:

- `mkdocs.yml` with site configuration
- `docs/` directory with content
- Optional: `docs/learning-graph/book-metrics.md` for statistics

## Customization

The skill will prompt for information if not found:

- Repository URL if not in .git/config
- GitHub Pages URL if not configured
- Project description if not in mkdocs.yml

## Output Sections

1. Title with badges
2. Overview/Description
3. Quick Links (site, repo, docs)
4. Key Metrics (if available)
5. Features list
6. Getting Started
7. Project Structure
8. Contributing
9. License
10. Contact/Author

## Integration

This skill is typically used when initializing a new project or before
major releases to ensure the repository has professional documentation.
