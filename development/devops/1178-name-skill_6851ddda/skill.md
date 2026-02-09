---
name: obsidian-deploy-integration
description: |
  Publish Obsidian plugins to the community plugin directory.
  Use when releasing your first plugin, updating existing plugins,
  or managing the community plugin submission process.
  Trigger with phrases like "publish obsidian plugin", "obsidian community plugins",
  "submit obsidian plugin", "obsidian plugin directory".
allowed-tools: Read, Write, Edit, Bash(git:*), Bash(gh:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Obsidian Deploy Integration

## Overview
Complete guide for publishing Obsidian plugins to the official community plugin directory.

## Prerequisites
- Completed and tested plugin
- GitHub account
- GitHub repository for your plugin
- `gh` CLI installed (optional but recommended)

## Community Plugin Directory Structure

### Your Repository Requirements
```
your-plugin-repo/
├── main.js              # Built plugin (required)
├── manifest.json        # Plugin metadata (required)
├── styles.css           # Optional styles
├── README.md            # Documentation
├── LICENSE              # License file
└── .github/
    └── workflows/       # CI/CD workflows
```

### obsidian-releases Repository
```
obsidian-releases/
├── community-plugins.json    # Plugin registry
└── community-css-themes.json # Theme registry
```

## Instructions

### Step 1: Prepare Your Repository
```bash
# Ensure all required files are in repo root
ls -la main.js manifest.json

# Validate manifest.json
cat manifest.json
# Should have: id, name, version, minAppVersion, description, author

# Create GitHub release
gh release create 1.0.0 \
  main.js \
  manifest.json \
  styles.css \
  --title "v1.0.0" \
  --notes "Initial release"
```

### Step 2: Fork obsidian-releases
```bash
# Fork the releases repository
gh repo fork obsidianmd/obsidian-releases --clone

# Navigate to fork
cd obsidian-releases

# Create branch for your submission
git checkout -b add-my-plugin
```

### Step 3: Add Plugin Entry
```bash
# Edit community-plugins.json
# Add your plugin in alphabetical order by id

# Example entry to add:
{
  "id": "my-plugin-id",
  "name": "My Plugin Name",
  "author": "Your Name",
  "description": "Brief description of what your plugin does",
  "repo": "your-username/your-repo-name"
}
```

### Step 4: Validate Entry
```bash
# Validate JSON syntax
cat community-plugins.json | jq empty

# Check your entry
cat community-plugins.json | jq '.[] | select(.id == "my-plugin-id")'

# Verify alphabetical ordering
cat community-plugins.json | jq '.[].id' | sort -c
```

### Step 5: Submit Pull Request
```bash
# Commit changes
git add community-plugins.json
git commit -m "Add my-plugin-id to community plugins"

# Push to your fork
git push origin add-my-plugin

# Create PR
gh pr create \
  --repo obsidianmd/obsidian-releases \
  --title "Add my-plugin-id" \
  --body "## Plugin Submission

**Plugin ID:** my-plugin-id
**Plugin Name:** My Plugin Name
**Repository:** https://github.com/your-username/your-repo-name

### Description
Brief description of what the plugin does and why it's useful.

### Checklist
- [x] Plugin works on desktop (Windows, macOS, Linux)
- [x] Plugin works on mobile (iOS, Android) or is marked desktop-only
- [x] README.md explains installation and usage
- [x] No console errors during normal operation
- [x] Follows Obsidian plugin guidelines

### Screenshots
[Include screenshots if visual]
"
```

### Step 6: Post-Submission Process
```markdown
## After Submitting

### Review Process
1. **Automated checks** run within minutes
2. **Manual review** by Obsidian team (1-7 days typically)
3. **Feedback** may be provided via PR comments
4. **Approval** leads to automatic publication

### Common Review Feedback
- Missing README documentation
- Console.log statements in production code
- Hardcoded API keys or sensitive data
- Missing error handling
- Incompatibility with mobile

### Responding to Feedback
1. Make requested changes in your plugin repo
2. Create new release with fixes
3. Comment on PR when ready for re-review
```

### Step 7: Post-Publication Maintenance
```bash
# After approval, your plugin is live!

# For updates, just create new GitHub releases
gh release create 1.0.1 \
  main.js manifest.json styles.css \
  --title "v1.0.1" \
  --notes "Bug fixes and improvements"

# Obsidian automatically picks up new releases
# Users will see update notifications

# For major updates, update community-plugins.json if needed
# (usually only for name/description changes)
```

## Output
- Plugin repository properly structured
- GitHub release created
- Pull request submitted to obsidian-releases
- Plugin published to community directory

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| PR validation fails | JSON syntax error | Validate with `jq` |
| Plugin not found | Wrong repo path | Check repo URL exactly matches |
| Release not picked up | No GitHub release | Create release with required files |
| Mobile crash | Desktop-only API | Set `isDesktopOnly: true` |

## Examples

### Updating Plugin Description
```bash
# If you need to update name/description after publication
# Edit community-plugins.json in your fork

# Create PR with changes
gh pr create \
  --repo obsidianmd/obsidian-releases \
  --title "Update my-plugin-id description" \
  --body "Updating plugin description to better reflect new features."
```

### Complete Submission Script
```bash
#!/bin/bash
# submit-plugin.sh

PLUGIN_ID="my-plugin"
PLUGIN_NAME="My Plugin"
AUTHOR="Your Name"
DESCRIPTION="Plugin description"
REPO="username/repo-name"

# Fork and clone obsidian-releases
gh repo fork obsidianmd/obsidian-releases --clone
cd obsidian-releases

# Create branch
git checkout -b add-$PLUGIN_ID

# Add entry (using jq to maintain proper format)
jq --arg id "$PLUGIN_ID" \
   --arg name "$PLUGIN_NAME" \
   --arg author "$AUTHOR" \
   --arg desc "$DESCRIPTION" \
   --arg repo "$REPO" \
   '. += [{"id": $id, "name": $name, "author": $author, "description": $desc, "repo": $repo}] | sort_by(.id)' \
   community-plugins.json > tmp.json && mv tmp.json community-plugins.json

# Commit and push
git add community-plugins.json
git commit -m "Add $PLUGIN_ID"
git push origin add-$PLUGIN_ID

# Create PR
gh pr create --repo obsidianmd/obsidian-releases \
  --title "Add $PLUGIN_ID" \
  --body "Plugin submission for $PLUGIN_NAME"

echo "PR submitted! Check: https://github.com/obsidianmd/obsidian-releases/pulls"
```

### BRAT Beta Distribution
```markdown
## Pre-Release Distribution via BRAT

For beta testing before community submission:

1. User installs BRAT plugin
2. BRAT Settings > Add Beta Plugin
3. Enter: `your-username/your-repo-name`
4. BRAT installs directly from your repo

No community-plugins.json entry needed for BRAT!
```

## Resources
- [Plugin Submission Guidelines](https://docs.obsidian.md/Plugins/Releasing/Submit+your+plugin)
- [obsidian-releases Repository](https://github.com/obsidianmd/obsidian-releases)
- [Plugin Guidelines](https://docs.obsidian.md/Plugins/Releasing/Plugin+guidelines)
- [BRAT Plugin](https://github.com/TfTHacker/obsidian42-brat)

## Next Steps
For event handling patterns, see `obsidian-webhooks-events`.
