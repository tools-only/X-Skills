---
name: Setup Instructions
source: https://raw.githubusercontent.com/Chat2AnyLLM/awesome-claude-skills/main/SETUP.md
original_path: SETUP.md
source_repo: Chat2AnyLLM/awesome-claude-skills
category: automation
subcategory: workflow
tags: ['automation']
collected_at: 2026-01-31T18:34:05.979346
file_hash: f916f107c2c70dfd2c5942133b3e644557e1d9b801124c4c6caa39bb04d0c0b9
---

# Setup Instructions

## Automated README Updates

This repository uses GitHub Actions to automatically update the README with the latest skills from configured marketplaces.

### Workflow Configuration

The workflow is configured to use `GITHUB_TOKEN` for all operations (checkout and push). For this to work in public repositories, you may need to adjust repository settings:

1. Go to Repository Settings → Actions → General
2. Under "Workflow permissions", select **"Read and write permissions"**
3. This allows `GITHUB_TOKEN` to push changes to the repository

### Workflow Behavior

The workflow runs hourly and will:
- Fetch the latest skills from configured marketplaces
- Update the README.md if changes are detected
- Commit and push the changes automatically

### Troubleshooting

If the workflow fails with permission errors:
- Ensure workflow permissions are set to "Read and write permissions" in repository settings
- Check that the repository allows GitHub Actions to create commits
- For public repositories, verify that branch protection rules don't block automated commits

### Manual Testing

You can trigger the workflow manually from the Actions tab or by dispatching it via the GitHub CLI:

```bash
gh workflow run "Update Skills README"
```

### Alternative: Personal Access Token (Optional)

If you prefer not to change repository permissions, you can still set up a PAT as a fallback:

1. Create a Personal Access Token with `repo` or `public_repo` scope
2. Add it as `SKILL_UPDATE_TOKEN` repository secret
3. The workflow will automatically use it if available

However, this should not be necessary with proper repository permission settings.