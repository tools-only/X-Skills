---
name: cover-image-generator
description: Guides users through generating book cover images using the generate-cover.sh script. Determines the best workflow based on user's available resources (OpenAI API billing, ChatGPT Pro subscription, macOS) and runs the appropriate command.
---

# Cover Image Generator

This guide helps you generate a professional book cover image using the automated `generate-cover.sh` script. The workflow adapts based on your available resources.

## Prerequisites

Before starting, ensure you have:

1. **Project Structure**: An MkDocs project with:
   - `mkdocs.yml` containing `site_name:` field
   - `docs/course-description.md` with book content description
   - `docs/img/` directory (will be created if missing)

2. **Script Location**: The generate-cover.sh script at:
   ```
   ~/.claude/skills/claude-skills/src/image-generation/generate-cover.sh
   ```
   Or the full path to the claude-skills repository.

## Workflow

### Step 1: Verify Project Requirements

First, check that the required files exist in the user's project:

```bash
# Check for mkdocs.yml
ls mkdocs.yml

# Check for course description
ls docs/course-description.md

# Check/create output directory
mkdir -p docs/img
```

If `docs/course-description.md` is missing, inform the user they need to create it first. This file provides the keywords for generating an appropriate cover image.

### Step 2: Determine User's Resources

Ask the user these questions to determine the best workflow:

**Question 1: OpenAI API Key**
```
Do you have an OpenAI API key set up?
(Check with: echo $OPENAI_API_KEY)

1. Yes, I have an API key with active billing
2. Yes, I have an API key but billing is not active
3. No, I don't have an API key
```

**Question 2: ChatGPT Subscription (if no active API billing)**
```
Do you have a ChatGPT Pro/Plus subscription ($20/month)?

1. Yes
2. No
```

**Question 3: Operating System (if has ChatGPT Pro)**
```
Are you running on macOS?

1. Yes (can use browser automation)
2. No (will use manual copy/paste)
```

### Step 3: Select and Run the Appropriate Command

Based on the answers, use one of these workflows:

#### Path A: Full Auto Mode (API key + active billing)

This is the fastest, fully automated option.

```bash
# Navigate to the user's project root
cd /path/to/project

# Run the script with no flags
~/.claude/skills/claude-skills/src/image-generation/generate-cover.sh
```

The script will:
1. Extract book title from mkdocs.yml
2. Generate an optimized image prompt
3. Call OpenAI Images API
4. Save the result to `docs/img/cover.png`

**Expected output:**
```
=== Cover Image Generator ===
Project directory: /path/to/project
Book Title: Your Book Title
...
Generating cover image...
Wrote: docs/img/cover.png
```

#### Path B: Browser Automation (ChatGPT Pro + macOS)

For users with ChatGPT Pro on macOS.

```bash
# Navigate to the user's project root
cd /path/to/project

# Run with --open-browser flag
~/.claude/skills/claude-skills/src/image-generation/generate-cover.sh --open-browser
```

The script will:
1. Extract book title from mkdocs.yml
2. Generate an optimized image prompt locally (no API call)
3. Open ChatGPT in the default browser
4. Automatically paste the prompt

**Note:** First run may require granting Accessibility permissions:
- System Settings > Privacy & Security > Accessibility
- Allow Terminal (or your terminal app) to control the computer

**After the script runs:**
1. Wait for ChatGPT to generate the image
2. Download the generated image
3. Save it to `docs/img/cover.png`
4. Resize to 1200x630 pixels if needed

#### Path C: Local Prompt (ChatGPT Pro, any OS)

For users with ChatGPT Pro on non-macOS systems.

```bash
# Navigate to the user's project root
cd /path/to/project

# Run with --local-prompt flag
~/.claude/skills/claude-skills/src/image-generation/generate-cover.sh --local-prompt
```

The script will:
1. Extract book title from mkdocs.yml
2. Generate an optimized image prompt locally (no API call)
3. Display the prompt for manual copying

**After the script runs:**
1. Copy the displayed IMAGE PROMPT
2. Go to https://chatgpt.com/
3. Paste the prompt
4. Download the generated image
5. Save it to `docs/img/cover.png`
6. Resize to 1200x630 pixels if needed

#### Path D: No Resources Available

If the user has neither API billing nor ChatGPT Pro:

**Option 1: Set up OpenAI API billing**
1. Go to https://platform.openai.com/account/billing
2. Add a payment method
3. Add credits ($5-10 is sufficient for many images)
4. Return and use Path A

**Option 2: Subscribe to ChatGPT Pro**
1. Go to https://openai.com/chatgpt/pricing
2. Subscribe to Plus ($20/month)
3. Return and use Path B or C

**Option 3: Use the prompt manually with free tier**
1. Run `--local-prompt` to get the prompt
2. Use the prompt with any free AI image generator:
   - Bing Image Creator (free)
   - Leonardo.ai (free tier)
   - Ideogram (free tier)

### Step 4: Verify the Cover Image

After generation, verify the image:

```bash
# Check file exists
ls -la docs/img/cover.png

# Check dimensions (if ImageMagick installed)
identify docs/img/cover.png
```

**Required specifications:**
- Format: PNG
- Size: 1200x630 pixels (1.91:1 aspect ratio)
- This is the Open Graph standard for social media previews

### Step 5: Update Home Page (Optional)

If the user wants to display the cover on their home page, add to `docs/index.md`:

```markdown
---
title: Your Book Title
description: Brief description
image: /img/cover.png
og:image: /img/cover.png
---

![Your Book Title](./img/cover.png){ width="100%" }
```

## Troubleshooting

### "billing_not_active" Error

Your OpenAI API key is valid but the account lacks billing.

**Solutions:**
- Add payment method at https://platform.openai.com/account/billing
- Or use `--local-prompt` with ChatGPT Pro

### "Could not extract site_name from mkdocs.yml"

The `mkdocs.yml` file is missing the `site_name:` field.

**Fix:**
```yaml
# Add to mkdocs.yml
site_name: Your Book Title
```

### "Course description not found"

The file `docs/course-description.md` doesn't exist.

**Fix:**
Create the file with a description of your book's content, topics, and themes. This provides keywords for the image generation.

### Auto-paste not working on macOS

First run requires Accessibility permissions.

**Fix:**
1. Open System Settings
2. Go to Privacy & Security > Accessibility
3. Enable your terminal app (Terminal, iTerm2, etc.)
4. Run the script again

### Image is wrong dimensions

ChatGPT may not generate exact dimensions.

**Fix:**
1. Open the image in an editor (Preview on macOS, GIMP, etc.)
2. Resize to 1200x630 pixels
3. Save as PNG

## Command Reference

| Command | Description | Requirements |
|---------|-------------|--------------|
| `generate-cover.sh` | Full auto via API | API key + billing |
| `generate-cover.sh --local-prompt` | Generate prompt only | None |
| `generate-cover.sh --open-browser` | Open ChatGPT + paste | macOS + ChatGPT Pro |
| `generate-cover.sh --prompt-only` | API prompt only | API key + billing |

## Related Resources

- [Cover Image Workflow Diagram](https://dmccreary.github.io/claude-skills/sims/cover-image-workflow/)
- [Home Page Template Guide](./home-page-template.md)
- [Image Generation README](https://github.com/dmccreary/claude-skills/tree/main/src/image-generation)
