---
name: readme-maker
description: >
  This skill should be used when the user asks to "create a README",
  "generate a README", "make a readme", "write a README for my project",
  "need a README", "add a README", "document my project", "set up project docs",
  "readme with badges".
allowed-tools: AskUserQuestion, Write, Read, Glob
---

# README Maker

Generate a professional `README.md` for any project — from 50-line minimal to
600-line comprehensive — with shields.io badges, styled headers, and structured
sections. Interview the user to determine preferences, then produce the full
README in one pass.

## Step 1 — Detect Project Context

Before asking, scan the working directory for context clues:

- Glob for `package.json`, `pyproject.toml`, `Cargo.toml`, `go.mod` — Read the
  first match to extract project name and version as pre-filled defaults
- Check if `README.md` already exists — warn the user before overwriting
- Note the primary language from the detected manifest file
- If no manifest is found, skip pre-filling — rely entirely on user input in Step 2

## Step 2 — Round 1 Interview (Always)

Ask 4 questions via `AskUserQuestion`. Pre-select options that match detected
context where possible.

| # | header | options |
|---|--------|---------|
| 1 | "Project type" | Library/Package · CLI Tool · Web App/API · Desktop App |
| 2 | "Language" | Python · JavaScript/Node · Go · Rust |
| 3 | "Depth" | Minimal · Standard · Comprehensive |
| 4 | "License" | MIT · Apache-2.0 · GPL-3.0 · No license |

Depth definitions (include as option descriptions):
- **Minimal** — small utilities, internal tools, or anything with one clear use case
- **Standard** — open-source projects expecting contributors and external users
- **Comprehensive** — mature projects with public APIs, active maintainers, or large communities

If depth = **Minimal**, skip Step 3 and proceed directly to Step 4.

## Step 3 — Round 2 Interview (Standard / Comprehensive Only)

Ask 4 questions via `AskUserQuestion`:

| # | header | multiSelect | options |
|---|--------|-------------|---------|
| 1 | "Sections" | true | Features list · Demo/Screenshots · Configuration docs · API Reference |
| 2 | "Extras" | true | Contributing guide · Roadmap · FAQ · Acknowledgments |
| 3 | "Badges" | true | License · CI/Build · Coverage · Version/Release · Downloads · Stars+Forks |
| 4 | "Style" | false | Simple · Centered · Styled |

Style definitions (include as option descriptions):
- **Simple** — plain `# Title` with `> tagline` on next line
- **Centered** — `<div align="center">` block, title + tagline + badge row centered
- **Styled** — centered + emoji prefix on every H2 (✨ Features, 🚀 Quick Start,
  ⚙️ Configuration, 📖 API, 🤝 Contributing, 🗺 Roadmap, ❓ FAQ, 📄 License)

## Step 4 — Generate README

### Minimal Structure

Five sections in order: title + tagline, description paragraph (what/who/why),
Installation with syntax-highlighted copy-paste command, Usage with a minimal
working example plus expected output as a comment, License with SPDX name and
link to `LICENSE` file.

### Standard / Comprehensive Structure

Assemble sections in this order, including only those selected in Round 2:

1. **Header block** — apply Style from Round 2 (see styles below)
2. **Table of Contents** — Comprehensive only; auto-generate anchored links for
   every H2 section present
3. **Description** — 2–3 sentences: problem solved, target user, key differentiator
4. **Features** (if selected) — bulleted list of capabilities, present tense, parallel form
5. **Demo/Screenshots** (if selected) — `![Demo](demo.gif)` placeholder; instruct
   user to replace with an actual GIF or screenshot path
6. **Installation** — prerequisites block, then numbered install steps
7. **Usage** — single example for Standard; multiple headed examples for Comprehensive
8. **Configuration** (if selected) — markdown table: Name · Type · Default · Description
9. **API Reference** (if selected) — each function/endpoint: signature, one-line
   description, minimal example
10. **Contributing** (if selected) — fork → clone → branch → commit → PR, 5 steps
11. **Roadmap** (if selected) — `- [ ]` checklist of planned features
12. **FAQ** (if selected) — H3 questions with short paragraph answers (2–3 entries)
13. **Acknowledgments** (if selected) — linked credits
14. **License** — SPDX name, link to LICENSE file

For Comprehensive: append `[⬆ back to top](#readme)` after each H2 section.

### Header Block Styles

**Simple:** `# Project Name` followed by `> One-line tagline.`

**Centered:** Wrap title, tagline, and badge row in `<div align="center">`,
using `<h1>` for title, `<p>` for tagline, and markdown badge syntax for badges.

**Styled:** Same as Centered; prepend the relevant emoji to every H2 per the
emoji list in Step 3.

### Badge URLs

Load `references/badges.md` when assembling the badge row — it contains the
complete URL patterns for all supported services with substitution instructions.
If Language = Other, omit the tech stack badge row unless the user provides
their language name.

## Step 5 — Write and Summarize

Write to `README.md` in the current working directory. If `README.md` already
exists, confirm overwrite before writing.

After writing, output a concise summary:
- File path written
- Sections included
- Placeholder values the user must replace: `USER/REPO`, `PKG`, workflow filename,
  screenshot/GIF paths, and any description body text
