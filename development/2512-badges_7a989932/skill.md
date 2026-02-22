# Badge URL Patterns

Replace `USER/REPO` with the GitHub username/repo slug, `PKG` with the package name.

## Status Badges (flat-square)

Use `flat-square` style for all status badges — consistent appearance across
the header block.

```markdown
![License](https://img.shields.io/badge/license-MIT-blue?style=flat-square)
![Build](https://img.shields.io/github/actions/workflow/status/USER/REPO/ci.yml?style=flat-square)
![Coverage](https://img.shields.io/codecov/c/github/USER/REPO?style=flat-square)
![npm](https://img.shields.io/npm/v/PKG?style=flat-square)
![PyPI](https://img.shields.io/pypi/v/PKG?style=flat-square)
![Downloads](https://img.shields.io/npm/dm/PKG?style=flat-square)
![Stars](https://img.shields.io/github/stars/USER/REPO?style=social)
```

## Tech Stack Row (for-the-badge)

Place after the header block, not inside `<div align="center">`. Use
`for-the-badge` exclusively for this row — mixing styles looks inconsistent.

```markdown
![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Node.js](https://img.shields.io/badge/Node.js-339933?style=for-the-badge&logo=node.js&logoColor=white)
![Go](https://img.shields.io/badge/Go-00ADD8?style=for-the-badge&logo=go&logoColor=white)
![Rust](https://img.shields.io/badge/Rust-000000?style=for-the-badge&logo=rust&logoColor=white)
![TypeScript](https://img.shields.io/badge/TypeScript-3178C6?style=for-the-badge&logo=typescript&logoColor=white)
```

## Substitution Rules

- **License badge:** replace `MIT` with the license display name. Shields.io badge
  path segments use `-` as a separator, so escape literal hyphens with `--`:
  `Apache-2.0` → `Apache--2.0`, `GPL-3.0` → `GPL--3.0`
- **CI badge:** replace `ci.yml` with the actual workflow filename from `.github/workflows/`
- **Version badge:** use `npm/v` for Node, `pypi/v` for Python,
  `github/v/release/USER/REPO?style=flat-square` for others
- **Stars badge:** use `?style=social` — this is the conventional choice; do not
  override with flat-square
