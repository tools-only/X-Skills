# Plugin Discussion: dasel

Date: 2026-02-19

## Scope Decisions

- Platforms: Linux (x86_64 + ARM64), Windows (native + WSL2). No macOS initially.
- Source patterns: Both official Anthropic skill examples AND this repo's existing skills
- Analyst model: Sonnet 4.6 (user-confirmed default)
- Explorer/guide agents: Haiku (user-specified)

## UX Preferences

- Install script: auto-download, auto-update, user-space bin dir
- Skills: adapt XML handling and data worker research patterns for dasel
- Agents: dedicated explorer (data finding/processing), guide (how-to), analyst (structural analysis)

## Technical Choices

- v3-only: no v2 backward compatibility (v3 is current stable, v2 syntax incompatible)
- Python PEP 723 for install script (repo convention)
- SHA256 verification via GitHub API digest field
