# Research Findings: dasel Plugin

Date: 2026-02-19

## 1. Dasel Tool (Researcher 1)

- **What**: Single-binary CLI for querying/modifying/converting structured data (JSON, YAML, TOML, XML, CSV, HCL, INI)
- **Current**: v3.2.2 (2026-02-13), major rewrite from v2 with breaking syntax changes
- **Key v3 changes**: `put`/`delete` subcommands removed; modifications via inline assignment with `--root`
- **Binary distribution**: GitHub Releases, no separate checksums file — SHA256 in API `digest` field
- **Asset naming**: `dasel_{os}_{arch}` (e.g., `dasel_linux_amd64`, `dasel_windows_amd64.exe`)
- **Version check**: `dasel --version`
- **Core syntax**: `dasel -f file.json 'selector'`, supports filter/map/each/search/conditionals/variables
- Full CLI reference: `./research-dasel-cli.md`

## 2. Existing Skill Patterns (Researcher 2)

- 23 plugins, 101 skills in repository
- Two archetypes for dasel: Tool/CLI Reference (like uv, gitlab-skill) and Data Processing (like summarizer)
- **Plugin skills MUST NOT have `name:` field** (Claude Code bug prevents slash command registration)
- Token-based complexity thresholds from plugin_validator.py
- Progressive disclosure: SKILL.md (~2000 tokens) + references/ for detail
- Agent frontmatter requires `name` + `description`; skills do not need `name`
- Study models: `plugins/python3-development/skills/uv/SKILL.md` (Tool/CLI pattern)

## 3. Schema Requirements (Researcher 3)

- Agent frontmatter: `name` (required), `description` (required), optional `tools`, `model`, `color`, `permissionMode`, `maxTurns`, `skills`
- Model aliases: `haiku`, `sonnet`, `opus`, `inherit`
- **All comma-separated strings** — tools, disallowedTools, skills (NOT YAML arrays)
- `agents` in plugin.json: MUST be array of individual file paths, NOT directory string
- All plugin paths start with `./`
- Full schema: `./research-schema.md`

## 4. Install Script Patterns (Researcher 4)

- **Approach**: Python PEP 723 with httpx + typer (repo convention: 27+ PEP 723 scripts)
- Platform detection: `platform.system()` + `platform.machine()` + `/proc/version` for WSL2
- Install dirs: `~/.local/bin` (Linux/WSL2), `%LOCALAPPDATA%\Programs\dasel` (Windows)
- SHA256 from GitHub API `assets[].digest` field (format: `sha256:<hex>`)
- Version comparison: tuple-based semver compare
- Full patterns: `./research-install-patterns.md`

## Synthesis

### Plugin Architecture

- **3 agents**: data-explorer (haiku), dasel-guide (haiku), data-analyst (sonnet)
- **3-4 skills**: dasel-reference (v3 syntax), data-exploration (query patterns), data-transformation (format conversion, mutations), tool-setup (install/update)
- **1 install script**: Python PEP 723, cross-platform (Linux x86_64/ARM64, Windows native/WSL2)
- **1 optional hook**: SessionStart to check dasel availability

### Key Design Decisions

1. v3-only (no v2 backward compatibility — v3 is current stable)
2. Sonnet for analyst (user confirmed), Haiku for explorer/guide agents
3. Skills use progressive disclosure with references/ subdirectories
4. Install script follows repo convention: PEP 723 with typer + httpx
5. No macOS support initially (user chose Linux + Windows)
