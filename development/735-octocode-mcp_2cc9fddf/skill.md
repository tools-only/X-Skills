# Octocode-MCP - Research Driven Development Platform for AI

**Research Date**: January 26, 2026
**Source URL**: <https://octocode.ai>
**GitHub Repository**: <https://github.com/bgauryy/octocode-mcp>
**npm Package**: <https://www.npmjs.com/package/octocode-mcp>
**Version at Research**: 12.0.0 (npm), GitHub release 9.1.1
**License**: MIT

---

## Overview

Octocode is an MCP server and methodology platform that transforms AI code assistants from pattern-matching systems into research-driven developers. It provides semantic code search across GitHub repositories, local filesystem analysis with LSP intelligence, and a GAN-inspired adversarial workflow for evidence-based software development.

**Core Value Proposition**: Shift AI from "guessing based on training data" to "knowing based on real-time code research" through the Research Driven Development (RDD) methodology.

---

## Problem Addressed

| Problem                                                      | How Octocode Solves It                                                                                |
| ------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------- |
| AI assistants hallucinate code patterns from training data   | Real-time GitHub search provides actual implementations as evidence                                   |
| Context windows get polluted with irrelevant information     | Minimal context principle: each session starts fresh with only task-relevant data                     |
| AI lacks semantic understanding of code relationships        | LSP integration provides go-to-definition, find-references, call hierarchy                            |
| Developers spend hours searching for implementation patterns | Natural language search across millions of GitHub repos finds real examples in seconds                |
| AI planning lacks validation                                 | GAN-inspired Planner/Verifier adversarial loop catches errors before implementation                   |
| Single-model blind spots cause repeated mistakes             | Cross-model validation ensures different models check each other's work                               |
| Code research tools require external services                | GitHub CLI authentication handles tokens automatically; local tools work offline                      |
| AI implementations lack evidence backing                     | RDD methodology requires every code decision to cite source evidence with file paths and line numbers |

---

## Key Statistics (as of January 26, 2026)

| Metric             | Value                |
| ------------------ | -------------------- |
| GitHub Stars       | 689                  |
| Forks              | 53                   |
| Contributors       | 8                    |
| Total Commits      | 358                  |
| Primary Language   | TypeScript (99.8%)   |
| Repository Created | June 5, 2025         |
| Latest npm Version | 12.0.0               |
| Latest Release     | 9.1.1 (Dec 15, 2025) |
| YouTube Tutorials  | Available            |

---

## Key Features

### GitHub Research Tools

- **githubSearchCode**: Semantic code search across millions of repositories
- **githubSearchRepositories**: Repository discovery by topic, language, stars, activity
- **githubViewRepoStructure**: Explore repository file/folder structure
- **githubGetFileContent**: Fetch file contents with token optimization
- **githubSearchPullRequests**: Find PRs by state, labels, authors
- **githubSearchCommits**: Search commit history with diffs (v6.0.0: removed for token reduction)

### Local Filesystem Tools (v8.0.0+)

- **local_ripgrep**: Fast local code search using ripgrep
- **local_fetch_content**: Read local file contents
- **Workspace Root Detection**: Automatic project boundary detection
- **Allowed Paths Configuration**: Security through path restriction

### LSP Intelligence

- **lspGotoDefinition**: Navigate to symbol definitions
- **lspFindReferences**: Find all usages of a symbol
- **lspCallHierarchy**: Understand function call relationships
- **Semantic Understanding**: Compiler-level code comprehension

### Research Skill (Agent Capability)

- **Auto-Prompting**: Injects RDD system prompts automatically
- **Advanced Planning**: Breaks problems into research questions
- **Deep Research Orchestration**: Coordinates tools in optimal order (Search → Go to Definition → Read)
- **Parallel Agents**: Spawns sub-agents for parallel research execution

### Slash Commands

| Command                | Description                                                          |
| ---------------------- | -------------------------------------------------------------------- |
| `/research`            | Deep code and product research with intelligent tool orchestration   |
| `/plan`                | Research-backed planning with step-by-step implementation guidance   |
| `/review_pull_request` | Comprehensive PR review with defects-first analysis                  |
| `/review_security`     | Security audit with vulnerability detection and remediation guidance |

### Security & Privacy

- **Automatic Secrets Detection**: Enterprise-grade data protection
- **Content Sanitization**: Removes sensitive data before processing
- **Telemetry Opt-Out**: `OCTOCODE_TELEMETRY_DISABLED=1`
- **Local-Only Mode**: Analyze code without external API calls

---

## Technical Architecture

### The RDD Adversarial Flow (GAN-Inspired)

```text
┌────────────────────────────────────────────────────────────────────────┐
│                        RDD ADVERSARIAL FLOW                            │
├────────────────────────────────────────────────────────────────────────┤
│                                                                        │
│  0. INIT RESEARCH     1. PLAN              2. VERIFY                   │
│  ┌───────────┐      ┌───────────┐       ┌───────────┐                  │
│  │ Researcher│─────►│ Planner   │─────►│ Verifier  │                   │
│  │ (Context) │      │ (Gen)     │       │ (Disc)    │                   │
│  └───────────┘      └───────────┘       └───────────┘                  │
│       │                  │                    │                         │
│       ▼                  ▼                    ▼                         │
│   [init-ctx]         [plan.md]           [plan.md']                     │
│                                                                        │
│  3. RESEARCH         4. VALIDATE                                       │
│  ┌───────────┐      ┌───────────┐                                      │
│  │ Researcher│─────►│ Verifier  │                                      │
│  │ (Gen)     │      │ (Disc)    │                                      │
│  └───────────┘      └───────────┘                                      │
│       │                  │                                             │
│       ▼                  ▼                                             │
│   [research.md]      [research.md']                                    │
│                                                                        │
│  5. IMPLEMENT        6. VALIDATE                                       │
│  ┌───────────┐      ┌───────────┐                                      │
│  │ Coder     │─────►│ Verifier  │                                      │
│  │ (Gen)     │      │ (Disc)    │                                      │
│  └───────────┘      └───────────┘                                      │
│       │                  │                                             │
│       ▼                  ▼                                             │
│   [code.ts]        [code+tests]                                        │
└────────────────────────────────────────────────────────────────────────┘
```

### Core Principles

| Principle                  | Implementation                                                            |
| -------------------------- | ------------------------------------------------------------------------- |
| Separate Sessions          | Each flow (Plan, Research, Implement) runs in isolated agent session      |
| Minimal Context per Action | Context window contains only task-relevant data                           |
| Output Bridges Actions     | `plan.md` → `research.md` → `code + tests`                                |
| Adversarial Checks         | Verifier (Discriminator) actively tries to find flaws in Generator output |
| Cross-Model Validation     | Different models check each other to eliminate shared blind spots         |

### The RDD Equation

$$Quality = \frac{Relevant\ Context}{Context\ Noise} \times Validation \times \epsilon$$

- **Static Context (Knowns)**: Local codebase truth via `octocode-local` tools
- **Dynamic Context (Unknowns)**: External patterns via GitHub search
- **Validation (Proof)**: Cross-referencing ensures map matches territory

### Monorepo Structure

| Package             | Description                               |
| ------------------- | ----------------------------------------- |
| `octocode-mcp`      | Core MCP server for GitHub, Local FS, LSP |
| `octocode-cli`      | Command-line interface for management     |
| `octocode-research` | Research Skill for autonomous RDD         |
| `octocode-vscode`   | VS Code extension for authentication      |
| `octocode-shared`   | Shared utilities and types                |

---

## Installation

### Quick Start (Recommended)

```bash
npx octocode-cli
```

Interactive wizard handles IDE detection, environment verification, and MCP configuration.

### Manual Configuration

```json
{
  "mcpServers": {
    "octocode": {
      "command": "npx",
      "args": ["octocode-mcp@latest"]
    }
  }
}
```

### Authentication

**Option 1: GitHub CLI (Recommended)**

```bash
gh auth login
```

**Option 2: Personal Access Token**

```json
{
  "mcpServers": {
    "octocode": {
      "command": "npx",
      "args": ["octocode-mcp@latest"],
      "env": {
        "GITHUB_TOKEN": "ghp_your_token_here"
      }
    }
  }
}
```

### Research Skill Installation

```bash
npx add-skill octocode-research
```

Or direct:

```bash
npx add-skill https://github.com/bgauryy/octocode-mcp/tree/main/skills/octocode-research
```

### Supported Clients

Cursor, VS Code, Claude Desktop, Claude Code, Codex, Cline, Windsurf, Warp, Goose, LM Studio, Amp, Qodo Gen, Kiro, Gemini CLI, Zed, opencode

---

## Configuration Reference

### Environment Variables

| Variable                      | Description                              |
| ----------------------------- | ---------------------------------------- |
| `GITHUB_TOKEN` / `GH_TOKEN`   | GitHub authentication                    |
| `GITLAB_TOKEN` / `GL_TOKEN`   | GitLab authentication                    |
| `GITHUB_API_URL`              | GitHub Enterprise URL                    |
| `GITLAB_HOST`                 | GitLab instance URL                      |
| `ENABLE_LOCAL`                | Enable local filesystem tools            |
| `ALLOWED_PATHS`               | Comma-separated allowed filesystem paths |
| `TOOLS_TO_RUN`                | Whitelist specific tools                 |
| `DISABLE_TOOLS`               | Blacklist specific tools                 |
| `REQUEST_TIMEOUT`             | Network timeout (ms, default: 30000)     |
| `OCTOCODE_TELEMETRY_DISABLED` | Disable telemetry                        |

### Configuration File (`~/.octocode/.octocoderc`)

```jsonc
{
  "version": 1,
  "github": {
    "apiUrl": "https://api.github.com",
    "defaultOrg": "my-company"
  },
  "local": {
    "enabled": true,
    "allowedPaths": ["/home/me/projects"],
    "excludePaths": ["node_modules", "dist", "vendor"]
  },
  "tools": {
    "disabled": ["githubSearchPullRequests"]
  },
  "lsp": {
    "enabled": true,
    "timeout": 10000
  },
  "research": {
    "defaultProvider": "github",
    "maxQueriesPerBatch": 3,
    "maxResultsPerQuery": 10
  }
}
```

---

## Relevance to Claude Code Development

### Direct Applications

1. **GitHub Research Integration**: Search real implementations instead of relying on training data patterns
2. **LSP Intelligence**: Compiler-level code understanding for semantic navigation
3. **Research Skill**: Drop-in capability for autonomous evidence-based development
4. **Cross-Model Validation**: Reduce hallucinations through adversarial checking

### Patterns Worth Adopting

1. **RDD Methodology**: Evidence-first development with plan → research → implement flow
2. **Adversarial Verification**: GAN-inspired Planner/Verifier separation
3. **Minimal Context Principle**: Fresh context windows per action, not mega-contexts
4. **Output Bridges Actions**: Artifacts (`plan.md`, `research.md`) connect isolated sessions
5. **Cross-Model Validation**: Different models validate each other to eliminate blind spots
6. **Vibe-Research UX**: Research flows that maintain developer "flow state"

### Integration Opportunities

1. **MCP Server**: Drop-in GitHub research capability for Claude Code
2. **Research Skill**: Autonomous research agent with RDD prompts
3. **Manifest Philosophy**: Could inform Claude Code skill design principles
4. **Adversarial Workflow**: Pattern for multi-agent verification in orchestration

### Comparison with Related Tools

| Feature              | octocode-mcp | narsil-mcp | GitHub MCP |
| -------------------- | ------------ | ---------- | ---------- |
| GitHub Search        | Yes          | Remote API | Native     |
| Local File Analysis  | Yes (v8+)    | Yes        | No         |
| LSP Integration      | Yes          | Yes        | No         |
| Semantic Code Search | Yes          | Yes        | No         |
| Call Graph Analysis  | No           | Yes        | No         |
| Security Scanning    | PR Review    | Yes (111)  | No         |
| RDD Methodology      | Yes          | No         | No         |
| Research Skill       | Yes          | No         | No         |
| GitLab Support       | Yes          | No         | No         |
| Cross-Model Valid.   | Yes          | No         | No         |

---

## Version History (Notable Releases)

| Version | Date     | Changes                                                                        |
| ------- | -------- | ------------------------------------------------------------------------------ |
| 12.0.0  | Jan 2026 | Current npm release                                                            |
| 9.1.1   | Dec 2025 | Latest GitHub release, added releases workflow                                 |
| 9.1.0   | Dec 2025 | Build system migration (Rollup → tsup/esbuild), standalone binary distribution |
| 8.0.0   | Nov 2025 | Local filesystem research via `octocode-mcp-local`, `local_ripgrep`            |
| 7.0.13  | Oct 2025 | Security Review Prompt (`/review_security`)                                    |
| 7.0.0   | Oct 2025 | Major architecture refactoring, centralized error handling                     |
| 6.0.0   | Sep 2025 | Token reduction: removed npm/package search and commit search tools            |

---

## References

1. **Official Website**: <https://octocode.ai> (accessed 2026-01-26)
2. **GitHub Repository**: <https://github.com/bgauryy/octocode-mcp> (accessed 2026-01-26)
3. **npm Package**: <https://www.npmjs.com/package/octocode-mcp> (accessed 2026-01-26)
4. **RDD Manifest**: <https://github.com/bgauryy/octocode-mcp/blob/main/MANIFEST.md> (accessed 2026-01-26)
5. **MCP Package README**: <https://github.com/bgauryy/octocode-mcp/blob/main/packages/octocode-mcp/README.md> (accessed 2026-01-26)
6. **Configuration Reference**: <https://github.com/bgauryy/octocode-mcp/blob/main/docs/CONFIGURATION_REFERENCE.md> (accessed 2026-01-26)
7. **Glama MCP Listing**: <https://glama.ai/mcp/servers/@bgauryy/octocode-mcp> (accessed 2026-01-26)
8. **MCP.so Listing**: <https://mcp.so/server/octocode/bgauryy> (accessed 2026-01-26)
9. **YouTube Channel**: <https://www.youtube.com/@Octocode-ai>
10. **GitHub Releases**: <https://github.com/bgauryy/octocode-mcp/releases> (accessed 2026-01-26)

---

## Freshness Tracking

| Field                        | Value                 |
| ---------------------------- | --------------------- |
| Last Verified                | 2026-01-26            |
| npm Version at Verification  | 12.0.0                |
| GitHub Release at Verif.     | 9.1.1                 |
| GitHub Stars at Verification | 689                   |
| Next Review Recommended      | 2026-04-26 (3 months) |

**Change Detection Indicators**:

- Monitor npm for version changes (active development: 12.x currently)
- Check GitHub releases for new features and breaking changes
- Review MANIFEST.md for methodology updates
- Track star growth as adoption indicator (689 stars in ~7 months)
- Verify tool additions/removals (token optimization strategy may affect tools)
- Monitor Research Skill updates for prompt improvements
