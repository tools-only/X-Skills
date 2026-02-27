# Everything Claude Code

**Research Date**: 2026-02-26
**Source URL**: <https://github.com/affaan-m/everything-claude-code>
**GitHub Repository**: <https://github.com/affaan-m/everything-claude-code>
**Version at Research**: v1.6.0 (Codex Edition)
**License**: MIT

---

## Overview

Everything Claude Code (ECC) is a complete, production-hardened Claude Code plugin delivering 13 specialized sub-agents, 48+ skills, 31+ slash commands, a trigger-based hooks pipeline, multi-language rules, and MCP server configurations assembled over 10+ months of daily use. Originating as an Anthropic hackathon winner (Cerebral Valley x Anthropic, Feb 2026), it has reached 52K+ GitHub stars and 6.4K+ forks, making it the largest community-maintained Claude Code configuration library by engagement. As of v1.6.0 it supports four AI coding tools in a single installation: Claude Code, Cursor, OpenCode, and OpenAI Codex CLI.

---

## Problem Addressed

| Problem | Solution |
|---------|----------|
| No standard starting point for Claude Code configuration | Ready-to-install plugin with 13 agents, 48 skills, 31 commands, and a hooks pipeline immediately available |
| Claude Code lacks memory between sessions | `hooks/memory-persistence/` hooks save and reload context at session start/end automatically |
| Configuring MCP servers for GitHub, Supabase, Vercel, Railway requires manual JSON editing | `mcp-configs/mcp-servers.json` ships pre-configured server definitions for common services |
| Code quality is advisory rather than enforced when using AI | Rules directory provides always-follow guidelines; PostToolUse hooks enforce linting and testing |
| Token budget wasted on language-agnostic rules in single-language projects | Rules restructured into `common/` + `typescript/` + `python/` + `golang/` — install only what applies |
| AI coding agents lack domain-specific knowledge (ClickHouse, Django, Spring Boot, Swift) | 48 skills covering analytics, backend, frontend, mobile, databases, deployment, and security |
| Claude Code configuration is not portable across teams | Plugin marketplace integration: one command installs the full config for any team member |
| AI sessions generate insights that are immediately lost | Continuous learning v2 (`/learn`, `/evolve`) extracts instincts from sessions into reusable skill files |

---

## Key Statistics

| Metric | Value | Date Gathered |
|--------|-------|---------------|
| GitHub Stars | 52,311 | 2026-02-26 |
| GitHub Forks | 6,473 | 2026-02-26 |
| Contributors | 50 (page 50 in paginated list) | 2026-02-26 |
| Watchers | 298 | 2026-02-26 |
| Open Issues | 22 | 2026-02-26 |
| Latest Release | v1.6.0 | 2026-02-24 |
| Repository Created | 2026-01-18 | 2026-02-26 |
| Last Push | 2026-02-25 | 2026-02-26 |
| Agents | 13 | 2026-02-24 (v1.6.0 release notes) |
| Skills | 48+ | 2026-02-24 (v1.6.0 release notes) |
| Commands | 31+ | 2026-02-24 (v1.6.0 release notes) |
| Internal Tests | 978 | 2026-02-24 (v1.6.0 release notes) |
| AgentShield Security Rules | 102 | 2026-02-24 (v1.6.0 release notes) |
| Community PRs Merged (v1.6.0 cycle) | 30 | 2026-02-24 (v1.6.0 release notes) |
| Primary Language | JavaScript | 2026-02-26 |
| Supported Languages | JavaScript, TypeScript, Shell, Python | 2026-02-26 |

---

## Key Features

### Agents (13 Specialized Sub-Agents)

Each agent is a Markdown file with YAML frontmatter declaring name, description, tools, and model. Claude Code delegates to these agents via the Task tool.

- **planner.md** — Feature implementation planning; decomposes tasks into ordered steps
- **architect.md** — System design decisions; evaluates trade-offs across architecture options
- **tdd-guide.md** — Test-driven development enforcement; RED → GREEN → REFACTOR protocol
- **code-reviewer.md** — Quality and security review with specific checklist categories
- **security-reviewer.md** — Vulnerability analysis; OWASP categories, secret scanning, injection risks
- **build-error-resolver.md** — Root cause analysis and fix generation for build failures
- **e2e-runner.md** — Playwright E2E test generation and execution coordination
- **refactor-cleaner.md** — Dead code removal and structural cleanup
- **doc-updater.md** — Documentation synchronization with code changes
- **go-reviewer.md / go-build-resolver.md** — Go-specific code review and build error resolution
- **python-reviewer.md** — Python-specific code review (added v1.4.0+)
- **database-reviewer.md** — Database/Supabase pattern review (added v1.4.0+)

### Skills (48+ Workflow Definitions)

Skills provide domain knowledge and workflow definitions. Coverage spans:

- **Language standards**: coding-standards, cpp-coding-standards, java-coding-standards, golang-patterns, python-patterns
- **Framework patterns**: frontend-patterns (React, Next.js), backend-patterns, django-patterns, springboot-patterns
- **Testing**: tdd-workflow, cpp-testing, django-tdd, golang-testing, python-testing, springboot-tdd, e2e-testing
- **Security**: security-review, django-security, springboot-security, security-scan (AgentShield integration)
- **Databases**: clickhouse-io, postgres-patterns, database-migrations, jpa-patterns
- **Infrastructure**: deployment-patterns, docker-patterns, api-design
- **LLM/AI cost**: cost-aware-llm-pipeline, regex-vs-llm-structured-text
- **Mobile/Apple**: swift-actor-persistence, swift-protocol-di-testing, liquid-glass-design, foundation-models-on-device, swift-concurrency-6-2
- **Session learning**: continuous-learning, continuous-learning-v2, iterative-retrieval, strategic-compact, eval-harness, verification-loop
- **Meta/tooling**: skill-stocktake, search-first, configure-ecc, project-guidelines-example

### Commands (31+ Slash Commands)

<eg>
/tdd            — Test-driven development workflow (RED → GREEN → REFACTOR)
/plan           — Feature implementation planning with clarifying questions
/e2e            — E2E test generation with Playwright
/code-review    — Quality and security code review
/build-fix      — Fix build errors via root cause analysis
/refactor-clean — Dead code and structural cleanup
/learn          — Extract patterns mid-session into reusable skills
/learn-eval     — Extract, evaluate, and save patterns with quality scoring
/checkpoint     — Save verification state for long-running tasks
/verify         — Run verification loop against saved checkpoint
/skill-create   — Generate SKILL.md files from git history analysis
/instinct-status  — View learned instincts with confidence scores
/instinct-import  — Import instincts from other sessions or contributors
/instinct-export  — Export instincts for team sharing
/evolve           — Cluster related instincts into higher-level skills
/pm2            — PM2 service lifecycle management
/multi-plan     — Multi-agent task decomposition
/multi-execute  — Orchestrated multi-agent workflow execution
/orchestrate    — Multi-agent coordination for complex workflows
/sessions       — Session history management
/security-scan  — AgentShield security auditor (102 rules, 1282 tests)
/codex-setup    — Generate codex.md for OpenAI Codex CLI compatibility
/setup-pm       — Configure package manager (npm/pnpm/yarn/bun detection)
/go-review / /go-test / /go-build  — Go language workflow commands
/python-review  — Python code review
/eval           — Evaluate implementation against criteria
/test-coverage  — Test coverage analysis
/update-docs    — Documentation sync with code
/update-codemaps — Auto-generate codebase maps
</eg>

### Hooks Pipeline

<eg>
hooks/hooks.json — JSON config with PreToolUse, PostToolUse, Stop event handlers
hooks/memory-persistence/
  session-start.js    — Load context from persistent storage at session start
  session-end.js      — Save state at session end
  pre-compact.js      — Snapshot active state before auto-compaction
  suggest-compact.js  — Strategic compaction timing suggestions
  evaluate-session.js — Auto-extract patterns from completed sessions
</eg>

### Multi-Language Rules Architecture

Rules are split into four directories, allowing selective installation:

<eg>
rules/
  common/           — Language-agnostic (coding-style, git-workflow, testing, security, patterns, agents, hooks)
  typescript/       — TypeScript/JavaScript specific
  python/           — Python specific
  golang/           — Go specific
</eg>

Install via `./install.sh typescript python golang` (the installer script handles merge/overwrite detection).

### Continuous Learning v2 (Instinct System)

The instinct system automatically extracts patterns from sessions:

- Instincts stored with confidence scores, evidence, and example sessions
- `/instinct-status` surfaces high-confidence patterns for review
- `/evolve` clusters related instincts into generalized skills via LLM analysis
- Import/export enables team-wide pattern sharing
- `continuous-learning-v2/` skill documents the full lifecycle

### AgentShield Security Auditor

Built at the Claude Code Hackathon (Cerebral Valley x Anthropic, Feb 2026). Scans Claude Code configuration for vulnerabilities.

```bash
npx ecc-agentshield scan         # Quick scan
npx ecc-agentshield scan --fix   # Auto-fix safe issues
npx ecc-agentshield scan --opus  # Deep analysis with three Opus 4.6 agents
```

Scans: CLAUDE.md, settings.json, MCP configs, hooks, agents, skills across 5 categories — secrets detection (14 patterns), permission auditing, hook injection analysis, MCP server risk profiling, and agent config review. Output formats: terminal (A-F grade), JSON, Markdown, HTML. Exit code 2 on critical findings for CI build gates.

### Cross-Platform Support

All hooks and scripts rewritten in Node.js. Supports Windows, macOS, Linux. Package manager detection priority: `CLAUDE_PACKAGE_MANAGER` env → `.claude/package-manager.json` → `package.json` `packageManager` field → lock file detection → global config → fallback to first available.

### Multi-Tool Support (v1.6.0)

Single config source for four AI coding tools:

- **Claude Code** — native via plugin marketplace
- **Cursor** — `/.cursor/` directory with rules and agents
- **OpenCode** — `/.opencode/` with plugin system, 20+ hook event types, 3 native custom tools
- **OpenAI Codex CLI** — `/.codex/` via `/codex-setup` command generating `codex.md`

### GitHub Marketplace App

Available at `github.com/marketplace/ecc-tools`. Free/Pro/Enterprise tiers. Provides `/skill-creator analyze` via issue comments and auto-triggers on push to default branch. Generates SKILL.md files, instinct collections, and pattern extraction from commit history.

---

## Technical Architecture

<eg>
everything-claude-code/
├── .claude-plugin/       — Plugin and marketplace manifests
│   ├── plugin.json         — Plugin metadata and component paths
│   └── marketplace.json    — Marketplace catalog
├── .claude/              — Claude Code native config (CLAUDE.md, settings)
├── .cursor/              — Cursor IDE config
├── .codex/               — OpenAI Codex CLI config
├── .opencode/            — OpenCode plugin config
├── agents/               — 13 agent Markdown files (YAML frontmatter)
├── skills/               — 48+ skill Markdown files organized by domain
├── commands/             — 31+ slash command Markdown files
├── rules/
│   ├── common/           — Language-agnostic guidelines (7 files)
│   ├── typescript/       — TS/JS rules
│   ├── python/           — Python rules
│   └── golang/           — Go rules
├── hooks/
│   ├── hooks.json        — All hooks config (auto-loaded by Claude Code v2.1+)
│   ├── memory-persistence/
│   └── strategic-compact/
├── scripts/
│   ├── lib/utils.js      — Cross-platform file/path/system utilities
│   ├── lib/package-manager.js — Package manager detection
│   └── hooks/            — Node.js hook implementations
├── mcp-configs/
│   └── mcp-servers.json  — GitHub, Supabase, Vercel, Railway, etc.
├── contexts/             — Dynamic system prompt injection (dev/review/research)
├── examples/             — Real-world CLAUDE.md configs (SaaS, Go microservice, Django, Rust)
├── tests/                — 978-test suite (lib, hooks, cross-platform)
├── install.sh            — Language-selective installer with merge detection
├── package.json          — npm package metadata
├── the-shortform-guide.md — Setup and foundations guide
├── the-longform-guide.md  — Token optimization, memory, evals, parallelization
├── the-security-guide.md  — Security configuration guide
└── the-openclaw-guide.md  — OpenCode integration guide
</eg>

**Installation mechanism**: Claude Code v2.1+ auto-loads `hooks/hooks.json` from any installed plugin by convention. Explicitly declaring it in `plugin.json` causes a duplicate detection error — the repo enforces this via a regression test.

**Plugin installation**:

```bash
/plugin marketplace add affaan-m/everything-claude-code
/plugin install everything-claude-code@everything-claude-code
```

Or manually add to `~/.claude/settings.json` under `extraKnownMarketplaces` and `enabledPlugins`.

**Rules require manual copy** (upstream Claude Code plugin limitation — rules cannot be distributed via plugins):

```bash
git clone https://github.com/affaan-m/everything-claude-code.git
./install.sh typescript    # or python or golang
```

---

## Installation and Usage

```bash
# Option 1: Install as Claude Code plugin (recommended)
/plugin marketplace add affaan-m/everything-claude-code
/plugin install everything-claude-code@everything-claude-code

# Option 2: Manual rules installation after plugin install
git clone https://github.com/affaan-m/everything-claude-code.git
cd everything-claude-code
./install.sh typescript          # single language
./install.sh typescript python golang  # multiple languages
./install.sh --target cursor typescript  # for Cursor

# Option 3: Direct settings.json (headless/CI environments)
# Add to ~/.claude/settings.json:
# "extraKnownMarketplaces": { "everything-claude-code": { "source": { "source": "github", "repo": "affaan-m/everything-claude-code" } } }
# "enabledPlugins": { "everything-claude-code@everything-claude-code": true }

# Run tests
node tests/run-all.js

# AgentShield security scan
npx ecc-agentshield scan
npx ecc-agentshield scan --fix
npx ecc-agentshield scan --opus --stream

# Configure OpenAI Codex CLI compatibility
/codex-setup

# Continuous learning workflow
/learn                   # Extract patterns from current session
/instinct-status         # View learned instincts with confidence scores
/evolve                  # Cluster instincts into skills

# Multi-agent orchestration
/multi-plan "Build user auth with OAuth and JWT"
/multi-execute
```

---

## Relevance to Claude Code Development

### Direct Applicability

This repository is the closest peer project to `claude_skills` in the Claude Code ecosystem. Both organize configuration as plugins with agents, skills, commands, hooks, and rules. Reviewing ECC's structure surfaces patterns directly applicable to `claude_skills` development:

- **Continuous learning v2 (instinct system)**: Automatic session pattern extraction with confidence scoring and cluster-to-skill promotion — a novel architecture not present in `claude_skills` that could address the gap between ephemeral session knowledge and persistent skills.
- **Language-selective rules**: The `common/` + per-language directory split reduces context waste for single-language projects. `claude_skills` uses a flat rules structure.
- **Skill creator command (`/skill-create`)**: Generates SKILL.md files from git history via local analysis — a workflow that could complement `claude_skills`' `/plugin-creator:skill-creator`.
- **`contexts/` directory**: Dynamic system prompt injection modes (dev/review/research) provide lightweight context switching without separate skills or commands.
- **Real-world CLAUDE.md examples**: `examples/saas-nextjs-CLAUDE.md`, `examples/go-microservice-CLAUDE.md` etc. are concrete references for users setting up their own projects.

### Patterns Worth Adopting

- **Blocking Stop hooks**: Prevent premature session end when a multi-phase workflow is in progress — applicable to `claude_skills` long-running agent workflows.
- **Pre/post-compact hooks**: `pre-compact.js` + session-start restore enables seamless continuation across context resets for any skill managing multi-phase state.
- **Regression test for plugin.json hooks field**: Auto-detecting a common misconfiguration (duplicate hooks declaration) via a test prevents repeated fix/revert cycles.
- **AgentShield scan in CI**: Running `npx ecc-agentshield scan` as a build gate (exit code 2 on critical findings) is a concrete security enforcement pattern for Claude Code plugin repositories.

### Competitive Analysis

| Dimension | ECC (affaan-m) | claude_skills (this repo) |
|-----------|----------------|--------------------------|
| Primary focus | Complete Claude Code config collection | Modular skills marketplace with plugin framework |
| Distribution | Plugin marketplace + git clone | Plugin marketplace (Claude Code native) |
| Rules distribution | Manual copy via install.sh (upstream limitation) | Same upstream limitation |
| Learning system | Continuous learning v2 (instinct-based, auto-extract) | No built-in equivalent |
| Security tooling | AgentShield (102 rules, npm package, CI-ready) | No built-in equivalent |
| Cross-tool support | Claude Code, Cursor, OpenCode, Codex CLI | Claude Code only |
| Testing | 978 tests (Node.js) | Python-based validation scripts |
| Languages supported | 6 (JS, TS, Python, Go, C++, Java, Swift, Rust, Django) | Language-agnostic |
| Stars at research | 52,311 | N/A (internal) |
| License | MIT | Internal/proprietary |

### Integration Opportunities

- The `mcp-configs/mcp-servers.json` pre-configured MCP definitions (GitHub, Supabase, Vercel, Railway) are directly reusable as references when documenting MCP server configurations in `claude_skills`.
- ECC's `examples/` CLAUDE.md files for real-world stacks (SaaS Next.js, Go microservice, Django API, Rust API) are valuable research references for `claude_skills` documentation.
- The `/security-scan` skill wrapping AgentShield could be evaluated as a security audit step in `claude_skills` CI/pre-commit pipeline.

---

## References

- [affaan-m/everything-claude-code GitHub repository](https://github.com/affaan-m/everything-claude-code) (accessed 2026-02-26)
- [v1.6.0 release — ECC Codex Edition + Github App](https://github.com/affaan-m/everything-claude-code/releases/tag/v1.6.0) (accessed 2026-02-26)
- [ECC GitHub Marketplace App](https://github.com/marketplace/ecc-tools) (accessed 2026-02-26)
- [AgentShield GitHub repository](https://github.com/affaan-m/agentshield) (accessed 2026-02-26)
- [ecc-agentshield npm package](https://www.npmjs.com/package/ecc-agentshield) (accessed 2026-02-26)
- [Repository README](https://github.com/affaan-m/everything-claude-code/blob/main/README.md) (accessed 2026-02-26)
- [Repository CLAUDE.md](https://github.com/affaan-m/everything-claude-code/blob/main/CLAUDE.md) (accessed 2026-02-26)

---

## Freshness Tracking

| Field | Value |
|-------|-------|
| Last Verified | 2026-02-26 |
| Version at Verification | v1.6.0 |
| Next Review Recommended | 2026-05-26 |
