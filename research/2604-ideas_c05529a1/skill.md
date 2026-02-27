# Ideas to Dig

> Research topics for future guide improvements. Curated and validated.

## Done

### MCP Security Hardening ✅
Unified security research covering MCP vulnerabilities, prompt injection, and secret detection.

**Completed**: [Security Hardening Guide](./guide/security-hardening.md) covers:
- CVE-2025-53109/53110, 54135, 54136 with mitigations
- MCP vetting workflow with 5-minute audit checklist
- MCP Safe List (community vetted)
- Prompt injection evasion techniques (Unicode, ANSI, null bytes)
- Secret detection tool comparison (Gitleaks, TruffleHog, GitGuardian)
- Incident response procedures (secret exposed, MCP compromised)
- 3 new hooks: `unicode-injection-scanner.sh`, `repo-integrity-scanner.sh`, `mcp-config-integrity.sh`

---

## High Priority

*(No items currently)*

---

## Medium Priority

### CI/CD Workflows Gallery
Concrete GitHub Actions examples for Claude Code integration.

**Topics:**
- Automated PR review workflows
- Test generation pipelines
- Cost optimization patterns for API calls in CI
- Pre-commit hooks integration

**Perplexity Query:**
```
GitHub Actions workflows for AI coding assistants 2024-2025:
- Automated code review with Claude/GPT
- Cost optimization for API calls in CI/CD
- Real examples from open source projects
```

### MCP Server Catalog
Exhaustive list of MCP servers with real-world use cases.

**Topics:**
- Available servers by category (dev tools, databases, APIs)
- Performance benchmarks vs native tools
- Security trust levels per server
- Custom server development patterns

**Perplexity Query:**
```
MCP Model Context Protocol servers catalog 2024-2025:
- Most useful servers for developers
- Performance comparison MCP vs native tools
- How to build custom MCP servers
```

---

## Lower Priority

### CLAUDE.md Patterns Library
Stack-specific templates for common project types.

**Topics:**
- React/Next.js optimized configurations
- Python/FastAPI patterns
- Go project conventions
- Monorepo configurations

**Perplexity Query:**
```
CLAUDE.md configuration examples by framework:
- React, Next.js, Vue patterns
- Python, FastAPI, Django patterns
- Best practices from GitHub repositories
```

---

## Watching (Waiting for Demand)

### Multi-LLM Consultation Patterns
Using external LLMs (Gemini, GPT-4) as "second opinion" from Claude Code.

**Status:** No proven demand. Add if 3+ reader requests.

**Research done (Jan 2026):**
- Simple approach: Bash script calling Gemini API
- Production approach: [Plano](https://github.com/katanemo/plano) (overkill for solo devs)
- Community adoption: Near zero in Claude Code users

**If implementing:**
- `examples/scripts/gemini-second-opinion.sh`

### Type-Driven API Design for AI Agent Efficiency
Schema-first development impact on Claude Code token consumption.

**Status:** Anecdotal only (no empirical data). Reevaluate if benchmarks emerge.

**Resource evaluated (Feb 2026):**
- [ShipTypes](https://shiptypes.com/) by Boris Tane (Cloudflare)
- **Score:** 2/5 (Marginal) — Claims "types → fewer tokens" unverified
- **Full evaluation:** `docs/resource-evaluations/shiptypes-evaluation.md`

**What's missing:**
- Benchmark comparing token consumption: typed APIs (tRPC/Zod) vs untyped (REST/docs)
- A/B test showing AI agent iterations with/sans types
- Case study with reproducible metrics

**Reevaluation triggers:**
- [ ] Academic paper/blog with empirical data (token consumption metrics)
- [ ] Anthropic official recommendation on schema-first for Claude Code
- [ ] 5+ community discussions/issues requesting this topic

**If validated (score upgrade to 4/5):**
- Add subsection in `guide/methodologies.md` (after CDD, line 172)
- Use micro-integration template: `docs/resource-evaluations/shiptypes-evaluation.md` (section "Integration Plan")

**Check again:** August 2026
- 3-line mention in "See Also" section
- No full guide (maintenance burden, scope creep)

**Source:** [daily.dev article](https://app.daily.dev/posts/make-claude-code-opus-talk-to-gemini-pro-b7pyiq394)

### Vibe Coding Discourse

Evolution of the "developer as architect" narrative in AI-assisted development.

**Reference:** [Craig Adam - "Agile is Out, Architecture is Back"](https://medium.com/@craig_32726/agile-is-out-architecture-is-back-7586910ab810)

**Status:** Watching. Term "vibe coding" now mainstream (Collins Word of the Year 2025).

---

## Discarded Ideas

| Idea | Reason Discarded |
|------|------------------|
| LLM Fine-tuning guide | Out of scope - users don't control model training |
| Model architecture internals | Too theoretical, not actionable |
| Token pricing optimization | Changes too frequently, use official docs |
| A2A Protocol (Agent-to-Agent) | Claude Code is single-agent with sub-agents, not true multi-agent |
| AgentOps Enterprise Dashboard | Infrastructure doesn't exist for CLI tool |
| LLM-as-a-Judge evaluation | Overkill for CLI, adds latency without proportional value |
| Decision Trajectory logging | No access to internal traces (black box) |
| 4 Pillars formal framework | Too academic - guide already covers symptoms pragmatically |
| Canary/Blue-Green deployments | Infrastructure patterns, not relevant for CLI |
| Memory Poisoning defenses | Theoretical risk requiring prior system compromise |
| Prompt Engineering for Code Gen | Already well covered (xml_prompting, prompt_formula) |
| Context Window Optimization | Already well covered (context_management, context_triage) |
| Task Decomposition Patterns | Covered via plan_mode, interaction_loop |
| Agent Architecture Comparisons | Out of scope - not multi-agent theory |
| Real-World Case Studies | Non-verifiable metrics, marketing-prone |
| Comparison with Other Tools | Out of scope, rapid obsolescence |

---

## Contributing

Found something interesting? Add it with:
1. Topic name and why it matters
2. Specific research questions
3. Perplexity query to start
