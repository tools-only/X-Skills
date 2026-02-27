# Reddit Post Draft — Differentiation Strategy

**Target**: r/ClaudeAI (primary)
**Posting Window**: Tuesday–Thursday, 9–11am EST
**Strategy**: Lead with WHY (educational depth), security-first (threat DB), complementary positioning

---

## Title (Recommended)

**Claude Code's docs don't teach you WHY. Here's the educational guide with threat intelligence (22 CVEs), 257-question quiz, and 18 production security hooks.**

### Alternative Titles

**Option B** (experience-driven):
> I spent 6 months with Claude Code daily. Here's the guide I wish existed: threat intelligence, methodologies, and production templates.

**Option C** (problem-solution):
> Claude Code's docs show you WHAT to do. This guide shows you WHY it works and WHEN it breaks.

---

## Post Body

```markdown
**The gap**: Claude Code's official docs are solid references. But they don't explain WHY patterns work, don't address security threats systematically, and don't provide structured learning paths from junior to power user.

**What I built** (6 months daily practice, zero marketing BS):

**Security-First** (unique in the ecosystem):
- 🛡️ Threat intelligence database: 22 CVEs, 341 malicious skills tracked
- 🔒 18 production security hooks (bash + PowerShell) — drop-in, no config
- ⚠️ MCP vetting workflow: how to audit third-party servers before they access your codebase

**Educational Depth** (not just configs):
- 📖 19K-line guide: explains concepts first, not just configs
- 🧠 257-question quiz (only comprehensive Claude Code assessment available)
- 📊 4 learning paths: Junior → Senior → Power User → PM/DevOps/Designer
- 🔬 TDD/SDD/BDD/GSD methodologies documented with real workflows

**Production Templates** (111 total, CC0 licensed):
- Agents: Technical Writer, Security Auditor, Code Reviewer (with behavioral mindsets)
- Commands: `/commit`, `/refactor`, `/security-audit`, `/release`
- Hooks: Pre-commit, post-merge, security scanning
- Skills: Git workflows, testing patterns, deployment automation

**Quality Assurance**:
- 55 external resource evaluations (5-point scoring system)
- Only guide with systematic competitive analysis
- Architecture doc: how Claude Code works under the hood

**Positioning** (important):
- **Complementary** to everything-claude-code (learn WHY here → leverage battle-tested configs there)
- **Educational** vs tooling (if you want multi-agent orchestration, check ruvnet/claude-flow)
- **Security-focused** vs generic tips

**Quick start** (no clone needed):
```bash
claude "Fetch https://raw.githubusercontent.com/FlorianBruniaux/claude-code-ultimate-guide/main/tools/onboarding-prompt.md"
```

Or dive in:
- [Cheatsheet (1 page)](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/guide/cheatsheet.md) — printable quick ref
- [Quiz](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/tree/main/quiz) — test your knowledge (4 profiles)
- [Threat DB](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/machine-readable/threat-db.yaml) — 22 CVEs, 341 malicious skills
- [Full Guide](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/guide/ultimate-guide.md) — 19K lines, concepts-first

**Licensing**:
- Guide: CC BY-SA 4.0 (attribution required)
- Templates: CC0 (copy-paste, zero attribution)

**Repo**: https://github.com/FlorianBruniaux/claude-code-ultimate-guide

If this helps your workflow, a ⭐ helps others discover it.

---

**Why I built this**: I kept hitting security issues and pattern failures that the docs didn't explain. Spent 6 months documenting WHY patterns work, WHEN they break, and HOW to recover. This is the resource I wish existed when I started.
```

---

## Post-Posting Strategy

### Critical Engagement Window (T+0 to T+3h)

**T+0 to T+1h** (highest priority):
- Reply to ALL comments within 15 minutes
- Use prepared responses from `claudedocs/reddit-responses-prepared.md`
- Tone: Bold Guy (direct, bienveillant, énergique) — not defensive

**T+1h to T+3h**:
- Active monitoring, 30-minute response cadence
- Identify emerging themes in comments
- Cross-link to specific guide sections (provide value in replies)

**T+3h to T+24h**:
- Passive monitoring, 1-hour response cadence
- Track metrics: upvotes, comments, stars gained

### Anticipated Critiques & Responses

| Critique | Response Template |
|----------|------------------|
| "Just use everything-claude-code" | "Exactly! Learn concepts here → apply their battle-tested configs. Complementary, not competitive." |
| "Too long, overwhelming" | "That's why we have the 1-page cheatsheet + quiz paths. Start there, expand when you need depth." |
| "Security is overkill" | "22 CVEs disagree. If you're running untrusted MCP servers or community skills, the threat DB is essential." |
| "Why not contribute to official docs?" | "Different goals. Official docs = reference. This = educational + threat intelligence + methodologies." |
| "Stars = marketing spam" | "Fair concern. Verifiable: 257 quiz questions, 22 CVEs tracked, 55 resource evals. Check the work." |

### Success Metrics (24h window)

| Metric | Target | Tool |
|--------|--------|------|
| Reddit upvotes | >50 | Reddit analytics |
| Comments | >20 | Reddit analytics |
| Stars gained | >20 | GitHub Insights |
| Traffic spike | >500 unique | GitHub Insights |
| Discussions created | >5 | GitHub Discussions |

**Success threshold**: 50+ upvotes, 20+ comments, 20+ stars in 24h.

---

## Pre-Flight Checklist

- [ ] Version badge accurate in README (currently 3.26.0)
- [ ] Templates count verified (111)
- [ ] Quiz questions count verified (257)
- [ ] Threat DB stats verified (22 CVEs, 341 malicious skills)
- [ ] Landing site synchronized (`./scripts/check-landing-sync.sh`)
- [ ] All links tested:
  - [ ] Onboarding prompt
  - [ ] Cheatsheet
  - [ ] Quiz
  - [ ] Threat DB
  - [ ] Full guide
- [ ] No typos in post body
- [ ] Prepared responses ready (`claudedocs/reddit-responses-prepared.md`)

---

## Red Flags to Avoid

| ❌ Never Say | ✅ Say Instead |
|--------------|----------------|
| "Best guide" | "Most comprehensive I've found" |
| "everything-claude-code is outdated" | "Complementary to everything-claude-code" |
| "Trust me" | "Verifiable via [specific link]" |
| "Sorry if..." | (No unnecessary apologies) |
| Aggressive self-promo | Factual claims with sources |
| "Blazingly fast" | Specific metrics (e.g., "22 CVEs tracked") |

---

## Differentiation Highlights (for replies)

When responding to comments, emphasize these unique angles:

**Security-First**:
- Only guide with threat intelligence database (22 CVEs, 341 malicious skills)
- Only guide with production security hooks (18 total, bash + PowerShell)
- MCP vetting workflow documented (how to audit before trusting)

**Educational Depth**:
- Explains WHY patterns work, not just configs
- 257-question quiz (only comprehensive assessment)
- 4 learning paths (skill-based progression)

**Methodologies**:
- TDD/SDD/BDD/GSD workflows documented
- Not just tools — systematic approaches to different project types

**Complementary Positioning**:
- Learn concepts → Apply everything-claude-code configs
- Educational layer vs production configs
- No competition with ecosystem leaders

---

## Timeline

| Time | Action | Owner |
|------|--------|-------|
| Pre-post | Verify checklist above | Florian |
| Post day, 9-11am EST | Publish post | Florian |
| T+0 to T+1h | Reply to ALL comments (<15min) | Florian |
| T+1h to T+3h | Active monitoring (30min cadence) | Florian |
| T+3h to T+24h | Passive monitoring (1h cadence) | Florian |
| T+24h | Analyze metrics, assess success | Florian |
| T+7d | Post-mortem, document learnings | Florian |

---

**Status**: ✅ Ready for review (pending checklist verification)
