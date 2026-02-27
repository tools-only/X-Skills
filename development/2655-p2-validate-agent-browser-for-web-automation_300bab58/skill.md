---
name: Validate agent-browser for Web Automation
description: "Test agent-browser (Playwright-based) on host with unrestricted network and Playwright browsers installed.\n**Validation steps**:\n- Install browsers: `npx playwright install`\n- Test: `npx agent-browser open https://code.claude.com/docs/en/skills`\n- Test: `npx agent-browser snapshot -i` (get element refs)\n- Test: `npx agent-browser get text body` (extract page text)\n- Verify snapshot/interact/re-snapshot workflow works\n- Document prerequisites for skill to function\n**Blocked on 2026-02-05**: Could not download Playwright browsers (DNS resolution failed, missing system libs)\n**Skill location**: `.claude/skills/agent-browser/SKILL.md`"
metadata:
  topic: validate-agent-browser-for-web-automation
  source: Session experimentation 2026-02-05
  added: '2026-02-05'
  priority: P2
  type: Feature
  status: open
  groomed: '2026-02-26'
  issue: '#128'
---

## Fact-Check

Claims checked: 5
VERIFIED: 4 | REFUTED: 1 | INCONCLUSIVE: 0

1. agent-browser is a real npm package — **VERIFIED** (v0.15.0, Vercel Labs, playwright-core dependency)
2. `npx agent-browser open <url>` works — **VERIFIED** (confirmed via `--help`)
3. `npx agent-browser snapshot -i` gets element refs — **VERIFIED** (documented with `-i, --interactive` flag)
4. `npx agent-browser get text body` extracts page text — **REFUTED** — documented syntax is `get text <selector>` where selector is `@e1` ref or CSS selector. `body` works as CSS selector but is not the documented pattern.
5. snapshot/interact/re-snapshot workflow — **VERIFIED** (core workflow in SKILL.md and `--help` examples)

Refuted claims: `get text body` syntax — misleading but functional (body is a CSS selector)
Citations: npm registry (v0.15.0), `npx agent-browser --help` output, GitHub repository docs (accessed 2026-02-26)

## RT-ICA

**Goal**: Verify agent-browser skill works end-to-end on a properly equipped host and document prerequisites.

| # | Condition | Status | Evidence |
|---|-----------|--------|----------|
| 1 | agent-browser CLI available | AVAILABLE | v0.15.0 via npx |
| 2 | Playwright browsers installed | AVAILABLE | Chromium 141.0 at /root/.cache/ms-playwright/ |
| 3 | Network access for test URLs | DERIVABLE | Environment has network; target URL accessibility unverified |
| 4 | Skill documentation exists | AVAILABLE | SKILL.md + 7 references + 3 templates |
| 5 | Validation step syntax correct | MISSING | `get text body` refuted — should use `get text <selector>` with @ref or CSS |
| 6 | System libs for rendering | DERIVABLE | Headless mode should work; headed may need X11 |

**Decision**: APPROVED — no blockers for headless validation
**Missing**: Update validation step to use `get text @e1` or `get text body` (as CSS selector) with clarifying note

## Groomed (2026-02-26)

### Reproducibility

1. Install agent-browser: `npm install -g agent-browser` (or use `npx`)
2. Install Playwright browsers: `agent-browser install` (or `npx playwright install chromium`)
3. Open a test URL: `agent-browser open https://example.com`
4. Take interactive snapshot: `agent-browser snapshot -i`
5. Extract text via element ref: `agent-browser get text @e1` (ref from snapshot output)
6. Chain open + wait + snapshot: `agent-browser open https://example.com && agent-browser wait --load networkidle && agent-browser snapshot -i`
7. Close browser: `agent-browser close`

### Dependencies

- Node.js 18+ with npx
- Playwright Chromium browser binary (installed via `agent-browser install`)
- Network access to target URLs
- No external backlog items blocking this work
