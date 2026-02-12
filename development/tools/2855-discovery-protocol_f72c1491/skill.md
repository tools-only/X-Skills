# Browser Discovery Protocol

Reusable protocol for experiencing a web application firsthand as a QA tester.
Any skill can reference this document for the EXPERIENCE phase.

## Purpose

Discover bugs, UX issues, and broken features by **interacting with the application
as a real user** — not by reading code. Code analysis supplements discovery but
cannot replace it. Many issues only surface through real interaction.

## Prerequisites

- Chrome MCP tools available (`claude --chrome`)
- Application running (default: `localhost:3000`)
- For authenticated apps: credentials available

## Discovery Flow

### 1. Reconnaissance

```
tabs_context_mcp → get browser context
tabs_create_mcp → create fresh tab
navigate → application root
computer(screenshot) → capture initial state
read_page(filter: "all") → accessibility tree overview
```

### 2. Route Discovery

Identify all navigable routes:
- Read navigation elements via `read_page(filter: "interactive")`
- Check for route config in code (grep for routes/pages)
- Visit each discoverable route

### 3. Systematic Exploration

For each route/page:

```
a. navigate → load the page
b. computer(screenshot) → capture visual state
c. read_page(filter: "interactive") → catalog interactive elements
d. read_console_messages(pattern: "error|warn|fail") → JS errors
e. read_network_requests → failed API calls (4xx, 5xx)
f. Exercise interactive elements:
   - Click buttons, submit forms, toggle states
   - Test edge cases: empty inputs, invalid data, rapid clicks
   - Test navigation flows: back/forward, deep links
g. Record any issue discovered
```

### 4. User Flow Testing

Beyond individual pages, test end-to-end flows:
- Authentication: login, logout, session persistence
- CRUD operations: create, read, update, delete
- State transitions: how data changes propagate
- Error recovery: what happens after failures

### 5. Discovery Completeness

The model determines when discovery is sufficient. Signals:
- All visible routes have been visited
- Core user flows have been exercised
- Interactive elements on key pages have been tested
- Console errors and network failures have been cataloged
- At least one issue has been found (if the app has bugs, you'll find them)

## Output: discovery-report.json

Write `.claude/discovery-report.json`:

```json
{
  "app_url": "http://localhost:3000",
  "routes_visited": ["/ ", "/dashboard", "/settings", "/profile"],
  "timestamp": "2026-02-12T00:00:00Z",
  "issues": [
    {
      "id": "issue_1",
      "type": "bug|ux|performance|accessibility",
      "severity": "critical|high|medium|low",
      "page": "/dashboard",
      "description": "What's broken and how to reproduce",
      "user_impact": "How this affects real users",
      "evidence": "Screenshot observation or console error",
      "reproducible": true
    }
  ],
  "console_errors": ["Uncaught TypeError at dashboard.js:42"],
  "network_failures": ["GET /api/workspaces returned 401"]
}
```

## Dual Trust Layers

| Layer | What | Verified By |
|-------|------|-------------|
| **System-verified** | Pages were visited, screenshots taken, tools used | tool-usage-log.json (temporal ordering) |
| **Model-reported** | Issue descriptions, severity, user impact | Model judgment (trusted) |

The stop-validator checks that browser discovery tools appear in the tool-usage-log
BEFORE code editing tools when workflow mode is active. This proves the EXPERIENCE
phase happened — it does not judge the quality of discoveries.

## Integration with Skills

| Skill | How It Uses Discovery |
|-------|----------------------|
| `/forge` | Full EXPERIENCE phase with heavy analysis |
| `/improve` | OBSERVE phase (evaluation of known dimensions) |
| `/appfix` | Phase 1 health check (diagnosis of known symptoms) |
| `/repair` | Routes to appropriate debugging skill |

Discovery is the broadest: it finds unknowns. Other skills handle knowns.
