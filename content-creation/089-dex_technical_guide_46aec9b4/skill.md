# Dex Technical Guide

**Version:** 1.11.0
**Last Updated:** February 19, 2026

**Audience:** People who want to understand how Dex works under the hood - whether to customize it, contribute to it, or learn from its design patterns.

---

## Table of Contents

1. [Architecture Philosophy](#architecture-philosophy)
2. [The Agent Skills Standard](#the-agent-skills-standard)
3. [Skill Capabilities](#skill-capabilities)
4. [MCP Servers Deep Dive](#mcp-servers-deep-dive)
5. [Context Management](#context-management)
6. [Memory Architecture](#memory-architecture)
7. [State Management & Syncing](#state-management--syncing)
8. [Planning Architecture](#planning-architecture)
9. [Integration Layer](#integration-layer)
10. [Self-Learning System](#self-learning-system)
11. [Design Constraints](#design-constraints)

---

## Architecture Philosophy

### Why File-Based PKM?

Dex is built entirely on plain text files (markdown, YAML, shell scripts). No database, no external services for storage. Why?

**Portability:** Your data is just files. Move them anywhere, sync with any tool (Git, Dropbox, iCloud), edit with any text editor.

**Version Control:** Every change is git-trackable. You get full history, branching, and collaboration for free.

**AI-Native:** LLMs excel at reading and writing structured text. Files are the perfect interface - no ORM, no serialization overhead.

**Longevity:** Markdown files will be readable in 20 years. Proprietary databases won't.

**Transparency:** You can open any file and see exactly what's stored. No hidden schema, no data lock-in.

### The PARA Method

Dex uses Tiago Forte's PARA system (Projects, Areas, Resources, Archives) with a numbered prefix:

```
00-Inbox/          # Capture zone (step zero)
01-Quarter_Goals/  # 3-month outcomes
02-Week_Priorities/# Weekly focus
03-Tasks/          # Task backlog
04-Projects/       # Time-bound initiatives
05-Areas/          # Ongoing responsibilities
06-Resources/      # Reference material
07-Archives/       # Historical records
System/            # Configuration
```

**Why numbered?** Three reasons:

1. **Sort order:** Folders appear in workflow order (capture → plan → execute)
2. **Intentional hierarchy:** Numbers signal "these are the core structure, don't mess with them"
3. **Reference stability:** Can say "03-Tasks" and everyone knows what you mean

**Why PARA?** It maps to how knowledge flows:

- **Inbox (00)** → Raw capture, zero friction
- **Planning (01-03)** → Strategic → Tactical hierarchy
- **PARA (04-07)** → Active work → Long-term storage → Historical record

### File vs. Database Trade-offs

**What files are good at:**
- Human readability
- Version control
- Portability
- Zero infrastructure

**What files are bad at:**
- Complex queries (no SQL)
- Relationships (no foreign keys)
- Atomic transactions
- High-frequency updates

**Dex's solution:** Use files for *state*, use MCP servers for *operations*. Files are "storage", MCP provides "business logic" (validation, deduplication, syncing).

---

## The Agent Skills Standard

### What Is It?

Agent Skills (from [agentskills.io](https://agentskills.io)) is a universal format for AI workflows. It's basically "markdown files with YAML frontmatter + structured instructions."

Every skill lives in `.claude/skills/[skill-name]/SKILL.md`:

```yaml
---
name: daily-plan
description: Generate context-aware daily plan with calendar, tasks, and priorities
---

## Purpose
...instructions for Claude...
```

### Why YAML Frontmatter?

**Discoverability:** Claude can list available skills by reading frontmatter without parsing the whole file.

**Metadata:** Name, description, version, dependencies all in structured format.

**Tooling:** Other systems can parse YAML to build UIs, validation, etc.

### Skills vs. Raw Prompts

**Without skills (just CLAUDE.md):**
- One giant prompt file (10K+ lines)
- Hard to navigate
- Can't version individual workflows
- No conditional loading

**With skills:**
- Modular workflows
- User chooses what to invoke (`/daily-plan` vs `/week-review`)
- Each skill is independently versioned
- Only loads when needed (saves tokens)

**Example from `.claude/skills/daily-plan/SKILL.md`:**

```yaml
---
name: daily-plan
description: Generate context-aware daily plan with calendar, tasks, and priorities
---

## Step 0: Demo Mode Check
...

## Step 1: Gather Calendar Context
...

## Step 2: Load Tasks
...
```

Claude reads this file when you type `/daily-plan`, follows the steps, and doesn't load it otherwise.

### Role-Specific Skills

Dex has 25 core skills, plus 27 role-specific skills (Product, Sales, Marketing, etc.) stored in `.claude/skills/_available/[role]/[skill-name]/`.

**Why separate?** Not everyone needs `/pipeline-health` or `/board-prep`. Skills are discovered via `/dex-level-up` based on your role and installed on demand.

**Implementation:** Skills in `_available/` aren't loaded into Cursor's context until you explicitly install them (by moving to `.claude/skills/`).

---

## Skill Capabilities

Skills declare capabilities through YAML frontmatter. These go beyond basic metadata (name, description) to control how the skill runs, what memory it has, and what automation it triggers.

### Context Fork

**Problem:** Heavy skills like `/daily-plan` or `/week-review` load thousands of tokens of working data — calendar events, task lists, meeting notes, person context. All of that stays in the main conversation's context window. If you run `/daily-plan` at 9am and keep chatting until 5pm, those thousands of tokens are still there, making the conversation slower and the context window smaller.

**Solution:** Skills declare `context: fork` in their frontmatter:

```yaml
---
name: daily-plan
description: Generate context-aware daily plan
context: fork
---
```

**What happens:**

1. User invokes `/daily-plan`
2. Claude spawns the skill in an **isolated context** (a "fork")
3. The skill reads files, calls MCP tools, generates the plan
4. Only the final result is returned to the main conversation
5. All the working data (raw calendar JSON, full task lists, intermediate reasoning) is discarded

**Effect:** The main conversation stays clean. You can run `/daily-plan`, `/week-review`, and `/process-meetings` in the same session without accumulating dead context.

**Skills using `context: fork` (7 total):**

| Skill | Why It Forks |
|-------|-------------|
| `/daily-plan` | Loads calendar, tasks, priorities, meetings, person context |
| `/daily-review` | Reads session transcript, learnings, task completions |
| `/week-plan` | Synthesizes quarter goals, tasks, last week's review |
| `/week-review` | Scans full week of activity, meetings, task progress |
| `/meeting-prep` | Loads attendee person pages, recent meetings, account context |
| `/process-meetings` | Processes multiple Granola transcripts with full content |
| `/career-coach` | Reads career evidence, ladder, goals, role descriptions |

**When NOT to fork:** Quick skills like `/triage` or `/save-insight` that touch minimal context. The overhead of forking isn't worth it for lightweight operations.

### Skill-Scoped Hooks

**Problem:** Global hooks (like the person context injector) fire on every relevant tool call across the entire session. Some automation should only run during a specific skill — updating person pages during meeting processing, capturing evidence during career coaching.

**Solution:** Skills declare their own hooks in frontmatter:

```yaml
---
name: process-meetings
description: Process synced Granola meetings
context: fork
hooks:
  PostToolUse:
    - matcher: Write
      type: command
      command: "node .claude/hooks/post-meeting-person-update.cjs"
  Stop:
    - type: command
      command: "node .claude/hooks/meeting-summary-generator.cjs"
---
```

**Two hook events are available:**

- **PostToolUse** — Fires after a specific tool completes. The `matcher` field specifies which tool (e.g., `Write`). The hook receives tool context via `CLAUDE_HOOK_CONTEXT` environment variable.
- **Stop** — Fires once when the skill finishes execution. Used for summary generation or cleanup.

**Currently active skill-scoped hooks:**

| Skill | Event | Hook | Purpose |
|-------|-------|------|---------|
| `/process-meetings` | PostToolUse (Write) | `post-meeting-person-update.cjs` | Updates person pages with meeting references |
| `/process-meetings` | Stop | `meeting-summary-generator.cjs` | Placeholder for future meeting summaries |
| `/career-coach` | PostToolUse (Write) | `career-evidence-capture.cjs` | Logs achievements to Evidence Log |
| `/daily-plan` | Stop | `daily-plan-quick-ref.cjs` | Generates condensed quickref file |

**How hook context works:**

When a PostToolUse hook fires, it receives the tool's input via environment variable:

```javascript
const input = JSON.parse(process.env.CLAUDE_HOOK_CONTEXT || '{}');
const filePath = input?.tool_input?.file_path || '';
```

The hook can then inspect what file was written and take action — extracting person mentions, detecting achievement markers, etc.

**Why this matters:** Automation stays scoped to the right context. The person page updater only fires during meeting processing, not every time any file is written anywhere in the system.

### Hook Files

Dex ships with these hook scripts in `.claude/hooks/`:

**Global hooks (always active):**

| File | Trigger | What It Does |
|------|---------|-------------|
| `person-context-injector.cjs` | PostToolUse (Read) | Injects person page summaries when files mention people |
| `company-context-injector.cjs` | PostToolUse (Read) | Injects company/account context when files mention organizations |
| `session-start.sh` | Session start | Shows strategic hierarchy, urgent tasks, active patterns |
| `ensure-mcp-user-scope.cjs` | PreToolUse | Validates MCP tool calls use correct vault scope |

**Skill-scoped hooks (only active during their parent skill):**

| File | Parent Skill | Trigger | What It Does |
|------|-------------|---------|-------------|
| `post-meeting-person-update.cjs` | `/process-meetings` | PostToolUse (Write) | Detects person mentions in meeting notes (WikiLinks and plain names), appends meeting reference to their person page under "Recent Meetings" section |
| `meeting-summary-generator.cjs` | `/process-meetings` | Stop | Stub for future meeting summary generation |
| `career-evidence-capture.cjs` | `/career-coach` | PostToolUse (Write) | Scans career files for achievement markers (percentages, dollar amounts, impact verbs), infers skill area, appends to `05-Areas/Career/Evidence_Log.md` |
| `daily-plan-quick-ref.cjs` | `/daily-plan` | Stop | Reads the generated daily plan, extracts top focus items, key meetings, and time blocks into a condensed `YYYY-MM-DD-quickref.md` file |

**Standalone utility:**

| File | How to Run | What It Does |
|------|-----------|-------------|
| `maintenance.cjs` | `node .claude/hooks/maintenance.cjs` | Vault health checks: stale inbox files (>30 days), broken WikiLinks (samples 100 files), orphaned person pages (no references in tasks or meetings), stale agent memory files (>90 days) |

### Fast Mode Hints

**Problem:** Some skills make quick, low-stakes decisions — routing a file to the right folder, scanning inbox items. These don't need deep reasoning. Running them at full deliberation speed wastes time.

**Solution:** Skills declare `model_hint: fast` in frontmatter:

```yaml
---
name: triage
description: Strategically route orphaned files and extract scattered tasks
model_hint: fast
---
```

**What it does:** Signals the AI client to use faster output settings for this skill. The same model runs, just optimized for speed over deep reasoning.

**Currently used by:** `/triage` (file and task routing decisions)

**When to use it:** Skills that process many items with simple decisions (route this file, classify this task, yes/no on inbox items). Not appropriate for skills that need synthesis or nuanced judgment (daily plan, career coaching, weekly review).

### Agent Memory

**Problem:** Background agents (like the ones powering morning intelligence briefings) scan your vault for attention items — stale deals, missed commitments, unlinked people. Without memory, they re-scan everything every time and report the same findings repeatedly.

**Solution:** Agents declare `memory: project` in their frontmatter. Each agent writes a JSON memory file to `.claude/memory/`:

```json
{
  "last_scan_date": "2026-02-19",
  "flagged_items": [
    { "id": "deal-acme", "first_flagged": "2026-02-15", "times_flagged": 3 }
  ],
  "resolved_items": [
    { "id": "deal-globex", "resolved_date": "2026-02-18" }
  ],
  "escalations": [
    { "id": "deal-acme", "reason": "Flagged 3 times with no action" }
  ]
}
```

**How it works:**

1. Agent starts and checks its memory file in `.claude/memory/`
2. Compares current scan results against previous findings
3. Only surfaces **new items**, **escalations** (flagged multiple times), and **resolved items** (previously flagged, now handled)
4. Writes updated memory back to the file

**Effect:**

- **First run:** Full scan, everything is new
- **Second run:** Only new or changed items surface. Resolved items quietly drop off.
- **Repeated flags:** Items flagged across multiple sessions get escalated — "I've flagged this three sessions running. Still no action."

**Memory ownership:** Agent memory is scoped per-agent. The deal attention agent can't read the commitment tracker's memory, and vice versa. See `06-Resources/Dex_System/Memory_Ownership.md` for the full boundary diagram across all four memory layers.

---

## MCP Servers Deep Dive

### What Is MCP?

**Model Context Protocol** (from Anthropic) is an open standard for connecting AI assistants to external data sources and tools.

Think of it like this:
- **Files** = static knowledge (your notes, tasks, meetings)
- **MCP Servers** = dynamic operations (create task, check calendar, sync data)

MCP servers expose **tools** (like API endpoints) that Claude can call. Each tool has:
- Schema (required/optional parameters)
- Validation rules
- Deterministic behavior

### Why MCP Instead of Just Files?

**Problem:** Claude can read/write files, but it can't enforce rules:
- No deduplication (might create duplicate tasks)
- No validation (might break task ID format)
- No syncing (updates task file but doesn't update person pages)

**Solution:** MCP servers provide structured operations with guardrails.

**Example:** Creating a task

**Without MCP (just file edits):**
```markdown
- [ ] Build new feature
```
Problems: No task ID, no pillar tag, no bidirectional link to person page.

**With MCP (`create_task` tool):**
```python
create_task(
    title="Build new feature",
    pillar="product",
    person="John_Doe"
)
```
Result:
- Task created with proper ID (`^task-20260128-001`)
- Pillar tag validated against `System/pillars.yaml`
- Task added to `03-Tasks/Tasks.md`
- Task added to `05-Areas/People/Internal/John_Doe.md` → Related Tasks section
- All atomic, all validated

### Dex's 7 MCP Servers

#### 1. **Work MCP** (`core/mcp/task_server.py`)

**Purpose:** Task, priority, and goal management with validation and syncing.

**Key features:**
- Task ID generation (`^task-YYYYMMDD-XXX`)
- Pillar validation (checks `System/pillars.yaml`)
- Deduplication (fuzzy matching to detect similar tasks)
- Bidirectional syncing (task ↔ person page ↔ meeting notes)

**Code snippet from `task_server.py`:**

```python
def generate_task_id() -> str:
    """Generate unique task ID: ^task-YYYYMMDD-XXX"""
    today = datetime.now().strftime('%Y%m%d')
    tasks_file = get_tasks_file()
    
    if not tasks_file.exists():
        return f"^task-{today}-001"
    
    # Find existing tasks for today
    content = tasks_file.read_text()
    pattern = rf'\^task-{today}-(\d{{3}})'
    matches = re.findall(pattern, content)
    
    if not matches:
        return f"^task-{today}-001"
    
    max_num = max(int(m) for m in matches)
    next_num = max_num + 1
    return f"^task-{today}-{next_num:03d}"
```

This ensures every task gets a unique, sortable ID that's stable across files.

**Why this matters:** Task IDs are how we maintain relationships. When a meeting note says "^task-20260128-001", Dex can find that task in `03-Tasks/Tasks.md` AND on the person page AND link back to the meeting.

#### 2. **Calendar MCP** (`user-dave-calendar-mcp`)

**Purpose:** Read-only access to Apple Calendar for meeting context.

**Why Apple Calendar?** It syncs with Google Calendar accounts locally, so it's a universal interface for macOS users.

**Key tools:**
- `calendar_list_calendars` - Show available calendars
- `calendar_list_events` - Get meetings for date range

**How it's used:** The `/daily-plan` skill calls `calendar_list_events(start_date="2026-01-28")` to show today's meetings. Claude then cross-references meeting attendees with person pages to inject context.

#### 3. **Granola MCP** (`core/mcp/granola_server.py`)

**Purpose:** Fetch meeting transcripts and notes from Granola.

**What's Granola?** AI meeting assistant that records, transcribes, and summarizes meetings.

**Architecture:** API-first with cache fallback (v2.0)
- **Primary:** Uses Granola's unofficial API for complete historical data (91% success rate)
- **Fallback:** Reads from local cache (`~/Library/Application Support/Granola/cache-v3.json`) if API fails
- **Protection:** Response caching (5 min TTL), exponential backoff, graceful degradation

**Key tools:**
- `granola_get_recent_meetings` - Get meetings within date range
- `granola_get_meeting_details` - Get full details + transcript
- `granola_search_meetings` - Search by title/attendee/content
- `granola_check_available` - Verify API + cache availability

**How it works:**
1. Reads auth token from `~/Library/Application Support/Granola/supabase.json`
2. Hits `https://api.granola.ai/v2/get-documents` with Bearer auth
3. Converts ProseMirror JSON to Markdown
4. Falls back to cache on rate limits (429), auth failures (401), or network errors
5. Caches responses (5 min) to avoid rate limits

**How it's used:** The `/process-meetings` skill calls `granola_get_recent_meetings(days_back=7)` to find recent meetings, then extracts:
- Action items
- Decisions made
- People mentioned
- Career development context (if it's a 1:1 with manager)

**Why API-first?** Granola's cache doesn't retain full content for older meetings. API provides 9x more complete historical data than cache-only approach.

#### 4. **Career MCP** (`core/mcp/career_server.py`)

**Purpose:** Manage career development artifacts (job descriptions, ladders, reviews, goals).

**Key tools:**
- `get_current_role` - Read job description
- `get_career_ladder` - Read promotion criteria
- `get_growth_goals` - Read active goals
- `add_evidence` - Save achievements for reviews

**How it's used:** The `/career-coach` skill uses these tools to provide personalized coaching based on your actual role and ladder.

#### 5. **Resume MCP** (`core/mcp/resume_server.py`)

**Purpose:** Build and update resume/LinkedIn profile based on evidence.

**Key tools:**
- `get_resume_sections` - Read current resume
- `update_resume_section` - Edit specific section
- `generate_linkedin_profile` - Convert resume to LinkedIn format

**How it's used:** The `/resume-builder` skill interviews you about your experience, then structures it into ATS-friendly resume format.

#### 6. **Dex Improvements MCP** (`core/mcp/dex_improvements_server.py`)

**Purpose:** Capture and rank ideas for improving Dex itself.

**Key tools:**
- `capture_idea` - Quick save improvement idea
- `get_backlog` - Read all ideas with rankings
- `update_idea_status` - Mark as implemented/rejected

**How it's used:** 
- You can call `capture_idea("Add weekly goal rollover")` from any context
- `/dex-backlog` ranks ideas by impact/alignment/token-efficiency
- `/dex-improve [idea]` workshops an idea into implementation plan

#### 7. **Onboarding MCP** (`core/mcp/onboarding_server.py`)

**Purpose:** Stateful onboarding with validation, dependency checking, and vault creation.

**Key features:**
- Session state management with resume capability
- Step-by-step validation enforcement (cannot skip required fields)
- Email domain validation (Step 4) with format checking (no @, must have dot)
- Dependency verification (Python packages, Calendar.app, Granola)
- Automatic MCP configuration with VAULT_PATH substitution
- PARA folder structure creation

**Key tools:**
- `start_onboarding_session()` - Initialize or resume from `System/.onboarding-session.json`
- `validate_and_save_step(step_number, step_data)` - Validate and save each step (1-6)
- `get_onboarding_status()` - Check completion status and missing steps
- `verify_dependencies()` - Check Python packages and system requirements
- `finalize_onboarding()` - Create vault structure, write configs, setup MCP

**Why it matters:** Email domain (Step 4) is critical for Internal/External person routing. Without it, the system can't automatically route people to the correct folder or create company pages for external organizations. The MCP enforces this validation - you cannot skip Step 4 or finalize without a valid email domain.

**Session state example:**
```json
{
  "version": "1.0",
  "completed_steps": [1, 2, 3, 4],
  "current_step": 5,
  "data": {
    "name": "Jane Doe",
    "role": "Product Manager",
    "email_domain": "acme.com",
    "pillars": []
  }
}
```

If interrupted, calling `start_onboarding_session()` resumes from the last completed step.

### External MCP Integrations

#### Pendo (Optional for Pendo Customers)

**Type:** Hosted external MCP server (not shipped with Dex)

**Purpose:** Product analytics for Pendo customers - track guide performance, feature adoption, visitor/account engagement.

**Setup:**
1. Admin must enable in Pendo: Settings → Subscription Settings → AI Features → Pendo MCP Server
2. Add to AI client config (Cursor example):
```json
{
  "mcpServers": {
    "pendo": {
      "url": "https://app.pendo.io/mcp/v0/shttp"
    }
  }
}
```
3. Authenticate with OAuth using Pendo login credentials

**Regional URLs:**
- US: `https://app.pendo.io/mcp/v0/shttp`
- US1: `https://us1.app.pendo.io/mcp/v0/shttp`
- EU: `https://app.eu.pendo.io/mcp/v0/shttp`
- Japan: `https://app.jpn.pendo.io/mcp/v0/shttp`
- Australia: `https://app.au.pendo.io/mcp/v0/shttp`

**Available tools:**
- Visitor and account metadata
- Page, Feature, and Track Event analytics
- Event-level aggregation queries
- Activity and engagement patterns

**Use cases:**
- "What's our top performing guide this month?"
- "Which accounts are most active in the last 30 days?"
- "How many users adopted the new dashboard feature?"

**Documentation:** https://support.pendo.io/hc/en-us/articles/41102236924955

### MCP Development Pattern

All Dex MCP servers follow this pattern:

```python
from mcp.server import Server
import mcp.server.stdio

server = Server("server-name")

@server.call_tool()
async def tool_name(arguments: dict) -> list:
    """Tool description"""
    # 1. Validate inputs
    # 2. Read/write files
    # 3. Return structured response
    
# Run the server
mcp.server.stdio.stdio_server()(server)
```

**Why Python?** It's the lingua franca of data work, has great YAML/markdown libraries, and the MCP SDK is well-maintained.

**Why async?** MCP servers run as background processes. Async ensures they don't block on file I/O.

---

## Context Management

### The Token Budget Problem

Claude Code has a context window (currently ~200K tokens). Every file you read, every skill you load, every message in the chat - all count against this budget.

**Naive approach:** Load everything at session start.
**Result:** 50K tokens gone before you type anything. Chat ends after 10 exchanges.

**Dex's approach:** Lazy loading + hooks + strategic injection.

### Hooks: Background Context Injection

Hooks are scripts that run automatically when Claude uses certain tools. They *augment* the tool output without Claude explicitly asking.

**Example: Person Context Hook**

**File:** `.claude/hooks/person-context-injector.cjs`

**When it runs:** Whenever Claude calls the `Read` tool and the file contains a person's name.

**What it does:**

```javascript
// 1. Detect person references in file
const content = fs.readFileSync(filePath, 'utf-8');
const personIndex = buildPersonIndex(); // All person page filenames

// 2. Find matches
const foundPeople = new Set();
for (const name of personIndex.keys()) {
  if (content.toLowerCase().includes(name)) {
    foundPeople.add(personIndex[name]);
  }
}

// 3. Inject person page summaries
for (const personFile of foundPeople) {
  const personContent = fs.readFileSync(personFile, 'utf-8');
  // Extract role, last interaction, open tasks
  console.log(`<person_context>${summary}</person_context>`);
}
```

**Why this is powerful:**

- **Automatic:** Claude doesn't need to "know" to check person pages
- **Targeted:** Only injects context for people actually mentioned
- **Token-efficient:** Summaries, not full files

**Similar hooks:**
- `company-context-injector.cjs` - Injects company/account context
- `session-start.sh` - Shows strategic hierarchy at session start

### Session Start Hook

**File:** `.claude/hooks/session-start.sh`

**When it runs:** Every time you open a chat with Claude.

**What it shows:**

```bash
=== Dex Session Context ===

--- Strategic Pillars ---
• Product — Ship features that delight users
• Growth — 10X user base in 2026

--- Quarter Goals ---
### 1. Launch mobile app (Q1)
**Progress:** Design complete, dev 40%

--- This Week's Top 3 ---
1. Finish onboarding flow
2. Partner API integration
3. Sprint planning

--- Urgent Tasks ---
- [ ] Review PR #245 (P0)

--- Working Preferences ---
• Writing: Terse, bullet points, no preamble

--- Active Mistake Patterns (2) ---
• Over-promising timelines without checking with eng

=== End Session Context ===
```

**Why this matters:**

- **Strategic alignment:** Claude sees your pillars/goals every session
- **Immediate context:** Knows what's urgent before you ask
- **Learning integration:** Shows past mistakes to avoid repeating

**Token cost:** ~500 tokens. Worth it for consistent context.

### Lazy Loading Pattern

Skills aren't loaded until invoked:

1. User types `/daily-plan`
2. Cursor finds `.claude/skills/daily-plan/SKILL.md`
3. Claude reads the skill file (~2K tokens)
4. Claude follows the instructions
5. Skill is unloaded after command completes

**Alternative (bad):** Load all 42 skills at session start = 84K tokens before chat begins.

---

## Memory Architecture

Dex has four distinct memory layers, each owning a different scope of information. They stack — they don't compete.

### Layer 1: Claude Auto-Memory (Native)

**Owns:** Preferences, style, communication patterns, formatting choices.

**Examples:** "User prefers bullet points", "Use neutral mermaid theme", "Direct communication style."

**How it works:** Automatically captured by Claude's built-in memory system. Persists across all sessions and AI clients. No Dex configuration needed.

**Implication for Dex:** Don't duplicate. If Claude already remembers a formatting preference, Dex doesn't need to store it in a file.

### Layer 2: Agent Memory (Per-Agent State)

**Owns:** Operational state for individual agents across sessions.

**Examples:** "Deal attention agent flagged Acme Corp 3 times", "Commitment tracker: pricing follow-up resolved."

**How it works:** Each agent reads/writes its own JSON file in `.claude/memory/`. Declared via `memory: project` in frontmatter. Scoped strictly to that agent — one agent cannot read another's memory.

**File format:**

```
.claude/memory/
├── deal-attention.json
├── commitment-tracker.json
├── people-radar.json
├── project-pulse.json
├── focus-guardian.json
└── pillar-balance.json
```

Each file contains:

```json
{
  "last_scan_date": "2026-02-19",
  "flagged_items": [],
  "resolved_items": [],
  "escalations": []
}
```

**Why this matters:** Without agent memory, the deal attention agent scans your entire vault every session and reports the same stale deal every morning. With memory, it knows what it already told you, tracks resolution, and escalates items you're ignoring.

### Layer 3: Dex Session Memory

**Owns:** Operational decisions, commitments, work patterns, system learnings.

**Examples:** "Agreed to deliver a deck by Friday", "Meeting-prep skill needs more account context."

**How it works:** Captured at session end (during `/daily-review`), stored in `System/Session_Learnings/`. Consolidated during weekly reviews.

### Layer 4: Vault Search (QMD)

**Owns:** Semantic search across all vault content.

**How it works:** If QMD (local semantic search) is configured, it indexes all markdown files and provides hybrid search (keyword + vector + LLM reranking). This is a read-only layer — it doesn't store new information, it makes existing information findable by meaning.

### Boundary Diagram

```
┌───────────────────────────────────────────────┐
│  Claude Auto-Memory (native)                  │
│  Preferences, style, tone, formatting         │
│  Scope: All sessions, all AI clients          │
├───────────────────────────────────────────────┤
│  Agent Memory (.claude/memory/*.json)         │
│  Per-agent operational state                  │
│  Scope: Per agent, across sessions            │
├───────────────────────────────────────────────┤
│  Session Memory (System/Session_Learnings/)   │
│  Decisions, commitments, patterns             │
│  Scope: Per session, consolidated weekly      │
├───────────────────────────────────────────────┤
│  Vault Search (QMD — optional)                │
│  Semantic index across all vault content      │
│  Scope: Read-only, full vault                 │
└───────────────────────────────────────────────┘
```

**Full reference:** `06-Resources/Dex_System/Memory_Ownership.md`

---

## State Management & Syncing

### The Canonical File Pattern

Dex has three "canonical" files that are the single source of truth:

1. `01-Quarter_Goals/Quarter_Goals.md` - Goals
2. `02-Week_Priorities/Week_Priorities.md` - Weekly priorities  
3. `03-Tasks/Tasks.md` - All tasks

**Why canonical?** Prevents sync conflicts. If the same task appears in 5 places, which one is "real"?

**Dex's rule:** The canonical file is truth. Other locations are *references* to it via task ID.

### Bidirectional Syncing

When you create a task via Work MCP:

```python
create_task(
    title="Review API design",
    pillar="product",
    person="John_Doe",
    meeting_source="00-Inbox/Meetings/2026-01-28/api-review.md"
)
```

**What happens:**

1. **Task is created in canonical file** (`03-Tasks/Tasks.md`):
   ```markdown
   - [ ] Review API design ^task-20260128-001 #product
   ```

2. **Task is added to person page** (`05-Areas/People/Internal/John_Doe.md`):
   ```markdown
   ## Related Tasks
   - [ ] Review API design ^task-20260128-001 #product
   ```

3. **Task is added to meeting note** (`00-Inbox/Meetings/2026-01-28/api-review.md`):
   ```markdown
   ## Action Items
   - [ ] Review API design ^task-20260128-001 #product
   ```

4. **If there's a project tag**, task is added to project file:
   ```markdown
   ## Next Actions
   - [ ] Review API design ^task-20260128-001 #product
   ```

**All four locations get updated atomically.** The task ID (`^task-20260128-001`) is how we maintain links.

### Task Completion Flow

When you say "I finished reviewing the API design":

1. Claude searches for the task (fuzzy match on title)
2. Finds task ID: `^task-20260128-001`
3. Calls Work MCP: `update_task_status(task_id="task-20260128-001", status="d")`
4. MCP updates **all four locations**:
   - Changes `- [ ]` to `- [x]` 
   - Adds completion timestamp
   - Archives if configured

**Code from `task_server.py`:**

```python
def update_task_status(task_id: str, status: str):
    """Update task status everywhere it appears"""
    # 1. Update canonical file
    tasks_file = get_tasks_file()
    content = tasks_file.read_text()
    content = re.sub(
        rf'- \[ \] (.*){task_id}',
        rf'- [x] \1{task_id} ✅ {datetime.now():%Y-%m-%d %H:%M}',
        content
    )
    tasks_file.write_text(content)
    
    # 2. Find all person pages that reference this task
    people_dir = get_people_dir()
    for person_file in people_dir.rglob('*.md'):
        if task_id in person_file.read_text():
            # Update person page too
            update_task_in_file(person_file, task_id, status)
    
    # 3. Find all meeting notes that reference this task
    # ... (similar logic)
    
    # 4. Find all project files that reference this task
    # ... (similar logic)
```

This is **deterministic business logic** that files alone can't provide.

### Why Not Just File Search-and-Replace?

**You could** do this with Claude directly editing files. Problems:

1. **Race conditions:** If two tasks are updated simultaneously, edits conflict
2. **Validation:** Claude might use wrong checkbox format, break task ID
3. **Discoverability:** Hard to know which files need updating without scanning everything
4. **Rollback:** If update fails halfway through, you have partial state

**MCP servers solve this:** They're single-threaded, validated, transactional (either all updates succeed or none do).

---

## Planning Architecture

### Why a Hierarchy?

Dex's planning structure is:

```
Strategic Pillars (System/pillars.yaml)
    ↓
Quarter Goals (01-Quarter_Goals/)
    ↓
Week Priorities (02-Week_Priorities/)
    ↓
Daily Plan (07-Archives/Plans/)
    ↓
Tasks (03-Tasks/)
```

**Why this matters:**

**Without hierarchy:**
- Tasks are disconnected
- No way to prioritize (everything feels urgent)
- No learning over time (just a treadmill)

**With hierarchy:**
- Every task ladders up to a goal
- Goals ladder up to pillars
- You can ask "Does this task advance my goals?" (strategic filter)
- Reviews compound knowledge (see patterns across quarters)

### Example Flow

**User's pillar:** "Product — Ship features that delight users"

**Quarter goal (Q1):** "Launch mobile app beta with 5 core features"

**Week priority:** "Finish onboarding flow (blockers: API auth, designs)"

**Daily plan:** 
- 9am: Review onboarding designs with Maya
- 10am: Pair with Jordan on OAuth flow
- 2pm: Test signup flow end-to-end

**Tasks:**
- `- [ ] Review onboarding mockups ^task-20260128-001 #product [Q1-1] [Week-3]`
- `- [ ] Implement OAuth with Google ^task-20260128-002 #product [Q1-1] [Week-3]`

**Why this works:**

1. **Tasks are contextualized:** You know *why* you're doing them (goal linkage)
2. **Priority is clear:** Week priorities = top 3 things advancing quarterly goals
3. **Reviews are meaningful:** "Did I advance Q1-1 this week?" (measurable)
4. **System learns:** Patterns emerge (e.g., "onboarding always takes 2x estimate")

### Tags as Metadata

Tasks use three tag types:

- **Pillar tags:** `#product`, `#growth`, `#operations`
- **Goal tags:** `[Q1-1]`, `[Q1-2]` (quarter-goal linkage)
- **Week tags:** `[Week-3]` (which week to focus on)

**Why tags not folders?**

- Tasks can relate to multiple contexts (folder = one location)
- Tags are grepable (`grep "#product" 03-Tasks/Tasks.md`)
- Tags are flexible (add new ones without restructuring)

---

## Integration Layer

### Calendar Integration

**Goal:** Surface today's meetings in daily plan with context about attendees.

**Tech stack:**
- **Calendar MCP** (`user-dave-calendar-mcp`)
- **Apple Calendar.app** (syncs Google Calendar accounts locally)

**Flow:**

1. User runs `/daily-plan`
2. Skill calls `calendar_list_events(start_date="2026-01-28")`
3. MCP returns meeting list:
   ```json
   [
     {
       "title": "API Review",
       "start": "2026-01-28T10:00:00",
       "attendees": ["john@company.com", "maya@company.com"]
     }
   ]
   ```
4. Skill cross-references attendees with person pages:
   - `john@company.com` → `05-Areas/People/Internal/John_Doe.md`
   - `maya@company.com` → `05-Areas/People/Internal/Maya_Patel.md`
5. Skill reads person pages, extracts context:
   - John: Tech Lead, last 1:1 was about API architecture
   - Maya: Designer, working on onboarding flow redesign
6. Skill injects context into daily plan:
   ```markdown
   **10:00 - API Review** (John, Maya)
   - John's context: Tech Lead, discussed API patterns in last 1:1
   - Maya's context: Designer, onboarding flow work
   - Prep: Review API design doc, bring questions about auth flow
   ```

**Why this matters:** You walk into meetings prepared, without manually digging through notes.

### Granola Integration

**Goal:** Process meeting transcripts to extract action items and update person pages.

**Tech stack:**
- **Granola MCP** (`user-granola`)
- **Granola app** (records + transcribes meetings)
- **Background automation** (`.scripts/meeting-intel/sync-from-granola.cjs`)

**Flow (automated, runs daily):**

1. **LaunchAgent triggers** (5pm daily)
2. **Script calls Granola MCP:** `search_notes(since="yesterday")`
3. **For each meeting:**
   - Extract action items
   - Detect decisions made
   - Identify people mentioned
   - Check for career development context (if manager 1:1)
4. **Save to Dex:**
   - Meeting note: `00-Inbox/Meetings/YYYY-MM-DD/meeting-slug.md`
   - Tasks: Create via Work MCP (`create_task`)
   - Person pages: Add meeting reference + action items
   - Career folder: If manager 1:1, save feedback to `05-Areas/Career/Evidence/`

**User experience:** Meetings auto-sync. Just run `/process-meetings` to review and triage.

### Why MCP for Integrations?

**Alternative:** Claude could call APIs directly (e.g., `curl https://granola-api.com/notes`).

**Problems:**
- Auth management (API keys in every chat)
- Rate limiting (no throttling logic)
- Response parsing (raw JSON, no validation)
- Error handling (API down = chat breaks)

**MCP solution:**
- Auth is configured once (in MCP server settings)
- MCP servers handle retries, rate limits, caching
- Responses are validated and structured
- Errors are graceful (MCP returns error message, chat continues)

---

## Self-Learning System

### Why Systems Should Learn

Most productivity tools are static: you configure them once, use them forever. They don't adapt.

**Dex's philosophy:** The system should learn from:
1. **Your mistakes** (patterns to avoid)
2. **Your preferences** (how you like to work)
3. **Claude's updates** (new capabilities)

### Three Learning Mechanisms

#### 1. Mistake Patterns

**File:** `06-Resources/Learnings/Mistake_Patterns.md`

**Captures:** Recurring mistakes with root causes and prevention strategies.

**Structure:**

```markdown
## Active Patterns

### Over-promising timelines without checking capacity
**Pattern:** Committing to dates in meetings without consulting eng
**Root cause:** Pressure to satisfy stakeholders + optimism bias
**How to prevent:** Always respond "Let me check with the team" in meetings
**Trigger:** When someone asks "Can you ship this by Friday?"
```

**How it's used:**

- **Session start hook** shows active patterns (reminder at start of day)
- **Daily review** prompts "Any mistakes today worth capturing?"
- **Weekly review** scans for patterns (e.g., "you over-promised 3 times this week")

**Implementation:** `.claude/hooks/session-start.sh` greps the file for `## Active Patterns` and shows top 3.

#### 2. Working Preferences

**File:** `06-Resources/Learnings/Working_Preferences.md`

**Captures:** How you prefer to work (communication style, collaboration preferences, focus times).

**Structure:**

```markdown
### Writing style
I prefer terse bullet points over long paragraphs. Get to the point.

### Meeting scheduling  
Block mornings (9-12) for deep work. Schedule meetings after lunch.

### Code reviews
Focus on architecture questions, not syntax nitpicks. I trust the team on details.
```

**How it's used:**

- **Session start hook** shows top preferences
- Claude adapts tone/style based on these (via `System/user-profile.yaml` too)
- `/daily-plan` respects focus times when suggesting schedule

#### 3. Claude Code Updates

**Goal:** Detect when Anthropic ships new Claude Code features.

**Tech stack:**
- **Background automation** (`.scripts/meeting-intel/daily-synthesis.cjs`)
- **Anthropic changelog monitoring** (every 6 hours)
- **LaunchAgent** (`com.dex.meeting-intel.plist`)

**Flow:**

1. **LaunchAgent runs script** (every 6 hours)
2. **Script checks Anthropic changelog** (via web scraping or API)
3. **If new updates detected:**
   - Save to `System/changelog-updates-pending.md`
   - Flag shows in session start hook
4. **User runs** `/dex-whats-new`
5. **Skill reads changelog, analyzes:**
   - Which updates apply to Dex?
   - Should we update docs?
   - New capabilities to leverage?
6. **Skill updates:**
   - `CLAUDE.md` (if behavior changes)
   - `06-Resources/Dex_System/Dex_System_Guide.md` (if features added)
   - `System/claude-code-state.json` (tracks last reviewed version)

**User experience:** Dex tells you "New Claude features available!" and explains what changed.

### Learning Capture Workflow

**During `/daily-review`:**

1. Claude scans session transcript
2. Asks: "Anything to capture?"
   - Mistakes or corrections
   - Preferences mentioned
   - Doc gaps discovered
   - Workflow inefficiencies
3. Writes to `System/Session_Learnings/YYYY-MM-DD.md`:

```markdown
## [14:32] - Folder structure confusion

**What happened:** User expected `00-Inbox/` and docs now consistently use `00-Inbox/`
**Why it matters:** Consistent paths prevent onboarding flow issues
**Suggested fix:** Audit FOLDER_STRUCTURE.md and CLAUDE.md for consistency
**Status:** pending
```

4. **Weekly review** consolidates session learnings:
   - Move to appropriate file (`Mistake_Patterns.md`, `Working_Preferences.md`, or fix immediately)
   - Mark as resolved in session learnings

### Background Automation Setup

Dex runs self-learning checks automatically through two mechanisms:

#### Inline Checks (Default - No Installation)

The system runs checks automatically during:
1. **Session start hook** (`.claude/hooks/session-start.sh`)
2. **During `/daily-plan`** command

**Smart throttling:**
- Changelog check: Only runs if 6+ hours since last check
- Learning review: Only runs once per day
- Both run in background (non-blocking, <1 second)
- Respect intervals even if triggered from multiple places

**State tracking:**
- `System/claude-code-state.json` - Tracks last changelog check, Claude version, features discovered
- `System/.last-learning-check` - Tracks last daily learning review

#### Optional: Launch Agent Installation (Background Optimization)

For faster execution without inline checks, install macOS Launch Agents:

```bash
# Install background automation (optional optimization)
bash .scripts/install-learning-automation.sh
```

**What it does:**
- Runs changelog check every 6 hours in background
- Runs learning review daily at 5pm
- Creates alert files ready before session start
- Reduces latency during session start and `/daily-plan`

**With Launch Agents:**
- Checks run continuously in background
- Alert files ready before you even start planning
- Lower latency during session start

**Without Launch Agents:**
- Checks run inline during session start and `/daily-plan`
- Still fast (<1 second) with interval throttling
- System works perfectly fine, just slightly more latency

**Uninstall:**
```bash
bash .scripts/install-learning-automation.sh --uninstall
```

**Manual testing:**
```bash
node .scripts/check-anthropic-changelog.cjs --force
bash .scripts/learning-review-prompt.sh
```

**Alert files created:**
- `System/changelog-updates-pending.md` - When new Claude features detected
- `System/learning-review-pending.md` - When 5+ pending learnings exist

### Dex System Improvement Backlog

**Purpose:** Systematically capture and prioritize improvements to Dex itself.

#### Workflow

1. **Capture** - Use `capture_idea` MCP tool from any context
   ```
   User: "This would be better if X"
   Claude: [calls capture_idea tool with description]
   ```

2. **Storage** - Ideas saved to `System/Dex_Backlog.md` with metadata:
   ```markdown
   ## Idea: [Title]
   **Status:** pending
   **Priority:** [High/Medium/Low]
   **Captured:** YYYY-MM-DD
   **Description:** [What user said]
   **Rationale:** [Why it matters]
   ```

3. **Ranking** - AI scores ideas on 5 dimensions (via `/dex-backlog`):
   - **Impact (35%)** - Daily workflow improvement potential
   - **Alignment (20%)** - Fits your usage patterns and needs
   - **Token Efficiency (20%)** - Reduces context/token usage
   - **Memory & Learning (15%)** - Enhances persistence, self-learning, compounding knowledge
   - **Proactivity (10%)** - Enables proactive concierge behavior

4. **Review** - Run `/dex-backlog` to see ranked priorities
   - High: 85+ (implement soon)
   - Medium: 60-84 (consider for next cycle)
   - Low: <60 (backlog)

5. **Workshop** - Run `/dex-improve [idea]` to plan implementation
   - Analyzes feasibility
   - Creates implementation plan
   - Suggests file changes

#### Cursor Feasibility Gate

Ideas must be implementable using Cursor's actual capabilities:
- ✅ File operations (read, write, search)
- ✅ MCP tools and servers
- ✅ Command/skill creation
- ✅ Hook scripts
- ❌ Edit tracking or change detection
- ❌ Internal event listeners
- ❌ Real-time UI modifications

Ideas requiring unavailable capabilities are rejected with explanation.

#### Automatic Integration

- **Weekly planning** checks for high-priority ideas
- **Quarterly reviews** assess implementation progress
- `/dex-level-up` mentions idea capture capability

**Note:** Effort is intentionally excluded from scoring. With AI coding, implementation is cheap. Focus on value and feasibility.

---

## Design Constraints

### What Cursor/Claude Code Can't Do

Understanding these constraints explains why Dex is designed the way it is.

#### 1. No Edit Tracking

**Constraint:** Claude can't "see" when you manually edit a file. If you change `03-Tasks/Tasks.md` directly, Claude doesn't know until it reads the file again.

**Implication:** Can't do "smart diffing" or "real-time sync." Must rely on explicit reads.

**Dex's approach:** MCP servers are "write APIs" - they handle edits deterministically. If you manually edit, just tell Claude "task X is done" and MCP syncs everywhere.

#### 2. No Internal Hooks

**Constraint:** Can't trigger logic "when file changes" or "when task completes." No event listeners.

**Implication:** Can't do automatic background syncing without external orchestration.

**Dex's approach:** Use macOS LaunchAgents for scheduled tasks (e.g., Granola sync every 24 hours).

#### 3. Limited Memory Between Sessions

**Constraint:** Each chat session starts fresh. Claude only remembers what's in:
- `CLAUDE.md` (loaded every session)
- Session start hook output
- Files it explicitly reads

**Implication:** Can't "accumulate knowledge" across sessions without writing to files.

**Dex's approach:**
- Session learnings → written to files
- Strategic context → injected via session start hook
- Preferences → stored in `System/user-profile.yaml`

#### 4. No Background Processes

**Constraint:** Claude isn't a daemon. It only runs when you're chatting.

**Implication:** Can't "auto-sync Granola every hour" from within Cursor.

**Dex's approach:** Use macOS LaunchAgents (`.scripts/meeting-intel/`) for background automation. Claude processes the results when you check in.

### What Cursor/Claude Code Is Good At

#### 1. File Operations

**Strength:** Reading, writing, searching files is fast and reliable.

**Dex's leverage:** Everything is files. Tasks, meetings, people, goals - all markdown. Claude navigates this effortlessly.

#### 2. Natural Language Understanding

**Strength:** Claude can parse messy input ("I finished the API thing with John") and map to structured operations.

**Dex's leverage:** User speaks naturally. Claude translates to MCP calls. No rigid forms or syntax.

#### 3. Contextual Reasoning

**Strength:** Claude can synthesize information across multiple files (meeting notes + person pages + calendar + tasks) to generate insights.

**Dex's leverage:** Daily plan isn't a template, it's a synthesis. "You have a 1:1 with John at 2pm. He mentioned being blocked on API auth in your last meeting. That task is still open. Want to prep a solution?"

#### 4. Tool Composition

**Strength:** Claude can chain tool calls (read file → call MCP → write file → read another file).

**Dex's leverage:** Skills orchestrate complex workflows. `/process-meetings` calls Granola MCP, parses output, creates tasks via Work MCP, updates person pages, all in one flow.

### Design Principles (Based on Constraints)

1. **Files as state, MCP as operations:** Files store data, MCP provides business logic
2. **Lazy loading:** Only load context when needed (skills, person pages)
3. **Idempotency:** MCP operations can be retried safely (e.g., `create_task` checks for duplicates)
4. **Explicit over implicit:** User says "mark task done" rather than system inferring completion
5. **External orchestration:** Use OS-level tools (LaunchAgents) for background work
6. **Progressive enhancement:** Core workflows work without integrations (Calendar, Granola optional)

---

## Appendix: Key Files Reference

### Core Configuration

- `CLAUDE.md` - Main AI behavior instructions
- `System/user-profile.yaml` - User preferences, company info, communication style
- `System/pillars.yaml` - Strategic pillars (focus areas)
- `.claude/settings.json` - Cursor settings (MCP server configs)

### Skills

- `.claude/skills/[skill-name]/SKILL.md` - All skills follow this structure
- `.claude/skills/_available/` - Role-specific skills (not loaded by default)

### Hooks

**Global (always active):**
- `.claude/hooks/person-context-injector.cjs` - Injects person context on file read
- `.claude/hooks/company-context-injector.cjs` - Injects company context on file read
- `.claude/hooks/session-start.sh` - Shows strategic context at session start
- `.claude/hooks/ensure-mcp-user-scope.cjs` - Validates MCP tool calls use correct vault scope

**Skill-scoped (active only during parent skill):**
- `.claude/hooks/post-meeting-person-update.cjs` - Updates person pages during `/process-meetings`
- `.claude/hooks/meeting-summary-generator.cjs` - Meeting summary stub (future)
- `.claude/hooks/career-evidence-capture.cjs` - Logs achievements during `/career-coach`
- `.claude/hooks/daily-plan-quick-ref.cjs` - Generates quickref after `/daily-plan`

**Standalone utility:**
- `.claude/hooks/maintenance.cjs` - Vault health checks (stale files, broken links, orphans)

### Agent Memory

- `.claude/memory/` - Per-agent state files (JSON)

### MCP Servers

- `core/mcp/task_server.py` - Work management (tasks, priorities, goals)
- `core/mcp/career_server.py` - Career development
- `core/mcp/resume_server.py` - Resume/LinkedIn building
- `core/mcp/dex_improvements_server.py` - System improvement ideas

### Background Automation

- `.scripts/meeting-intel/sync-from-granola.cjs` - Daily Granola sync
- `.scripts/meeting-intel/com.dex.meeting-intel.plist` - LaunchAgent config

### Templates

- `System/Templates/Person_Page.md` - Person page structure
- `System/Templates/Company.md` - Company page structure
- `System/Templates/Career_Evidence_*.md` - Career evidence templates

---

## Contributing & Customization

### Forking Dex

Want to adapt Dex for your use case? Key customization points:

1. **Pillars:** Edit `System/pillars.yaml` to match your focus areas
2. **Folder structure:** Keep PARA (it's battle-tested), but adjust subdirs (e.g., add `05-Areas/Customers/`)
3. **Skills:** Create role-specific skills in `.claude/skills/` (see `anthropic-skill-creator` skill)
4. **MCP servers:** Add integrations for your tools (Notion, Linear, etc.) - see `.claude/reference/mcp-servers.md`

### Best Practices

- **Don't edit canonical files manually:** Use MCP tools to ensure syncing works
- **Backup regularly:** `git commit` daily (or use `.claude/hooks/session-end.sh` to auto-commit)
- **Test in demo mode:** Use `System/Demo/` to experiment without touching real data
- **Document learnings:** Use `System/Session_Learnings/` to capture improvements

---

## Questions or Issues?

This guide assumes you're comfortable with:
- Python/JavaScript basics
- Command line tools
- Git version control
- MCP protocol concepts

If you get stuck, check:
- `06-Resources/Dex_System/Dex_System_Guide.md` (user-facing guide)
- `06-Resources/Dex_System/Dex_Jobs_to_Be_Done.md` (why each piece exists)
- `.claude/reference/mcp-servers.md` (MCP setup details)

---

*This guide is a living document. As Dex evolves and Claude Code gains new capabilities, this will be updated to reflect best practices.*
