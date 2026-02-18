# Virtual MCP Server - How It Works

This document explains how Virtual MCP Servers work using diagrams and examples. For detailed implementation specifics, see [virtual-mcp-server.md](virtual-mcp-server.md).

---

## What Problem Are We Solving?

Consider a typical development setup: you have separate MCP servers for GitHub (code search, PRs), Slack (messaging), and Jira (issue tracking). Your AI agent needs tools from all three, which means:

- Managing three separate connections
- Handling three different sessions
- Dealing with tool name conflicts (both GitHub and Jira have a `search` tool)

A Virtual MCP Server solves this by providing a **single endpoint** that aggregates tools from multiple backends. Your agent connects once and gets access to all the tools it needs.

```
WITHOUT Virtual Server:              WITH Virtual Server:

  You                                  You
   |                                    |
   +---> GitHub Server                  |
   |        |-> search                  v
   |        |-> create_pr         +------------+
   |                              |  Virtual   |
   +---> Slack Server             |  Server    |
   |        |-> send_message      +-----+------+
   |        |-> list_channels           |
   |                              +-----+-----+-----+
   +---> Jira Server              |           |     |
            |-> create_issue      v           v     v
            |-> search_issues   GitHub     Slack   Jira
                                Server     Server  Server
```

**Benefits:**
- Your app only connects to ONE server instead of many
- You can pick exactly which tools you want from each backend
- You can rename tools to avoid confusion (like "github_search" vs "jira_search")
- You can control who has access to which tools

---

## The Big Picture

Request flow when a client connects to a Virtual MCP Server:

```
+----------------+                    +------------------+
|   Your App     |                    |   MCP Gateway    |
|                |                    |                  |
|  "I want to    | ---(1) Request --> |  Nginx receives  |
|   search on    |                    |  your request    |
|   GitHub"      |                    +--------+---------+
|                |                             |
|                |                             v
|                |                    +------------------+
|                |                    |  Lua Router      |
|                |                    |  (the brain)     |
|                |                    |                  |
|                |                    |  "Ah, this tool  |
|                |                    |   belongs to     |
|                |                    |   GitHub backend"|
|                |                    +--------+---------+
|                |                             |
|                |                             v
|                |                    +------------------+
|                |                    |  GitHub Backend  |
|                |                    |                  |
|                | <--(4) Response -- |  (does the       |
|                |                    |   actual work)   |
+----------------+                    +------------------+
```

Each component is described below.

---

## The Three Key Players

### 1. Nginx (Reverse Proxy)

Nginx receives incoming requests and handles:

- JWT authentication via `auth_request` subrequest
- Path-based routing to determine which virtual server
- Invoking the Lua content handler for MCP protocol processing

```
Request arrives at /virtual/dev-tools
                |
                v
        +---------------+
        |    Nginx      |
        |               |
        |  1. Check JWT |  <-- "Is this token valid?"
        |  2. Read path |  <-- "Which virtual server?"
        |  3. Call Lua  |  <-- "Hand off to the router"
        +---------------+
```

### 2. Lua Router (Content Handler)

The Lua router (`virtual_router.lua`) runs as an nginx content handler. It:

- Reads tool-to-backend mappings from JSON config files
- Translates tool aliases back to original names
- Manages session multiplexing across backends
- Issues concurrent subrequests for aggregation methods

### 3. Backend Servers

The actual MCP servers (GitHub, Slack, Jira, etc.) that execute tool calls. The virtual server coordinates requests but delegates all execution to backends.

---

## How Tool Mapping Works

This is the core mechanism. The process works as follows:

### Step 1: Configuration is Created

When someone creates a virtual server, they specify which tools to include:

```
Virtual Server: "dev-tools"
Path: /virtual/dev-tools

Tool Mappings:
  +------------------+------------------+------------------+
  | Tool Name        | Backend Server   | Alias            |
  +------------------+------------------+------------------+
  | search           | /github          | github_search    |
  | search           | /jira            | jira_search      |
  | send_message     | /slack           | (none - use as-is)|
  +------------------+------------------+------------------+
```

Both GitHub and Jira have a tool called "search". Aliases resolve this naming conflict.

### Step 2: Mapping File is Generated

The system writes a JSON file that the Lua router will read:

```
File: /etc/nginx/lua/virtual_mappings/dev-tools.json

{
  "tool_backend_map": {
    "github_search": {
      "original_name": "search",
      "backend_location": "/_backend/github"
    },
    "jira_search": {
      "original_name": "search",
      "backend_location": "/_backend/jira"
    },
    "send_message": {
      "original_name": "send_message",
      "backend_location": "/_backend/slack"
    }
  }
}
```

This file is a lookup table mapping tool names to their backend locations.

### Step 3: Request Comes In

When your app calls a tool:

```
Your app sends:
{
  "method": "tools/call",
  "params": {
    "name": "github_search",      <-- The alias you see
    "arguments": { "query": "bug fixes" }
  }
}
```

### Step 4: Lua Router Translates

The Lua router:
1. Reads the mapping file
2. Looks up "github_search"
3. Finds: backend is "/_backend/github", original name is "search"
4. Rewrites the request:

```
Forwarded to /_backend/github:
{
  "method": "tools/call",
  "params": {
    "name": "search",             <-- Original name the backend knows
    "arguments": { "query": "bug fixes" }
  }
}
```

### Step 5: Response Goes Back

The GitHub backend responds. The Lua router passes it back to your app unchanged.

---

## The Complete Request Flow (Sequence Diagram)

Sequence diagram for a `tools/call` request:

```
Your App          Nginx           Lua Router        Backend
   |                |                  |               |
   |  POST /virtual/dev-tools          |               |
   |  tools/call: github_search        |               |
   |--------------->|                  |               |
   |                |                  |               |
   |                |  auth_request    |               |
   |                |  (check JWT)     |               |
   |                |----------------->|               |
   |                |  OK + scopes     |               |
   |                |<-----------------|               |
   |                |                  |               |
   |                |  content_by_lua  |               |
   |                |----------------->|               |
   |                |                  |               |
   |                |      Read mapping file           |
   |                |      "github_search" ->          |
   |                |        backend: /_backend/github |
   |                |        original: search          |
   |                |                  |               |
   |                |      Check session cache         |
   |                |      (do we have a session       |
   |                |       with this backend?)        |
   |                |                  |               |
   |                |      Rewrite tool name           |
   |                |      github_search -> search     |
   |                |                  |               |
   |                |                  | POST to       |
   |                |                  | /_backend/github
   |                |                  |-------------->|
   |                |                  |               |
   |                |                  |   Response    |
   |                |                  |<--------------|
   |                |                  |               |
   |<----------------------------------|               |
   |      Response                     |               |
```

---

## Session Management (The Tricky Part)

Each backend server requires its own session.

When your app connects to the virtual server, it gets ONE session ID:

```
Your app <---> Virtual Server (session: vs-abc123)
```

But behind the scenes, the virtual server maintains SEPARATE sessions with each backend:

```
Virtual Server:
  +-- Session with GitHub: sess-gh-001
  +-- Session with Slack:  sess-sl-002
  +-- Session with Jira:   sess-jr-003
```

The Lua router keeps track of this mapping so you don't have to.

### The Two-Tier Cache

Looking up sessions from the database on every request would be slow. So we use two levels of caching:

```
Request: "What's the GitHub session for client vs-abc123?"

        +-------------------+
        |  Level 1 Cache    |  <-- Super fast (in nginx memory)
        |  (Shared Dict)    |      TTL: 30 seconds
        |                   |
        |  "Do I have it?"  |
        +--------+----------+
                 |
        MISS     |
                 v
        +-------------------+
        |  Level 2 Cache    |  <-- Fast (MongoDB lookup)
        |  (MongoDB)        |      TTL: 1 hour
        |                   |
        |  "Check database" |
        +--------+----------+
                 |
        MISS     |
                 v
        +-------------------+
        |  Create New       |  <-- Send "initialize" to backend
        |  Session          |      Store in both caches
        +-------------------+
```

**Why two levels?**
- Level 1 is in memory - no network call, extremely fast
- Level 2 is in MongoDB - survives server restarts
- If the server restarts, we lose Level 1 but Level 2 still has sessions

---

## Listing Tools (Aggregation)

When your app asks "what tools do you have?", the virtual server needs to ask ALL backends:

```
Your app asks: tools/list

Lua Router:
  +-- Ask GitHub: "What tools do you have?"
  |     Response: [search, create_pr, list_repos]
  |
  +-- Ask Slack: "What tools do you have?"
  |     Response: [send_message, list_channels]
  |
  +-- Ask Jira: "What tools do you have?"
        Response: [create_issue, search_issues]

Lua Router combines them:
  [github_search, create_pr, list_repos,     <-- Applied aliases
   send_message, list_channels,
   jira_search, create_issue]                 <-- Renamed "search" to "jira_search"
```

**Important optimization:** These backend calls happen IN PARALLEL, not one after another. This makes the aggregation fast.

```
Time ------>

Sequential (slow):
  [GitHub call]---[Slack call]---[Jira call]---Done

Parallel (fast):
  [GitHub call]-----
  [Slack call]------+---Done
  [Jira call]-------
```

---

## What the Nginx Config Looks Like

When a virtual server is enabled, the system generates two things:

### 1. A Location Block (for routing)

```nginx
location /virtual/dev-tools {
    # Tell Lua which virtual server this is
    set $virtual_server_id "dev-tools";

    # Check authentication first
    auth_request /validate;

    # Run the Lua router
    content_by_lua_file /etc/nginx/lua/virtual_router.lua;
}
```

### 2. Internal Backend Locations

```nginx
# These are marked "internal" - only Lua can use them
# Regular users can't access them directly

location /_backend/github {
    internal;
    proxy_pass https://github-mcp.example.com/mcp;
}

location /_backend/slack {
    internal;
    proxy_pass https://slack-mcp.example.com/mcp;
}
```

The Lua router uses these internal locations to talk to backends.

---

## Error Handling

What happens when things go wrong?

### Backend is Down

```
Lua Router tries to call /_backend/github
  |
  v
Connection fails or returns error
  |
  v
Lua Router returns error to your app:
{
  "error": {
    "code": -32000,
    "message": "Backend server unreachable: /github"
  }
}
```

### Session Expired

```
Lua Router uses cached session sess-gh-001
  |
  v
GitHub returns: "400 Bad Request - Invalid session"
  |
  v
Lua Router:
  1. Delete sess-gh-001 from both caches
  2. Send new "initialize" to GitHub
  3. Get new session: sess-gh-002
  4. Cache it in both levels
  5. Retry the original request with new session
```

### User Lacks Permission

```
User has scopes: ["mcp-access"]
Tool "create_pr" requires: ["github-write"]
  |
  v
Lua Router checks scopes... DENIED
  |
  v
Response: 403 Forbidden
{
  "error": "Missing required scope: github-write"
}
```

---

## Access Control in Simple Terms

Access control works at two levels:

### Level 1: Server Access

To use the virtual server at all, you need certain scopes:

```
Virtual Server: /virtual/dev-tools
Required Scopes: ["mcp-access"]

User with scopes ["mcp-access"] -> Allowed in
User with scopes ["other-stuff"] -> Blocked at the door
```

### Level 2: Tool Access

Individual tools can require additional scopes:

```
Tool: create_pr
Required Scopes: ["github-write"]

User with ["mcp-access", "github-read"] -> Can't use this tool
User with ["mcp-access", "github-write"] -> Can use this tool
```

When listing tools, the Lua router hides tools the user can't access:

```
Full tool list:    [search, create_pr, delete_repo]
User scopes:       ["mcp-access", "github-read"]

Filtered list:     [search]  <-- Only shows tools user can actually use
```

---

## How Changes Are Applied

When you create or update a virtual server:

```
1. You call the API: POST /api/virtual-servers
          |
          v
2. Service validates the configuration
   - Does each backend server exist?
   - Does each tool exist on its backend?
   - Are all alias names unique?
          |
          v
3. Configuration saved to MongoDB
          |
          v
4. Nginx config regenerated
   - New location block written
   - New mapping JSON file written
          |
          v
5. Nginx reloaded
   - nginx -s reload
   - New config takes effect immediately
          |
          v
6. Virtual server is live!
```

---

## Quick Reference

### Files You Should Know

| File | What It Does |
|------|-------------|
| `virtual_router.lua` | The Lua brain that routes requests |
| `nginx_service.py` | Generates nginx config + mapping files |
| `virtual_server_service.py` | Business logic and validation |
| `virtual_server_routes.py` | REST API endpoints |
| `/etc/nginx/lua/virtual_mappings/*.json` | Tool mapping files read by Lua |

### Key Concepts

| Term | Plain English |
|------|--------------|
| Virtual Server | A fake server that coordinates real servers |
| Tool Mapping | "This tool comes from that backend" |
| Alias | A renamed tool to avoid confusion |
| Backend Location | Where to forward requests (internal nginx path) |
| Session Multiplexing | One client session, many backend sessions |
| Scope | A permission string that controls access |

### Common Operations

| What You Want | What Happens |
|---------------|--------------|
| List tools | Asks all backends in parallel, combines results |
| Call a tool | Looks up backend, translates name, forwards request |
| Initialize | Creates client session, backend sessions are lazy |
| Ping | Responds immediately, no backend calls |

---

## Summary

1. **Virtual servers aggregate tools** from multiple backends into one endpoint
2. **Nginx routes requests** to the Lua router based on path
3. **Lua router reads mapping files** to know which tool goes where
4. **Aliases solve naming conflicts** when two backends have same tool names
5. **Sessions are cached in two levels** for speed and reliability
6. **Access control works at server and tool level** using scopes
7. **Backend calls happen in parallel** when listing tools

The virtual server acts as a coordinator - all tool execution happens on the backend servers. The virtual server's role is to present a unified endpoint to clients.
