# Debugging Failing MCP Servers in Claude Code

Research date: 2026-02-26

## Executive Summary

Claude Code stores per-session JSONL log files for each MCP server in `~/.cache/claude-cli-nodejs/<project-dir>/mcp-logs-<server-name>/`. The primary diagnostic tools are: reading these JSONL logs, running `claude mcp list` and `/mcp` inside Claude Code, running `/doctor`, and launching with `claude --debug`. The most common startup failures are connection timeouts (default 30s), stderr output misinterpreted as errors, Python/Node dependency crashes, and invalid JSON configuration.

## Source Inventory

| # | URL | Description |
|---|-----|-------------|
| 1 | <https://code.claude.com/docs/en/troubleshooting> | Official Claude Code troubleshooting docs (accessed 2026-02-26) |
| 2 | <https://code.claude.com/docs/en/mcp> | Official Claude Code MCP configuration docs (accessed 2026-02-26) |
| 3 | <https://code.claude.com/docs/en/settings> | Official Claude Code settings and environment variables (accessed 2026-02-26) |
| 4 | <https://github.com/anthropics/claude-code/issues/72> | GitHub issue #72: "How to Debug MCP Server that fails to connect?" (accessed 2026-02-26) |
| 5 | <https://github.com/anthropics/claude-code/issues/17653> | GitHub issue #17653: "MCP servers incorrectly shown as failed when stderr contains non-error output" (accessed 2026-02-26 via search excerpt) |
| 6 | <https://github.com/anthropics/claude-code/blob/main/plugins/plugin-dev/skills/mcp-integration/SKILL.md?plain=1#L443#debugging> | Official claude-code repo: MCP integration skill debugging section (accessed 2026-02-26) |
| 7 | <https://github.com/anthropics/claude-code/blob/main/plugins/plugin-dev/skills/mcp-integration/references/tool-usage.md?plain=1#L502#troubleshooting> | Official claude-code repo: MCP tool usage troubleshooting (accessed 2026-02-26) |
| 8 | <https://institute.sfeir.com/en/claude-code/claude-code-mcp-model-context-protocol/troubleshooting/> | SFEIR Institute MCP troubleshooting guide (accessed 2026-02-26) |
| 9 | Local filesystem observation | Direct observation of `~/.cache/claude-cli-nodejs/` on this Linux machine (2026-02-26) |

## 1. MCP Server Log File Locations (Linux)

### Primary log directory

```text
~/.cache/claude-cli-nodejs/<project-dir-encoded>/mcp-logs-<server-name>/
```

**Path encoding**: The project directory path is encoded by replacing `/` with `-` and stripping the leading dash. For example, launching Claude Code from `/home/ubuntulinuxqa2/repos/claude_skills` creates the cache directory:

```text
~/.cache/claude-cli-nodejs/-home-ubuntulinuxqa2-repos-claude-skills/
```

Within that directory, each MCP server gets its own subdirectory named `mcp-logs-<server-name>`:

```text
mcp-logs-exa/
mcp-logs-context7/
mcp-logs-git-xray/
mcp-logs-claude-ai-Ref/
mcp-logs-plugin-episodic-memory-episodic-memory/
mcp-logs-mcp-server-docker/
```

SOURCE: Confirmed by Anthropic engineer @ashwin-ant in GitHub issue #72 [4]: "MCP errors should appear within `~/Library/Caches/claude-cli-nodejs/` (if you're on a Mac) or `~/.cache/claude-cli` (if you're on Linux)." Validated on this system by direct observation [9] -- the actual path on this Linux system is `~/.cache/claude-cli-nodejs/` (with the `-nodejs` suffix).

### Platform-specific log paths

| Platform | Log directory |
|----------|---------------|
| Linux | `~/.cache/claude-cli-nodejs/<project-dir>/mcp-logs-<server>/` |
| macOS | `~/Library/Caches/claude-cli-nodejs/<project-dir>/mcp-logs-<server>/` |
| Windows | UNVERIFIED -- not documented in sources fetched |

### Log file format

Each session creates a JSONL file named with the UTC timestamp of session start:

```text
2026-02-26T15-07-20-951Z.jsonl
```

Each line is a JSON object with these fields:

```json
{
  "debug": "message text (or absent if error)",
  "error": "error message text (or absent if debug)",
  "timestamp": "2026-02-26T15:07:24.312Z",
  "sessionId": "ecf72873-baea-46b8-950c-19a6af6acf98",
  "cwd": "/home/ubuntulinuxqa2/repos/claude_skills"
}
```

### MCP configuration file locations

| File | Purpose | Source |
|------|---------|--------|
| `~/.claude.json` | Global state including user-scope MCP servers | [1] |
| `.mcp.json` (project root) | Project-scoped MCP servers (checked into source control) | [2] |
| `~/.claude/settings.json` | User settings including MCP permissions | [1] |
| `.claude/settings.json` | Project settings | [1] |
| `.claude/settings.local.json` | Local project settings (not committed) | [1] |
| `/etc/claude-code/managed-mcp.json` | Managed MCP servers (Linux system-wide, admin-deployed) | [2] |

SOURCE: Official troubleshooting docs [1], official MCP docs [2].

## 2. How to Read MCP Server Startup Errors from Logs

### Successful startup sequence (4 log lines)

```json
{"debug":"Starting connection with timeout of 30000ms","timestamp":"...","sessionId":"...","cwd":"..."}
{"error":"Server stderr: [informational output from server process]\n","timestamp":"...","sessionId":"...","cwd":"..."}
{"debug":"Successfully connected to stdio server in 2272ms","timestamp":"...","sessionId":"...","cwd":"..."}
{"debug":"Connection established with capabilities: {\"hasTools\":true,...}","timestamp":"...","sessionId":"...","cwd":"..."}
```

**Key observation**: The `"error"` field labeled "Server stderr" does NOT necessarily mean failure. Claude Code logs ALL stderr output from the MCP server process under the `"error"` key, even informational messages. This is a known bug documented in GitHub issue #17653 [5].

### Failed startup sequence (connection closed)

```json
{"debug":"Starting connection with timeout of 30000ms","timestamp":"...","sessionId":"...","cwd":"..."}
{"error":"Server stderr: Traceback (most recent call last):\n  File ... \nNameError: name 'Annotated' is not defined\n","timestamp":"...","sessionId":"...","cwd":"..."}
{"debug":"Connection failed after 2584ms: MCP error -32000: Connection closed","timestamp":"...","sessionId":"...","cwd":"..."}
{"error":"Connection failed: MCP error -32000: Connection closed","timestamp":"...","sessionId":"...","cwd":"..."}
```

### Failed startup sequence (timeout)

```json
{"debug":"Starting connection with timeout of 30000ms","timestamp":"...","sessionId":"...","cwd":"..."}
{"error":"Connection failed: Connection to MCP server \"<name>\" timed out after 5000ms","timestamp":"...","sessionId":"...","cwd":"..."}
```

SOURCE: Direct observation of log files on this system [9], GitHub issue #72 [4].

### How to read the logs

```bash
# Find the most recent log file for a specific server
ls -lt ~/.cache/claude-cli-nodejs/<project-dir>/mcp-logs-<server-name>/ | head -5

# Read the most recent log (JSONL format, one JSON object per line)
cat ~/.cache/claude-cli-nodejs/<project-dir>/mcp-logs-<server-name>/<latest-file>.jsonl

# Pretty-print for readability
cat <log-file>.jsonl | python3 -c "import sys,json; [print(json.dumps(json.loads(l),indent=2)) for l in sys.stdin]"

# Filter for errors only
cat <log-file>.jsonl | python3 -c "
import sys, json
for line in sys.stdin:
    obj = json.loads(line)
    if 'error' in obj:
        print(json.dumps(obj, indent=2))
"

# Search all log files for connection failures
grep -r 'Connection failed' ~/.cache/claude-cli-nodejs/<project-dir>/mcp-logs-*/
```

## 3. Commands and Files That Show MCP Server Status/Errors

### In-session commands (interactive REPL)

| Command | What it shows | Source |
|---------|---------------|--------|
| `/mcp` | Lists all MCP servers, their connection status, available tools, and authentication status | [2] |
| `/doctor` | Checks MCP server configuration errors, malformed JSON, unreachable permission rules, high MCP token usage, plugin/agent loading errors | [1] |

SOURCE: Official troubleshooting docs [1]: "/doctor... checks: MCP server configuration errors."

### CLI commands (outside Claude Code)

```bash
# List all configured MCP servers
claude mcp list

# Get details for a specific server
claude mcp get <server-name>
```

SOURCE: Official MCP docs [2].

### Startup flags for debugging

```bash
# Full debug mode (replaces deprecated --mcp-debug)
claude --debug

# Debug with category filtering
claude --debug "mcp"

# Write debug logs to a file
claude --debug-file /tmp/claude-debug.log

# Deprecated but still functional
claude --mcp-debug
```

SOURCE: `claude --help` output observed on this system; GitHub issue #72 [4] confirming `--mcp-debug` was added in March 2025 and later deprecated in favor of `--debug`.

### Environment variables for MCP debugging

| Variable | Purpose | Default | Source |
|----------|---------|---------|--------|
| `MCP_TIMEOUT` | Timeout in milliseconds for MCP server startup | 30000 (30s) | [3] |
| `MCP_TOOL_TIMEOUT` | Timeout in milliseconds for MCP tool execution | Not documented | [3] |
| `MAX_MCP_OUTPUT_TOKENS` | Maximum tokens allowed in MCP tool responses (warning at 10,000) | 25000 | [3] |
| `MCP_CLIENT_SECRET` | OAuth client secret for MCP servers requiring pre-configured credentials | None | [3] |

```bash
# Increase startup timeout to 60 seconds
MCP_TIMEOUT=60000 claude

# Increase output token limit
MAX_MCP_OUTPUT_TOKENS=50000 claude
```

SOURCE: Official settings docs [3].

## 4. Common Causes of MCP Server Startup Failures

### 4.1 Connection timeout (most common)

The MCP server process takes longer than `MCP_TIMEOUT` (default 30s) to complete initialization and respond to the MCP `initialize` handshake.

**Symptoms in logs**:

```text
Connection failed: Connection to MCP server "<name>" timed out after 30000ms
```

**Common causes**:

- `npx` downloading packages on first run (observed in issue #72 [4])
- Python `uv run` installing dependencies at startup (observed on this system [9])
- Slow network connectivity for remote package registries
- Heavy server initialization logic

**Fix**: Increase timeout:

```bash
MCP_TIMEOUT=60000 claude
```

### 4.2 Server process crash (stderr traceback)

The MCP server process exits before completing the MCP handshake.

**Symptoms in logs**:

<console>
Server stderr: Traceback (most recent call last): ...
Connection failed after 2584ms: MCP error -32000: Connection closed
</console>

**Common causes**:

- Python import errors / missing dependencies
- Node.js module resolution failures
- Incompatible dependency versions (e.g., Pydantic version mismatch)
- Type annotation errors (e.g., `NameError: name 'Annotated' is not defined`)

**Fix**: Run the server command manually outside Claude Code to see the full error:

```bash
# For a stdio server, run its command directly:
npx -y @some/mcp-server
# or
uv run python server.py
```

SOURCE: Observed in kaizen-analysis server logs on this system [9]; issue #72 [4].

### 4.3 Stderr output misinterpreted as errors (false positive)

Claude Code reports "X MCP servers failed" even though the servers are running correctly, because ANY output to stderr is flagged as an error in the UI.

**Symptoms**: UI shows failure count, but `/mcp` shows servers are connected. Log files show `"error": "Server stderr: ..."` entries containing informational messages like version banners, deprecation warnings, or progress indicators.

SOURCE: GitHub issue #17653 [5] (open bug as of 2026-02-12).

### 4.4 Invalid JSON configuration

Malformed JSON in `.mcp.json`, `~/.claude.json`, or other config files prevents servers from loading.

**Symptoms**: `/doctor` reports configuration errors. Servers do not appear in `/mcp` listing.

**Fix**:

```bash
# Validate project MCP config
python3 -m json.tool .mcp.json

# Validate global config
python3 -m json.tool ~/.claude.json
```

SOURCE: Official troubleshooting docs [1]; SFEIR guide [8].

### 4.5 Incorrect binary path or missing executable

The `command` field in MCP server configuration points to a binary that does not exist or is not executable.

**Symptoms in logs**:

```text
Connection failed: spawn <command> ENOENT
```

**Fix**:

```bash
# Verify the command exists
which npx
which node
which uv

# Check execute permissions
ls -la /path/to/server
chmod +x /path/to/server
```

SOURCE: SFEIR guide [8]; official MCP docs [2].

### 4.6 Environment variable misconfiguration

MCP servers requiring API keys or tokens fail when environment variables are not set or not passed through.

**Symptoms**: Server starts but returns authentication errors, or fails to connect to external services.

**Fix**: Use `--env` flag when adding servers:

```bash
claude mcp add --transport stdio --env API_KEY=xxx my-server -- npx my-server
```

Or set in config JSON:

```json
{
  "mcpServers": {
    "my-server": {
      "command": "npx",
      "args": ["-y", "my-server"],
      "env": {
        "API_KEY": "xxx"
      }
    }
  }
}
```

SOURCE: Official MCP docs [2].

### 4.7 Port conflicts (SSE/HTTP transport)

For servers using SSE or HTTP transport, the configured port may already be in use.

**Fix**:

```bash
# Check what is using the port
lsof -i :<port>
# or
ss -tlnp | grep <port>
```

SOURCE: SFEIR guide [8].

## 5. Diagnostic Procedure: Identifying Which MCP Servers Are Failing and Why

### Step-by-step procedure

```text
Step 1: Run /doctor inside Claude Code
    - Reports MCP server configuration errors
    - Reports malformed JSON, incorrect types
    - Reports plugin and agent loading errors
    SOURCE: [1]

Step 2: Run /mcp inside Claude Code
    - Shows each server with connection status
    - Shows available tools per server
    - Plugin servers appear with plugin indicators
    SOURCE: [2]

Step 3: Check the JSONL log files
    - Navigate to ~/.cache/claude-cli-nodejs/<project-dir>/
    - List mcp-logs-* directories
    - Read the most recent .jsonl file in each failing server's directory
    - Look for:
      * "Connection failed" entries (definitive failure)
      * "Server stderr" entries containing tracebacks or error messages
      * "timed out" messages
    SOURCE: [4], [9]

Step 4: Run the server command manually
    - Extract the command from the config (claude mcp get <name>)
    - Run it in a terminal to see raw stdout/stderr
    SOURCE: [6], [8]

Step 5: Enable debug mode for detailed logs
    - claude --debug
    - Or: claude --debug-file /tmp/mcp-debug.log
    SOURCE: [4], CLI help output
```

### Quick diagnostic commands

```bash
# List all configured servers
claude mcp list

# Get config details for a specific server
claude mcp get <server-name>

# Find all recent error logs across all servers for current project
grep -r "Connection failed" ~/.cache/claude-cli-nodejs/-$(pwd | tr '/' '-')/mcp-logs-*/ 2>/dev/null

# Find servers that timed out
grep -r "timed out" ~/.cache/claude-cli-nodejs/-$(pwd | tr '/' '-')/mcp-logs-*/ 2>/dev/null

# List log directories to see which servers have logs
ls -d ~/.cache/claude-cli-nodejs/-$(pwd | tr '/' '-')/mcp-logs-*/

# Read the most recent log for a specific server
ls -t ~/.cache/claude-cli-nodejs/-$(pwd | tr '/' '-')/mcp-logs-<server>/*.jsonl | head -1 | xargs cat
```

### Lifecycle events visible in logs

The JSONL logs capture the full MCP server lifecycle:

1. **Startup**: `"Starting connection with timeout of 30000ms"`
2. **Stderr capture**: `"Server stderr: ..."` (all stderr output, not only errors)
3. **Connection success**: `"Successfully connected to stdio server in Xms"`
4. **Capabilities**: `"Connection established with capabilities: {...}"`
5. **Shutdown**: `"Sending SIGINT to MCP server process"` then `"SIGINT failed, sending SIGTERM"` then `"MCP server process exited cleanly"`
6. **Failure**: `"Connection failed after Xms: MCP error -32000: Connection closed"` or `"Connection to MCP server ... timed out after Xms"`

### Known bug: false failure reporting

GitHub issue #17653 [5] documents that Claude Code interprets ANY stderr output as an error. Servers that write informational messages to stderr (version banners, deprecation warnings, progress bars) trigger the "X MCP servers failed" UI message even when the servers are fully functional. The workaround is to check `/mcp` in-session or read the JSONL logs to verify whether the server actually connected successfully (look for the `"Successfully connected"` debug line).

## Appendix: Observed Log Directory Structure on This System

```text
~/.cache/claude-cli-nodejs/
  -home-ubuntulinuxqa2-repos-claude-skills/
    mcp-logs-claude-ai-Atlassian/
    mcp-logs-claude-ai-Canva/
    mcp-logs-claude-ai-Github-mcp/
    mcp-logs-claude-ai-Gmail/
    mcp-logs-claude-ai-Hex/
    mcp-logs-claude-ai-Hugging-Face/
    mcp-logs-claude-ai-Jam/
    mcp-logs-claude-ai-Mermaid-Chart/
    mcp-logs-claude-ai-Notion/
    mcp-logs-claude-ai-Ref/
    mcp-logs-claude-in-chrome/
    mcp-logs-context7/
    mcp-logs-exa/
    mcp-logs-git-forensics/
    mcp-logs-GitLab-Public-API/
    mcp-logs-git-xray/
    mcp-logs-ide/
    mcp-logs-kernel/
    mcp-logs-mcp-json-yaml-toml/
    mcp-logs-mcp-server-docker/
    mcp-logs-plugin-agentskill-kaizen-kaizen-analysis/
    mcp-logs-plugin-agentskill-kaizen-kaizen-duckdb/
    mcp-logs-plugin-context7-context7/
    mcp-logs-plugin-episodic-memory-episodic-memory/
    mcp-logs-plugin-Notion-notion/
    mcp-logs-plugin-serena-serena/
```

Each directory contains timestamped `.jsonl` files (one per session), e.g.:

```text
2026-02-26T15-07-20-951Z.jsonl
2026-02-26T14-54-43-316Z.jsonl
```
