# FastMCP Community Practices

Best practices for FastMCP 3.x (2025-2026), including packaging, security, observability, testing, performance, and ecosystem patterns from the MCP developer community.

## Overview

Since 2025-2026, the MCP developer community (spanning OpenAI, Anthropic, Hugging Face, and independent developers) has converged on best practices that make MCP servers more reliable, secure, and easier to integrate. FastMCP 3.x builds on these patterns with native support for many community-requested features.

## One-Click Deployment & Packaging

### Claude Desktop Extensions (.mcpb Format)

Anthropic introduced **Claude Desktop Extensions** - packaged MCP servers with all dependencies included.

#### The .mcpb (MCP Bundle) Format

A `.mcpb` file is a ZIP containing:

- Server code (Python/Node)
- Manifest JSON with metadata
- All required libraries and dependencies

#### Installation Experience

Users can install by:

1. Double-clicking the `.mcpb` file
2. No Python or Node setup required
3. Non-technical users can use MCP servers easily

### Manifest Specification

The manifest JSON specifies:

```json
{
  "name": "my-mcp-server",
  "version": "1.0.0",
  "description": "Short description",
  "runtime": "python",
  "command": "python server.py",
  "user_config": {
    "api_key": {
      "type": "string",
      "description": "API key for service",
      "required": true
    },
    "allowed_directories": {
      "type": "array",
      "description": "Allowed file system paths",
      "default": []
    }
  },
  "compatibility": {
    "min_version": "1.0.0"
  },
  "metadata": {
    "author": "Your Name",
    "homepage": "https://github.com/...",
    "screenshots": ["screenshot1.png"]
  }
}
```

### Manifest Features

- **Semver versioning** - Version compatibility checks
- **Compatibility requirements** - "requires Claude Desktop >=1.0"
- **User configuration prompts** - For API keys, paths, etc.
- **Screenshots and docs** - For better discovery
- **Runtime specification** - Python, Node, or other

### Creating .mcpb Bundles

1. Create manifest.json with required fields
2. Bundle server code + manifest + dependencies
3. ZIP everything together
4. Rename .zip to .mcpb

**Community adoption**: Many open-source servers now provide `.mcpb` downloads for easy installation in Claude Desktop or VS Code.

### Best Practice

If targeting desktop AI apps, offer your server in `.mcpb` format. Makes integration as smooth as installing a browser extension.

## Secure by Design

MCP's philosophy: **Expose only what's intended**.

### Path Restrictions

If server interacts with file system:

1. **Require user configuration** of allowed paths
2. **Pass via environment variable** or manifest field
3. **Enforce in code** - Reject file access outside allowed paths

**Example manifest configuration:**

```json
{
  "user_config": {
    "allowed_directories": {
      "type": "array",
      "description": "File paths this server can access",
      "required": true
    }
  }
}
```

**Example server enforcement:**

```python
import os
from fastmcp.exceptions import ToolError

ALLOWED_DIRS = os.getenv("ALLOWED_DIRECTORIES", "").split(",")

@mcp.tool()
def read_file(path: str) -> dict:
    """Read a file from allowed directories."""
    # Normalize and check path
    abs_path = os.path.abspath(path)
    if not any(abs_path.startswith(os.path.abspath(d)) for d in ALLOWED_DIRS):
        raise ToolError(
            f"Access denied: {path} is not in allowed directories. "
            f"Configure ALLOWED_DIRECTORIES environment variable."
        )

    # Safe to read
    with open(abs_path) as f:
        return {"content": f.read()}
```

### Destructive Hint Annotations

Mark tools that modify state:

```python
@mcp.tool(
    annotations={
        "destructiveHint": True,  # Requires user confirmation
        "openWorldHint": True    # Accesses external systems
    }
)
def delete_file(path: str) -> dict:
    """Delete a file."""
    # ...
```

Clients may prompt: "Are you sure?" before the AI deletes files.

### Capability Whitelisting

- **Scope API keys** to least privilege (read-only if possible)
- **Handle access denials gracefully**
- **Fail safely** when permissions missing

### Safety Checks Within Tools

Assume AI might misuse a tool - build in guardrails:

```python
@mcp.tool()
def send_email(
    to: str,
    subject: str,
    body: str,
    confirm: Annotated[bool, Field(description="Set to true to confirm sending")] = False
) -> dict:
    """Send an email (requires confirmation)."""
    if not confirm:
        raise ToolError(
            "Email not sent. Set confirm=true to actually send the email."
        )

    # Additional validation
    if len(to.split(",")) > 10:
        raise ToolError(
            "Cannot send to more than 10 recipients at once. "
            "This prevents accidental mass emails."
        )

    # Safe to send
    send_email_internal(to, subject, body)
    return {"status": "sent", "to": to}
```

### FastMCP Validation as Defense

Use Pydantic Field constraints as first line of defense:

```python
# Prevent suspiciously long inputs
command: Annotated[str, Field(max_length=1000)]

# Whitelist allowed values
action: Annotated[Literal["read", "list", "search"], Field(description="Allowed action")]

# Validate email format
email: Annotated[str, Field(pattern=r'^[\w\.-]+@[\w\.-]+\.\w+$')]
```

### Advanced Security Patterns

PATTERN: Confirmation Flags PURPOSE: Prevent accidental destructive operations IMPLEMENTATION:

```python
@mcp.tool()
def delete_data(
    id: str,
    confirm: Annotated[bool, Field(description="Set true to confirm deletion")] = False
) -> dict:
    """Delete data (requires confirmation)."""
    if not confirm:
        raise ToolError("Set confirm=true to proceed with deletion")
    perform_deletion(id)
    return {"status": "deleted", "id": id}
```

PATTERN: Rate Limiting PURPOSE: Prevent abuse and resource exhaustion IMPLEMENTATION:

```python
from functools import wraps
from time import time

def rate_limit(calls_per_minute: int):
    """Decorator to rate limit tool calls."""
    calls = []
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            now = time()
            # Remove calls older than 1 minute
            calls[:] = [t for t in calls if now - t < 60]
            if len(calls) >= calls_per_minute:
                raise ToolError(f"Rate limit exceeded: {calls_per_minute}/min")
            calls.append(now)
            return func(*args, **kwargs)
        return wrapper
    return decorator

@mcp.tool()
@rate_limit(calls_per_minute=10)
def expensive_operation(query: str) -> dict:
    """Rate-limited expensive operation."""
    return process_query(query)
```

PATTERN: Suspicious Parameter Detection PURPOSE: Block SQL injection, path traversal, command injection IMPLEMENTATION:

```python
def validate_safe_path(path: str) -> str:
    """Validate path does not contain traversal sequences."""
    if ".." in path or path.startswith("/"):
        raise ToolError("Invalid path: contains traversal or absolute path")
    return path

def validate_safe_query(query: str) -> str:
    """Validate SQL query does not contain suspicious patterns."""
    forbidden = ["DROP", "DELETE", "TRUNCATE", "--", ";"]
    if any(pattern in query.upper() for pattern in forbidden):
        raise ToolError("Query contains forbidden pattern")
    return query

@mcp.tool()
def read_file(path: Annotated[str, Field(max_length=500)]) -> dict:
    """Read file with path validation."""
    safe_path = validate_safe_path(path)
    # Additional allowed directory check
    if not safe_path.startswith(ALLOWED_PREFIX):
        raise ToolError(f"Path must start with {ALLOWED_PREFIX}")
    return {"content": open(safe_path).read()}
```

CONSTRAINTS:

- The model must validate all user-provided paths for traversal sequences
- The model must rate limit expensive operations (API calls, database queries)
- The model must require confirmation flags for destructive operations
- The model must whitelist allowed values using Literal types when possible
- The model must reject suspiciously long inputs using Field(max_length)

## Observability & Testing

### Testing Utilities

TOOL: FastMCP CLI COMMAND: `fastmcp cli` PURPOSE: Interactive testing without AI client USAGE:

1. Start server in one terminal: `python server.py`
2. Run CLI in another: `fastmcp cli`
3. Send test tool calls interactively
4. Verify responses and error handling

TOOL: mcp-cli Project REPOSITORY: <https://github.com/wong2/mcp-cli> PURPOSE: Automated integration testing USAGE:

```bash
npm install -g mcp-cli
mcp-cli call server.py tool_name '{"param": "value"}'
```

TESTING_CHECKLIST:

- The model must verify tool calls with invalid parameters raise ToolError
- The model must verify tool calls with missing required parameters are rejected
- The model must verify async tools complete without blocking
- The model must verify error masking works when mask_error_details=True
- The model must verify Context.info() messages appear in client logs
- The model must verify streaming tools yield incremental results
- The model must verify excluded parameters do not appear in tool schema

### Logging and Notifications

NOTIFICATION: tools/list_changed TRIGGER: Tools added or removed at runtime PURPOSE: Notify clients of dynamic tool availability SUPPORTED: FastMCP v2.5+

```python
# Dynamic tool registration based on runtime conditions
@mcp.tool()
def base_tool():
    """Always available."""
    pass

# Conditionally add tools
if user_has_permission("admin"):
    @mcp.tool()
    def admin_tool():
        """Only available to admin users."""
        pass
    # Client receives tools/list_changed notification

# Remove tools dynamically
mcp.remove_tool("deprecated_tool")
# Client receives tools/list_changed notification
```

CONSTRAINTS:

- The model must emit tools/list_changed only when tool set changes
- The model must not rely on notifications for initial tool discovery (use tools/list)
- The model must handle clients that ignore notifications gracefully

USE_CASES:

- Role-based tool availability
- Feature flags enabling/disabling tools
- A/B testing tool variants
- Progressive tool loading for large servers

### Debugging Practices

1. **Verbose mode** - Watch JSON messages being passed
2. **Capture MCP traffic** - Cursor dev tools or Claude inspector show raw messages
3. **Structured logging** - Use Context.info() for AI-visible progress

**Testing checklist**:

- Can it handle multiple simultaneous tool calls?
- Does it recover if a tool crashes?
- Are tool descriptions clear enough that the model chooses the right tool?
- Does it work with different AI models (Claude, GPT-4)?

### Early Feedback Loop

Community feedback often points out if an AI is confused or does something odd. Address by:

- Tweaking descriptions
- Adding constraints
- Improving error messages

## Performance Tuning

### Caching Strategies

Use in-memory caches for repeated queries:

```python
from functools import lru_cache

@lru_cache(maxsize=128)
def get_file_index(directory: str):
    """Cache file listings to avoid re-reading on each call."""
    return list_all_files(directory)

@mcp.tool()
def search_files(pattern: str, directory: str = ".") -> dict:
    """Search for files matching pattern."""
    file_list = get_file_index(directory)  # Uses cache
    matches = [f for f in file_list if pattern in f]
    return {"matches": matches}
```

### Concurrency

FastMCP is thread-safe and async-friendly:

```python
import asyncio

@mcp.tool()
async def parallel_requests(urls: list[str]) -> dict:
    """Fetch multiple URLs concurrently."""
    async with httpx.AsyncClient() as client:
        tasks = [client.get(url) for url in urls]
        responses = await asyncio.gather(*tasks)
        return {"results": [r.json() for r in responses]}
```

### Streaming Large Outputs

For big text or files, yield in chunks:

```python
@mcp.tool()
async def stream_large_file(path: str):
    """Stream file content in chunks."""
    with open(path) as f:
        while chunk := f.read(4096):
            # FastMCP supports streaming in HTTP mode
            yield {"chunk": chunk}
```

**Updated for FastMCP 3.x**: Streaming supported in HTTP mode with automatic threadpool for sync tools.

### Performance Best Practices

1. **Cache file indexes** in memory using lru_cache
2. **Use asyncio.gather** for independent operations
3. **Stream partial results** to avoid timeout limits
4. **Use background tasks** for operations >30 seconds
5. **Configure OpenTelemetry** for performance monitoring
6. **Leverage automatic threadpool** for CPU-bound sync operations
7. **Use component versioning** instead of breaking changes

## Ecosystem and Compatibility

### Multi-Platform Support

PLATFORM: Anthropic Claude Desktop STATUS: Native MCP support INTEGRATION: Built-in, no adapter needed CONFIG: `~/Library/Application Support/Claude/claude_desktop_config.json`

PLATFORM: OpenAI GPT-4/ChatGPT STATUS: Community adapters INTEGRATION: Via bridge libraries EXAMPLES:

- Glama Chat with MCP support
- LobeChat MCP integration
- Custom OpenAI → MCP bridges

IMPLEMENTATION_PATTERN:

```python
# OpenAI function calling → MCP tool adapter
def convert_mcp_to_openai_function(tool_schema: dict) -> dict:
    """Convert MCP tool schema to OpenAI function format."""
    return {
        "name": tool_schema["name"],
        "description": tool_schema["description"],
        "parameters": tool_schema["inputSchema"]
    }

# Use FastMCP server with OpenAI client
openai_functions = [
    convert_mcp_to_openai_function(tool)
    for tool in mcp.list_tools()
]
```

PLATFORM: VS Code STATUS: Native MCP support (2025+) INTEGRATION: Extension API for auto-registering MCP servers CONFIG: `.vscode/mcp.json`

```json
{
  "mcpServers": {
    "my-server": {
      "command": "python",
      "args": ["server.py"],
      "env": {
        "API_KEY": "${env:API_KEY}"
      }
    }
  }
}
```

CONSTRAINTS:

- The model must provide OpenAI function schemas when targeting OpenAI clients
- The model must use VS Code MCP extension API for VS Code integration
- The model must configure transport appropriately for each platform

### Cross-Platform Clients

Community projects enabling MCP across different AI platforms:

- **Swarm** (OpenAI/Swarm) - Bridging GPT to MCP
- **CrewAI** - Multi-agent framework with MCP support
- Various chat interfaces with MCP integration

### MCP as De Facto Standard

Community consensus: MCP is becoming the standard for AI tool integration.

**Evidence**:

- Hugging Face blog coverage
- GitHub star growth charts
- Major platform adoptions
- Active community development

### Staying Updated

**Resources**:

- **Awesome MCP** lists for clients and servers
- **Official MCP website** (modelcontextprotocol.io)
- **Community Discord** channels
- **GitHub trending** MCP repositories

**Compatibility notes**: Lists often mention if new versions of Claude or Cursor changed connection formats or configuration.

## Agent Orchestration Patterns

### Multi-Model Orchestration

Use "planner" LLM that breaks requests into steps, and specialist MCP servers to execute.

**Architecture**:

```text
User Request
    ↓
Planner LLM (GPT-4 / Claude)
    ↓
[Step 1] → MCP Server A
[Step 2] → MCP Server B
[Step 3] → MCP Server C
    ↓
Planner composes results
    ↓
Response to User
```

### Design for Composability

**Server design principles**:

1. **Focused scope** - Do one thing well
2. **Consistent naming** - Prefix related tools
3. **Clear metadata** - Help orchestrators discover capabilities
4. **Tag tools** - For categorization and routing

**Example tool tagging**:

```python
@mcp.tool(
    tags=["calendar", "scheduling"],
    annotations={"category": "productivity"}
)
def schedule_meeting(...):
    ...
```

### Performance Metadata

PATTERN: Performance Annotations for Orchestration PURPOSE: Enable intelligent tool routing by agent orchestrators STATUS: Emerging practice (mid-2025+)

```python
@mcp.tool(
    annotations={
        "estimatedLatencyMs": 150,
        "costPerCall": 0.001,
        "preferredModel": "fast",
        "cacheable": True,
        "idempotent": True
    }
)
def search_database(query: str) -> dict:
    """Fast, cacheable database search."""
    return execute_query(query)

@mcp.tool(
    annotations={
        "estimatedLatencyMs": 5000,
        "costPerCall": 0.05,
        "preferredModel": "accurate",
        "cacheable": False,
        "idempotent": False
    }
)
def ml_inference(data: dict) -> dict:
    """Expensive ML model inference."""
    return run_model(data)
```

ANNOTATION_SCHEMA:

- estimatedLatencyMs: Average response time (int)
- costPerCall: Dollar cost per invocation (float)
- preferredModel: Routing hint ("fast", "accurate", "balanced")
- cacheable: Safe to cache results (bool)
- idempotent: Multiple calls with same params produce same result (bool)

ORCHESTRATION_USAGE:

- Route cheap queries to fast models
- Route expensive operations to accurate models
- Cache idempotent + cacheable tool results
- Parallelize independent, low-latency tools
- Serialize expensive, high-latency tools

CONSTRAINTS:

- The model must estimate latency conservatively (prefer overestimate)
- The model must mark non-idempotent tools explicitly
- The model must not cache non-cacheable tool results

### Projects Demonstrating Orchestration

- **WayStation.ai** - Multi-server orchestration
- **MetaMCP** - Agent that discovers and installs servers
- **Ultimate MCP Server** - Model delegation patterns

## FastMCP 3.x Community Patterns (2025-2026)

### Provider-Based Server Organization

Community pattern: Use FileSystemProvider for modular tool organization:

```python
from fastmcp import FastMCP
from fastmcp.server.providers import FileSystemProvider

# Organize tools in separate files
mcp = FastMCP("organized-server", providers=[
    FileSystemProvider("tools/", reload=True)
])

# tools/search.py
from fastmcp.tools import tool

@tool
def search(query: str) -> dict:
    """Search operation."""
    ...

# tools/fetch.py
from fastmcp.tools import tool

@tool
def fetch(url: str) -> dict:
    """Fetch operation."""
    ...
```

Benefits:
- No import coupling between tools
- Hot reload during development
- Clear separation of concerns
- Easy to test individual tools

### Component Versioning Strategy

Community best practice: Use semantic versioning for API evolution:

```python
# Maintain v1 for backward compatibility
@mcp.tool(version="1.0")
def process_data(data: str) -> dict:
    """Legacy v1 implementation."""
    return {"result": simple_process(data)}

# Add v2 with enhanced features
@mcp.tool(version="2.0")
def process_data(data: str, options: dict = None) -> dict:
    """Enhanced v2 with options support."""
    return {
        "result": enhanced_process(data, options),
        "version": "2.0"
    }
```

Pattern adopted by:
- Large-scale API servers
- Enterprise MCP deployments
- Multi-tenant platforms

### Session-Scoped Feature Gating

Community pattern: Use session visibility for premium features:

```python
@mcp.tool(tags={"premium"})
def premium_analysis(data: str) -> dict:
    """Premium-only analysis tool."""
    ...

@mcp.tool()
async def upgrade_session(ctx: Context, license_key: str) -> str:
    """Unlock premium features with license."""
    if validate_license(license_key):
        await ctx.enable_components(tags={"premium"})
        await ctx.set_state("license", license_key)
        return "Premium features unlocked"
    raise ToolError("Invalid license key")

# Globally disable premium - users unlock per-session
mcp.disable(tags={"premium"})
```

Use cases:
- SaaS MCP servers with tiered access
- Trial/demo modes
- Role-based feature access

### Background Task Patterns

Community best practice: Use background tasks for long operations:

```python
from fastmcp.server.tasks import TaskConfig

@mcp.tool(task=TaskConfig(mode="required"))
async def process_dataset(dataset_id: str) -> dict:
    """Process large dataset (background only)."""
    # Long-running computation
    results = await analyze_dataset(dataset_id)
    return {
        "dataset_id": dataset_id,
        "status": "complete",
        "summary": results
    }

@mcp.tool()
async def check_job_status(ctx: Context, job_id: str) -> dict:
    """Check background job status."""
    status = await ctx.get_state(f"job_{job_id}")
    return {"job_id": job_id, "status": status}
```

Adopted by:
- Data processing servers
- Report generation tools
- Batch operation platforms

### OpenTelemetry Integration Pattern

Community standard: Configure OTEL for production observability:

```python
from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.trace.export import BatchSpanProcessor
from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter

# Configure once at startup
provider = TracerProvider()
provider.add_span_processor(
    BatchSpanProcessor(OTLPSpanExporter(endpoint="http://otel-collector:4317"))
)
trace.set_tracer_provider(provider)

# FastMCP automatically traces all operations
mcp = FastMCP("traced-server")

@mcp.tool()
async def traced_operation(data: str) -> dict:
    # Automatic span creation with context propagation
    result = await process(data)
    return {"result": result}
```

Benefits:
- Distributed tracing across services
- Performance monitoring
- Error tracking
- Request correlation

### Transform-Based Tool Curation

Community pattern: Use ToolTransform to optimize auto-generated tools:

```python
from fastmcp.server.providers.openapi import OpenAPIProvider
from fastmcp.server.transforms import ToolTransform
from fastmcp.tools.tool_transform import ToolTransformConfig

# Start with OpenAPI spec
api_provider = OpenAPIProvider(openapi_spec=spec, client=httpx_client)

# Curate auto-generated tools for agents
api_provider.add_transform(ToolTransform({
    "get_api_v1_users_user_id_": ToolTransformConfig(
        name="get_user",
        description="Get user by ID - returns profile, settings, and activity",
        tags={"users", "read"}
    ),
    "post_api_v1_users_": ToolTransformConfig(
        name="create_user",
        description="Create new user account",
        tags={"users", "write"}
    )
}))

mcp = FastMCP("curated-api", providers=[api_provider])
```

Solves: Auto-generated tools often have poor names/descriptions for LLM consumption.

## Production Deployment Practices

### Standard Production Setup

1. **Containerization** - Docker for consistency
2. **Process supervision** - Uvicorn, Gunicorn, supervisord
3. **Environment variables** - Secrets management
4. **Logging** - Appropriate levels (INFO for production)
5. **Health checks** - `/health` endpoints
6. **Monitoring** - Track usage, errors, latency

### FastMCP Cloud

One-command deployment:

```bash
fastmcp deploy server.py
```

Provides:

- Automatic HTTPS
- Built-in authentication
- Auto-scaling
- Managed infrastructure

### Custom Production Deployment

Embed in Starlette/FastAPI app for custom routes:

```python
from starlette.applications import Starlette
from starlette.routing import Route

async def health(request):
    return JSONResponse({"status": "healthy"})

async def mcp_endpoint(request):
    payload = await request.json()
    response = await mcp.handle_http(payload)
    return JSONResponse(response)

app = Starlette(routes=[
    Route("/health", health),
    Route("/mcp", mcp_endpoint, methods=["POST"])
])
```

Deploy with:

```bash
uvicorn app:app --host 0.0.0.0 --port 8080
```

## Community Resources

### Discovery and Learning

- **Awesome MCP Servers** (GitHub)
- **Hugging Face Turing Post** articles
- **MCP Discord** communities
- **Reddit** r/ClaudeAI, r/LocalLLaMA
- **GitHub Discussions** on MCP repos

### Contributing Back

- Many projects accept PRs
- Share your servers on Awesome MCP lists
- Write blog posts about patterns
- Help others in Discord/Reddit

### Staying Current

MCP is rapidly evolving. Stay updated:

1. Watch official MCP repository
2. Follow Anthropic engineering blog
3. Join community Discord servers
4. Check Awesome MCP lists monthly
5. Read Hugging Face AI news

## Summary

MANDATORY_PRACTICES:

1. INTEGRATION:

   - The model must package servers as .mcpb for Claude Desktop
   - The model must provide manifest.json with user_config fields
   - The model must support environment variable configuration

2. SECURITY:

   - The model must validate file paths against allowed directories
   - The model must use confirmation flags for destructive operations
   - The model must set destructiveHint annotation for state-changing tools
   - The model must rate limit expensive operations
   - The model must use authorization system for granular access control

3. PERFORMANCE:

   - The model must cache repeated queries using lru_cache or similar
   - The model must use async for I/O-bound operations
   - The model must stream large responses in HTTP mode
   - The model must add performance annotations for orchestration
   - The model must use background tasks for long-running operations

4. OBSERVABILITY:

   - The model must test servers using fastmcp cli or mcp-cli
   - The model must emit tools/list_changed for dynamic tool updates
   - The model must use Context.info() for progress logging
   - The model must verify tool schemas before deployment
   - The model must configure OpenTelemetry for production tracing

5. COMPATIBILITY:

   - The model must provide OpenAI function schemas when targeting GPT-4
   - The model must configure transport per platform requirements
   - The model must follow platform-specific config formats

6. FASTMCP 3.x PATTERNS:

   - The model must use FileSystemProvider for modular organization
   - The model must version components when making breaking changes
   - The model must use session-scoped visibility for feature gating
   - The model must use TaskConfig for background operations
   - The model must use ToolTransform to curate auto-generated tools
   - The model must configure distributed state storage for multi-instance deployments

## Sources

- [FastMCP GitHub Repository (PrefectHQ/fastmcp)](https://github.com/prefecthq/fastmcp) (accessed 2026-02-22)
- [FastMCP Official Documentation](https://gofastmcp.com/) (accessed 2026-02-22)
- [FastMCP on PyPI](https://pypi.org/project/fastmcp/) (accessed 2026-02-22)
- [MCP Protocol Specification](https://spec.modelcontextprotocol.io/) (accessed 2026-02-22)
- Claude Desktop Extensions (.mcpb format) — Anthropic official announcement and FastMCP packaging docs via <https://gofastmcp.com/> (accessed 2026-02-22)
