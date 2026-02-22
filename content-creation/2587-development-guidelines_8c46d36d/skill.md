# FastMCP 3.x Development Guidelines

Complete guide to building MCP servers with FastMCP 3.x, covering decorators, provider/transform architecture, component versioning, authorization, session state, validation, async patterns, error handling, background tasks, and deployment.

**VERSION**: This guide covers FastMCP 3.0+. For FastMCP 2.x, see legacy documentation.

## Architecture Overview (FastMCP 3.x)

FastMCP 3.x introduces a composable architecture built on three core concepts:

**Components** - The atoms of MCP: tools, resources, prompts. These are what clients interact with. Components have names, schemas, metadata, versions, and behavior.

**Providers** - Answer "where do components come from?" Anything that can list and retrieve components:
- Your decorated functions (LocalProvider)
- A directory of files (FileSystemProvider)
- A remote MCP server (ProxyProvider)
- Another FastMCP server (FastMCPProvider)
- An OpenAPI spec (OpenAPIProvider)
- Your custom sources (subclass Provider)

**Transforms** - Middleware for the component pipeline. Intercept and modify components as they flow from providers to clients:
- Namespace: Add prefixes to avoid collisions
- ToolTransform: Reshape tool names, descriptions, schemas
- VersionFilter: Expose only certain versions
- Visibility: Enable/disable components dynamically
- Custom: Write your own by subclassing Transform

This composability means features emerge from combining primitives rather than special-case code.

## Core Abstractions

RULE: The model must instantiate FastMCP object before registering tools
RULE: The model must use decorators to expose Python functions as MCP capabilities
RULE: The model must not implement low-level MCP protocol details (FastMCP handles this)

**FastMCP 3.x Basic Pattern:**

```python
from fastmcp import FastMCP

mcp = FastMCP("server-name")

@mcp.tool()
def tool_name(...):
    ...

@mcp.prompt()
def prompt_name(...):
    ...

@mcp.resource("uri://pattern")
def resource_name(...):
    ...

if __name__ == "__main__":
    mcp.run()  # STDIO by default
```

**Important v3 Change**: Decorated functions now return callable functions (not objects). This enables direct testing:

```python
@mcp.tool()
def greet(name: str) -> str:
    return f"Hello, {name}!"

# Now you can call it directly for testing
result = greet("World")  # Works in v3!
```

For v2 compatibility, set environment variable: `FASTMCP_DECORATOR_MODE=object`

## Define Tools with Decorators

Annotate functions with `@mcp.tool` to turn them into callable tools for the AI client.

### Tool Name and Description

AUTO_GENERATION:

- Tool name ← function name (override: `@mcp.tool(name="custom")`)
- Description ← first line of docstring (override: `@mcp.tool(description="...")`)
- JSON Schema ← type hints

### Type Annotations Required

RULE: The model must annotate all parameters and return values RULE: The model must use type hints to enable automatic validation RULE: The model must return dict or Pydantic models for structured output

DECISION_TREE:

```text
1. IF parameter needs validation THEN
   - Use Annotated[type, Field(...)]
   - Add constraints (ge, le, pattern, etc)
2. ELSE IF parameter is optional THEN
   - Use type | None with default = None
3. ELSE
   - Use plain type hint

4. IF return value is structured THEN
   - Return dict or Pydantic model
5. ELSE IF return value is simple THEN
   - Return str, int, bool (discouraged for AI consumption)
```

CONSTRAINTS:

- The model must not omit type hints (causes schema generation failure)
- The model must not return unstructured strings (AI cannot parse reliably)
- The model must use structured dict for all non-trivial outputs

### Example Tool Definition

```python
from fastmcp import FastMCP
from pydantic import Field
from typing import Annotated

mcp = FastMCP("product-search")

@mcp.tool()
def search_products(
    query: str,
    category: Annotated[str | None, Field(description="Filter by product category")] = None,
    max_results: Annotated[int, Field(ge=1, le=100, description="Maximum results to return")] = 10
) -> dict:
    """Search the product catalog with optional category filtering.

    Returns product information including name, price, availability, and ratings.
    Useful for finding products based on keywords or browsing by category.
    """
    # Implementation
    results = database.search(query, category, limit=max_results)
    return {"results": results, "count": len(results)}
```

## Parameter Metadata & Validation

FastMCP integrates with Pydantic for robust parameter handling.

### Field Constraints

Attach `typing.Annotated` with `Field(...)` to impose constraints and descriptions:

```python
from pydantic import Field
from typing import Annotated, Literal

# Numeric ranges
width: Annotated[int, Field(ge=1, le=2000, description="Image width in pixels")]
price: Annotated[float, Field(ge=0.01, description="Price must be positive")]

# String patterns
email: Annotated[str, Field(pattern=r'^[\w\.-]+@[\w\.-]+\.\w+$')]
code: Annotated[str, Field(min_length=3, max_length=10)]

# Literal types for enums
action: Annotated[Literal["read", "write", "delete"], Field(description="Allowed actions")]
format: Annotated[Literal["json", "xml", "yaml"], Field(description="Output format")]

# Lists with item constraints
tags: Annotated[list[str], Field(min_length=1, max_length=10, description="1-10 tags")]

# Exclude from schema (runtime injection)
api_key: Annotated[str, Field(exclude=True)] = None

# Optional parameters
category: Annotated[str | None, Field(description="Optional category filter")] = None
```

### Parameter Validation Benefits

- **Safety**: AI cannot call tools with invalid or dangerous inputs
- **Documentation**: Constraints appear in tool schema for the AI to see
- **Automatic checking**: FastMCP validates before calling your function

### Excluding Parameters from Schema

Use `Field(exclude=True)` for parameters that should not be visible to the LLM (e.g., secrets, runtime-injected values):

```python
@mcp.tool()
def api_call(
    query: str,
    api_key: Annotated[str, Field(exclude=True)] = None  # Injected at runtime
) -> dict:
    # api_key never appears in the tool schema the LLM sees
    headers = {"Authorization": f"Bearer {api_key}"}
    # ...
```

**Note**: Only optional arguments can be excluded. Required arguments must be visible to the AI.

## Asynchronous Tools for Performance

FastMCP supports both sync and async tool functions.

### When to Use Async

Use `async def` for:

- I/O-bound operations (database queries, API calls)
- High-latency operations (network requests)
- Any operation that would block the event loop

### Example Async Tool

```python
import httpx

@mcp.tool()
async def fetch_weather(city: str) -> dict:
    """Fetch current weather for a city."""
    async with httpx.AsyncClient() as client:
        response = await client.get(f"https://api.weather.com/{city}")
        return response.json()
```

### Performance Benefits

The framework handles running async tools concurrently, boosting throughput. In a real server, mix:

- **Synchronous tools** - For quick in-memory tasks
- **Asynchronous tools** - For network or disk operations

This maximizes performance without blocking.

### Streaming Responses

FastMCP v2.10+ supports streaming large responses in HTTP mode:

```python
@mcp.tool()
async def stream_large_data(query: str):
    """Stream results incrementally to avoid timeout."""
    for chunk in process_query_streaming(query):
        yield {"chunk": chunk, "progress": chunk.index}
```

CONSTRAINTS:

- The model must use HTTP transport for streaming (STDIO does not support streaming)
- The model must yield dict objects (not strings) for structured streaming
- The model must include progress indicators when streaming long operations

Use streaming for:

- Large file processing
- Database queries returning many rows
- Long-running computations with incremental results
- API responses exceeding context window limits

## Error Handling & Robustness

Proper error handling is baked into FastMCP.

### Automatic Error Handling

If a tool function raises any Python exception, FastMCP:

1. Intercepts the exception
2. Returns an MCP error response to the AI client
3. Allows the LLM to react (apologize, try different approach)

### Error Masking for Production

By default, exception details are included in error responses. For production, use `mask_error_details=True` to replace error traces with generic messages:

```python
mcp = FastMCP("my-server", mask_error_details=True)
```

This prevents leaking internal implementation details.

### ToolError for Business Logic Errors

Use `FastMCP.exceptions.ToolError` for expected errors that the AI should know about:

```python
from fastmcp.exceptions import ToolError

@mcp.tool()
def get_user(user_id: str) -> dict:
    user = database.find_user(user_id)
    if not user:
        raise ToolError(f"User {user_id} not found")
    return user
```

**Key behavior**: Messages from `ToolError` are **always sent to the client**, even if masking is enabled.

### Error Handling Pattern

Common pattern for production servers:

```python
@mcp.tool()
def process_order(order_id: str) -> dict:
    # Business logic errors - AI should see these
    if not order_id.startswith("ORD-"):
        raise ToolError("Invalid order ID format. Must start with 'ORD-'")

    order = database.get_order(order_id)
    if not order:
        raise ToolError(f"Order {order_id} not found")

    # Let unexpected exceptions bubble up as masked generic errors
    # (e.g., database connection failures)
    result = order.process()
    return {"status": "processed", "order": result}
```

This dual approach ensures:

- **User-friendly errors** for business logic issues
- **Security** by masking internal implementation errors

## Component Versioning (FastMCP 3.x)

RULE: The model must use semantic versioning (PEP 440) for component versions
RULE: The model must register multiple versions when making breaking changes to tools
RULE: The model must let FastMCP automatically expose the highest version

### Declaring Versions

Register multiple versions of the same component:

```python
@mcp.tool(version="1.0")
def add(x: int, y: int) -> int:
    """Add two numbers (v1)."""
    return x + y

@mcp.tool(version="2.0")
def add(x: int, y: int, z: int = 0) -> int:
    """Add two or three numbers (v2)."""
    return x + y + z

# Only v2.0 exposed via list_tools()
# Calling "add" invokes v2.0 implementation
# v1.0 still available for compatibility
```

### Version Selection

FastMCP automatically exposes the highest version. Clients can request specific versions:

```python
# FastMCP client supports direct version selection
from fastmcp import Client

async with Client(server) as client:
    # Call latest version (default)
    result = await client.call_tool("add", {"x": 1, "y": 2})

    # Call specific version
    result = await client.call_tool("add", {"x": 1, "y": 2}, version="1.0")
```

Generic MCP clients use `_meta` in arguments:

```json
{
  "x": 1,
  "y": 2,
  "_meta": {
    "fastmcp": {
      "version": "1.0"
    }
  }
}
```

### Version Metadata

All available versions exposed in component metadata:

```python
tools = await client.list_tools()
# Each tool includes:
# - meta["fastmcp"]["version"]: "2.0" (current)
# - meta["fastmcp"]["versions"]: ["2.0", "1.0"] (all available)
```

### When to Version

Use versioning when:
- Changing parameter names or types (breaking change)
- Removing required parameters (breaking change)
- Changing return value structure (breaking change)
- Adding optional parameters (non-breaking, version bump optional)

## Session-Scoped State (FastMCP 3.x)

RULE: The model must use async methods for session state (v3 requirement)
RULE: The model must understand state persists across tool calls in a session
RULE: The model must configure storage backend for distributed deployments

### Basic Session State

State persists across tool calls within a session:

```python
from fastmcp import Context

@mcp.tool()
async def increment_counter(ctx: Context) -> int:
    """Increment session counter."""
    count = await ctx.get_state("counter") or 0  # async in v3
    await ctx.set_state("counter", count + 1)     # async in v3
    return count + 1

@mcp.tool()
async def reset_counter(ctx: Context) -> str:
    """Reset session counter."""
    await ctx.delete_state("counter")
    return "Counter reset"
```

**Key Changes from v2**:
- Methods are now async: `await ctx.get_state()`, `await ctx.set_state()`, `await ctx.delete_state()`
- State automatically keyed by session ID (isolation between clients)
- State expires after 1 day (TTL) to prevent unbounded growth

### Distributed Storage

For production deployments, configure a distributed backend:

```python
from key_value.aio.stores.redis import RedisStore

mcp = FastMCP(
    "server",
    session_state_store=RedisStore(host="localhost", port=6379)
)
```

FastMCP uses [pykeyvalue](https://github.com/strawgate/py-key-value) for pluggable storage backends.

### Stateless HTTP Sessions

For stateless HTTP deployments, FastMCP respects the `mcp-session-id` header that most clients send. If you've configured a storage backend, virtual sessions are created automatically.

## Authorization (FastMCP 3.x)

RULE: The model must use component-level auth for granular access control
RULE: The model must use server-wide auth for global requirements
RULE: The model must understand STDIO transport bypasses all auth

### Component-Level Authorization

Protect individual components with auth decorators:

```python
from fastmcp.server.auth import require_scopes, restrict_tag

@mcp.tool(auth=require_scopes("write"))
def protected_tool():
    """Requires 'write' scope."""
    ...

@mcp.resource("data://secret", auth=require_scopes("read"))
def secret_data():
    """Requires 'read' scope."""
    ...

@mcp.prompt(auth=require_scopes("admin"))
def admin_prompt():
    """Requires 'admin' scope."""
    ...
```

Built-in auth checks:
- `require_scopes(*scopes)`: Requires specific OAuth scopes
- `restrict_tag(tag, scopes)`: Requires scopes only for tagged components

### Server-Wide Authorization

Apply authorization globally via AuthMiddleware:

```python
from fastmcp.server.middleware import AuthMiddleware
from fastmcp.server.auth import require_scopes, restrict_tag

# Require specific scopes for all components
mcp = FastMCP(middleware=[AuthMiddleware(auth=require_scopes("read"))])

# Tag-based restrictions
mcp = FastMCP(middleware=[
    AuthMiddleware(auth=restrict_tag("admin", scopes=["admin"]))
])
```

### Custom Auth Checks

Custom checks receive `AuthContext` with token and component:

```python
from fastmcp.server.auth import AuthContext

def custom_check(ctx: AuthContext) -> bool:
    """Custom authorization logic."""
    return ctx.token is not None and "admin" in ctx.token.scopes

@mcp.tool(auth=custom_check)
def custom_protected():
    ...
```

**Note**: STDIO transport bypasses all auth checks (no OAuth concept in local subprocess execution).

## Visibility System (FastMCP 3.x)

RULE: The model must use visibility system to dynamically enable/disable components
RULE: The model must understand server-level vs session-level visibility
RULE: The model must send notifications when visibility changes

### Server-Level Visibility

Control which components are exposed globally:

```python
mcp = FastMCP("Server")

# Disable by name
mcp.disable(names={"dangerous_tool"}, components=["tool"])

# Disable by tag
mcp.disable(tags={"admin"})

# Allowlist mode - only show components with these tags
mcp.enable(tags={"public"}, only=True)

# Enable overrides earlier disable (later transform wins)
mcp.disable(tags={"internal"})
mcp.enable(names={"safe_tool"})  # safe_tool visible despite internal tag
```

### Session-Level Visibility

Control visibility per-session for feature gating:

```python
@mcp.tool(tags={"premium"})
def premium_analysis(data: str) -> str:
    """Premium feature - disabled by default."""
    return f"Premium analysis of: {data}"

@mcp.tool()
async def unlock_premium(ctx: Context) -> str:
    """Unlock premium features for this session."""
    await ctx.enable_components(tags={"premium"})
    return "Premium features unlocked"

@mcp.tool()
async def reset_features(ctx: Context) -> str:
    """Reset to default feature set."""
    await ctx.reset_visibility()
    return "Features reset to defaults"

# Globally disabled - individual sessions unlock
mcp.disable(tags={"premium"})
```

Session visibility methods:
- `await ctx.enable_components(...)`: Enable for this session
- `await ctx.disable_components(...)`: Disable for this session
- `await ctx.reset_visibility()`: Clear session rules, return to global defaults

FastMCP automatically sends `ToolListChangedNotification` (and resource/prompt equivalents) when visibility changes.

### Blocklist vs Allowlist

- **Blocklist mode** (default): All components visible except explicitly disabled
- **Allowlist mode** (`only=True`): Only explicitly enabled components visible

## Provider/Transform Architecture (FastMCP 3.x)

RULE: The model must understand providers source components
RULE: The model must understand transforms modify components in the pipeline
RULE: The model must use composability to build complex server behaviors

### Using FileSystemProvider

Organize tools as self-contained files instead of imports:

```python
from fastmcp import FastMCP
from fastmcp.server.providers import FileSystemProvider

# Create provider pointing at directory
provider = FileSystemProvider("mcp/", reload=True)

# Attach to server
mcp = FastMCP("server", providers=[provider])
```

Tool files use standalone imports:

```python
# mcp/greet.py
from fastmcp.tools import tool

@tool
def greet(name: str) -> str:
    """Greet someone by name."""
    return f"Hello, {name}!"
```

Benefits:
- No server import coupling
- Hot reload with `reload=True`
- Modular tool organization

### Mounting Servers

Compose servers together with namespace isolation:

```python
from fastmcp import FastMCP

main = FastMCP("Main")
sub = FastMCP("Sub")

@sub.tool()
def greet(name: str) -> str:
    return f"Hello, {name}!"

# Mount with prefix - greet becomes "sub_greet"
main.mount(sub, prefix="sub")
```

Under the hood: FastMCPProvider + Namespace transform.

### Using Transforms

Apply transforms at provider or server level:

```python
from fastmcp.server.transforms import Namespace, ToolTransform, VersionFilter
from fastmcp.tools.tool_transform import ToolTransformConfig

# Provider-level transform (affects only this provider)
provider.add_transform(Namespace("api"))

# Server-level transform (affects all providers)
mcp.add_transform(ToolTransform({
    "verbose_name": ToolTransformConfig(
        name="short_name",
        description="Better description for agents",
        tags={"category"}
    )
}))

# Version filtering for API versioning
api_v1 = FastMCP("API v1", providers=[components])
api_v1.add_transform(VersionFilter(version_lt="2.0"))

api_v2 = FastMCP("API v2", providers=[components])
api_v2.add_transform(VersionFilter(version_gte="2.0"))
```

## Context and Annotations

### Context Parameter

For advanced use, tools can accept a `Context` parameter to access the runtime MCP context:

```python
from fastmcp import Context

@mcp.tool()
def long_operation(data: str, context: Context) -> dict:
    """Perform a long-running operation with progress updates."""
    context.info("Starting operation...")

    # Process data
    for i in range(10):
        context.info(f"Processing chunk {i+1}/10")
        process_chunk(data, i)

    context.info("Operation complete")
    return {"status": "success"}
```

Context object allows:

- **Logging** - Emit info/warning messages the client can display
- **Progress callbacks** - Update the AI on long-running operations
- **Reading resources** - Access other resources exposed by the server

### MCP Annotations

FastMCP supports MCP Annotations on tools - metadata not in the prompt but informing clients about tool behavior:

```python
@mcp.tool(
    annotations={
        "readOnlyHint": True,      # Tool only reads data, doesn't modify state
        "openWorldHint": False,    # Tool doesn't access external systems
        "destructiveHint": False,  # Tool doesn't delete/modify data
        "title": "Search Products" # Display name in UI
    }
)
def search_products(query: str) -> dict:
    # ...
```

### Annotation Usage

AI-centric IDEs like Cursor and Claude Desktop use these hints to:

- **Require user approval** for destructive actions
- **Label tools** nicely in their UI
- **Apply safety checks** based on annotations

**Always set annotations accurately** (flag tools that write to disk or call external APIs) so client applications can apply proper safety checks.

## Resources & Data Access

Besides tools (which perform actions), FastMCP lets you expose **Resources** - essentially read-only data endpoints.

### Basic Resource

```python
@mcp.resource("data://config")
def get_config() -> dict:
    """Provide server configuration data."""
    return {
        "api_version": "1.0",
        "max_requests_per_minute": 60,
        "supported_formats": ["json", "xml"]
    }
```

### Resource Characteristics

- **Fetched by clients** via `resources/read` call (not tool invocation)
- **Read-only** - For providing reference data, documents, or images to the AI
- **URI-based** - Use custom URI schemes like `data://`, `config://`, `file://`

### Parameterized Resource Templates

URI patterns with placeholders for dynamic content:

```python
@mcp.resource("user://{user_id}/profile")
def get_user_profile(user_id: str) -> dict:
    """Get user profile by ID."""
    user = database.get_user(user_id)
    return {
        "id": user_id,
        "name": user.name,
        "email": user.email,
        "preferences": user.preferences
    }
```

### When to Use Resources

Use resources for:

- **Large data** that doesn't fit in prompts
- **Frequently-needed reference data** (schemas, configs)
- **Dynamic content generation** (user profiles, reports)

This keeps interactions efficient and contextual.

## Transport & Deployment Choices

DECISION_TREE: Transport Selection

```text
1. IF single-user desktop integration (Claude Desktop, Cursor) THEN
   - Use STDIO transport
   - Client launches server as subprocess
   - Communication via stdin/stdout
   - GOTO setup_stdio

2. ELSE IF multi-user or remote access THEN
   - Use HTTP transport
   - Server runs as network service
   - Communication via HTTP endpoint
   - GOTO setup_http

3. ELSE IF production deployment THEN
   - Use HTTP with FastMCP Cloud OR custom deployment
   - GOTO production_deployment

setup_stdio:
  mcp.run()  # Defaults to STDIO

setup_http:
  mcp.run(transport="http", host="0.0.0.0", port=8000)

production_deployment:
  Option A: fastmcp deploy server.py  # FastMCP Cloud
  Option B: uvicorn app:app --host 0.0.0.0 --port 8080  # Custom
```

CONSTRAINTS:

- The model must use STDIO for Claude Desktop/Cursor integration
- The model must use HTTP for remote or multi-client scenarios
- The model must not use STDIO for network-accessible servers

### STDIO Transport (Default)

Calling `mcp.run()` with no arguments starts the server in **STDIO transport mode**:

```python
if __name__ == "__main__":
    mcp.run()  # Defaults to STDIO
```

#### STDIO Characteristics

- **Local integration** - Ideal for Claude Desktop or IDE plugins
- **Subprocess model** - AI client launches server as subprocess
- **Pipe communication** - Uses stdin/stdout for messages
- **Isolated sessions** - One server instance per client session
- **No networking** - Great for desktop apps or CLI tools

### HTTP Transport

For broader usage, run the server as a network service over HTTP:

```python
if __name__ == "__main__":
    mcp.run(transport="http", host="0.0.0.0", port=8000)
```

#### HTTP Characteristics

- **Network service** - Can handle multiple clients
- **Streamable protocol** - MCP's preferred network transport
- **Endpoint** - Accessible at `http://localhost:8000/mcp`
- **Bidirectional** - Full streaming and communication support
- **Remote/multi-user** - Suitable for cloud deployments

### Legacy Transports

There is also legacy SSE (Server-Sent Events) transport and WebSocket support via community add-ons, but HTTP has essentially replaced SSE in the latest MCP spec.

### Production Deployment

Treat an MCP server like any web service:

#### Containerization

- **Docker** - Containerize for consistent deployment
- **Process managers** - Use Uvicorn, Gunicorn, or supervisord
- **Cloud services** - Deploy to AWS, GCP, Azure, or serverless platforms

#### FastMCP Cloud

The FastAPI team launched **FastMCP Cloud** - a hosting platform for one-command deployment:

```bash
fastmcp deploy server.py
```

Features:

- Automatic HTTPS
- Built-in authentication
- Auto-scaling
- "Remote MCP that just works"

#### Manual HTTP Server

FastMCP's built-in CLI can launch a production-ready HTTP server:

```bash
fastmcp run server.py --transport http --port 8080 --log-level INFO
```

#### Custom Integration

For advanced use, embed FastMCP in a Starlette app:

```python
from starlette.applications import Starlette
from starlette.routing import Route
from starlette.responses import JSONResponse

async def health_check(request):
    return JSONResponse({"status": "healthy"})

async def mcp_handler(request):
    payload = await request.json()
    response = await mcp.handle_http(payload)
    return JSONResponse(response)

app = Starlette(routes=[
    Route("/health", health_check),
    Route("/mcp", mcp_handler, methods=["POST"])
])
```

### Production Best Practices

1. **Environment variables** - Secure secrets via env vars (FastMCP passes them through)
2. **Logging** - Enable appropriate log levels (`--log-level INFO`)
3. **Process supervision** - Use container orchestration or supervisord
4. **Health checks** - Add `/health` endpoints for monitoring
5. **Secret management** - Never hardcode API keys or credentials

### Example Production Configuration

```python
import os
from fastmcp import FastMCP

mcp = FastMCP(
    "production-server",
    mask_error_details=True  # Hide internal errors
)

# Read configuration from environment
API_KEY = os.getenv("API_KEY")
ALLOWED_DIRS = os.getenv("ALLOWED_DIRECTORIES", "").split(",")

@mcp.tool()
def secure_operation(
    query: str,
    api_key: Annotated[str, Field(exclude=True)] = API_KEY
) -> dict:
    # Use injected API key
    # ...
```

## Background Tasks (FastMCP 3.x)

RULE: The model must use TaskConfig for long-running operations
RULE: The model must install fastmcp[tasks] extra for Docket integration
RULE: The model must understand task execution modes (forbidden/optional/required)

FastMCP 3.x implements SEP-1686 for persistent background tasks with Docket integration.

### Installation

```bash
pip install "fastmcp[tasks]>=3.0.0"
# or
uv add "fastmcp[tasks]>=3.0.0"
```

### Task Configuration

```python
from fastmcp.server.tasks import TaskConfig

# Required mode - must execute as background task
@mcp.tool(task=TaskConfig(mode="required"))
async def long_running_task(data: str) -> dict:
    """Process large dataset (background only)."""
    # Long-running computation
    result = await process_large_dataset(data)
    return {"status": "complete", "result": result}

# Optional mode - supports both sync and task execution
@mcp.tool(task=TaskConfig(mode="optional"))
async def flexible_task(query: str) -> dict:
    """Flexible execution mode."""
    results = await perform_query(query)
    return {"results": results}

# Shorthand for optional mode
@mcp.tool(task=True)
async def simple_task(data: str) -> dict:
    """Simple background-capable task."""
    return await process(data)
```

### Task Modes

- **`"forbidden"`** (default): Does not support task execution
- **`"optional"`**: Supports both synchronous and task execution
- **`"required"`**: Must be executed as background task

### Task Lifecycle

Background tasks run in persistent queue with:
- Progress tracking
- Error handling
- Cancellation support
- Result storage

Backend options:
- SQLite (development)
- PostgreSQL (production)
- Horizontal scaling support

## Production Features (FastMCP 3.x)

### OpenTelemetry Tracing

RULE: The model must configure OpenTelemetry for production observability
RULE: The model must understand FastMCP auto-instruments all operations

Native OpenTelemetry support traces all operations:

```python
from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.trace.export import BatchSpanProcessor
from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter

# Configure OpenTelemetry
provider = TracerProvider()
provider.add_span_processor(BatchSpanProcessor(OTLPSpanExporter()))
trace.set_tracer_provider(provider)

# Use FastMCP normally - automatic tracing
mcp = FastMCP("traced-server")

@mcp.tool()
async def traced_operation(data: str) -> dict:
    # Automatically traced with context propagation
    return await process(data)
```

Automatic span attributes:
- Component key (tool/resource/prompt name)
- Provider type
- Session ID
- Auth context
- W3C trace context propagation

### Tool Timeouts

Limit foreground execution time:

```python
import httpx
from urllib.parse import urlparse
from fastmcp.exceptions import ToolError

# Define allowed destinations to prevent SSRF attacks
ALLOWED_SCHEMES = {"https"}
ALLOWED_HOSTS = {
    "api.example.com",
    "services.internal.example",
}

@mcp.tool(timeout=30.0)
async def fetch_data(url: str) -> dict:
    """Fetch JSON from an allowed URL with a 30-second timeout.

    The URL is validated against an allowlist to prevent SSRF attacks.
    """
    parsed = urlparse(url)

    if parsed.scheme not in ALLOWED_SCHEMES or parsed.hostname not in ALLOWED_HOSTS:
        raise ToolError("URL not allowed; must use https and an approved host")

    async with httpx.AsyncClient() as client:
        response = await client.get(url)
        response.raise_for_status()
        return response.json()
```

When timeout is exceeded, clients receive MCP error code `-32000`.

**Security Note**: Always validate URLs against an allowlist to prevent SSRF (Server-Side Request Forgery) attacks where malicious users could access internal services, cloud metadata endpoints, or protected resources.

**Note**: Timeouts don't apply to background tasks (those use Docket's lifecycle).

### Pagination

For servers with many components:

```python
server = FastMCP("ComponentRegistry", list_page_size=50)
```

When set, list operations paginate with `nextCursor`:

```python
from fastmcp import Client

async with Client(server) as client:
    # Auto-fetches all pages
    tools = await client.list_tools()

    # Manual pagination
    result = await client.list_tools_mcp()
    while result.nextCursor:
        result = await client.list_tools_mcp(cursor=result.nextCursor)
```

### Hot Reload Development

```bash
# Development mode with auto-reload
fastmcp dev server.py

# Explicit reload with custom directories
fastmcp run server.py --reload --reload-dir ./src --reload-dir ./lib
```

Changes trigger automatic restart without manual intervention.

### Composable Lifespans

Combine setup/teardown with `|` operator:

```python
from fastmcp import FastMCP
from fastmcp.server.lifespan import lifespan

@lifespan
async def db_lifespan(server):
    db = await connect_db()
    try:
        yield {"db": db}
    finally:
        await db.close()

@lifespan
async def cache_lifespan(server):
    cache = await connect_cache()
    try:
        yield {"cache": cache}
    finally:
        await cache.close()

# Compose lifespans - both enter in order, exit in reverse (LIFO)
mcp = FastMCP("server", lifespan=db_lifespan | cache_lifespan)
```

### PingMiddleware

Keep long-lived connections alive:

```python
from fastmcp.server.middleware import PingMiddleware

mcp = FastMCP("server")
mcp.add_middleware(PingMiddleware(interval_ms=5000))
```

Sends periodic pings to prevent connection timeouts.

### Context Transport Detection

Tools can detect transport type:

```python
@mcp.tool()
def adaptive_tool(ctx: Context) -> str:
    if ctx.transport == "stdio":
        return "short response"  # Local desktop client
    return "detailed response with more context"  # HTTP client

# Returns: "stdio", "sse", or "streamable-http"
```

### Automatic Threadpool

Synchronous tools automatically run in threadpool:

```python
import time

@mcp.tool()
def slow_tool():
    time.sleep(10)  # No longer blocks other requests
    return "done"

# Three concurrent calls execute in parallel (~10s)
# Not sequentially (30s) like in synchronous execution
```

## Summary

FastMCP 3.x encourages a composable, production-ready architecture:

**You focus on**:

- Writing Python functions (tools/prompts/resources)
- Proper types, docs, and safety checks
- Business logic and workflows
- Component versioning for API evolution
- Authorization rules for security

**Framework handles**:

- Protocol details (MCP 1.25+)
- Validation (Pydantic)
- Concurrency (async + threadpool)
- Transport layer (STDIO/HTTP/SSE)
- Provider/Transform composability
- Session state management
- Background task orchestration
- OpenTelemetry tracing
- Visibility and feature gating
- Component versioning and metadata

**Key FastMCP 3.x Improvements**:

- Composable architecture (Providers + Transforms)
- Native component versioning
- Session-scoped state with distributed backends
- Granular authorization (component + server level)
- Background tasks (SEP-1686 with Docket)
- Built-in OpenTelemetry instrumentation
- Dynamic visibility system
- Hot reload for development
- Automatic threadpool for sync tools
- Rich result classes for structured responses

By following these guidelines - clear function schemas, thorough validation, async for I/O, proper versioning, and careful auth/visibility handling - you build a reliable, scalable, production-ready MCP server that leverages the full power of FastMCP 3.x.

## Sources

- [FastMCP GitHub Repository (PrefectHQ/fastmcp)](https://github.com/prefecthq/fastmcp) (accessed 2026-02-22)
- [FastMCP Official Documentation](https://gofastmcp.com/) (accessed 2026-02-22)
- [FastMCP on PyPI](https://pypi.org/project/fastmcp/) (accessed 2026-02-22)
- [FastMCP 3.0 GA launch post](https://www.jlowin.dev/blog/fastmcp-3-launch) (accessed 2026-02-22)
- [MCP Protocol Specification](https://spec.modelcontextprotocol.io/) (accessed 2026-02-22)
- Authorization API — `fastmcp/server/auth/__init__.py` lines 8–13, `fastmcp/server/auth/authorization.py` lines 48, 78, 106 (verified against FastMCP 3.0.0rc2, 2026-02-21)
- Authorization middleware — `fastmcp/server/middleware/authorization.py` line 51 (verified against FastMCP 3.0.0rc2, 2026-02-21)
- Key–value storage backend — [py-key-value (strawgate/py-key-value)](https://github.com/strawgate/py-key-value) (accessed 2026-02-22)
