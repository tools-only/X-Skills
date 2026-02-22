# Console Toolset

The console toolset provides ready-to-use pydantic-ai tools for file operations and shell execution.

!!! info "Requires pydantic-ai"
    ```bash
    pip install pydantic-ai-backend[console]
    ```

## Quick Start

```python
from dataclasses import dataclass
from pydantic_ai import Agent
from pydantic_ai_backends import LocalBackend, create_console_toolset

@dataclass
class Deps:
    backend: LocalBackend

# Create toolset
toolset = create_console_toolset()

# Create agent with tools
agent = Agent("openai:gpt-4o", deps_type=Deps)
agent = agent.with_toolset(toolset)

# Run
backend = LocalBackend(root_dir="/workspace")
result = agent.run_sync("List all Python files", deps=Deps(backend=backend))
```

## Available Tools

| Tool | Description |
|------|-------------|
| `ls` | List files in a directory |
| `read_file` | Read file content with line numbers |
| `write_file` | Create or overwrite a file |
| `edit_file` | Replace strings in a file |
| `glob` | Find files matching a pattern |
| `grep` | Search for patterns in files |
| `execute` | Run shell commands (optional) |

## Configuration

```python
# Default: execute enabled, requires approval
toolset = create_console_toolset()

# Without shell execution
toolset = create_console_toolset(include_execute=False)

# Auto-approve writes
toolset = create_console_toolset(
    require_write_approval=False,
    require_execute_approval=False,
)

# Custom toolset ID
toolset = create_console_toolset(id="my-console")

# Include hidden files by default for grep
toolset = create_console_toolset(default_ignore_hidden=False)
```

## Permission-based Configuration

For fine-grained control, use the permission system:

```python
from pydantic_ai_backends import create_console_toolset
from pydantic_ai_backends.permissions import (
    DEFAULT_RULESET,
    READONLY_RULESET,
    PermissionRuleset,
    OperationPermissions,
)

# Use pre-configured presets
toolset = create_console_toolset(permissions=DEFAULT_RULESET)

# Read-only toolset
toolset = create_console_toolset(permissions=READONLY_RULESET)

# Custom permissions
custom = PermissionRuleset(
    write=OperationPermissions(default="allow"),  # No approval needed
    execute=OperationPermissions(default="ask"),   # Requires approval
)
toolset = create_console_toolset(permissions=custom)
```

When `permissions` is provided, it overrides the legacy `require_write_approval` and `require_execute_approval` flags.

See [Permissions](permissions.md) for full documentation.

## Image Support

When working with multimodal models (e.g., GPT-4o, Claude), you can enable image support so that `read_file` returns image data that the model can see, instead of garbled binary text.

```python
# Enable image support
toolset = create_console_toolset(image_support=True)
```

When `image_support=True`, reading a file with a recognized image extension returns a `BinaryContent` object that pydantic-ai sends to the model as an inline image. For all other file types, `read_file` behaves normally and returns text.

### Recognized Image Types

| Extension | Media Type |
|-----------|------------|
| `.png` | `image/png` |
| `.jpg` | `image/jpeg` |
| `.jpeg` | `image/jpeg` |
| `.gif` | `image/gif` |
| `.webp` | `image/webp` |

### Size Limits

Large images are rejected to avoid excessive token usage. The default limit is 50 MB:

```python
# Default: 50 MB max
toolset = create_console_toolset(image_support=True)

# Custom limit: 5 MB max
toolset = create_console_toolset(
    image_support=True,
    max_image_bytes=5 * 1024 * 1024,
)
```

Images exceeding the limit return an error message like: `Error: Image 'photo.png' too large (12.3MB, max 5.0MB)`.

### Example: Visual Analysis Agent

```python
from dataclasses import dataclass
from pydantic_ai import Agent
from pydantic_ai_backends import LocalBackend, create_console_toolset

@dataclass
class Deps:
    backend: LocalBackend

# Enable image support for multimodal model
toolset = create_console_toolset(image_support=True)

agent = Agent(
    "openai:gpt-4o",  # Multimodal model
    system_prompt="You can read and analyze images using read_file.",
    deps_type=Deps,
)
agent = agent.with_toolset(toolset)

result = agent.run_sync(
    "Read screenshot.png and describe what you see",
    deps=Deps(backend=LocalBackend(root_dir="/workspace")),
)
```

!!! tip "When to enable image support"
    Only enable `image_support` when using a multimodal model that can process images.
    With text-only models, image data will be wasted tokens. When disabled (the default),
    reading an image file returns the raw text representation, which is not useful
    but avoids unexpected behavior.

## ConsoleDeps Protocol

Your dependencies class must have a `backend` property:

```python
from pydantic_ai_backends import BackendProtocol

class ConsoleDeps(Protocol):
    @property
    def backend(self) -> BackendProtocol:
        ...
```

Any class with a `backend` attribute works:

```python
from dataclasses import dataclass

@dataclass
class MyDeps:
    backend: LocalBackend
    user_id: str  # Additional fields are fine
```

## System Prompt

Include the console system prompt for better tool usage:

```python
from pydantic_ai_backends import get_console_system_prompt

system_prompt = f"""You are a helpful coding assistant.

{get_console_system_prompt()}
"""

agent = Agent(
    "openai:gpt-4o",
    system_prompt=system_prompt,
    deps_type=Deps,
)
```

## Tool Details

### ls

```python
async def ls(ctx, path: str = ".") -> str:
    """List files and directories at the given path."""
```

### read_file

```python
async def read_file(ctx, path: str, offset: int = 0, limit: int = 2000) -> str:
    """Read file content with line numbers."""
```

### write_file

```python
async def write_file(ctx, path: str, content: str) -> str:
    """Write content to a file (creates or overwrites)."""
```

### edit_file

```python
async def edit_file(
    ctx, path: str, old_string: str, new_string: str, replace_all: bool = False
) -> str:
    """Edit a file by replacing strings."""
```

### glob

```python
async def glob(ctx, pattern: str, path: str = ".") -> str:
    """Find files matching a glob pattern."""
```

### grep

```python
async def grep(
    ctx,
    pattern: str,
    path: str | None = None,
    glob_pattern: str | None = None,
    output_mode: str = "files_with_matches",
    ignore_hidden: bool = True,
) -> str:
    """Search for a regex pattern in files.

    Args:
        pattern: Regex pattern to search for.
        path: Optional file/directory scope.
        glob_pattern: Glob filter applied before searching.
        output_mode: "content", "files_with_matches", or "count".
        ignore_hidden: Whether to skip hidden files (defaults to the toolset setting).
    """
```

### execute

```python
async def execute(ctx, command: str, timeout: int | None = 120) -> str:
    """Execute a shell command."""
```

## With Different Backends

The toolset works with any backend:

=== "LocalBackend"

    ```python
    backend = LocalBackend(root_dir="/workspace")
    ```

=== "StateBackend"

    ```python
    backend = StateBackend()
    ```

=== "DockerSandbox"

    ```python
    backend = DockerSandbox(image="python:3.12-slim")
    ```

## Next Steps

- [Permissions](permissions.md) - Fine-grained access control
- [CLI Agent Example](../examples/cli-agent.md) - Build a CLI coding assistant
- [API Reference](../api/toolsets.md) - Complete API
