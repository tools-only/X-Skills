# API Reference

Complete API documentation for pydantic-ai-backend.

## Quick Example

```python
from dataclasses import dataclass
from pydantic_ai import Agent
from pydantic_ai_backends import LocalBackend, create_console_toolset

@dataclass
class Deps:
    backend: LocalBackend

backend = LocalBackend(root_dir=".")
toolset = create_console_toolset()
agent = Agent("openai:gpt-4o", deps_type=Deps).with_toolset(toolset)

result = agent.run_sync("Create hello.py and run it", deps=Deps(backend=backend))
```

## Modules

| Module | Description |
|--------|-------------|
| [Backends](backends.md) | LocalBackend, StateBackend, CompositeBackend |
| [Docker](docker.md) | DockerSandbox, SessionManager, RuntimeConfig |
| [Toolsets](toolsets.md) | Console toolset for pydantic-ai |
| [Types](types.md) | Type definitions |

## Import Reference

```python
# Toolset for pydantic-ai agents (requires [console] extra)
from pydantic_ai_backends import (
    create_console_toolset,
    get_console_system_prompt,
    ConsoleDeps,
)

# Backends
from pydantic_ai_backends import (
    LocalBackend,
    StateBackend,
    CompositeBackend,
)

# Docker (requires [docker] extra)
from pydantic_ai_backends import (
    DockerSandbox,
    SessionManager,
    RuntimeConfig,
    BUILTIN_RUNTIMES,
)

# Types
from pydantic_ai_backends import (
    FileInfo,
    WriteResult,
    EditResult,
    ExecuteResponse,
    GrepMatch,
)

# Protocols
from pydantic_ai_backends import (
    BackendProtocol,
    SandboxProtocol,
)
```

## Protocols

### BackendProtocol

All backends implement this interface:

```python
class BackendProtocol(Protocol):
    def ls_info(self, path: str) -> list[FileInfo]: ...
    def read(self, path: str, offset: int = 0, limit: int = 2000) -> str: ...
    def write(self, path: str, content: str | bytes) -> WriteResult: ...
    def edit(self, path: str, old: str, new: str, replace_all: bool = False) -> EditResult: ...
    def glob_info(self, pattern: str, path: str = ".") -> list[FileInfo]: ...
    def grep_raw(
        self,
        pattern: str,
        path: str | None = None,
        glob: str | None = None,
        ignore_hidden: bool = True,
    ) -> list[GrepMatch] | str: ...
```

### SandboxProtocol

Extends BackendProtocol with execution:

```python
class SandboxProtocol(BackendProtocol, Protocol):
    def execute(self, command: str, timeout: int | None = None) -> ExecuteResponse: ...
    @property
    def id(self) -> str: ...
```
