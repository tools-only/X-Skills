# CLAUDE.md

Guidance for Claude Code when working on this repository.

## What This Project Is

**pydantic-ai-backend** provides file storage and sandbox backends for AI agents. It's designed to work with pydantic-ai and pydantic-deep.

Key pattern: **Protocol-based backends** - all backends implement `BackendProtocol` for consistent file operations.

## Commands

```bash
uv sync --all-extras --group dev  # Install all dependencies
uv run pytest                      # Run tests
uv run coverage run -m pytest && uv run coverage report  # Test with coverage
uv run ruff check .                # Lint
uv run ruff format .               # Format
uv run pyright                     # Type check
uv run mypy src/pydantic_ai_backends  # MyPy check
```

## Structure

```
src/pydantic_ai_backends/
├── types.py          # FileData, FileInfo, WriteResult, EditResult, etc.
├── protocol.py       # BackendProtocol, SandboxProtocol
├── state.py          # StateBackend (in-memory)
├── filesystem.py     # FilesystemBackend (real filesystem)
├── composite.py      # CompositeBackend (routing)
├── sandbox.py        # BaseSandbox, DockerSandbox, LocalSandbox
├── session.py        # SessionManager
├── runtimes.py       # RuntimeConfig, BUILTIN_RUNTIMES
└── __init__.py       # Public API with lazy loading
```

## Core Pattern

```python
class BackendProtocol(Protocol):
    def ls_info(self, path: str) -> list[FileInfo]: ...
    def read(self, path: str, offset: int = 0, limit: int = 2000) -> str: ...
    def write(self, path: str, content: str | bytes) -> WriteResult: ...
    def edit(self, path: str, old: str, new: str, replace_all: bool = False) -> EditResult: ...
    def glob_info(self, pattern: str, path: str = "/") -> list[FileInfo]: ...
    def grep_raw(self, pattern: str, path: str | None = None, glob: str | None = None) -> list[GrepMatch] | str: ...

class SandboxProtocol(BackendProtocol, Protocol):
    def execute(self, command: str, timeout: int | None = None) -> ExecuteResponse: ...
    @property
    def id(self) -> str: ...
```

## Requirements

- **100% test coverage** - every PR must maintain this
- **Type annotations** - pyright and mypy strict mode
- **Lazy loading** - optional deps (docker, pypdf, chardet) loaded on-demand

## Testing

```bash
# Run specific test
uv run pytest tests/test_backends.py::TestStateBackend -v

# Debug mode
uv run pytest -v -s
```

## Integration

This library is used by [pydantic-deep](https://github.com/vstorm-co/pydantic-deep) which re-exports its API. Changes here affect pydantic-deep users.
