# Toolsets API

## create_console_toolset

::: pydantic_ai_backends.toolsets.console.create_console_toolset
    options:
      show_root_heading: true

## get_console_system_prompt

::: pydantic_ai_backends.toolsets.console.get_console_system_prompt
    options:
      show_root_heading: true

## ConsoleDeps

::: pydantic_ai_backends.toolsets.console.ConsoleDeps
    options:
      show_root_heading: true

## Console Tools

The toolset provides these tools:

### ls

```python
async def ls(ctx: RunContext[ConsoleDeps], path: str = ".") -> str:
    """List files and directories at the given path.

    Args:
        path: Directory path to list. Defaults to current directory.
    """
```

### read_file

```python
async def read_file(
    ctx: RunContext[ConsoleDeps],
    path: str,
    offset: int = 0,
    limit: int = 2000,
) -> str:
    """Read file content with line numbers.

    Args:
        path: Path to the file to read.
        offset: Line number to start reading from (0-indexed).
        limit: Maximum number of lines to read.
    """
```

### write_file

```python
async def write_file(
    ctx: RunContext[ConsoleDeps],
    path: str,
    content: str,
) -> str:
    """Write content to a file (creates or overwrites).

    Args:
        path: Path to the file to write.
        content: Content to write to the file.
    """
```

### edit_file

```python
async def edit_file(
    ctx: RunContext[ConsoleDeps],
    path: str,
    old_string: str,
    new_string: str,
    replace_all: bool = False,
) -> str:
    """Edit a file by replacing strings.

    Args:
        path: Path to the file to edit.
        old_string: String to find and replace.
        new_string: Replacement string.
        replace_all: If True, replace all occurrences.
    """
```

### glob

```python
async def glob(
    ctx: RunContext[ConsoleDeps],
    pattern: str,
    path: str = ".",
) -> str:
    """Find files matching a glob pattern.

    Args:
        pattern: Glob pattern to match (e.g., "**/*.py").
        path: Base directory to search from.
    """
```

### grep

```python
async def grep(
    ctx: RunContext[ConsoleDeps],
    pattern: str,
    path: str | None = None,
    glob_pattern: str | None = None,
    output_mode: Literal["content", "files_with_matches", "count"] = "files_with_matches",
    ignore_hidden: bool = True,
) -> str:
    """Search for a regex pattern in files.

    Args:
        pattern: Regex pattern to search for.
        path: Specific file or directory to search.
        glob_pattern: Glob pattern to filter files.
        output_mode: Output format.
        ignore_hidden: Whether to skip hidden files (defaults to the toolset setting).
    """
```

### execute

```python
async def execute(
    ctx: RunContext[ConsoleDeps],
    command: str,
    timeout: int | None = 120,
) -> str:
    """Execute a shell command.

    Args:
        command: The shell command to execute.
        timeout: Maximum execution time in seconds.
    """
```
