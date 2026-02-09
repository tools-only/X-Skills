# IoTHackBot Tool Development Guide

This guide documents the standard project structure and development patterns used for iothackbot tools.

## Project Structure Overview

All iothackbot tools follow a consistent architecture separating CLI, core functionality, and shared interfaces.

### Directory Structure
```
tools/iothackbot/
├── __init__.py                 # Package initialization
├── tool_name.py               # CLI entry point with argparse and colorama
├── core/
│   ├── tool_name_core.py      # Core tool logic implementing ToolInterface
│   └── interfaces.py          # Shared interfaces (ToolInterface, ToolConfig, etc.)
└── bin/
    └── tool_name             # Executable binary (imports from tools/iothackbot/)
```

## Development Patterns

### 1. Core Tool Implementation (`core/tool_name_core.py`)

```python
"""
Core tool_name functionality - Description.
Separated from CLI logic for automation and chaining.
"""

from .interfaces import ToolInterface, ToolConfig, ToolResult

class ToolNameTool(ToolInterface):
    """Tool implementation."""

    @property
    def name(self) -> str:
        return "tool_name"

    @property
    def description(self) -> str:
        return "What the tool does"

    def run(self, config: ToolConfig) -> ToolResult:
        """Execute the tool."""
        # Implementation here
        pass
```

### 2. CLI Implementation (`tool_name.py`)

```python
#!/usr/bin/env python3
import argparse
from colorama import init, Fore, Style
from .core.tool_name_core import ToolNameTool
from .core.interfaces import ConfigBuilder, OutputFormatter

class ToolNameOutputFormatter(OutputFormatter):
    """Custom output formatter for tool results."""

    def _format_text(self, result: 'ToolResult') -> str:
        """Format results as human-readable text."""
        if not result.success:
            return "\n".join(result.errors)
        # Custom formatting logic
        return formatted_output

def tool_name():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(description="Tool description")
    parser.add_argument("input", help="Input description")
    parser.add_argument("-o", "--output", help="Output option")
    parser.add_argument("--format", choices=['text', 'json', 'quiet'], default='text')
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    init()  # Initialize colorama

    config = ConfigBuilder.from_args(args, 'tool_name')
    tool = ToolNameTool()
    result = tool.run(config)

    formatter = ToolNameOutputFormatter()
    output = formatter.format_result(result, config.output_format)
    if output:
        print(output)

    return 0 if result.success else 1
```

### 3. Binary Executable (`bin/tool_name`)

```python
#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from iothackbot.tool_name import tool_name
if __name__ == "__main__":
    sys.exit(tool_name())
```

### 4. ConfigBuilder Integration

Add custom argument handling to `interfaces.py`:

```python
# In ConfigBuilder.from_args()
if hasattr(args, 'custom_arg'):
    custom_args['custom_arg'] = args.custom_arg
```

## Output Formatting Standards

### Color Scheme
- **Green**: Success messages, found items
- **Yellow**: Warnings, directory paths
- **Cyan**: Detailed information, file listings
- **Red**: Errors, failures

### Emoji Usage
- **NO EMOJIS** in tool output - use text labels instead
- Use descriptive text like `[HIGH RISK]`, `FAILED`, `SUCCESS` instead of emoji symbols
- Maintain professional, parseable output that works in all environments

### Format Support
All tools must support:
- `text`: Human-readable colored output (default)
- `json`: Structured JSON output
- `quiet`: Minimal/no output

## Development Workflow

1. **Create core functionality** in `core/tool_name_core.py`
2. **Implement CLI wrapper** in `tool_name.py`
3. **Create binary** in `bin/tool_name`
4. **Update ConfigBuilder** for custom arguments
5. **Test integration** with existing tools
6. **Add to registry** if needed for chaining

## Key Interfaces

### ToolInterface
- `name`: Tool identifier
- `description`: Tool purpose
- `run(config: ToolConfig) -> ToolResult`: Main execution

### ToolConfig
- `input_path`: Primary input file/directory
- `output_format`: text/json/quiet
- `verbose`: Enable verbose output
- `custom_args`: Tool-specific arguments

### ToolResult
- `success`: Boolean success indicator
- `data`: Tool-specific results
- `errors`: List of error messages
- `metadata`: Additional execution info
- `execution_time`: Performance metric

## Testing

- Test with various input types
- Verify all output formats work
- Check error handling
- Validate color output
- Test binary execution
- Ensure chaining compatibility

## Example Tool Template

Use this as a starting point for new tools:

```bash
# Create files
touch tools/iothackbot/new_tool.py
touch tools/iothackbot/core/new_tool_core.py
touch bin/new_tool
chmod +x bin/new_tool

# Implement following the patterns above
# Test integration with existing tools
```

## Chaining Support

Tools should be designed for chaining:
- Accept ToolResult from previous tools
- Return standardized ToolResult
- Support pipeline operations
- Maintain compatibility with iothackbot workflow
