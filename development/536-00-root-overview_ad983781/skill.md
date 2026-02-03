---
title: TunaCode Root Module
path: tunacode
type: directory
depth: 0
description: Core application package with constants, exceptions, and entry points
seams: [constants, exceptions, types]
---

# TunaCode Root Module (`src/tunacode`)

## Where
Located at the root of the `src/tunacode` package, serving as the top-level application module.

## What
The root module contains foundational elements that support the entire TunaCode CLI application:
- `constants.py` - Global constants, configuration values, UI themes, and magic numbers
- `exceptions.py` - Hierarchical exception system with enhanced error reporting
- `py.typed` - Type hint marker for PEP 561 compliance
- `__init__.py` - Package initialization

## How
The root module follows a **centralized foundation pattern**:

### Constants System (`constants.py`)
- **Numeric Constants**: File size limits, timeouts, viewport dimensions
- **String Constants**: File paths, environment variables, configuration keys
- **Enumerations**: Tool names for type-safe tool references
- **Theme Definitions**: Two complete UI themes (default "TunaCode" and "NeXTSTEP")
- **Error Messages**: Centralized error message templates
- **UI Constants**: Color palettes, styling directives, message formats

### Exception Hierarchy (`exceptions.py`)
```
TunaCodeError (base)
├── ConfigurationError
│   └── ModelConfigurationError
├── ValidationError
│   └── SetupValidationError
├── ToolExecutionError
│   ├── TooBroadPatternError
│   └── ToolBatchingJSONError
├── AgentError
├── StateError
├── ServiceError
│   └── GitOperationError
├── FileOperationError
├── UserAbortError
├── GlobalRequestTimeoutError
└── AggregateToolError
```

Each exception includes:
- **Contextual Attributes**: `original_error`, `suggested_fix`, `troubleshooting_steps`
- **Enhanced Messages**: Formatted with emojis and actionable guidance
- **Recovery Commands**: Suggested shell commands for error resolution

## Why
**Centralization Strategy**: By collocating constants and exceptions at the root level:
- **Maintainability**: Single source of truth prevents duplication
- **Consistency**: All modules reference identical error messages and constants
- **Type Safety**: Enumerations prevent typos in tool names and configuration keys
- **Developer Experience**: Import errors show clear paths to root module
- **Theme Consistency**: UI styling is centralized for coherent visual design

**Design Principles**:
1. **Fail Fast, Fail Loud**: Exceptions carry rich context and recovery guidance
2. **No Magic Numbers**: All numeric constants have descriptive names
3. **Symbolic Constants**: String literals are replaced with named constants
4. **Modular Enhancement**: Exception attributes enable progressive disclosure of error details

## Integration Points
- **All submodules** import from `tunacode.constants` for configuration values
- **Tool implementations** raise `ToolExecutionError` for consistent error handling
- **Configuration system** raises `ConfigurationError` for setup issues
- **UI components** use `UI_COLORS` and theme builders for visual consistency

## Code Quality Notes
- **Type Hints**: Full type annotation coverage throughout
- **Enum Usage**: `ToolName` enum prevents string-based tool reference errors
- **Documentation**: Comprehensive docstrings explain constant purposes
- **Theme Builders**: Lazy import of `textual.theme` prevents import cycles
