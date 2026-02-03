---
title: Configuration Module
path: src/tunacode/configuration
type: directory
depth: 1
description: User settings, model registry, and pricing management
exports: [load_user_config, ModelRegistry, get_pricing]
seams: [M]
---

# Configuration Module

## Purpose
Manages application configuration including user settings, model registry, pricing information, and default values.

## Key Components

### settings.py
**load_user_config()**
- Loads user configuration from config directory
- Merges with defaults
- Validates configuration structure
- Returns UserConfig dictionary

**Config Locations:**
- `~/.config/tunacode/config.json`
- `.tunacode/config.json` (project-specific)

### models.py
**ModelRegistry Class**
- Loads model definitions from models_registry.json
- Provides model lookup by name
- Validates model configurations
- Caches registry for performance

**Model Registry Format:**
```json
{
  "models": {
    "claude-opus-4-5": {
      "api": "anthropic",
      "max_tokens": 200000,
      "supports_tool_use": true
    }
  }
}
```

### defaults.py
**DEFAULT_USER_CONFIG**
Default configuration values:
- **default_model** - Default LLM to use
- **theme** - UI theme preference
- **max_iterations** - Agent loop limit
- And more...

### pricing.py
**get_pricing()**
- Retrieves pricing information for models
- Calculates token costs
- Supports multiple providers (Anthropic, OpenAI, etc.)
- Tracks session costs

**Pricing Data:**
- Input token costs (per 1M tokens)
- Output token costs (per 1M tokens)
- Cached in pricing.json

## Configuration Schema

**UserConfig Structure:**
```python
{
  "default_model": str,
  "theme": str,
  "max_iterations": int,
  "timeout_seconds": int,
  # ... additional settings
}
```

## Integration Points

- **core/state.py** - Session configuration loading
- **core/agents/** - Model selection and agent creation
- **ui/** - Theme and preference application
- **types/** - UserConfig type definition

## Seams (M)

**Modification Points:**
- Add new configuration options
- Customize default values
- Extend model registry format
- Add new pricing models
