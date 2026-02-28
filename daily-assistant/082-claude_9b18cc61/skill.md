# config/ — Application Configuration & Prompt Templates

All application configuration defined as Pydantic BaseModel schemas with hardcoded defaults. No YAML/TOML/JSON config files — the system boots with sane defaults and secrets come from Vault at runtime. The global `config` singleton is the canonical access point for all settings.

## Files

- **`__init__.py`** — Re-exports `AppConfig` and the `config` singleton. Also initializes `tools.registry` first (import order matters for circular dep avoidance)
- **`config.py`** — Pure schema definitions (~800 lines of Pydantic BaseModel classes). Every config section lives here: `ApiConfig`, `EmbeddingsConfig`, `LTMemoryConfig`, `PeanutGalleryConfig`, etc. No logic, no side effects — just data shapes with defaults and descriptions
- **`config_manager.py`** — Runtime orchestrator. Defines `AppConfig` (root config aggregating all sections), handles system prompt loading, provides `config.get("api.max_tokens")` dot-notation access, dynamic tool config via `__getattr__`, and creates the singleton via `config = initialize_config()` at module load
- **`system_prompt.txt`** — Mira's core system prompt (identity, behavioral directives, output format). Loaded once at startup by `AppConfig._load_system_prompt()`. Template variables: `{first_name}`, `{user_context}`, `{relative time since account creation}` — substituted by `working_memory/core.py` at runtime
- **`announcement.py`** — Announcement system module. Reads `announcement.json` once at startup into a module-level cache. `load_announcement()` called during lifespan, `get_cached_announcement()` used by `cns/api/data.py` to serve to frontend
- **`announcement.json`** — Live announcement config. Set `id` + `message` to show a banner, set both to `null` to hide. Requires app restart
- **`announcement.sample.json`** — Reference example for announcement format
- **`vault.hcl`** — Local dev HashiCorp Vault config (file storage, TCP on 127.0.0.1:8200, no TLS). Consumed by Vault binary, not Python

## Subdirectories

- **`prompts/`** — LLM prompt templates for all subsystems (extraction, feedback, synthesis, subcortical, etc.). See `prompts/CLAUDE.md`

## Key Entry Points

| What | Import | Used by |
|------|--------|---------|
| `config` singleton | `from config import config` | ~60+ files — the universal config access pattern |
| Schema classes | `from config.config import ExtractionConfig` | Services accepting config sections as constructor params |
| `config.system_prompt` | Property on singleton | Prompt composition layer |
| `config.<tool_name>` | Dynamic `__getattr__` | Auto-resolves tool configs via registry |

## Config Sections

| Access path | Schema class | Scope |
|-------------|-------------|-------|
| `config.api` | `ApiConfig` | Anthropic API (default model, max_tokens, temperature, timeout), emergency fallback, thinking budgets, subcortical toggle |
| `config.api_server` | `ApiServerConfig` | FastAPI host/port, CORS, rate limiting |
| `config.tools` | `ToolConfig` | Essential tools list, idle thresholds |
| `config.system` | `SystemConfig` | Logging, timezone, streaming, segments, manifest, session cache |
| `config.embeddings` | `EmbeddingsConfig` | Fast model settings |
| `config.lt_memory` | `LTMemoryConfig` | All LT_Memory subsystems (extraction, batching, linking, refinement, proactive, search, entity GC, scheduled jobs). `ScheduledJobsConfig` uses `*_use_days` fields — intervals measured in user activity days via modular arithmetic, not calendar days |
| `config.context` | `ContextConfig` | Context window management, topic drift |
| `config.peanutgallery` | `PeanutGalleryConfig` | Metacognitive observer settings |
| `config.lattice` | `LatticeConfig` | Federation service URL/timeout |

## Patterns to Follow

### Accessing Config
Always use the singleton: `from config import config`. Then dot-access: `config.api.max_tokens`, `config.lt_memory.extraction.fuzzy_dedup_threshold`. For dynamic key access: `config.get("api.max_tokens")` or `config.require("api.max_tokens")` (raises on None).

### Adding New Config
1. Define a Pydantic `BaseModel` in `config.py` with `Field(default=..., description=...)`
2. Add it as a field on `AppConfig` in `config_manager.py`
3. Every field must have a default — no required fields (the system must boot without external config)

### Secrets vs Config
- **Config values** (thresholds, model names, timeouts): Hardcoded defaults in Pydantic schemas
- **Secrets** (API keys, DB URLs, passwords): Lazy Vault lookups via properties on `AppConfig` (e.g., `config.api_key`). These use deferred imports to avoid circular deps

### Dynamic Tool Config
Accessing `config.<tool_name>_tool` triggers `__getattr__` → registry lookup. If the tool registered a custom Pydantic config class, it's used; otherwise a default `{ToolName}Config(enabled=True)` is auto-generated.

### System Prompt Template Variables
`system_prompt.txt` uses `{first_name}`, `{user_context}`, and `{relative time since account creation}` — substituted by the working memory composition layer (`working_memory/core.py:_handle_compose_prompt`), not by config itself. LoRA behavioral directives are injected as a separate system prompt section by `LoraTrinket`, not inline.

### Prompt Loading Convention
Each service loads its own prompts from `config/prompts/` via `Path("config/prompts")`. There is no centralized loader. See `prompts/CLAUDE.md` for the full catalog.
