# clients/ — Infrastructure Clients & External Service Adapters

Thin clients wrapping external infrastructure (PostgreSQL, Valkey, Vault, LLMs, embeddings, federation). Each client owns its connection lifecycle and credential sourcing. Vault is the foundation — all other clients get secrets from it. Clients follow fail-fast semantics: infrastructure failures propagate, never return silent defaults.

## Files

- **`__init__.py`** — Re-exports the public API: `HybridEmbeddingsProvider`, `LLMProvider`, `PostgresClient`, `SQLiteClient`, `ValkeyClient` + their factory functions. Note: `LatticeClient` and `FilesManager` are NOT exported here — import them directly
- **`vault_client.py`** — HashiCorp Vault client (AppRole auth). Foundational — every other client depends on this for secrets. Key functions: `get_database_url()`, `get_api_key()`, `get_auth_secret()`, `get_service_config(field)`, `preload_secrets()`. Exports `VaultHealthCheck` TypedDict. `preload_secrets()` is fail-fast — raises `RuntimeError` if any secret group fails to load
- **`postgres_client.py`** — `PostgresClient` with connection pooling and automatic RLS user isolation via `SET app.current_user_id`. Dict-style results via raw SQL (`execute_query`, `execute_single`, `execute_update`, `execute_insert`, `execute_returning`, `execute_scalar`, `execute_transaction`). `get_pool_stats()` returns `Dict[str, PoolStats]` TypedDict. **Note:** All execute_* methods are monkey-patched by `utils/perf.py` at startup when the `mira.perf` logger is at INFO or DEBUG — see `install_db_instrumentation()`
- **`valkey_client.py`** — `ValkeyClient` for caching, sessions, rate limiting. Dual sync/async clients plus a binary client for raw bytes (numpy arrays). Factory: `get_valkey()` / `get_valkey_client()`
- **`sqlite_client.py`** — `SQLiteClient` for per-user tool data. Manual user_id filtering (no RLS). Factory: `get_sqlite_client(db_path, user_id)`. Methods: `execute_query`, `execute_insert`, `execute_update`, `execute_delete`, `create_table`
- **`llm_provider.py`** — `LLMProvider` — universal LLM entry point. Anthropic SDK primary, emergency failover to generic OpenAI-compatible providers. Streaming events, parallel tool execution with circuit breaker, prompt caching, billing integration. NOT a singleton — instantiated per-use. ALWAYS use `LLMProvider.generate_response()` for ALL LLM calls — never instantiate GenericOpenAIClient directly. **`internal_llm=` param** is the single source of truth for all LLM-tuning params — resolves endpoint, model, API key, max_tokens, and effort from `InternalLLMConfig` (DB-backed). Callers just pass `internal_llm='purpose'`. Caller-provided params override DB values (explicit > implicit). Mutually exclusive with `endpoint_url`/`model_override`/`api_key_override` (which remain for one-off overrides). Thinking/effort control: auto-resolved from DB config's `effort` field, or overridden via `effort` (named level) or `thinking_tokens` (exact budget) params. Resolved into a `ThinkingConfig` frozen dataclass at entry, translated per-backend by `_anthropic_thinking_params()` (adaptive for 4.6, budget for older) and `_generic_thinking_params()` (effort passthrough). **Two public entry points**: `generate_response()` (convenience wrapper adding `internal_llm` DB resolution + non-streaming consumption) and `stream_events()` (self-sufficient streaming entry point with explicit keyword params for all LLM-tuning: `effort`, `thinking_tokens`, `thinking`, `temperature`, `max_tokens`, `container_id`, `endpoint_url`, `model_override`, `model_preference`, `api_key_override`, `system_override`, `allow_negative`). `stream_events()` resolves `effort`/`thinking_tokens` → `ThinkingConfig` internally — direct callers (orchestrator) get the same resolution as `generate_response()` callers. **No `**kwargs` anywhere in the call chain** — all 5 methods (`generate_response`, `stream_events`, `_generate_non_streaming`, `_stream_response`, `_execute_with_tools`) use explicit keyword-only params so Python raises `TypeError` on unknown arguments at the call site. `EFFORT_BUDGET_MAP` maps named levels to token counts. **`build_batch_params(purpose, system_prompt, messages, *, cache_ttl=None)`** — standard way to construct Anthropic Batch API param dicts. All LLM-tuning params (model, max_tokens, effort) resolved from `InternalLLMConfig` — callers just pass the purpose key. Two states: effort set → thinking/effort params; effort NULL → vanilla API call (no thinking, no temperature). Wraps system prompt in content blocks with `cache_control`, optional `cache_ttl="1h"` for paid 1-hour caching. All 4 batch sites (extraction, relationship, consolidation, entity_gc) use this. Exports `ThinkingConfig` dataclass, `EFFORT_BUDGET_MAP` dict, `build_batch_params` function, `ToolCall` TypedDict, `ToolExecution` NamedTuple, `ToolExecutionResult` NamedTuple. Default model from `config.api.model`
- **`hybrid_embeddings_provider.py`** — `HybridEmbeddingsProvider` — local asymmetric embeddings via `MongoDB/mdbr-leaf-ir-asym` (768-dim). Separate `encode_realtime()` (query) and `encode_deep()` (document) methods. Valkey-cached. Singleton via `get_hybrid_embeddings_provider()`
- **`lattice_client.py`** — `LatticeClient` — thin HTTP client for Lattice federation service. Only used by pager_tool. Methods: `send_message() -> SendMessageResponse`, `get_identity() -> LatticeIdentity`. Singleton via `get_lattice_client()`
- **`files_manager.py`** — `FilesManager` — Anthropic Files API operations (upload/delete) with segment-scoped cleanup. NOT a singleton — requires an initialized Anthropic client

## Patterns to Follow

### Credential Sourcing
1. **VaultClient itself** — configured via env vars (`VAULT_ADDR`, `VAULT_ROLE_ID`, `VAULT_SECRET_ID`). The only place env vars are the primary source
2. **All clients** — source secrets from Vault: `get_database_url()`, `get_api_key()`, `get_service_config(field)`. No env var fallbacks — Vault failure propagates

### Singleton Factories
Most clients use a module-level singleton pattern with a `get_*()` factory function. Use these factories instead of constructing clients directly:
- `get_valkey()` / `get_valkey_client()` — ValkeyClient
- `get_hybrid_embeddings_provider()` — HybridEmbeddingsProvider
- `get_sqlite_client(db_path, user_id)` — SQLiteClient (cached per user_id:path)
- `get_lattice_client()` — LatticeClient
- **Exceptions**: `LLMProvider` (instantiated per-use), `FilesManager` (requires Anthropic client), `PostgresClient` (shared pool but not singleton instances)

### Connection Pooling
- **PostgresClient**: Class-level `ConnectionPool` dict keyed by database name. Thread-safe via `RLock`. Max 20 connections, with automatic connection recycling (1 hour lifetime) and idle timeout (5 minutes)
- **ValkeyClient**: Module-level global `ConnectionPool`, max 50 connections. Dual sync + async + binary clients
- **SQLiteClient**: No pooling (appropriate for SQLite — fresh connection per request)
- **LLMProvider**: Anthropic SDK manages its own pooling internally

### Billing: Auto-Resolved pricing_key
The billing hooks in `generate_response()` auto-resolve `pricing_key` by matching `(model, endpoint_url)` against config caches via `billing.pricing.resolve_pricing_key()`. Callsites do NOT need to pass `pricing_key`. Unknown model+endpoint in user context raises `BillingConfigurationError` (fails closed). No user context (system/startup calls) → billing skipped.

### Fail-Fast Infrastructure
All clients propagate infrastructure failures. Never catch and return `None`/`[]`/defaults. Specific patterns:
- `ValkeyClient.__init__()` pings on construction — fails immediately if unreachable
- `ValkeyClient._load_config()` sources URL from Vault — fails if Vault is unavailable
- `VaultClient` raises typed exceptions: `FileNotFoundError`, `PermissionError`, `RuntimeError`, `KeyError`
- `PostgresClient` lets pool creation and connection exhaustion errors propagate
- **Only exception**: `FilesManager.delete_file()` logs warnings on 404 (intentional — cleanup is best-effort)

### RLS User Context
`PostgresClient.get_connection()` automatically sets `app.current_user_id` via SQL `SET` on checkout and `RESET` on release. This is the core mechanism for PostgreSQL Row Level Security. No manual user filtering needed in queries.

### Deferred Imports
Clients use inline imports to break circular dependencies (e.g., `PostgresClient.__init__()` imports `vault_client.get_database_url` inline). Follow this pattern when adding new cross-references.

### Vault Preloading
Call `vault_client.preload_secrets()` at startup to bulk-load all secrets into memory cache. This avoids individual lookups and token expiration issues during long-running service operation. Raises `RuntimeError` if any secret group (API keys, database, auth) fails to load — this is intentional fail-fast behavior.

### LLMProvider User Context Propagation
`LLMProvider._execute_with_tools()` explicitly propagates contextvars to `ThreadPoolExecutor` worker threads. If you add new threaded execution paths, use `contextvars.copy_context()` to preserve user context for RLS enforcement.
