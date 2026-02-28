# CNS (Central Nervous System)

Conversation orchestration layer. Manages the Continuum aggregate (immutable conversation state),
LLM interaction, tool execution, memory surfacing, and segment lifecycle.

## Subdirectories

- `core/` — Immutable domain model: Continuum aggregate, Message value objects, domain events, stream events
- `services/` — Stateless service layer: orchestrator, subcortical processing, summary generation, segment collapse, feedback loop
- `infrastructure/` — Persistence and caching: PostgreSQL repository, Valkey message cache, Unit of Work, feedback storage
- `api/` — FastAPI endpoints: chat (HTTP + WebSocket), actions (domain routing), data queries, health, tool config
- `integration/` — Event bus (pub/sub) and factory (dependency wiring for the entire CNS graph)

## Cross-Cutting Patterns

### Event-Driven Coordination
All inter-component communication goes through `EventBus.publish(event)`. Subscribers register
by event class name string. Callbacks execute synchronously. New events belong to one of four
categories defined in `core/events.py`: MessageEvent, ToolEvent, WorkingMemoryEvent, ContinuumCheckpointEvent.
Always use the `.create()` classmethod on event subclasses — it auto-generates event_id, occurred_at,
and pulls user_id from contextvar.

### Immutability
`Message` and `ContinuumState` are frozen dataclasses. `ContinuumEvent` subclasses are frozen.
Mutations produce new objects. Segment sentinel metadata dicts are the one exception (mutated in place
for atomic DB updates via `jsonb_set`).

### User Context
All user-scoped operations rely on `utils.user_context.get_current_user_id()` contextvar.
Set once at the API boundary, flows through services → repositories → PostgreSQL RLS.
When spawning threads, use `contextvars.copy_context()` (see PeanutGalleryService for the pattern).

### LLM Configuration
Never hardcode model names or API keys. For `generate_response()` calls, use `internal_llm='purpose'`
param — it resolves endpoint, model, and API key from the internal_llm database schema. Purpose
keys include `'summary'`, `'analysis'`, `'domaindoc_summary'`, etc. For batch API paths that build
raw param dicts, use `get_internal_llm(purpose)` directly to access `.model` and `.effort`.

### XML Tag Parsing
All LLM output uses `<mira:tag>content</mira:tag>` format. Parse with regex allowing flexible
whitespace. Use `utils.tag_parser.TagParser` for standard tags. For custom tags, follow the
regex patterns in `subcortical.py` and `assessment_extractor.py`.

### Singleton Services
Infrastructure services use module-level `_instance` + `get_*()` / `initialize_*()` pairs.
The factory in `integration/factory.py` calls the initializers; everyone else calls the getters.
Don't create service instances directly outside the factory.

### Fail-Fast
Required infrastructure (DB, Valkey, embeddings, LLM) failures propagate. Never catch and return
None/[]/defaults for required services. Non-blocking exceptions are only acceptable for genuinely
optional features (feedback loop, retrieval logging, peanut gallery).

### Memory ID Convention
Memories use 8-char ID prefixes (first 8 of UUID, stripped of hyphens) formatted as `mem_XXXXXXXX`
via `utils.tag_parser.format_memory_id()`. All memory matching (pinning, retention)
uses this prefix, case-insensitive.

### Embedding Dimension
All embeddings are 768-dimensional via `get_hybrid_embeddings_provider()`. Use `encode_deep()`
for persistence-quality embeddings (segment summaries, memories). Generate once, pass to all consumers.
