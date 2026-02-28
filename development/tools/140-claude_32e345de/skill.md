# cns/services/

Stateless service layer. Services never touch the database directly — they go through
repositories in `cns/infrastructure/` or receive data from the orchestrator.

## Files

- `orchestrator.py` — Main flow coordinator (~1300 lines). Owns the message→LLM→response cycle: memory surfacing pipeline, system prompt composition, LLM provider routing, context overflow remediation (3-tier), tool execution loop. Singleton via `get_orchestrator()`. Exports `TurnMetadata` TypedDict (process_message return metadata), `LLMKwargs` TypedDict (kwargs passed to stream_events). `MemorySurfacingResult.surfaced_memories` typed as `list[SurfacedMemory]` (from subcortical.py).
- `subcortical.py` — Pre-LLM processing. `generate()` returns `SubcorticalResult` (never None — raises on failure, orchestrator catches). Generates retrieval query expansion (replaces original query), evaluates memory retention (pinned IDs), extracts named entities, assesses complexity (straightforward/complex → thinking budget). Uses `<query_expansion>`, `<entities>`, `<complexity>`, `<passage>` tags.
- `summary_generator.py` — Segment summary generation via LLM with thinking enabled. Handles context overflow via hierarchical chunking. Returns `SummaryResult(synopsis, display_title, complexity)` dataclass. Previous summaries passed for narrative continuity.
- `segment_collapse_handler.py` — Handles `SegmentTimeoutEvent`. Orchestrates: sentinel lookup → `_generate_summary()` → embedding → sentinel collapse → downstream extraction → feedback loop. `lt_memory_factory` is required (not optional). Synthesis trigger is a direct `should_synthesize()` call (modular arithmetic on `cumulative_activity_days`). Singleton via `get_segment_collapse_handler()`.
- `segment_helpers.py` — Pure utility functions for segment sentinels: creation, collapse, display formatting, marker creation. The canonical place for sentinel operations.
- `segment_timeout_service.py` — Scheduled job (every 5 min). Admin-level query across all users for active segments past timeout threshold. Publishes `SegmentTimeoutEvent`. `check_timeouts()` returns `TimeoutCheckResult` TypedDict.
- `assessment_extractor.py` — User model pipeline step 1: evaluates conversation against system prompt sections, producing alignment/misalignment/contextual_pass signals with evidence. Sees anonymized system prompt and user model context. Uses `get_internal_llm('assessment')`.
- `user_model_synthesizer.py` — User model pipeline step 2: synthesizes section-anchored observations from assessment signals. Includes critic validation loop (Haiku, up to 3 attempts) returning `CriticResult(passed, feedback)`. Exports `UserObservation`, `CriticResult`, `SynthesisResult`. Uses `get_internal_llm('synthesis')` for synthesis, `get_internal_llm('critic')` for validation.
- `system_prompt_parser.py` — Pure utility: extracts section IDs and content from system prompt XML. Provides anonymization and blocklisting. Used by assessment_extractor and user_model_synthesizer.
- `memory_relevance_service.py` — Thin wrapper around `lt_memory.ProactiveService`. Validates 768d embedding, delegates search. Returns `list[MemoryDict]`.
- `peanutgallery_model.py` — Two-stage LLM pipeline: fast prerunner filters seed memories, Sonnet observer evaluates conversation. Returns noop/compaction/concern/coaching.
- `peanutgallery_service.py` — Fire-and-forget async service. Subscribes to `TurnCompletedEvent`, runs observation every N turns via ThreadPoolExecutor with `contextvars.copy_context()`.
- `domaindoc_summary_service.py` — Generates one-sentence summaries for DomainDoc sections via Gemini Flash. Module-level singletons with lazy init.
- `manifest_query_service.py` — Retrieves segment data for manifest display. Valkey-cached with event-driven invalidation on `ManifestUpdatedEvent`. Singleton via `initialize_manifest_query_service()` / `get_manifest_query_service()` split. Exports `ManifestSegment` TypedDict for segment data shape.

## Patterns to Follow

### Service Structure
- Services receive dependencies via constructor injection. Use `internal_llm='purpose'` param on `generate_response()` for LLM routing — it resolves endpoint, model, and API key internally. For one-off overrides, use explicit `endpoint_url`/`model_override`/`api_key_override` params.
- Prompts live in `config/prompts/` as `.txt` files. Load them in `__init__` or a `_load_prompts()` method. Never inline prompts as string constants.
- Module-level singletons follow `_instance = None` + `get_*()` / `initialize_*()` pattern.

### Billing: Auto-Resolved
The billing hook in `LLMProvider.generate_response()` auto-resolves `pricing_key` from config caches. Callsites do NOT pass `pricing_key` — it's matched by `(model, endpoint_url)` at the billing hook layer.

### LLM Output Parsing
- All LLM responses use `<mira:tag>` XML format with regex parsing.
- Always allow flexible whitespace in regex: `<mira:tag\s*type\s*=\s*"..."`.
- Return empty list / log warning on parse failure for non-critical paths.
- Use `json_repair` library for JSON-in-LLM responses (see subcortical.py).

### Async / Threading
- Fire-and-forget async uses `ThreadPoolExecutor` with `contextvars.copy_context().run()`.
- Non-critical services (peanut gallery, feedback loop, logging) swallow exceptions and log.
- Critical path (orchestrator, summary generation, collapse handler) propagates failures.

### Segment Sentinel Operations
- Always use `segment_helpers.*` functions for sentinel creation/collapse.
- Never construct sentinel Messages manually or parse sentinel metadata directly.

### Memory Surfacing Pipeline (in orchestrator)
Order: subcortical (query expansion + retention + entities + pressure alert) → hard pinned cap (`max_pinned_memories=15`) → sliding fresh budget (`max(min_fresh, max_surfaced - pinned)`) → embedding generation → fresh memory retrieval with dynamic limit → pinned+fresh merge → UpdateTrinketEvent to ProactiveMemoryTrinket. Total surfaced memories bounded by `max_surfaced_memories=20`. Subcortical uses graduated two-tier pressure alerts (warning at `max_pinned-4`, critical at `max_pinned`) derived from `max_pinned_memories`. Linked memories capped at `max_linked_per_primary=2` per primary and rendered as compact `<context>` annotations instead of full XML blocks.

### Context Overflow Remediation (in orchestrator)
Two tiers: (1) embedding-based topic drift pruning (drops old conversation before the drift boundary), (2) oldest-first fallback. Memories are untouched — conversation history is where the token bulk lives. Remediation events logged via standard `logger.info()`.
