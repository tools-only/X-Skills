# cns/core/

Immutable domain model. No database access, no HTTP, no external dependencies.
Downstream layers (services, infrastructure) consume these types.

## Files

- `continuum.py` — Continuum aggregate root. Holds state + message cache. `get_messages_for_api()` handles all API formatting (timestamps, multimodal, tool calls, cache control). Messages loaded externally via `apply_cache()`.
- `message.py` — Frozen Message dataclass. Value object with UUID identity, UTC timestamps, role validation. `to_db_tuple()` for DB insertion, `with_metadata()` for immutable updates. Exports `MessageMetadata` TypedDict (all known metadata keys, `total=False`), `ContentBlock` union type (`TextBlock | ImageBlock | DocumentBlock | ContainerUploadBlock`).
- `state.py` — Minimal ContinuumState (id, user_id, metadata). Exports `ContinuumStateDict` TypedDict for `to_dict()`/`from_dict()` round-trips. Intentionally lightweight — real state lives in the message cache.
- `events.py` — Domain event hierarchy. Four base categories: MessageEvent, ToolEvent, WorkingMemoryEvent, ContinuumCheckpointEvent. All frozen, all have `.create()` classmethods. `TurnCompletedEvent.continuum` is typed as `Continuum` via `TYPE_CHECKING` guard (circular import with continuum.py).
- `stream_events.py` — Mutable dataclasses for LLM streaming: TextEvent, ThinkingEvent, ToolDetectedEvent/Executing/Completed/Error, CompleteEvent, ErrorEvent, CircuitBreakerEvent, RetryEvent, FileArtifactEvent (code execution file downloads). Type-discriminated via `type` field. `CompleteEvent.response` is typed as `anthropic.types.Message`.
- `segment_cache_loader.py` — Reconstructs message history on cache miss. Loads collapsed summaries (complexity-scored selection), continuity messages, active segment. Uses `segment_helpers` for markers.

## Patterns to Follow

- **New events**: Subclass one of the four base categories. Add a `.create()` classmethod that auto-generates event_id/occurred_at and pulls user_id from contextvar. Use `frozen=True, kw_only=True`.
- **New stream events**: Subclass `StreamEvent`, set a unique `type` string default. These are mutable (streaming perf).
- **Message creation**: Always use `Message(content=..., role=...)` — id and created_at auto-generate. Never construct raw dicts when a Message would do.
- **Serialization**: Use `to_dict()` / `from_dict()` round-trip methods. UUIDs become strings, datetimes become ISO format at the boundary.
- **Metadata keys**: All known keys are documented in `MessageMetadata` TypedDict in `message.py`. Common keys: `is_segment_boundary`, `status`, `segment_id`, `system_notification`, `has_tool_calls`, `tool_calls`, `tool_call_id`, `complexity_score`, `embedding_value`.
