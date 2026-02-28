# cns/infrastructure/

Persistence and caching. All database access for CNS flows through this layer.

## Files

- `continuum_repository.py` — PostgreSQL repository for continuums and messages. Per-user DB client caching. Handles segment sentinel creation (`_ensure_active_segment`), message persistence with UPSERT + `::vector(768)` cast, JSONB metadata operations (`metadata->>'key'`, `jsonb_set()`), segment queries, history pagination, full-text search. Admin queries use `AdminSession` to bypass RLS. Session resume: `find_resumable_segment()` finds most-recently-collapsed non-tombstoned segment (when no active exists), `reactivate_collapsed_segment()` atomically clears collapse artifacts and resets to active, `delete_segment_memories()` removes memories extracted from a specific segment. Exports `HistoryResult`, `FailedSegment`, `ActiveSegmentRow` TypedDicts for return type contracts.
- `continuum_pool.py` — Valkey-backed session pool. `get_or_create()` returns cached Continuum or creates new session with segment boundary. `UnitOfWork` accumulates messages + metadata changes, commits atomically. Cache miss = new session. Event-driven invalidation (not TTL). Dead methods `get_by_id()` and `get_session_info()` were removed (zero callers).
- `valkey_message_cache.py` — Message serialization/deserialization for Valkey. Converts Message objects to/from JSON. Event-driven invalidation via segment timeout, no TTL expiry. No constructor DI — uses `get_valkey_client()` singleton directly.
- `feedback_repository.py` — CRUD for `feedback_signals` table. Saves assessment signals (with section_id, strength, evidence), retrieves unsynthesized signals, marks as synthesized. Exports `FeedbackSignalRow` TypedDict for `get_unsynthesized_signals()` return type.
- `feedback_tracker.py` — Synthesis lifecycle for the user model pipeline. Uses modular arithmetic on `users.cumulative_activity_days` (not its own counter) to determine synthesis timing. Stores `activity_days_at_last_synthesis` snapshot, user model XML, and needs_checkin flag. Exports `LoraContent` TypedDict for `get_lora_content()` return type, `TrackingStatus` TypedDict for `get_tracking_status()` return type.

## Patterns to Follow

### Database Access
- Use `PostgresClient("mira_service", user_id=user_id)` for user-scoped queries (RLS enforced).
- Use `AdminSession` from `utils.database_session_manager` for cross-user admin queries only.
- Repository caches per-user DB clients in `_db_cache` dict — don't create new PostgresClient per call.
- JSONB filtering: `metadata->>'key' = %s` for string comparison, `metadata ? 'key'` for existence.
- Atomic metadata updates: `jsonb_set(metadata, '{key}', to_jsonb(%s))` with `execute_returning()`.
- UUIDs stay native for PostgreSQL queries. Convert to string only at serialization boundaries.

### Message Persistence
- Use `save_message()` or `save_messages_batch()` — they handle segment boundary creation automatically via `_ensure_active_segment()`.
- The UPSERT pattern (`ON CONFLICT (id) DO UPDATE`) handles replay safety.
- `_parse_message_rows()` is the canonical DB-row-to-Message converter. Don't reimplement.

### Caching
- Valkey cache is event-driven (no TTL). Invalidated when segment times out.
- Cache miss in `continuum_pool.get_or_create()` triggers `SegmentCacheLoader` to reconstruct from DB.
- `ValkeyMessageCache` handles Message ↔ JSON serialization — use it, don't roll your own.

### Unit of Work
- `continuum_pool.begin_work(continuum)` returns a `UnitOfWork`.
- Call `uow.add_messages(msg1, msg2)` during processing, then `uow.commit()` once at the end.
- Commit persists to both PostgreSQL and Valkey cache atomically.

### User Model Storage
- `FeedbackRepository` for assessment signal CRUD (section_id, strength, evidence), `FeedbackTracker` for synthesis lifecycle.
- Synthesis timing uses modular arithmetic: `cumulative_activity_days % threshold == 0` with dedup guard (`activity_days > activity_days_at_last_synthesis`). No counters to reset — completely stateless.
- `get_tracking_status()` computes `use_days_since_synthesis` on the fly as `cumulative_activity_days - activity_days_at_last_synthesis` to preserve API contract.
- `initialize_user()` uses `ON CONFLICT DO NOTHING` for race-safe account setup.
- `get_lora_content()` returns `{synthesis_xml, needs_checkin}` in a single query for LoraTrinket.
