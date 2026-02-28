# lt_memory/

Long-term memory system. Batch extraction from conversation segments, entity linking, consolidation, and proactive retrieval. All services created and managed by `LTMemoryFactory`.

## Files

- `models.py` — Pydantic models and shared TypedDicts. Central type registry for all cross-file structures. Defines `Memory`, `ExtractedMemory`, `MemoryLink`, `ExtractionBatch`, `PostProcessingBatch`, `ProcessingChunk`, `ConsolidationCluster`, `PendingManualMemory`. Literal aliases: `RelationshipType`, `BatchStatus`, `BatchKind`. Shared TypedDicts (18): see TypedDict Registry below.
- `db_access.py` — Database gateway. All SQL operations, memory CRUD, batch tracking, entity operations, scoring formula. Batch operations use `BatchKind` discriminator (`"extraction"` | `"post_processing"`) with 6 generic methods + 2 create methods. `get_all_memories()` returns `List[Memory]`; `get_memories_paginated()` returns `MemoryPageResult`. Uses `LTMemorySessionManager` for RLS-scoped connections.
- `factory.py` — Creates and manages all service instances via `LTMemoryFactory`. Layered initialization with reverse-order cleanup. Singleton via `get_lt_memory_factory()`.
- `vector_ops.py` — Embedding operations: store memories with embeddings, find similar memories, hybrid search integration. Uses `HybridEmbeddingsProvider` (768d).
- `hybrid_search.py` — BM25 + vector hybrid search with reciprocal rank fusion. `HybridSearcher.hybrid_search()` returns `List[Memory]`.
- `linking.py` — Relationship classification between memories. `build_classification_payload() -> ClassificationPayload`, `traverse_related() -> List[TraversalResult]`, bidirectional link creation. Dead code: `classify_relationship_sync()` and `_parse_classification_response()` have no active callers.
- `proactive.py` — Proactive memory retrieval for conversation context. Merges similarity-based and hub-based memory pools, reranks with link traversal. Returns `List[MemoryDict]`.
- `hub_discovery.py` — Entity-driven memory retrieval. Matches extracted entities to DB entities via pg_trgm fuzzy matching, collects linked memories, ranks by expansion embedding similarity. Fail-fast: DB errors propagate (no broad exception catching).
- `refinement.py` — Consolidation cluster identification and payload building. `build_consolidation_payload() -> ConsolidationPayload`.
- `entity_extraction.py` — spaCy NER entity extraction. `extract_entities_with_types() -> List[NamedEntity]`.
- `entity_gc.py` — Entity garbage collection: BFS connected-components grouping, LLM review (batch or synchronous), merge/delete/keep execution. File-local types: `EntityNode`, `EntityGCDecision`, `ReviewPrompt`.
- `memory_formatter.py` — XML formatting for memory output in prompts. Uses `AnnotationEntry` TypedDict.
- `batch_result_handlers.py` — Result handlers for extraction, relationship, consolidation, and entity GC batch types. Implements `BatchResultProcessor` ABC. `PostProcessingBatchDispatcher` routes by batch type.
- `__init__.py` — Public exports for all shared types and models.
- `processing/` — Extraction and post-processing pipeline (see `processing/CLAUDE.md`).

## TypedDict Registry (models.py)

Shared types that cross file boundaries. File-local types stay in their files.

| TypedDict | Key consumers |
|-----------|--------------|
| `MemoryLinkEntry` | Memory fields, db_access, memory_formatter, linking, consolidation_handler |
| `EntityLinkEntry` | Memory.entity_links |
| `AnnotationEntry` | Memory.annotations, memory_formatter, consolidation_handler |
| `LinkMetadata` | Memory.link_metadata |
| `TraversalResult` | linking.traverse_related; memory_tool.py, peanutgallery_model.py |
| `ClassificationPayload` | linking.build_classification_payload |
| `ClassificationResult` | linking._parse_classification_response |
| `ClassificationPair` | execution_strategy, batch_result_handlers |
| `LinkingPair` | execution_strategy, batch_result_handlers |
| `EntityPairRow` | db_access, entity_gc |
| `GCStats` | entity_gc (3 sites) |
| `UserMemorySettings` | db_access, orchestrator |
| `MemoryPageResult` | db_access.get_memories_paginated; cns/api/data.py |
| `NamedEntity` | entity_extraction, entity linking |
| `MemoryContext` | extraction_engine, memory_processor |
| `MemoryContextSnapshot` | orchestrator, ProcessingChunk |
| `ChunkMetadata` | ExtractionBatch |
| `MemoryDict` | proactive; peanutgallery_service, memory_relevance_service |
| `ConsolidationPayload` | refinement, post_processing_orchestrator |

## Patterns to Follow

### Strong Typing
- `RelationshipType` Literal replaces runtime `@field_validator` on `MemoryLink.link_type`
- `BatchStatus` Literal replaces runtime validators on `ExtractionBatch.status` and `PostProcessingBatch.status`
- `BatchKind` Literal (`"extraction"` | `"post_processing"`) discriminates batch DB operations — 6 generic methods parameterized by kind replace 12 duplicated methods
- Use TypedDicts from `models.py` for well-defined dict structures instead of `Dict[str, Any]`
- File-local types (used in only one file) are defined in that file, not `models.py`
- `ProcessingChunk.messages` stays `List[Any]` (Pydantic can't resolve TYPE_CHECKING-only imports at runtime)

### Batch Submission Contract
- `BatchCoordinator.submit_batch()` is the **only** path for submitting batches to Anthropic — no direct `anthropic_client.beta.messages.batches.create()` calls elsewhere
- `BatchCoordinator.submit_batch()` takes `(requests, batch_type, user_id)` — no `input_data` param
- Callers store `input_data` in their own batch records via `create_extraction_batch()` or `create_post_processing_batch()`

### Execution Contract
- `ExecutionStrategy.execute_extraction()` returns `str` (always) or raises `ValueError` — never `None`

### Fail-Fast
- Hub discovery: DB failures propagate from `_match_entities()` — no broad exception catching
- `update_access_stats()`: raises `RuntimeError` if memory disappears between UPDATE and SELECT
- Required infrastructure failures propagate; only catch exceptions for genuinely optional features
