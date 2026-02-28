# lt_memory/processing/

Memory extraction and post-processing pipeline. Strategy pattern separates business logic (shared) from transport (batch vs immediate). All LLM configs resolved via `get_internal_llm(name)`.

## Files

- `orchestrator.py` — High-level extraction workflow. Two entry points: `submit_segment_extraction(user_id, boundary_message_id)` for on-demand extraction, `extract_unprocessed_segments()` for safety-net sweep (6-hour schedule). Holds a single `ExecutionStrategy` reference. Loads messages from ContinuumRepository, builds chunks, delegates to strategy.
- `execution_strategy.py` — Strategy pattern for extraction execution. `ExecutionStrategy` ABC with `execute_extraction(user_id, chunks) -> str` (always returns batch ID or raises `ValueError`). `BatchExecutionStrategy` submits via `BatchCoordinator.submit_batch()` (no direct Anthropic client). `ImmediateExecutionStrategy` calls `LLMProvider.generate_response()` directly (failover mode). `TierAwareExecutionStrategy` composite resolves billing tier → internal LLM config → routes to batch (Anthropic endpoint) or immediate (non-Anthropic). Shared business logic in `_process_and_store_memories()` and `_persist_llm_entities()`. File-local types: `ExtractionMessage(TypedDict)` for message format, `MemoryContext(TypedDict)` for context snapshots.
- `extraction_engine.py` — Builds extraction payloads: prompt loading from `config/prompts/`, UUID shortening/mapping, memory context retrieval, message formatting (conversation → XML turns), Anthropic message batch building. Produces `ExtractionPayload` consumed by execution strategies. File-local type: `ExtractionMessage(TypedDict)`.
- `memory_processor.py` — Pure data processing: parse LLM JSON responses with `json_repair` fallback, validate structure, remap short→full UUIDs, sanitize memory fields, fuzzy+vector duplicate detection, index remapping, linking pair construction. No side effects. File-local types: `DuplicateCheckResult(NamedTuple)`, `RawMemoryDict(TypedDict)`.
- `batch_coordinator.py` — Generic Anthropic Batch API infrastructure. `BatchCoordinator.submit_batch(requests, batch_type, user_id)` is the **single submission point** for all batches. `poll_extraction_batches()` and `poll_post_processing_batches()` are convenience wrappers over the generic `poll_batches()`. `BatchResultProcessor` ABC for pluggable result handling.
- `consolidation_handler.py` — Memory consolidation with link transfer. Merges multiple memories into one, rewrites outbound links from source memories, archives consolidated originals. Pure business logic, no orchestration.
- `post_processing_orchestrator.py` — Submits relationship classification, consolidation, and verbose refinement batches after extraction completes. Thin coordinator delegating to `RefinementService`, `BatchCoordinator`, `ConsolidationHandler`. Uses `get_internal_llm('relationship')`, `get_internal_llm('consolidation')`, `get_internal_llm('refinement')`.

## Patterns to Follow

### Tier-Aware Execution (Extraction)
Billing tier determines which internal LLM config is used for extraction. The `TierAwareExecutionStrategy` composite:
1. Resolves billing tier via `get_billing_backend().get_billing_tier(user_id)` → `(tier_name, multiplier)`
2. Resolves config via `get_internal_llm(f'extraction_{tier_name}')` — no fallbacks, missing config = `KeyError` = configuration bug
3. Routes: `api.anthropic.com` in `endpoint_url` → `BatchExecutionStrategy`, else → `ImmediateExecutionStrategy`
4. Passes resolved `InternalLLMConfig` explicitly to the concrete strategy — strategies never call `get_internal_llm()` internally

Config naming convention: `{operation}_{billing_tier}` (e.g., `extraction_free`, `extraction_friend`, `extraction_standard`).

### LLM Config Resolution
`InternalLLMConfig` is the single source of truth for all LLM-tuning params (model, endpoint, API key, temperature, max_tokens, thinking_effort). Both batch and immediate paths auto-resolve everything — callers just pass the purpose key.

For `generate_response()` calls (immediate/failover paths), use `internal_llm='purpose'` — it resolves endpoint, model, API key, temperature, max_tokens, and thinking_effort internally:
```python
response = self.llm_provider.generate_response(
    messages=[...],
    internal_llm='relationship',
    allow_negative=True
)
```
For batch API paths, use `build_batch_params(purpose, ...)` from `clients.llm_provider` — all LLM-tuning params resolved from `InternalLLMConfig`:
```python
from clients.llm_provider import build_batch_params
params = build_batch_params(
    'extraction',
    system_prompt=payload.system_prompt,
    messages=payload.messages,
    cache_ttl="1h",  # optional — omit for default 5-minute
)
```
All 4 batch sites (extraction, relationship, consolidation, entity_gc) use this helper. Never construct batch param dicts inline — use `build_batch_params()` to ensure consistent system prompt wrapping, cache_control, and thinking/temperature resolution.

### Execution Strategy Selection
`LTMemoryFactory` creates both concrete strategies and wraps them in `TierAwareExecutionStrategy`. The composite holds references to both `BatchExecutionStrategy` and `ImmediateExecutionStrategy` and routes per-call based on the resolved config's endpoint. The `ExtractionOrchestrator` is unaware of tier routing — it sees a single `ExecutionStrategy`.

### Shared Business Logic
`_process_and_store_memories()` on the base `ExecutionStrategy` is the single path for: response parsing → memory storage with embeddings → entity persistence → extraction-time link creation. Both batch result handling (`batch_result_handlers.py`) and immediate execution use this method. Never duplicate this logic.

### source_segment_id Pipeline
`BatchExecutionStrategy.execute_extraction()` stores `chunk.segment_id` in `ExtractionBatch.chunk_metadata["segment_id"]`. `ExtractionBatchResultHandler.process_result()` reads it back and sets `memory.source_segment_id` on each `ExtractedMemory` before storage. This enables segment-scoped memory cleanup on session resume.

### Billing: Auto-Resolved
The billing hook in `LLMProvider.generate_response()` auto-resolves `pricing_key` from `(model, endpoint_url)`. Strategies do NOT pass `pricing_key` — it matches automatically. Each tiered config needs a corresponding `usage_pricing` row (auto-inserted by `build_config_lookup()` at startup, prices auto-resolved from OpenRouter).
