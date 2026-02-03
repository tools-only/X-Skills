# Research – Pydantic-AI to Alchemy-RS Migration

**Date:** 2026-01-27
**Owner:** Claude (Opus 4.5)
**Phase:** Research

## Goal

Map out what needs to change to migrate from pydantic-ai to alchemy-rs for LLM API calls, while keeping Python for state management and agent loop.

## Architecture Summary

**alchemy-rs** provides:
- Unified LLM API abstraction (8+ providers)
- Streaming-first async design
- Tool/function calling support
- Cross-provider message transformation

**Current tunacode** uses pydantic-ai for:
1. **LLM API calls** (what alchemy-rs will replace)
2. **Message types** (ModelRequest, ModelResponse, parts)
3. **Agent iteration loop** (agent.iter() streaming)
4. **Tool registration & execution** (Tool class, ModelRetry)

---

## Migration Difficulty Ranking

### Tier 1: EASY - Direct Replacements

| Component | Files | What Changes | Difficulty |
|-----------|-------|--------------|------------|
| **LLM API Client** | `agent_config.py` | Replace `AnthropicModel`/`OpenAIChatModel` with alchemy-rs `get_model()` | Easy |
| **HTTP Transport** | `agent_config.py:381-404` | alchemy-rs handles retries internally | Easy |
| **Model Factory** | `agent_config.py:267-316` | Single `get_model(provider, model_name)` call | Easy |

**Why Easy:** alchemy-rs is a drop-in replacement for the model layer. Same inputs (model name, API key), same outputs (streaming tokens).

---

### Tier 2: MEDIUM - Need Adapter Layer

| Component | Files | What Changes | Difficulty |
|-----------|-------|--------------|------------|
| **Message Types** | `types/pydantic_ai.py`, `types/canonical.py` | New wire format from Rust, need Python wrapper types | Medium |
| **Streaming Protocol** | `streaming.py:122-384` | Replace `PartDeltaEvent` with alchemy-rs stream events | Medium |
| **Tool Call Format** | `tool_dispatcher.py:219-352` | Parse alchemy-rs tool call format | Medium |
| **ModelRetry Exception** | 8 tool files | Replace pydantic-ai's `ModelRetry` with custom exception | Medium |

**Why Medium:** The message format will change. You already have a canonical message layer (`types/canonical.py`, `utils/messaging/adapter.py`) that isolates format changes. Extend this pattern for alchemy-rs.

---

### Tier 3: KEEP IN PYTHON - No Migration Needed

| Component | Files | Reason to Keep |
|-----------|-------|----------------|
| **StateManager** | `core/state.py` | Pure Python state, no LLM dependency |
| **SessionState** | `core/types/state_structures.py` | Dataclasses, Python-native |
| **ToolCallRegistry** | `core/types/tool_registry.py` | In-memory tracking, Python dict |
| **Session Persistence** | `state.py:247-353` | JSON serialization, filesystem I/O |
| **Message Sanitization** | `resume/sanitize.py` | Logic is format-agnostic via canonical types |
| **Token Counting** | `state.py:63-81` | `tiktoken` library, pure Python |

**Why Keep:** State management has zero coupling to pydantic-ai's LLM features. It uses pydantic-ai message types for storage, but those are just dataclasses.

---

### Tier 4: KEEP - Agent Loop (Python Orchestration)

| Component | Files | Reason to Keep |
|-----------|-------|----------------|
| **RequestOrchestrator** | `main.py:216-465` | Orchestration logic, timeout handling |
| **Node Iteration** | `main.py:346-425` | Control flow stays in Python |
| **Tool Executor** | `tool_executor.py:54-141` | Python async, parallel execution |
| **Response State Machine** | `response_state.py` | Python state machine |
| **Cleanup Loop** | `sanitize.py:459-503` | Python iteration |

**Why Keep:** The agent loop is control flow, not LLM protocol. Replace `agent.iter()` call with alchemy-rs streaming, keep the surrounding Python orchestration.

---

## Concrete Changes Required

### Phase 1: Create Alchemy-RS Bridge (Python bindings)

```
src/tunacode/core/llm/
├── __init__.py
├── alchemy_client.py      # PyO3 bindings to alchemy-rs
├── message_types.py       # Python wrappers for Rust types
└── streaming.py           # Async stream adapter
```

**What this does:**
1. Compile alchemy-rs as Python extension (PyO3/maturin)
2. Expose `get_model()`, `stream()`, `Context` to Python
3. Convert Rust stream events to Python async generator

### Phase 2: Replace Model Layer

**File:** `src/tunacode/core/agents/agent_components/agent_config.py`

**Before:**
```python
from pydantic_ai import Agent
from pydantic_ai.providers.anthropic import AnthropicModel, AnthropicProvider
from pydantic_ai.providers.openai import OpenAIChatModel, OpenAIProvider

def create_model_instance(...):
    if is_anthropic:
        return AnthropicModel(model, provider=AnthropicProvider(...))
    else:
        return OpenAIChatModel(model, provider=OpenAIProvider(...))
```

**After:**
```python
from tunacode.core.llm.alchemy_client import get_model

def create_model_instance(...):
    return get_model(provider, model_name, api_key=api_key, base_url=base_url)
```

### Phase 3: Update Streaming

**File:** `src/tunacode/core/agents/agent_components/streaming.py`

**Before:**
```python
async with node.stream(agent_run_ctx) as request_stream:
    async for event in request_stream:
        if isinstance(event, PartDeltaEvent):
            delta_text = event.delta.content_delta
```

**After:**
```python
from tunacode.core.llm.streaming import stream_completion

async for chunk in stream_completion(model, context):
    delta_text = chunk.content
```

### Phase 4: Update Message Adapter

**File:** `src/tunacode/utils/messaging/adapter.py`

Add conversion functions:
```python
def alchemy_to_canonical(msg: AlchemyMessage) -> CanonicalMessage:
    ...

def canonical_to_alchemy(msg: CanonicalMessage) -> AlchemyMessage:
    ...
```

### Phase 5: Replace ModelRetry

**Before (8 files):**
```python
from pydantic_ai.exceptions import ModelRetry
raise ModelRetry("error message")
```

**After:**
```python
from tunacode.exceptions import ToolRetryError
raise ToolRetryError("error message")
```

---

## Files Affected by Migration

### MUST CHANGE (11 files)

| File | Changes |
|------|---------|
| `types/pydantic_ai.py` | Remove, replace with alchemy types |
| `agent_config.py` | New model factory |
| `streaming.py` | New streaming protocol |
| `tool_dispatcher.py` | Parse alchemy tool calls |
| `message_handler.py` | Remove pydantic_ai import |
| `bash.py` | Replace ModelRetry |
| `web_fetch.py` | Replace ModelRetry |
| `grep.py` | Replace ModelRetry |
| `update_file.py` | Replace ModelRetry |
| `list_dir.py` | Replace ModelRetry |
| `glob.py` | Replace ModelRetry |
| `write_file.py` | Replace ModelRetry |
| `decorators.py` | Replace ModelRetry |

### MIGHT CHANGE (Adapter Layer)

| File | Changes |
|------|---------|
| `utils/messaging/adapter.py` | Add alchemy conversion |
| `types/canonical.py` | Possibly extend for alchemy parts |
| `resume/sanitize.py` | Only if alchemy format differs significantly |

### NO CHANGES NEEDED (16+ files)

| File | Reason |
|------|--------|
| `core/state.py` | Pure Python state |
| `core/types/state_structures.py` | Dataclasses |
| `core/types/tool_registry.py` | In-memory tracking |
| `main.py` | Keep orchestration, just call alchemy instead |
| `tool_executor.py` | Python async stays |
| `response_state.py` | State machine |
| `configuration/` | Config loading |
| `ui/` | Unchanged |
| `tests/` | Update mocks only |

---

## Estimated Complexity

| Phase | Effort | Dependencies |
|-------|--------|--------------|
| **Phase 1: PyO3 Bridge** | 2-3 days | Rust build toolchain, maturin |
| **Phase 2: Model Layer** | 1 day | Phase 1 complete |
| **Phase 3: Streaming** | 1-2 days | Phase 2 complete |
| **Phase 4: Message Adapter** | 1 day | Phase 3 complete |
| **Phase 5: ModelRetry** | 2-3 hours | Independent |

**Total: ~1 week of focused work**

---

## Key Patterns That Help

1. **Canonical Message Layer** - Already isolates format changes. Just add alchemy converters.

2. **Protocol-Based Types** - `StateManagerProtocol` means state code doesn't import concrete LLM types.

3. **Wrapper Pattern** - `types/pydantic_ai.py` provides single import point. Replace with alchemy types.

4. **Callback Architecture** - Streaming and tool callbacks decouple UI from agent internals.

---

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| PyO3 async compatibility | Use `pyo3-asyncio` crate, proven pattern |
| Message format mismatch | Canonical layer absorbs differences |
| Tool call format changes | Tool dispatcher is already modular |
| Serialization changes | Session persistence uses dicts, not pydantic-ai types directly |

---

## Recommendation

**Yes, this migration is achievable and relatively easy because:**

1. **Clear boundaries** - LLM API is isolated in `agent_config.py` and `streaming.py`
2. **Adapter pattern exists** - `utils/messaging/adapter.py` already abstracts message formats
3. **State is decoupled** - StateManager has no LLM dependencies beyond type annotations
4. **Agent loop is Python** - Only replace the `agent.iter()` call, keep the orchestration

**Suggested approach:**

1. Build alchemy-rs Python bindings first (PyO3)
2. Run both pydantic-ai and alchemy-rs in parallel during transition
3. Switch one model at a time (start with Anthropic)
4. Verify message round-trips through canonical layer
5. Remove pydantic-ai after all providers migrated

---

## References

- **alchemy-rs repo:** https://github.com/tunahorse/alchemy-rs
- **pydantic-ai usage:** `/home/tuna/tunacode/src/tunacode/types/pydantic_ai.py`
- **Message adapter:** `/home/tuna/tunacode/src/tunacode/utils/messaging/adapter.py`
- **Agent config:** `/home/tuna/tunacode/src/tunacode/core/agents/agent_components/agent_config.py`
- **Streaming:** `/home/tuna/tunacode/src/tunacode/core/agents/agent_components/streaming.py`
- **State manager:** `/home/tuna/tunacode/src/tunacode/core/state.py`
