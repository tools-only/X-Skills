# Research – API Error Logging Gaps

**Date:** 2026-01-28
**Owner:** claude-opus
**Phase:** Research

## Goal

Understand why the current debug logs for API/connection errors are insufficient for diagnosing failures, and identify what context is available but not being captured.

## Trigger

User experienced a `ModelAPIError: Connection error` at `2026-01-28T00:23:39.729952+00:00` during a streaming session. The debug logs only showed:

```
[DEBUG  ] [LIFECYCLE] Stream failed: ModelAPIError: Connection error.
```

This provides no actionable information about:
- Which API endpoint failed
- What HTTP status code (if any)
- What model/provider was being used
- Which retry attempt failed
- Request ID for correlation

## Findings

### The Core Problem

**No `ModelAPIError` class exists in tunacode.** The exception comes from pydantic-ai's internal handling. When it reaches tunacode, it's caught and logged with only:

```python
# src/tunacode/core/agents/agent_components/streaming.py:381-384
except Exception as e:
    logger.lifecycle(f"Stream failed: {type(e).__name__}: {e}")
    raise
```

### What IS Captured

| Information | Available | Logged on Error |
|------------|-----------|-----------------|
| Exception type | ✅ | ✅ |
| Exception message | ✅ | ✅ |

### What is NOT Captured (But Available)

| Information | Available At | Not Logged |
|------------|--------------|------------|
| HTTP status code | `HTTPStatusError.response.status_code` | ❌ |
| Request URL | `HTTPStatusError.request.url` | ❌ |
| Model name | `self.model` in main.py:153 | ❌ |
| Provider name | Parsed from model string | ❌ |
| Request ID | `ctx.request_id` in main.py:243 | ❌ |
| Iteration number | `iteration_index` in streaming.py:128 | ❌ |
| Retry attempt | Tenacity internal state | ❌ |
| Response headers | `HTTPStatusError.response.headers` | ❌ |
| Timeout threshold | `timeout` parameter | ❌ |
| Actual duration | Available via timing | ❌ |
| Base URL | Provider config | ❌ |

### Key Files & Line Numbers

**Exception Definitions:**
- `src/tunacode/exceptions.py:81-314` - Custom exceptions (no ModelAPIError)
- `GlobalRequestTimeoutError`: lines 268-278

**Stream Error Handler:**
- `src/tunacode/core/agents/agent_components/streaming.py:381-384` - Minimal logging

**Request Context (has request_id but unused on error):**
- `src/tunacode/core/agents/main.py:172-181` - Creates request_id
- `src/tunacode/core/agents/main.py:243-244` - Logs request_id only on SUCCESS

**HTTP Retry Config:**
- `src/tunacode/core/agents/agent_components/agent_config.py:385-393` - Tenacity transport

**Logging Manager:**
- `src/tunacode/core/logging/manager.py:119-124` - `lifecycle()` only works in debug_mode

### Critical Architecture Gaps

1. **lifecycle() is debug-only**: Stream failures use `logger.lifecycle()` which is gated by `debug_mode`. Production errors get NO log entry.

2. **No HTTPStatusError enrichment**: The retry transport catches `HTTPStatusError` for retries, but when retries are exhausted, the exception bubbles up with no context extraction.

3. **Request ID only logged on success**: `main.py:243-244` logs the request ID, but only in the success path. Error paths never log it.

4. **Silent exception handlers**: `streaming.py` has 8+ `except Exception: pass` blocks for instrumentation that could capture diagnostic info but discard it instead.

5. **No exception wrapping for API errors**: Tool errors get wrapped in `ToolExecutionError` with context. API errors bubble as raw pydantic-ai/httpx exceptions.

## Key Patterns / Solutions Found

### Pattern: Error Context at Each Layer

The codebase follows "fail fast, fail loud" but doesn't enrich errors as they propagate:

```
HTTP Layer (httpx) → pydantic-ai → streaming.py → main.py → UI
     ↑                                    ↑
     Has full context              Logs only type+message
```

### Pattern: Attribute-Based Error Rendering

`ui/renderers/errors.py` extracts context via `getattr()`:

```python
EXCEPTION_CONTEXT_ATTRS = {
    "tool_name": "Tool",
    "path": "Path",
    "model": "Model",  # ← exists but never set on API errors
    ...
}
```

If API errors carried `model`, `status_code`, `request_id` attributes, the renderer would display them automatically.

## Recommended Fix

### Immediate: Enrich Exception Before Logging

At `streaming.py:381`:

```python
except Exception as e:
    # Extract context before logging
    context = {
        "request_id": request_id,
        "iteration": iteration_index,
        "model": model_name,
    }
    if hasattr(e, 'response'):
        context["status_code"] = e.response.status_code
        context["url"] = str(e.request.url)

    logger.error(f"Stream failed: {type(e).__name__}: {e}", extra=context)
    raise
```

### Medium-term: Create ServiceAPIError

Add to `exceptions.py`:

```python
class ServiceAPIError(ServiceError):
    def __init__(
        self,
        message: str,
        *,
        provider: str | None = None,
        model: str | None = None,
        status_code: int | None = None,
        request_id: str | None = None,
        suggested_fix: str | None = None,
    ):
        self.provider = provider
        self.model = model
        self.status_code = status_code
        self.request_id = request_id
        ...
```

Then wrap at the boundary where pydantic-ai errors surface.

## Knowledge Gaps

- What does pydantic-ai's `ModelAPIError` actually contain? Need to inspect the exception object.
- Does tenacity expose retry attempt number in a callback hook?
- Should error logging be unconditional (not gated by debug_mode)?

## References

- `src/tunacode/core/agents/agent_components/streaming.py`
- `src/tunacode/core/agents/main.py`
- `src/tunacode/core/agents/agent_components/agent_config.py`
- `src/tunacode/exceptions.py`
- `src/tunacode/ui/renderers/errors.py`
- `src/tunacode/core/logging/manager.py`
