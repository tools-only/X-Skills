# Pydantic Runtime Removal Attempt (PR #334, closed)

Date: 2026-01-31

Context
- We attempted to decouple runtime usage of pydantic-ai by normalizing messages to canonical wire dicts and removing rehydration.
- The PR was closed after review due to context loss in agent history.

What we did (record only)
- Normalized conversation messages to canonical wire dicts in memory and persistence (StateManager serialization/deserialization, RequestOrchestrator message persistence).
- Removed pydantic message rehydration on load and eliminated the message_handler shim.
- Removed truncation checker from the response flow.
- Switched UI/history filtering to kind/part accessors (_get_attr, _get_parts) instead of pydantic-ai classes.
- Switched tool retry flow to ToolRetryError (no ModelRetry conversion).
- Reduced agent cache to session-only and fixed request_delay default source.
- Added pydantic usage ratchet test + baseline update to prevent creep.

Outcome observed
- pydantic-ai’s agent loop expects ModelRequest/ModelResponse objects; when history is wire dicts, `_clean_message_history` drops prior turns, losing context.
- This showed the runtime removal cannot be staged in multiple steps: message history, adapters, and the agent loop need to change together.
