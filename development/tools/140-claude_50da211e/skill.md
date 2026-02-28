# cns/api/

FastAPI endpoints. Thin routing layer — validates input, sets user context, delegates to services.

## Files

- `chat.py` — `POST /chat`. Non-streaming chat. Handles image compression (two tiers: 1200px inference, 512px storage), document processing (PDF/DOCX/XLSX/TXT/CSV/JSON), per-user distributed request lock (60s TTL). Supports `include_thinking` flag to include thinking trace in response. Delegates to `orchestrator.process_message()`.
- `websocket_chat.py` — `WebSocket /ws/chat`. Real-time bidirectional chat with custom auth protocol (first message = auth token). Queue-based streaming: orchestrator pushes events via callback, async loop consumes and sends. Thinking events gated behind per-message `include_thinking` flag (default false). Thread context propagation via `run_in_threadpool()`.
- `actions.py` — `POST /actions`. Domain-routed state mutations (~2200 lines). Domains: REMINDER, MEMORY, CONTACTS, USER, DOMAIN_KNOWLEDGE, CONTINUUM, LORA. Schema-based validation per action. Also `GET /tools/{tool_name}/query` for direct tool polling.
- `data.py` — `GET /data?type=...`. Unified read endpoint. Types: HISTORY, MEMORIES, DASHBOARD, USER, DOMAINDOCS, WORKING_MEMORY, LORA. Pagination via offset/limit.
- `health.py` — `GET /health`, `/health/threads`, `/health/thread-dump`. No auth required.
- `tool_config.py` — CRUD for per-user tool configuration. Lists configurable tools, gets/sets config via `UserCredentialService`, validates against Pydantic schemas.
- `update.py` — `GET /check_update`. Version comparison against `/VERSION` file. No auth required.
- `demo.py` — `POST /demo/session`, `POST /demo/chat`. Ephemeral demo sessions in Valkey (15 min TTL, 5 msg/min rate limit). Supports account conversion flow.
- `files.py` — `GET /files/{file_id}`. Serves code execution file artifacts from `data/users/{user_id}/tmp/`. Security: auth gate, file_id regex validation, path traversal guard, forced-download headers (`Content-Disposition: attachment`, `application/octet-stream`).
- `federation.py` — `POST /federation/deliver`. Server-to-server webhook for Lattice federation. Receives inbound federated messages and delivers to local users via `PagerTool`. No user auth (Lattice handles verification).
- `base.py` — Response types (`SuccessResponse`, `ErrorResponse` frozen dataclasses), `ErrorDetail`/`ResponseMeta` TypedDicts, `APIResponse = SuccessResponse | ErrorResponse` union alias, error hierarchy (`APIError`, `ValidationError`, `NotFoundError`, `ServiceUnavailableError`), `BaseHandler` template method pattern.

## Patterns to Follow

### Endpoint Structure
- Protected routes use `Depends(get_current_user)` from `auth.api`.
- Set user context immediately: `set_current_user_id(user_id)` at the top of the handler.
- Use `BaseHandler` for request processing with automatic error handling and request IDs.
- Return responses via `create_success_response()` / `create_error_response()`.

### Response Types
`base.py` defines two frozen dataclasses for API responses:
- `SuccessResponse(data: dict[str, Any], meta: ResponseMeta)` — `.to_dict()` hardcodes `"success": True`
- `ErrorResponse(error: ErrorDetail, meta: ResponseMeta)` — `.to_dict()` hardcodes `"success": False`
- `APIResponse = SuccessResponse | ErrorResponse` — union alias (cannot be used with `isinstance`; use the concrete types)
- `ErrorDetail(TypedDict)` — `{code: str, message: str, details: dict[str, Any]}`
- `ResponseMeta(TypedDict, total=False)` — `{timestamp: str, request_id: str}`
- Both are frozen; use `dataclasses.replace()` for mutation (see `auth/api.py` for pattern)
- Construct via `create_success_response(data, meta)` / `create_error_response(exception, request_id)`

Wire format (unchanged):
```json
{"success": true, "data": {...}, "meta": {"request_id": "..."}}
{"success": false, "error": {"code": "...", "message": "...", "details": {...}}, "meta": {"request_id": "..."}}
```

### Actions Domain Routing
- Each domain has a handler class extending `BaseDomainHandler`.
- Actions define validation schemas: `required`, `optional`, `types` fields.
- Domain handlers validate, then delegate to tools or repositories.
- New domains: add handler class, add to domain router in `ActionsEndpoint`.

### Concurrency
- `UserRequestLock` ensures one active request per user (distributed via Valkey, 60s TTL).
- WebSocket uses `run_in_threadpool()` with context propagation for sync orchestrator code.

### Error Handling
- All errors extend `APIError(Exception)` with `message: str`, `code: str`, `details: dict[str, Any]`
- `create_error_response(exception, request_id)` maps any `APIError` subclass to `ErrorResponse` with the correct `ErrorDetail`
- HTTP 400: `ValidationError` (bad input, schema violations)
- HTTP 402: `InsufficientBalanceError` (billing)
- HTTP 404: `NotFoundError` (resource not found)
- HTTP 500: Unhandled exceptions (logged, generic message to client)
- HTTP 503: `ServiceUnavailableError` (infrastructure down)
- Infrastructure failures propagate (no catch-and-return-None). Only `except` for adding context before re-raising.
