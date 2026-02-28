# working_memory/

Event-driven system prompt composition. Coordinates trinkets (modular content providers) through CNS events to build the final system prompt delivered to the LLM.

## Files

- `types.py` — TypedDicts for structures produced by composer and core: `ComposedPrompt`, `TrinketState`, `AllTrinketStates`
- `composer.py` — `SystemPromptComposer`: collects sections from trinkets, routes by `SECTION_LAYOUT` (system/post_history/notification/conversation_prefix), returns `ComposedPrompt`
- `core.py` — `WorkingMemory`: event-driven coordinator. Manages trinket registration, routes `UpdateTrinketEvent` to trinkets, collects `TrinketContentEvent`, composes and publishes `SystemPromptComposedEvent`

## Subdirectory: trinkets/

- `base.py` — `EventAwareTrinket` (ABC): base class with `@abstractmethod generate_content()`. `StatefulTrinket`: adds turn tracking with `@abstractmethod _expire_items()` and `_clear_all_state()`
- `domaindoc_trinket.py` — Injects enabled domain knowledge documents with section collapse/expand state
- `lora_trinket.py` — Injects user model observations from feedback synthesis pipeline
- `manifest_trinket.py` — Displays conversation segment manifest in notification center
- `peanutgallery_trinket.py` — Displays metacognitive guidance with TTL expiry (StatefulTrinket)
- `proactive_memory_trinket.py` — Displays surfaced long-term memories in notification center
- `reminder_manager.py` — Fetches and formats active reminders for notification center
- `getcontext_trinket.py` — Displays async context search results (StatefulTrinket)
- `time_manager.py` — Injects current datetime into notification center
- `HOW_TO_BUILD_A_TRINKET.md` — Guide for creating new trinkets

## Patterns to Follow

### Trinket Base Classes
All trinkets extend `EventAwareTrinket` (ABC). Must implement `generate_content(context: Dict[str, Any]) -> str`. Trinkets with turn-scoped state extend `StatefulTrinket` and implement `_expire_items()` and `_clear_all_state()`.

### Event Flow
`ComposeSystemPromptEvent` → WorkingMemory broadcasts `UpdateTrinketEvent` to all trinkets → each trinket calls `generate_content()` → publishes `TrinketContentEvent` → composer assembles `ComposedPrompt` → `SystemPromptComposedEvent`

### Type Contracts
- Composer returns `ComposedPrompt` (TypedDict from `types.py`)
- `get_trinket_state()` returns `TrinketState | None`
- `get_all_trinket_states()` returns `AllTrinketStates`
- Trinkets are stored as `Dict[str, EventAwareTrinket]`, not `object`
- Event handler parameters are typed via `TYPE_CHECKING` imports

### Infrastructure Failures
Trinket failures are isolated in `core.py:_handle_update_trinket()` — individual trinket exceptions don't crash the prompt composition. But trinkets themselves must NOT catch infrastructure failures (no try/except around DB/Valkey calls) — let them propagate to the isolation layer.
