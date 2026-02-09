# Notification Center Implementation

## Problem
Trinkets in the system prompt feel "cold/received" to MIRA - like reading a biography, not recalling memories.

## Solution
Move dynamic content into a sliding assistant message positioned AFTER the cached conversation but BEFORE the current user message.

## Context Window Structure
```
SYSTEM PROMPT [CACHED]
├── base_prompt
├── domaindoc
├── tool_guidance
└── tool_hints

CONVERSATION HISTORY [PROGRESSIVELY CACHED]
└── ... [cache marker on last assistant msg]

════════════════════════════════════════════════════════════  ← user msg (delimiter)
NOTIFICATION CENTER                                            ← assistant msg
This section moves to the front of your context each turn
to keep important information front-of-mind.
════════════════════════════════════════════════════════════
├── datetime_section
├── conversation_manifest
├── active_reminders
├── context_search_results
└── relevant_memories
════════════════════════════════════════════════════════════

CURRENT USER MESSAGE
```

## Key Design Decisions

### 1. Placement Registry in Base Class
Placement is 100% inside `EventAwareTrinket` - trinkets don't import or declare placement:
```python
# base.py
_NOTIFICATION_CENTER_TRINKETS = frozenset({
    'TimeManager', 'ManifestTrinket',
    'ReminderManager', 'GetContextTrinket', 'ProactiveMemoryTrinket',
})

class EventAwareTrinket:
    @property
    def placement(self) -> TrinketPlacement:
        if self.__class__.__name__ in _NOTIFICATION_CENTER_TRINKETS:
            return TrinketPlacement.NOTIFICATION_CENTER
        return TrinketPlacement.SYSTEM
```

### 2. Single Dict with NamedTuple in Composer
Avoid parallel dicts that can drift:
```python
class SectionData(NamedTuple):
    content: str
    cache_policy: bool
    placement: str

self._sections: Dict[str, SectionData] = {}
```

### 3. Constants Over Magic Strings
```python
PLACEMENT_SYSTEM = "system"
PLACEMENT_NOTIFICATION = "notification"
```

### 4. Delimiter-as-Marker
The user message between conversation and notification center uses the same `═` delimiter, making the role boundary invisible. The delimiter is the opening frame; the notification center content starts with the header. This prevents the model from thinking it spoke the notification center aloud while keeping the infrastructure visually cohesive.

## Files Changed

| File | Change |
|------|--------|
| `working_memory/trinkets/base.py` | `TrinketPlacement` enum + `_NOTIFICATION_CENTER_TRINKETS` registry + `placement` property |
| `working_memory/composer.py` | `SectionData` namedtuple, routing by `section.placement`, `_build_notification_center()` |
| `cns/core/events.py` | `placement` field on `TrinketContentEvent`, `notification_center` field on `SystemPromptComposedEvent` |
| `working_memory/core.py` | Pass `event.placement` to `composer.add_section()` |
| `cns/services/orchestrator.py` | Inject delimiter + notification center after conversation history |
| `cns/integration/factory.py` | Remove `PunchclockTrinket` |
| `working_memory/trinkets/punchclock_trinket.py` | **DELETED** |

## Data Flow
```
EventAwareTrinket.placement (property)
    ↓ .value
TrinketContentEvent.placement ("notification")
    ↓
composer.add_section(placement=...)
    ↓
SectionData.placement
    ↓
compose() routes to notification_parts
    ↓
_build_notification_center() formats with ═══ delimiters
    ↓
SystemPromptComposedEvent.notification_center
    ↓
orchestrator injects after [context refresh]
```

## Cache Behavior (Empirically Verified)
- System prompt: cached
- Conversation history: progressively cached
- Notification center: NOT cached (changes each turn without invalidating conversation cache)
- `[context refresh]` marker: NOT cached
