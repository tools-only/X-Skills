# How to Build a Tool

*Technical guide based on successful tool development patterns*

## 🚀 START HERE

Tools in MIRA follow consistent patterns. This guide shows you where to find each pattern in existing, battle-tested implementations. Use the **Pattern Index** below to jump directly to the code that demonstrates what you need.

**Learning approach:**
1. Scan the Pattern Index to see what's available
2. Find the pattern you need in an existing tool
3. Read that specific implementation with line numbers provided
4. Copy the pattern and adapt it to your use case

## Core Concepts

Effective tool building leverages two different types of pattern recognition:
- **Human perspective**: Real-world usage patterns, social dynamics, UX intuition
- **AI perspective**: System architecture, technical constraints, code consistency

Neither perspective alone produces optimal tools. The best solutions emerge from the intersection.

## 📋 Pattern Index

### Essential Patterns

Every tool needs these core patterns:

| Pattern | Where to Find | What It Shows |
|---------|---------------|---------------|
| **Tool Base Class** | tools/repo.py:40-143 | Base Tool class, properties, abstract methods |
| **Tool Metadata** | reminder_tool.py:38-112 | name, simple_description, anthropic_schema |
| **Configuration** | reminder_tool.py:29-35 | Pydantic config with registry.register() |
| | contacts_tool.py:28-41 | Config with Field descriptions |
| **Deferred Initialization** | reminder_tool.py:114-145 | has_user_context() check, table creation |
| **Database Schema** | reminder_tool.py:124-145 | Tables with encrypted__ fields, indexes |
| | contacts_tool.py:264-279 | Schema with UUID primary keys |
| **Operation Routing** | reminder_tool.py:179-252 | JSON kwargs parsing, operation dispatch |
| | contacts_tool.py:173-220 | Clean operation routing pattern |
| **Input Validation** | contacts_tool.py:242-258 | Validation with clear error messages |
| | reminder_tool.py:282-294 | Required field validation |
| **CRUD Operations** | contacts_tool.py:222-298 | Complete CRUD with UUID generation |
| | reminder_tool.py:254-354 | CRUD with cross-tool linking |
| **Timezone Handling** | reminder_tool.py:788-852 | Natural language parsing, UTC storage |
| | reminder_tool.py:335-338 | Display conversion to user timezone |
| **Encryption** | contacts_tool.py:264-279 | encrypted__ prefix for automatic encryption |
| | userdata_manager.py:298-330 | How encryption/decryption works internally |
| **Fuzzy Name Matching** | contacts_tool.py:115-171 | Score-based matching with ambiguity handling |
| **File Operations** | tools/repo.py:118-130 | make_dir, get_file_path, open_file, file_exists |
| **Response Formatting** | reminder_tool.py:340-354 | Consistent response format with message |
| | contacts_tool.py:282-298 | Formatted responses with success flag |
| **Error Handling** | reminder_tool.py:501-531 | Helpful errors with available options |
| | contacts_tool.py:363-374 | Error responses with suggestions |

### Specialized Patterns

For advanced functionality, reference these specific implementations:

| Pattern | Tool/File | Lines | When You Need It |
|---------|-----------|-------|------------------|
| **Batch Operations** | contacts_tool.py | 300-351 | Bulk imports with error tracking |
| **Gated Tools** | tools/repo.py | 281-298 | Self-determining availability via is_available() |
| **Tool Dependencies** | tools/repo.py | 212-213, 522-555 | Tools that depend on other tools |
| **Dependency Injection** | tools/repo.py | 341-382 | LLMProvider, WorkingMemory injection |
| **Working Memory Integration** | tools/repo.py | 635-662 | Publishing tool hints via trinkets |
| **Tool Discovery** | tools/repo.py | 557-595 | Auto-discovery from tools/implementations/ |
| **Natural Language Dates** | reminder_tool.py | 788-852 | Parsing "tomorrow", "in 3 weeks", etc. |
| **UUID Cross-Tool Linking** | reminder_tool.py | 306-312, 777-784 | Linking reminders to contacts |
| **Partial Name Matching** | contacts_tool.py | 139-170 | starts-with, then contains matching |
| **Duplicate Detection** | contacts_tool.py | 246-258 | Case-insensitive duplicate checking |
| **Config Validation** | tools/repo.py | 57-75 | validate_config() classmethod for connection tests |
| **Manager Caching** | userdata_manager.py | 451-465 | Per-user UserDataManager caching |

## Architecture Deep Dive

### Tool Base Class (tools/repo.py:40-143)

The `Tool` base class provides automatic user scoping and file operations:

```python
class Tool(ABC):
    name = "base_tool"
    description = "Base class for all tools"

    @property
    def user_id(self) -> str:
        """Current user from context."""
        return get_current_user_id()

    @property
    def user_data_path(self) -> Path:
        """User-specific directory for this tool."""
        user_data = get_user_data_manager(self.user_id)
        return user_data.get_tool_data_dir(self.name)

    @property
    def db(self):
        """Lazy-loaded, user-scoped UserDataManager."""
        current_user_id = self.user_id
        if not self._db or self._db.user_id != current_user_id:
            self._db = get_user_data_manager(current_user_id)
        return self._db

    # File operations - user-scoped automatically
    def make_dir(self, path: str) -> Path: ...
    def get_file_path(self, filename: str) -> Path: ...
    def open_file(self, filename: str, mode: str = 'r'): ...
    def file_exists(self, filename: str) -> bool: ...
```

**Key Points:**
- **Automatic User Scoping**: `self.db` is always scoped to the current user - no manual filtering needed
- **Lazy Initialization**: Database connection created on first access
- **File Isolation**: `self.user_data_path` returns `data/users/{user_id}/tools/{tool_name}/`

### Database Operations (utils/userdata_manager.py)

UserDataManager provides a simple API for SQLite operations with automatic encryption:

```python
# Create table with schema
self.db.create_table('my_items', """
    id TEXT PRIMARY KEY,
    encrypted__title TEXT NOT NULL,
    encrypted__notes TEXT,
    created_at TEXT NOT NULL
""")

# Create indexes
self.db.execute("CREATE INDEX IF NOT EXISTS idx_items_date ON my_items(created_at)")

# CRUD operations - encryption is automatic
self.db.insert('my_items', {
    'id': item_id,
    'encrypted__title': 'Secret Title',  # Will be encrypted
    'created_at': format_utc_iso(utc_now())
})

items = self.db.select('my_items')  # Returns decrypted data
# Note: On read, encrypted__ prefix is KEPT in the field name but value is decrypted

self.db.update('my_items',
    {'encrypted__title': 'New Title'},  # Will be encrypted
    'id = :id',
    {'id': item_id}
)

self.db.delete('my_items', 'id = :id', {'id': item_id})
```

**Critical Encryption Details (userdata_manager.py:298-330):**
- Fields prefixed with `encrypted__` are automatically encrypted on write
- On read, values are automatically decrypted but **the prefix is kept in the field name**
- Access decrypted values as `item['encrypted__title']`, not `item['title']`
- Encryption key is derived deterministically from user_id (persistent across sessions)

### Connection Management (userdata_manager.py:45-66)

```python
@property
def connection(self) -> sqlite3.Connection:
    """Lazy persistent connection (thread-safe for cross-thread reuse)."""
    if self._conn is None:
        self._conn = sqlite3.connect(
            str(self.db_path),
            check_same_thread=False  # WAL mode handles concurrency
        )
        self._conn.row_factory = sqlite3.Row
    return self._conn
```

- **Lazy Creation**: Connection created on first database access
- **Thread-Safe**: `check_same_thread=False` allows cross-ThreadPoolExecutor usage
- **Cached Per-User**: `get_user_data_manager()` returns cached instances
- **Automatic Cleanup**: Connections closed on session collapse via event subscription

## Development Process

### Phase 1: Requirements Discovery

Initial descriptions often use metaphors or analogies. Extract concrete requirements:

```
Example: "Like 90s pagers"
Extracted requirements:
- High urgency messaging only
- Minimal UI complexity
- Respects user attention
- No feature creep
```

Essential questions:
- "Can you walk me through a typical usage scenario?"
- "What problem does this solve that existing tools don't?"
- "What should this tool explicitly NOT do?"
- "What would indicate success for users?"

### Phase 2: Specification Analysis

Detailed specifications often contain both explicit features and implicit design philosophy. Minor details frequently encode critical constraints.

Example: If a spec mentions "no notification fatigue," this implies rate limiting, priority systems, or other attention-management features.

### Phase 3: Codebase Pattern Study

**Use the Pattern Index above** to find exactly what you need. Each entry shows the tool and line numbers where the pattern is implemented.

**Recommended reading order:**
1. **Start simple**: `reminder_tool.py` - Clean CRUD operations, timezone handling, basic patterns
2. **Add complexity**: `contacts_tool.py` - Fuzzy search, encryption, batch operations
3. **Learn infrastructure**: `tools/repo.py` - Base class, dependency injection, discovery
4. **Understand data**: `utils/userdata_manager.py` - Database API, encryption internals

**Infrastructure references:**
```python
tools/repo.py                   # Base Tool class, ToolRepository, dependency injection
utils/userdata_manager.py       # Database API (self.db operations), encryption
utils/timezone_utils.py         # utc_now(), format_utc_iso(), convert_from_utc()
utils/user_context.py           # get_current_user_id(), get_user_preferences(), has_user_context()
```

**Critical:** Deviating from established patterns causes integration issues and maintenance debt. Always check the Pattern Index first.

### Phase 4: Incremental Implementation

Build order matters:
1. **Tool structure** - Define class, metadata, schema (see reminder_tool.py:38-112)
2. **Deferred init** - Check has_user_context() before table creation (see reminder_tool.py:114-122)
3. **Database schema** - Create tables with indexes (see reminder_tool.py:124-145)
4. **Basic CRUD** - Add/get/update/delete operations (see contacts_tool.py:222-589)
5. **Advanced features** - Search, export, etc. as needed

**Key implementation principles:**

- **Track progress**: Use TodoWrite to break down complex implementations
- **Deferred table creation**: Only create tables when user context exists (prevents startup failures)
- **Validate inputs**: Collect ALL errors before raising (see contacts_tool.py:242-258)
- **Log before raising**: Always log errors before propagating (see reminder_tool.py:244-252)
- **Consistent responses**: Use `{"success": bool, "message": str, ...}` format
- **Timezone everywhere**: Use `utc_now()` and `format_utc_iso()` for all timestamps
- **Encrypted fields**: Prefix sensitive data with `encrypted__` - access them the same way on read
- **Helpful errors**: Include suggestions and available options (see reminder_tool.py:501-531)

### Phase 5: Collaborative Solution Design

Real workshopping between skilled practitioners - direct but not theatrical:

```
H: "Need a high-urgency notification tool. Think 90s pagers - when it goes off,
    it matters. No spam, no noise."

AI: "So restricted senders? Like an allowlist of who can page you?"

H: "Allowlists are a pain to manage upfront. What if we make people solve a
    captcha before they can send?"

AI: "That won't work - you'd need an external captcha service, and if someone has
    an actual emergency they shouldn't be solving puzzles."

H: "Fair point. What about a time delay? Message sits for 30 seconds before
    sending, gives them time to reconsider?"

AI: "Same problem - real emergencies need immediate delivery. And a delay doesn't
    actually filter importance, it just annoys people."

H: "Okay, back to allowlist but somehow less annoying?"

AI: "What if we flip it - first message from a new sender goes to a review queue.
    You see it, approve or block them once, then future messages are automatic."

H: "So it's still an allowlist, just built retroactively instead of upfront?"

AI: "Right. You're not maintaining a list ahead of time. Someone sends, you review
    that specific message, make a decision, done."

H: "I like that better. But what happens when someone gets a new phone? Their
    device changed, so..."

AI: "Yeah we'd need device fingerprinting. Store a hash of their browser/device,
    check it when they send. If the hash changes, route back to review queue."

H: "That's gonna annoy people though. People upgrade phones all the time."

AI: "True. We could show context when it happens - 'New device detected for Mom,
    last device was used 2 weeks ago' so you can decide if that timing makes sense."

H: "Okay yeah, that way it's not just silently blocking them. What about message
    length? I want them short."

AI: "How short are we talking?"

H: "I don't know, maybe 300 characters? That's long enough for 'stuck in traffic,
    start the meeting without me' but short enough to force brevity."

AI: "A hard 300 char limit is pretty restrictive. What if someone has a real
    reason to send more detail?"

H: "Then it's not urgent enough for this tool, use regular messaging."

AI: "Sure, but people will push back on that. What if we auto-summarize messages
    over 300 chars instead of blocking them?"

H: "Using the LLM?"

AI: "Yeah, run it through summarization, show the summary to you by default, but
    the full message is still there if you click to expand."

H: "Does that defeat the whole brevity thing though?"

AI: "Not really - you still see brevity on your end. The sender just doesn't get
    blocked from expressing themselves. We're translating verbose to concise
    instead of punishing verbose."

H: "Okay I can see that working. So schema-wise we need operations for sending,
    approving, and blocking. Device hash field. Message field and summary field."

AI: "Standard user-scoped table. Status enum for pending, approved, blocked.
    Timestamps for created_at and maybe last_device_seen."

H: "How long should messages stick around? I don't want this to become a chat app."

AI: "24 hours? Keeps the ephemeral pager feel."

H: "Do we need a scheduled cleanup job for that?"

AI: "Could just check on read - delete anything older than 24h whenever someone
    queries their messages. Simpler than scheduling."

H: "Alright, let's build it."
```

Key observations from this exchange:
- Initial requirement was too vague, needed clarification
- Human's first two solutions (captcha, delay) were technically flawed
- AI had to directly reject bad ideas rather than being diplomatic
- Solution emerged through iteration, not a single insight
- Neither party had the answer immediately - they built it together
- Conversation included uncertainty ("I don't know", "Can we?")
- Real decision-making sounds tentative ("Let's try it") not triumphant
- Technical patterns emerged naturally from problem-solving, not name-dropped

### Phase 6: Handling Mid-Implementation Feedback

Interruptions during tool use contain valuable course corrections:

```
[Request interrupted by user for tool use]
"You're setting the default expiry to 48 hours but that's too long
for a pager metaphor. These should be ephemeral - 24 hours max."
```

This feedback indicates design misalignment. Parse for:
- Specific parameter corrections
- Underlying philosophy mismatches
- Missing requirements

## Technical Requirements

### User Scoping

**Every tool automatically gets user-scoped access via `self.db`** - no manual filtering needed.

**See reminder_tool.py:124-145** for complete database patterns including:
- Table creation with proper schema
- Encrypted fields (use `encrypted__` prefix)
- Indexes for performance
- CRUD operations

**Key Points:**

1. **Automatic User Scoping**: All `self.db` operations are scoped to the current user
2. **Automatic Encryption**: Fields prefixed with `encrypted__` are encrypted on write, decrypted on read
3. **Prefix Retained on Read**: Access decrypted fields as `item['encrypted__title']` (prefix is kept)

**Example:**
```python
# Creating a table with encrypted fields
schema = """
    id TEXT PRIMARY KEY,
    encrypted__title TEXT NOT NULL,
    encrypted__notes TEXT,
    created_at TEXT NOT NULL
"""
self.db.create_table('my_items', schema)

# Insert - encryption happens automatically
self.db.insert('my_items', {
    'id': item_id,
    'encrypted__title': 'Secret Meeting',  # Will be encrypted
    'encrypted__notes': 'Confidential',     # Will be encrypted
    'created_at': timestamp
})

# Select - decryption happens automatically, prefix is KEPT
items = self.db.select('my_items')
# Returns: [{'id': '...', 'encrypted__title': 'Secret Meeting', 'encrypted__notes': 'Confidential', ...}]
#          ^^^^ Note: 'encrypted__title' NOT 'title'
```

### Deferred Table Creation

Tools should only create tables when user context exists. This prevents startup failures during tool discovery:

```python
def __init__(self):
    super().__init__()
    self.logger = logging.getLogger(__name__)

    # Only create tables if user context is available (not during startup/discovery)
    from utils.user_context import has_user_context
    if has_user_context():
        self._ensure_tables()

def _ensure_tables(self):
    """Create tables if they don't exist."""
    schema = """
        id TEXT PRIMARY KEY,
        encrypted__data TEXT,
        created_at TEXT NOT NULL
    """
    self.db.create_table('my_data', schema)
```

Alternatively, ensure tables exist on first use in the `run()` method:

```python
def run(self, operation: str, **kwargs) -> Dict[str, Any]:
    # Ensure tables exist on first use
    self._ensure_tables()
    # ... rest of operation routing
```

### Credential Storage

For API keys and sensitive credentials, use `UserCredentialService`:

```python
from utils.user_credentials import UserCredentialService

cred_service = UserCredentialService()
api_key = cred_service.get_credential(user_id, 'api_key', 'service_name')
cred_service.store_credential(user_id, 'api_key', 'service_name', key_value)
```

### File Operations

**See tools/repo.py:118-130** for file operation methods.

Tools get automatic file methods that are user-scoped:
- `self.open_file(filename, mode)` - Open file in tool's user directory
- `self.get_file_path(filename)` - Get full path to file
- `self.file_exists(filename)` - Check if file exists
- `self.make_dir(path)` - Create subdirectory

**Example:**
```python
# Export data to JSON
filename = f"export_{utc_now().strftime('%Y%m%d_%H%M%S')}.json"
with self.open_file(filename, 'w') as f:
    json.dump(data, f, indent=2)

full_path = self.get_file_path(filename)
return {"success": True, "file_path": str(full_path)}
```

### Timezone Handling

**See reminder_tool.py:788-852** for natural language date parsing and **reminder_tool.py:335-338** for display conversion.

**CRITICAL: Always use UTC internally, convert only for display**

```python
from utils.timezone_utils import (
    utc_now, format_utc_iso, parse_utc_time_string, convert_from_utc
)
from utils.user_context import get_user_preferences

# Store as UTC ISO strings (ALWAYS)
timestamp = format_utc_iso(utc_now())
self.db.insert('items', {'created_at': timestamp})

# Parse stored UTC strings
stored_dt = parse_utc_time_string(item['created_at'])

# Convert to user's timezone ONLY for display
user_tz = get_user_preferences().timezone
local_dt = convert_from_utc(stored_dt, user_tz)
display_string = format_datetime(local_dt, "date_time_short", include_timezone=True)
```

**Available timezone utilities (utils/timezone_utils.py):**
- `utc_now()` - Get current UTC datetime (lines 170-179)
- `format_utc_iso(dt)` - Format as ISO 8601 string (lines 307-320)
- `parse_utc_time_string(time_str)` - Parse ISO string to UTC datetime (lines 491-502)
- `ensure_utc(dt)` - Ensure datetime is UTC-aware (lines 145-167)
- `convert_from_utc(dt, to_tz)` - Convert UTC to user timezone (lines 225-242)
- `convert_to_timezone(dt, target_tz)` - Convert between timezones (lines 182-206)
- `format_relative_time(dt)` - "5 hours ago" format (lines 323-398)

### Code Organization

Tools don't require specific section markers, but consistency helps. Look at existing tools for organization patterns:

```python
# Standard tool structure
import statements
logging setup
configuration class (if needed)
tool class with:
    - metadata (name, descriptions)
    - anthropic_schema
    - __init__
    - run() method
    - operation handlers (_add_item, _get_items, etc.)
    - helper methods
```

### Tool Registration and Auto-Discovery

**See reminder_tool.py:29-35** for configuration examples.

**Auto-Discovery**: Place your tool in `tools/implementations/` and restart MIRA. The ToolRepository automatically scans this package using `pkgutil.iter_modules()` (tools/repo.py:557-595).

**Configuration (Optional)**:
```python
from pydantic import BaseModel, Field
from tools.registry import registry

class MyToolConfig(BaseModel):
    enabled: bool = Field(default=True, description="Whether enabled")
    max_items: int = Field(default=10, description="Max items to return")

# Register if you have custom config beyond 'enabled'
registry.register("my_tool", MyToolConfig)
```

If you don't register a config, a default one with just `enabled: bool = True` is auto-created (tools/repo.py:81-97).

### Tool Descriptions

**See reminder_tool.py:48-53** for complete metadata examples.

Two required fields:
- `simple_description`: Ultra-concise action phrase (used by invokeother_tool for discovery)
- `anthropic_schema`: Full schema with `name`, `description`, and `input_schema`

### Anthropic Schema

**See reminder_tool.py:51-112** for complete schema with operations.

**Critical points:**
- Set `"additionalProperties": false` - prevents unexpected params
- Use `enum` for fixed options
- Clear descriptions - Claude uses these to understand how to call the tool
- Mark required fields in `"required"` array

### Database Operations Quick Reference

**See utils/userdata_manager.py:332-405** for full implementation.

```python
# Create table
schema = """
    id TEXT PRIMARY KEY,
    encrypted__title TEXT NOT NULL,
    created_at TEXT NOT NULL
"""
self.db.create_table('items', schema)

# Create indexes
self.db.execute("CREATE INDEX IF NOT EXISTS idx_items_created ON items(created_at)")

# CRUD operations
self.db.insert('items', {'id': id, 'encrypted__title': title, 'created_at': timestamp})
items = self.db.select('items', 'status = :status', {'status': 'active'})
self.db.update('items', {'encrypted__title': new_title}, 'id = :id', {'id': id})
self.db.delete('items', 'id = :id', {'id': id})

# Raw SQL for complex queries
results = self.db.execute("SELECT * FROM items WHERE created_at > :date", {'date': cutoff})
single = self.db.fetchone("SELECT * FROM items WHERE id = :id", {'id': id})
```

### Quick Start Template

Start with this minimal structure, then add patterns from the index as needed:

```python
import logging
import uuid
from typing import Dict, Any, Optional
from pydantic import BaseModel, Field

from tools.repo import Tool
from tools.registry import registry
from utils.timezone_utils import utc_now, format_utc_iso

logger = logging.getLogger(__name__)


class MyToolConfig(BaseModel):
    enabled: bool = Field(default=True, description="Whether enabled")

registry.register("my_tool", MyToolConfig)


class MyTool(Tool):
    name = "my_tool"
    simple_description = "does something useful"

    anthropic_schema = {
        "name": "my_tool",
        "description": "What this tool does and when to use it",
        "input_schema": {
            "type": "object",
            "properties": {
                "operation": {
                    "type": "string",
                    "enum": ["add", "get", "delete"],
                    "description": "The operation to perform"
                },
                "title": {
                    "type": "string",
                    "description": "Title for the item (required for add)"
                },
                "item_id": {
                    "type": "string",
                    "description": "ID of item (required for get, delete)"
                }
            },
            "required": ["operation"],
            "additionalProperties": False
        }
    }

    def __init__(self):
        super().__init__()
        self.logger = logging.getLogger(__name__)

        # Only create tables if user context is available
        from utils.user_context import has_user_context
        if has_user_context():
            self._ensure_tables()

    def _ensure_tables(self):
        schema = """
            id TEXT PRIMARY KEY,
            encrypted__title TEXT NOT NULL,
            encrypted__notes TEXT,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        """
        self.db.create_table('my_items', schema)
        self.db.execute("CREATE INDEX IF NOT EXISTS idx_my_items_created ON my_items(created_at)")

    def run(self, operation: str, **kwargs) -> Dict[str, Any]:
        try:
            # Ensure tables exist on first use
            self._ensure_tables()

            if operation == "add":
                return self._add(**kwargs)
            elif operation == "get":
                return self._get(**kwargs)
            elif operation == "delete":
                return self._delete(**kwargs)
            else:
                raise ValueError(f"Unknown operation: {operation}. Valid: add, get, delete")
        except Exception as e:
            self.logger.error(f"Error in {operation}: {e}")
            raise

    def _add(self, title: str, notes: Optional[str] = None) -> Dict[str, Any]:
        if not title:
            raise ValueError("Title is required")

        item_id = f"item_{uuid.uuid4().hex[:8]}"
        timestamp = format_utc_iso(utc_now())

        self.db.insert('my_items', {
            'id': item_id,
            'encrypted__title': title,
            'encrypted__notes': notes,
            'created_at': timestamp,
            'updated_at': timestamp
        })

        return {
            "success": True,
            "item": {"id": item_id, "encrypted__title": title},
            "message": f"Created item: {title}"
        }

    def _get(self, item_id: Optional[str] = None) -> Dict[str, Any]:
        if item_id:
            items = self.db.select('my_items', 'id = :id', {'id': item_id})
            if not items:
                return {"success": False, "message": f"Item {item_id} not found"}
            return {"success": True, "item": items[0]}
        else:
            items = self.db.select('my_items')
            return {"success": True, "items": items, "count": len(items)}

    def _delete(self, item_id: str) -> Dict[str, Any]:
        if not item_id:
            raise ValueError("item_id is required for delete")

        items = self.db.select('my_items', 'id = :id', {'id': item_id})
        if not items:
            return {"success": False, "message": f"Item {item_id} not found"}

        self.db.delete('my_items', 'id = :id', {'id': item_id})
        return {
            "success": True,
            "message": f"Deleted item: {items[0]['encrypted__title']}"
        }
```

### Error Handling

**See reminder_tool.py:501-531** for helpful error patterns.

**Key principles:**
- Always log before raising: `self.logger.error(f"Error: {e}")` then `raise`
- Use `ValueError` for invalid input, `RuntimeError` for operation failures
- Provide helpful messages with suggestions
- Show available options when something isn't found

```python
def _get_item_not_found_error(self, item_id: str) -> str:
    """Generate helpful error message with available items."""
    items = self.db.select('my_items')

    error_msg = f"Item '{item_id}' not found."

    if not items:
        error_msg += " No items available."
    else:
        available = [f"  - {i['id']}: {i['encrypted__title']}" for i in items[:5]]
        error_msg += f"\n\nAvailable items ({len(items)} total):\n"
        error_msg += "\n".join(available)
        if len(items) > 5:
            error_msg += f"\n  ... and {len(items) - 5} more"

    return error_msg
```

### Dependency Injection

**See tools/repo.py:341-382** for how injection works.

Tools can request dependencies in `__init__`:
```python
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from working_memory.core import WorkingMemory

class MyTool(Tool):
    def __init__(self, working_memory: Optional["WorkingMemory"] = None):
        super().__init__()
        self._working_memory = working_memory
```

ToolRepository automatically injects these when creating instances:
- `LLMProvider` / `LLMBridge` - For LLM calls
- `ToolRepository` - For invoking other tools
- `WorkingMemory` - For publishing trinket updates

### Gated Tools

**See tools/repo.py:281-298** for gated tool registration.

Gated tools self-determine their availability via `is_available()` method:

```python
class MyGatedTool(Tool):
    name = "my_gated_tool"

    def is_available(self) -> bool:
        """Check if tool should appear in tool list."""
        # Example: only available if config file exists
        return self.file_exists("config.json")
```

Register as gated in ToolRepository:
```python
tool_repo.register_gated_tool("my_gated_tool")
```

Unlike regular enabled tools, gated tools:
- Cannot be enabled/disabled via `enable_tool()`/`disable_tool()`
- Are checked at invocation time via `is_available()`
- Automatically appear/disappear from tool list based on state

### Config Validation (Optional)

**See tools/repo.py:57-75** for the base pattern.

Tools that need custom validation (connection tests, auto-discovery) can override the `validate_config` classmethod. This is called by the `/actions/tools/{tool}/validate` API endpoint.

```python
@classmethod
def validate_config(cls, config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate configuration and return discovered data.

    Args:
        config: The configuration dict to validate

    Returns:
        Dict with discovered data (e.g., {"folders": [...], "connection_test": "success"})

    Raises:
        ValueError: If validation fails (e.g., bad credentials, unreachable server)
    """
    server = config.get("server")
    password = config.get("password")

    if not server or not password:
        raise ValueError("Missing required fields: server, password")

    try:
        connection = connect_to_server(server, password)
        discovered_settings = connection.discover()
        return {
            "connection_test": "success",
            "discovered": discovered_settings
        }
    except ConnectionError as e:
        raise ValueError(f"Connection failed: {e}")
```

**Key points:**
- Use `@classmethod` - validation happens before tool instantiation
- Return discovered data that the frontend can use (folders, capabilities, etc.)
- Raise `ValueError` with helpful messages on failure - these are shown to users
- The base `Tool` class returns `{}` by default (no validation needed)

## Common Failure Patterns

| Pattern | Indicator | Resolution |
|---------|-----------|------------|
| Building without understanding | No clarifying questions asked | Stop coding, explore use case |
| Feature creep | Adding unrequested functionality | Return to core requirements |
| Over-engineering | Complex solution to simple problem | Discuss simpler alternatives |
| Under-communication | Assumptions instead of questions | Increase dialogue frequency |
| Misaligned mental models | "That's not quite right" feedback | Deep dive into user's vision |
| Table creation at startup | Errors during tool discovery | Use `has_user_context()` check |
| Assuming prefix stripping | Using `item['title']` for encrypted fields | Use `item['encrypted__title']` |

## Best Practices

1. **Question assumptions**: Initial requirements are rarely complete
2. **Propose alternatives**: When you see technical issues, suggest better approaches
3. **Embrace iteration**: First implementations reveal hidden requirements
4. **Respect both perspectives**: Technical elegance without usability fails; user-friendly but broken also fails
5. **Document decisions**: When choosing between approaches, note why
6. **Defer table creation**: Always check `has_user_context()` before creating tables in `__init__`
7. **Keep encrypted prefix**: Access decrypted data as `item['encrypted__field']`, not `item['field']`

## Summary

Tool building is a discovery process. The human's initial vision and the final implementation often differ significantly - not due to miscommunication, but because dialogue reveals better solutions. Success requires engaged technical discussion where both perspectives challenge and refine each other.

The extended pager example above isn't just documentation - it's a template for how these discoveries happen. Study the pattern: question, propose, identify issues, iterate, converge on solution neither party initially envisioned.
