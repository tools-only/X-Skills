# Sessions

Sessions manage conversation state, track context usage, and extract long-term memories.

## API Reference

### client.session()

Create a new session or load an existing one.

**Signature**

```python
def session(self, session_id: Optional[str] = None) -> Session
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| session_id | str | No | None | Session ID. Creates new session with auto-generated ID if None |

**Returns**

| Type | Description |
|------|-------------|
| Session | Session object |

**Example: Create New Session**

```python
import openviking as ov

client = ov.OpenViking(path="./data", user="alice")
client.initialize()

# Create new session (auto-generated ID)
session = client.session()
print(f"Session URI: {session.uri}")

client.close()
```

**Example: Load Existing Session**

```python
import openviking as ov

client = ov.OpenViking(path="./data", user="alice")
client.initialize()

# Load existing session
session = client.session(session_id="abc123")
session.load()
print(f"Loaded {len(session.messages)} messages")

client.close()
```

---

### Session.add_message()

Add a message to the session.

**Signature**

```python
def add_message(
    self,
    role: str,
    parts: List[Part],
) -> Message
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| role | str | Yes | - | Message role: "user" or "assistant" |
| parts | List[Part] | Yes | - | List of message parts (TextPart, ContextPart, ToolPart) |

**Returns**

| Type | Description |
|------|-------------|
| Message | Created message object |

**Part Types**

```python
from openviking.message import TextPart, ContextPart, ToolPart

# Text content
TextPart(text="Hello, how can I help?")

# Context reference
ContextPart(
    uri="viking://resources/docs/auth/",
    context_type="resource",  # "resource", "memory", or "skill"
    abstract="Authentication guide..."
)

# Tool call
ToolPart(
    tool_id="call_123",
    tool_name="search_web",
    skill_uri="viking://skills/search-web/",
    tool_input={"query": "OAuth best practices"},
    tool_output="",
    tool_status="pending"  # "pending", "running", "completed", "error"
)
```

**Example: Text Message**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# Add user message
session.add_message("user", [
    TextPart(text="How do I authenticate users?")
])

# Add assistant response
session.add_message("assistant", [
    TextPart(text="You can use OAuth 2.0 for authentication...")
])

client.close()
```

**Example: With Context Reference**

```python
import openviking as ov
from openviking.message import TextPart, ContextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

session.add_message("assistant", [
    TextPart(text="Based on the documentation..."),
    ContextPart(
        uri="viking://resources/docs/auth/",
        context_type="resource",
        abstract="Authentication guide covering OAuth 2.0..."
    )
])

client.close()
```

**Example: With Tool Call**

```python
import openviking as ov
from openviking.message import TextPart, ToolPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# Add message with tool call
msg = session.add_message("assistant", [
    TextPart(text="Let me search for that..."),
    ToolPart(
        tool_id="call_123",
        tool_name="search_web",
        skill_uri="viking://skills/search-web/",
        tool_input={"query": "OAuth best practices"},
        tool_status="pending"
    )
])

# Later, update tool result
session.update_tool_part(
    message_id=msg.id,
    tool_id="call_123",
    output="Found 5 relevant articles...",
    status="completed"
)

client.close()
```

---

### Session.used()

Track which contexts and skills were actually used in the conversation.

**Signature**

```python
def used(
    self,
    contexts: Optional[List[str]] = None,
    skill: Optional[Dict[str, Any]] = None,
) -> None
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| contexts | List[str] | No | None | List of context URIs that were used |
| skill | Dict | No | None | Skill usage info with uri, input, output, success |

**Skill Dict Structure**

```python
{
    "uri": "viking://skills/search-web/",
    "input": "search query",
    "output": "search results...",
    "success": True  # default True
}
```

**Example**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# Search for relevant contexts
results = client.find("authentication")

# Use the contexts in your response
session.add_message("assistant", [
    TextPart(text="Based on the documentation...")
])

# Track which contexts were actually helpful
session.used(contexts=[
    "viking://resources/auth-docs/"
])

# Track skill usage
session.used(skill={
    "uri": "viking://skills/code-search/",
    "input": "search for auth examples",
    "output": "Found 3 example files",
    "success": True
})

session.commit()

client.close()
```

---

### Session.update_tool_part()

Update a tool call's output and status.

**Signature**

```python
def update_tool_part(
    self,
    message_id: str,
    tool_id: str,
    output: str,
    status: str = "completed",
) -> None
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| message_id | str | Yes | - | ID of the message containing the tool call |
| tool_id | str | Yes | - | ID of the tool call to update |
| output | str | Yes | - | Tool execution output |
| status | str | No | "completed" | Tool status: "completed" or "error" |

**Example**

```python
import openviking as ov
from openviking.message import ToolPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# Add tool call
msg = session.add_message("assistant", [
    ToolPart(
        tool_id="call_456",
        tool_name="execute_code",
        skill_uri="viking://skills/code-runner/",
        tool_input={"code": "print('hello')"},
        tool_status="pending"
    )
])

# Execute tool and update result
session.update_tool_part(
    message_id=msg.id,
    tool_id="call_456",
    output="hello",
    status="completed"
)

client.close()
```

---

### Session.commit()

Commit the session, archiving messages and extracting long-term memories.

**Signature**

```python
def commit(self) -> Dict[str, Any]
```

**Returns**

| Type | Description |
|------|-------------|
| Dict | Commit result with status and statistics |

**Return Structure**

```python
{
    "session_id": "abc123",
    "status": "committed",
    "memories_extracted": 3,
    "active_count_updated": 5,
    "archived": True,
    "stats": {
        "total_turns": 10,
        "contexts_used": 4,
        "skills_used": 2,
        "memories_extracted": 3
    }
}
```

**What Happens on Commit**

1. **Archive**: Current messages are archived to `history/archive_N/`
2. **Memory Extraction**: Long-term memories are extracted using LLM
3. **Deduplication**: New memories are deduplicated against existing ones
4. **Relations**: Links are created between memories and used contexts
5. **Statistics**: Usage statistics are updated

**Memory Categories**

| Category | Location | Description |
|----------|----------|-------------|
| profile | `user/memories/.overview.md` | User profile information |
| preferences | `user/memories/preferences/` | User preferences by topic |
| entities | `user/memories/entities/` | Important entities (people, projects) |
| events | `user/memories/events/` | Significant events |
| cases | `agent/memories/cases/` | Problem-solution cases |
| patterns | `agent/memories/patterns/` | Interaction patterns |

**Example**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# Add conversation
session.add_message("user", [
    TextPart(text="I prefer dark mode and vim keybindings")
])
session.add_message("assistant", [
    TextPart(text="I've noted your preferences for dark mode and vim keybindings.")
])

# Commit session
result = session.commit()
print(f"Status: {result['status']}")
print(f"Memories extracted: {result['memories_extracted']}")
print(f"Stats: {result['stats']}")

client.close()
```

---

### Session.load()

Load session data from storage.

**Signature**

```python
def load(self) -> None
```

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Load existing session
session = client.session(session_id="existing-session-id")
session.load()

print(f"Loaded {len(session.messages)} messages")
for msg in session.messages:
    print(f"  [{msg.role}]: {msg.parts[0].text[:50]}...")

client.close()
```

---

### Session.get_context_for_search()

Get session context for search query expansion.

**Signature**

```python
def get_context_for_search(
    self,
    query: str,
    max_archives: int = 3,
    max_messages: int = 20
) -> Dict[str, Any]
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| query | str | Yes | - | Query to match relevant archives |
| max_archives | int | No | 3 | Maximum number of archives to retrieve |
| max_messages | int | No | 20 | Maximum number of recent messages |

**Returns**

| Type | Description |
|------|-------------|
| Dict | Context with summaries and recent messages |

**Return Structure**

```python
{
    "summaries": ["Archive 1 overview...", "Archive 2 overview...", ...],
    "recent_messages": [Message, Message, ...]
}
```

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session(session_id="existing-session")
session.load()

context = session.get_context_for_search(
    query="authentication",
    max_archives=3,
    max_messages=10
)

print(f"Summaries count: {len(context['summaries'])}")
print(f"Recent messages count: {len(context['recent_messages'])}")

client.close()
```

---

## Session Properties

| Property | Type | Description |
|----------|------|-------------|
| uri | str | Session Viking URI (`viking://session/{session_id}/`) |
| messages | List[Message] | Current messages in the session |
| stats | SessionStats | Session statistics |
| summary | str | Compression summary |
| usage_records | List[Usage] | Context and skill usage records |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# Access properties
print(f"URI: {session.uri}")
print(f"Messages: {len(session.messages)}")
print(f"Stats: {session.stats}")

client.close()
```

---

## Session Storage Structure

```
viking://session/{session_id}/
├── .abstract.md              # L0: Session overview
├── .overview.md              # L1: Key decisions
├── messages.jsonl            # Current messages
├── tools/                    # Tool executions
│   └── {tool_id}/
│       └── tool.json
├── .meta.json                # Metadata
├── .relations.json           # Related contexts
└── history/                  # Archived history
    ├── archive_001/
    │   ├── messages.jsonl
    │   ├── .abstract.md
    │   └── .overview.md
    └── archive_002/
```

---

## Full Example

```python
import openviking as ov
from openviking.message import TextPart, ContextPart, ToolPart

# Initialize client
client = ov.OpenViking(path="./my_data")
client.initialize()

# Create new session
session = client.session()

# Add user message
session.add_message("user", [
    TextPart(text="How do I configure embedding?")
])

# Search with session context
results = client.search("embedding configuration", session=session)

# Add assistant response with context reference
session.add_message("assistant", [
    TextPart(text="Based on the documentation, you can configure embedding..."),
    ContextPart(
        uri=results.resources[0].uri,
        context_type="resource",
        abstract=results.resources[0].abstract
    )
])

# Track actually used contexts
session.used(contexts=[results.resources[0].uri])

# Commit session (archive messages, extract memories)
result = session.commit()
print(f"Memories extracted: {result['memories_extracted']}")

client.close()
```

## Best Practices

### Commit Regularly

```python
# Commit after significant interactions
if len(session.messages) > 10:
    session.commit()
```

### Track What's Actually Used

```python
# Only mark contexts that were actually helpful
if context_was_useful:
    session.used(contexts=[ctx.uri])
```

### Use Session Context for Search

```python
# Better search results with conversation context
results = client.search(query, session=session)
```

### Load Before Continuing

```python
# Always load when resuming an existing session
session = client.session(session_id="existing-id")
session.load()
```

---

## Related Documentation

- [Context Types](../concepts/context-types.md) - Memory types
- [Retrieval](./05-retrieval.md) - Search with session
- [Client](./01-client.md) - Creating sessions
