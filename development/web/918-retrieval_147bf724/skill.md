# Retrieval

OpenViking provides two search methods: `find` for simple semantic search and `search` for complex retrieval with session context.

## find vs search

| Aspect | find | search |
|--------|------|--------|
| Intent Analysis | No | Yes |
| Session Context | No | Yes |
| Query Expansion | No | Yes |
| Default Limit | 10 | 10 |
| Use Case | Simple queries | Conversational search |

## API Reference

### find()

Basic vector similarity search.

**Signature**

```python
def find(
    self,
    query: str,
    target_uri: str = "",
    limit: int = 10,
    score_threshold: Optional[float] = None,
    filter: Optional[Dict] = None,
) -> FindResult
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| query | str | Yes | - | Search query string |
| target_uri | str | No | "" | Limit search to specific URI prefix |
| limit | int | No | 10 | Maximum number of results |
| score_threshold | float | No | None | Minimum relevance score threshold |
| filter | Dict | No | None | Metadata filters |

**Returns**

| Type | Description |
|------|-------------|
| FindResult | Search results containing contexts |

**FindResult Structure**

```python
class FindResult:
    memories: List[MatchedContext]   # Memory contexts
    resources: List[MatchedContext]  # Resource contexts
    skills: List[MatchedContext]     # Skill contexts
    query_plan: Optional[QueryPlan]  # Query plan (search only)
    query_results: Optional[List[QueryResult]]  # Detailed results
    total: int                       # Total count (auto-calculated)
```

**MatchedContext Structure**

```python
class MatchedContext:
    uri: str                         # Viking URI
    context_type: ContextType        # "resource", "memory", or "skill"
    is_leaf: bool                    # Whether it's a leaf node
    abstract: str                    # L0 content
    category: str                    # Category
    score: float                     # Relevance score (0-1)
    match_reason: str                # Why this matched
    relations: List[RelatedContext]  # Related contexts
```

**Example: Basic Search**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.find("how to authenticate users")

for ctx in results.resources:
    print(f"URI: {ctx.uri}")
    print(f"Score: {ctx.score:.3f}")
    print(f"Type: {ctx.context_type}")
    print(f"Abstract: {ctx.abstract[:100]}...")
    print("---")

client.close()
```

**Example: Search with Target URI**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Search only in resources
results = client.find(
    "authentication",
    target_uri="viking://resources/"
)

# Search only in user memories
results = client.find(
    "preferences",
    target_uri="viking://user/memories/"
)

# Search only in skills
results = client.find(
    "web search",
    target_uri="viking://skills/"
)

# Search in specific project
results = client.find(
    "API endpoints",
    target_uri="viking://resources/my-project/"
)

client.close()
```

---

### search()

Search with session context and intent analysis.

**Signature**

```python
def search(
    self,
    query: str,
    target_uri: str = "",
    session: Optional[Session] = None,
    limit: int = 3,
    score_threshold: Optional[float] = None,
    filter: Optional[Dict] = None,
) -> FindResult
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| query | str | Yes | - | Search query string |
| target_uri | str | No | "" | Limit search to specific URI prefix |
| session | Session | No | None | Session for context-aware search |
| limit | int | No | 3 | Maximum number of results |
| score_threshold | float | No | None | Minimum relevance score threshold |
| filter | Dict | No | None | Metadata filters |

**Returns**

| Type | Description |
|------|-------------|
| FindResult | Search results with query plan and contexts |

**Example: Session-Aware Search**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

# Create session with conversation context
session = client.session()
session.add_message("user", [
    TextPart(text="I'm building a login page with OAuth")
])
session.add_message("assistant", [
    TextPart(text="I can help you with OAuth implementation.")
])

# Search understands the conversation context
results = client.search(
    "best practices",
    session=session
)

for ctx in results.resources:
    print(f"Found: {ctx.uri}")
    print(f"Abstract: {ctx.abstract[:200]}...")

client.close()
```

**Example: Search Without Session**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# search can also be used without session
# It still performs intent analysis on the query
results = client.search(
    "how to implement OAuth 2.0 authorization code flow",
)

for ctx in results.resources:
    print(f"Found: {ctx.uri} (score: {ctx.score:.3f})")

client.close()
```

---

## Retrieval Pipeline

```
Query → Intent Analysis → Vector Search (L0) → Rerank (L1) → Results
```

1. **Intent Analysis** (search only): Understand query intent, expand queries
2. **Vector Search**: Find candidates using Embedding
3. **Rerank**: Re-score using content for accuracy
4. **Results**: Return top-k contexts

## Working with Results

### Read Content Progressively

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.find("authentication")

for ctx in results.resources:
    # Start with L0 (abstract) - already in ctx.abstract
    print(f"Abstract: {ctx.abstract}")

    if not ctx.is_leaf:
        # Get L1 (overview)
        overview = client.overview(ctx.uri)
        print(f"Overview: {overview[:500]}...")
    else:
        # Load L2 (content)
        content = client.read(ctx.uri)
        print(f"File content: {content}")

client.close()
```

### Get Related Resources

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.find("OAuth implementation")

for ctx in results.resources:
    print(f"Found: {ctx.uri}")

    # Get related resources
    relations = client.relations(ctx.uri)
    for rel in relations:
        print(f"  Related: {rel['uri']} - {rel['reason']}")

client.close()
```

## Best Practices

### Use Specific Queries

```python
# Good - specific query
results = client.find("OAuth 2.0 authorization code flow implementation")

# Less effective - too broad
results = client.find("auth")
```

### Scope Your Searches

```python
# Search in relevant scope for better results
results = client.find(
    "error handling",
    target_uri="viking://resources/my-project/"
)
```

### Use Session Context for Conversations

```python
# For conversational search, use session
from openviking.message import TextPart

session = client.session()
session.add_message("user", [
    TextPart(text="I'm building a login page")
])

# Search understands the context
results = client.search("best practices", session=session)
```

### Related Documentation

- [Resources](resources.md) - Resource management
- [Sessions](sessions.md) - Session context
- [Context Layers](../concepts/context-layers.md) - L0/L1/L2
