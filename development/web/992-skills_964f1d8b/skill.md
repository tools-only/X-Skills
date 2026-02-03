# Skills

Skills are callable capabilities that agents can invoke. This guide covers how to add and manage skills.

## API Reference

### add_skill()

Add a skill to the knowledge base.

**Signature**

```python
def add_skill(
    self,
    data: Any,
    wait: bool = False,
    timeout: float = None,
) -> Dict[str, Any]
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| data | Any | Yes | - | Skill data (dict, string, or path) |
| wait | bool | No | False | Wait for vectorization to complete |
| timeout | float | No | None | Timeout in seconds |

**Supported Data Formats**

1. **Dict (Skill format)**:
```python
{
    "name": "skill-name",
    "description": "Skill description",
    "content": "Full markdown content",
    "allowed_tools": ["Tool1", "Tool2"],  # optional
    "tags": ["tag1", "tag2"]  # optional
}
```

2. **Dict (MCP Tool format)** - Auto-detected and converted:
```python
{
    "name": "tool_name",
    "description": "Tool description",
    "inputSchema": {
        "type": "object",
        "properties": {...},
        "required": [...]
    }
}
```

3. **String (SKILL.md content)**:
```python
"""---
name: skill-name
description: Skill description
---

# Skill Content
"""
```

4. **Path (file or directory)**:
   - Single file: Path to `SKILL.md` file
   - Directory: Path to directory containing `SKILL.md` (auxiliary files included)

**Returns**

| Type | Description |
|------|-------------|
| Dict | Result containing status and skill URI |

**Return Structure**

```python
{
    "status": "success",
    "uri": "viking://agent/skills/skill-name/",
    "name": "skill-name",
    "auxiliary_files": 0
}
```

**Example: Add Skill from Dict**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

skill = {
    "name": "search-web",
    "description": "Search the web for current information",
    "content": """
# search-web

Search the web for current information.

## Parameters
- **query** (string, required): Search query
- **limit** (integer, optional): Max results, default 10

## Usage
Use when the user needs current information.
"""
}

result = client.add_skill(skill)
print(f"Added: {result['uri']}")

client.close()
```

**Example: Add from MCP Tool**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# MCP tool format is auto-detected and converted
mcp_tool = {
    "name": "calculator",
    "description": "Perform mathematical calculations",
    "inputSchema": {
        "type": "object",
        "properties": {
            "expression": {
                "type": "string",
                "description": "Mathematical expression to evaluate"
            }
        },
        "required": ["expression"]
    }
}

result = client.add_skill(mcp_tool)
print(f"Added: {result['uri']}")

client.close()
```

**Example: Add from SKILL.md File**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Add from file path
result = client.add_skill("./skills/search-web/SKILL.md")
print(f"Added: {result['uri']}")

# Add from directory (includes auxiliary files)
result = client.add_skill("./skills/code-runner/")
print(f"Added: {result['uri']}")
print(f"Auxiliary files: {result['auxiliary_files']}")

client.close()
```

---

## SKILL.md Format

Skills can be defined using SKILL.md files with YAML frontmatter.

**Structure**

```markdown
---
name: skill-name
description: Brief description of the skill
allowed-tools:
  - Tool1
  - Tool2
tags:
  - tag1
  - tag2
---

# Skill Name

Full skill documentation in Markdown format.

## Parameters
- **param1** (type, required): Description
- **param2** (type, optional): Description

## Usage
When and how to use this skill.

## Examples
Concrete examples of skill invocation.
```

**Required Fields**

| Field | Type | Description |
|-------|------|-------------|
| name | str | Skill name (kebab-case recommended) |
| description | str | Brief description |

**Optional Fields**

| Field | Type | Description |
|-------|------|-------------|
| allowed-tools | List[str] | Tools this skill can use |
| tags | List[str] | Tags for categorization |

---

## Managing Skills

### List Skills

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# List all skills
skills = client.ls("viking://agent/skills/")
for skill in skills:
    print(f"{skill['name']}")

# Simple list (names only)
names = client.ls("viking://agent/skills/", simple=True)
print(names)

client.close()
```

### Read Skill Content

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

uri = "viking://agent/skills/search-web/"

# L0: Brief description
abstract = client.abstract(uri)
print(f"Abstract: {abstract}")

# L1: Parameters and usage overview
overview = client.overview(uri)
print(f"Overview: {overview}")

# L2: Full skill documentation
content = client.read(uri)
print(f"Content: {content}")

client.close()
```

### Search Skills

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Semantic search for skills
results = client.find(
    "search the internet",
    target_uri="viking://agent/skills/",
    limit=5
)

for ctx in results.skills:
    print(f"Skill: {ctx.uri}")
    print(f"Score: {ctx.score:.3f}")
    print(f"Description: {ctx.abstract}")
    print("---")

client.close()
```

### Remove Skills

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Remove a skill
client.rm("viking://agent/skills/old-skill/", recursive=True)

client.close()
```

---

## MCP Conversion

OpenViking automatically detects and converts MCP tool definitions to skill format.

**Detection**

A dict is treated as MCP format if it contains an `inputSchema` field:

```python
if "inputSchema" in data:
    # Convert to skill format
    skill = mcp_to_skill(data)
```

**Conversion Process**

1. Name is converted to kebab-case
2. Description is preserved
3. Parameters are extracted from `inputSchema.properties`
4. Required fields are marked from `inputSchema.required`
5. Markdown content is generated

**Example Conversion**

Input (MCP format):
```python
{
    "name": "search_web",
    "description": "Search the web",
    "inputSchema": {
        "type": "object",
        "properties": {
            "query": {
                "type": "string",
                "description": "Search query"
            },
            "limit": {
                "type": "integer",
                "description": "Max results"
            }
        },
        "required": ["query"]
    }
}
```

Output (Skill format):
```python
{
    "name": "search-web",
    "description": "Search the web",
    "content": """---
name: search-web
description: Search the web
---

# search-web

Search the web

## Parameters

- **query** (string) (required): Search query
- **limit** (integer) (optional): Max results

## Usage

This tool wraps the MCP tool `search-web`. Call this when the user needs functionality matching the description above.
"""
}
```

---

## Skill Storage Structure

Skills are stored at `viking://agent/skills/`:

```
viking://agent/skills/
├── search-web/
│   ├── .abstract.md      # L0: Brief description
│   ├── .overview.md      # L1: Parameters and usage
│   ├── SKILL.md          # L2: Full documentation
│   └── [auxiliary files] # Any additional files
├── calculator/
│   ├── .abstract.md
│   ├── .overview.md
│   └── SKILL.md
└── ...
```

---

## Best Practices

### Clear Descriptions

```python
# Good - specific and actionable
skill = {
    "name": "search-web",
    "description": "Search the web for current information using Google",
    ...
}

# Less helpful - too vague
skill = {
    "name": "search",
    "description": "Search",
    ...
}
```

### Comprehensive Content

Include in your skill content:
- Clear parameter descriptions with types
- When to use the skill
- Concrete examples
- Edge cases and limitations

```python
skill = {
    "name": "search-web",
    "description": "Search the web for current information",
    "content": """
# search-web

Search the web for current information using Google.

## Parameters
- **query** (string, required): Search query. Be specific for better results.
- **limit** (integer, optional): Maximum number of results. Default: 10, Max: 100.

## Usage
Use this skill when:
- User asks about current events
- Information is not in the knowledge base
- User explicitly asks to search the web

Do NOT use when:
- Information is already available in resources
- Query is about historical facts

## Examples
- "What's the weather today?" → search-web(query="weather today")
- "Latest news about AI" → search-web(query="AI news 2024", limit=5)

## Limitations
- Rate limited to 100 requests per hour
- Results may not include paywalled content
"""
}
```

### Consistent Naming

Use kebab-case for skill names:
- `search-web` (good)
- `searchWeb` (avoid)
- `search_web` (avoid)

---

## Related Documentation

- [Context Types](../concepts/context-types.md) - Skill concept
- [Retrieval](./05-retrieval.md) - Finding skills
- [Sessions](./04-sessions.md) - Tracking skill usage
