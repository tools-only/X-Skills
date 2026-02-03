# File System

OpenViking provides Unix-like file system operations for managing context.

## API Reference

### abstract()

Read L0 abstract (~100 tokens summary).

**Signature**

```python
def abstract(self, uri: str) -> str
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI (must be a directory) |

**Returns**

| Type | Description |
|------|-------------|
| str | L0 abstract content (.abstract.md) |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

abstract = client.abstract("viking://resources/docs/")
print(f"Abstract: {abstract}")
# Output: "Documentation for the project API, covering authentication, endpoints..."

client.close()
```

---

### overview()

Read L1 overview, applies to directories.

**Signature**

```python
def overview(self, uri: str) -> str
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI (must be a directory) |

**Returns**

| Type | Description |
|------|-------------|
| str | L1 overview content (.overview.md) |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

overview = client.overview("viking://resources/docs/")
print(f"Overview:\n{overview}")

client.close()
```

---

### read()

Read L2 full content.

**Signature**

```python
def read(self, uri: str) -> str
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI |

**Returns**

| Type | Description |
|------|-------------|
| str | Full file content |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

content = client.read("viking://resources/docs/api.md")
print(f"Content:\n{content}")

client.close()
```

---

### ls()

List directory contents.

**Signature**

```python
def ls(self, uri: str, **kwargs) -> List[Any]
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI |
| simple | bool | No | False | Return only relative paths |
| recursive | bool | No | False | List all subdirectories recursively |

**Returns**

| Type | Description |
|------|-------------|
| List[Dict] | List of entries (when simple=False) |
| List[str] | List of paths (when simple=True) |

**Entry Structure**

```python
{
    "name": "docs",           # File/directory name
    "size": 4096,             # Size in bytes
    "mode": 16877,            # File mode
    "modTime": "2024-01-01T00:00:00Z",  # ISO timestamp
    "isDir": True,            # True if directory
    "uri": "viking://resources/docs/",  # Viking URI
    "meta": {}                # Optional metadata
}
```

**Example: Basic Listing**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

entries = client.ls("viking://resources/")
for entry in entries:
    type_str = "dir" if entry['isDir'] else "file"
    print(f"{entry['name']} - {type_str}")

client.close()
```

---

### tree()

Get directory tree structure.

**Signature**

```python
def tree(self, uri: str) -> List[Dict]
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI |

**Returns**

| Type | Description |
|------|-------------|
| List[Dict] | Flat list of entries with rel_path |

**Entry Structure**

```python
[
    {
        "name": "docs",
        "size": 4096,
        "mode": 16877,
        "modTime": "2024-01-01T00:00:00Z",
        "isDir": True,
        "rel_path": "docs/",      # Relative path from base URI
        "uri": "viking://resources/docs/"
    },
    ...
]
```

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

entries = client.tree("viking://resources/")
for entry in entries:
    type_str = "dir" if entry['isDir'] else "file"
    print(f"{entry['rel_path']} - {type_str}")

client.close()
```

---

### rm()

Remove file or directory.

**Signature**

```python
def rm(self, uri: str, recursive: bool = False) -> None
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI to remove |
| recursive | bool | No | False | Remove directory recursively |

**Returns**

| Type | Description |
|------|-------------|
| None | - |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Remove single file
client.rm("viking://resources/docs/old.md")

# Remove directory recursively
client.rm("viking://resources/old-project/", recursive=True)

client.close()
```

---

### mv()

Move file or directory.

**Signature**

```python
def mv(self, from_uri: str, to_uri: str) -> None
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| from_uri | str | Yes | - | Source Viking URI |
| to_uri | str | Yes | - | Destination Viking URI |

**Returns**

| Type | Description |
|------|-------------|
| None | - |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

client.mv(
    "viking://resources/old-name/",
    "viking://resources/new-name/"
)

client.close()
```

---

### grep()

Search content by pattern.

**Signature**

```python
def grep(
    self,
    uri: str,
    pattern: str,
    case_insensitive: bool = False
) -> Dict
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI to search in |
| pattern | str | Yes | - | Search pattern (regex) |
| case_insensitive | bool | No | False | Ignore case |

**Returns**

| Type | Description |
|------|-------------|
| Dict | Search results with matches |

**Return Structure**

```python
{
    "matches": [
        {
            "uri": "viking://resources/docs/auth.md",
            "line": 15,
            "content": "User authentication is handled by..."
        }
    ],
    "count": 1
}
```

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.grep(
    "viking://resources/",
    "authentication",
    case_insensitive=True
)

print(f"Found {results['count']} matches")
for match in results['matches']:
    print(f"  {match['uri']}:{match['line']}")
    print(f"    {match['content']}")

client.close()
```

---

### glob()

Match files by pattern.

**Signature**

```python
def glob(self, pattern: str, uri: str = "viking://") -> Dict
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| pattern | str | Yes | - | Glob pattern (e.g., `**/*.md`) |
| uri | str | No | "viking://" | Starting URI |

**Returns**

| Type | Description |
|------|-------------|
| Dict | Matching URIs |

**Return Structure**

```python
{
    "matches": [
        "viking://resources/docs/api.md",
        "viking://resources/docs/guide.md"
    ],
    "count": 2
}
```

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Find all markdown files
results = client.glob("**/*.md", "viking://resources/")
print(f"Found {results['count']} markdown files:")
for uri in results['matches']:
    print(f"  {uri}")

# Find all Python files
results = client.glob("**/*.py", "viking://resources/")
print(f"Found {results['count']} Python files")

client.close()
```

---

### link()

Create relations between resources.

**Signature**

```python
def link(
    self,
    from_uri: str,
    uris: Any,
    reason: str = ""
) -> None
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| from_uri | str | Yes | - | Source URI |
| uris | str or List[str] | Yes | - | Target URI(s) |
| reason | str | No | "" | Reason for the link |

**Returns**

| Type | Description |
|------|-------------|
| None | - |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Single link
client.link(
    "viking://resources/docs/auth/",
    "viking://resources/docs/security/",
    reason="Security best practices for authentication"
)

# Multiple links
client.link(
    "viking://resources/docs/api/",
    [
        "viking://resources/docs/auth/",
        "viking://resources/docs/errors/"
    ],
    reason="Related documentation"
)

client.close()
```

---

### relations()

Get relations for a resource.

**Signature**

```python
def relations(self, uri: str) -> List[Dict[str, Any]]
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI |

**Returns**

| Type | Description |
|------|-------------|
| List[Dict] | List of related resources |

**Return Structure**

```python
[
    {"uri": "viking://resources/docs/security/", "reason": "Security best practices"},
    {"uri": "viking://resources/docs/errors/", "reason": "Error handling"}
]
```

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

relations = client.relations("viking://resources/docs/auth/")
for rel in relations:
    print(f"Related: {rel['uri']}")
    print(f"  Reason: {rel['reason']}")

client.close()
```

---

### unlink()

Remove a relation.

**Signature**

```python
def unlink(self, from_uri: str, uri: str) -> None
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| from_uri | str | Yes | - | Source URI |
| uri | str | Yes | - | Target URI to unlink |

**Returns**

| Type | Description |
|------|-------------|
| None | - |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

client.unlink(
    "viking://resources/docs/auth/",
    "viking://resources/docs/security/"
)

client.close()
```

---

## Related Documentation

- [Viking URI](../concepts/viking-uri.md) - URI specification
- [Context Layers](../concepts/context-layers.md) - L0/L1/L2
- [Resources](resources.md) - Resource management
