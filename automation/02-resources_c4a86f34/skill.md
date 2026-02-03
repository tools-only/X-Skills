# Resources

Resources are external knowledge that agents can reference. This guide covers how to add, manage, and retrieve resources.

## Supported Formats

| Format | Extensions | Processing |
|--------|------------|------------|
| PDF | `.pdf` | Text and image extraction |
| Markdown | `.md` | Native support |
| HTML | `.html`, `.htm` | Cleaned text extraction |
| Plain Text | `.txt` | Direct import |
| JSON/YAML | `.json`, `.yaml`, `.yml` | Structured parsing |
| Code | `.py`, `.js`, `.ts`, `.go`, `.java`, etc. | Syntax-aware parsing |
| Images | `.png`, `.jpg`, `.jpeg`, `.gif`, `.webp` | VLM description |
| Video | `.mp4`, `.mov`, `.avi` | Frame extraction + VLM |
| Audio | `.mp3`, `.wav`, `.m4a` | Transcription |
| Documents | `.docx` | Text extraction |

## Processing Pipeline

```
Input → Parser → TreeBuilder → AGFS → SemanticQueue → Vector Index
```

1. **Parser**: Extracts content based on file type
2. **TreeBuilder**: Creates directory structure
3. **AGFS**: Stores files in virtual file system
4. **SemanticQueue**: Generates L0/L1 asynchronously
5. **Vector Index**: Indexes for semantic search

## API Reference

### add_resource()

Add a resource to the knowledge base.

**Signature**

```python
def add_resource(
    self,
    path: str,
    target: Optional[str] = None,
    reason: str = "",
    instruction: str = "",
    wait: bool = False,
    timeout: float = None,
) -> Dict[str, Any]
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| path | str | Yes | - | Local file path, directory path, or URL |
| target | str | No | None | Target Viking URI (must be in `resources` scope) |
| reason | str | No | "" | Why this resource is being added (improves search relevance) |
| instruction | str | No | "" | Special processing instructions |
| wait | bool | No | False | Wait for semantic processing to complete |

**Returns**

| Type | Description |
|------|-------------|
| Dict[str, Any] | Result containing status and resource information |

**Return Structure**

```python
{
    "status": "success",           # "success" or "error"
    "root_uri": "viking://resources/docs/",  # Root resource URI
    "source_path": "./docs/",      # Original source path
    "errors": [],                  # List of errors (if any)
    "queue_status": {...}          # Queue status (only when wait=True)
}
```

**Example: Add Single File**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

result = client.add_resource(
    "./documents/guide.md",
    reason="User guide documentation"
)
print(f"Added: {result['root_uri']}")

client.wait_processed()
client.close()
```

**Example: Add from URL**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

result = client.add_resource(
    "https://example.com/api-docs.md",
    target="viking://resources/external/",
    reason="External API documentation"
)

# Wait for processing
client.wait_processed()
client.close()
```

**Example: Wait for Processing**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Option 1: Wait inline
result = client.add_resource(
    "./documents/guide.md",
    wait=True
)
print(f"Queue status: {result['queue_status']}")

# Option 2: Wait separately (for batch processing)
client.add_resource("./file1.md")
client.add_resource("./file2.md")
client.add_resource("./file3.md")

status = client.wait_processed()
print(f"All processed: {status}")

client.close()
```

---

### export_ovpack()

Export a resource tree as a `.ovpack` file.

**Signature**

```python
def export_ovpack(self, uri: str, to: str) -> str
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| uri | str | Yes | - | Viking URI to export |
| to | str | Yes | - | Target file path |

**Returns**

| Type | Description |
|------|-------------|
| str | Path to the exported file |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Export a project
path = client.export_ovpack(
    "viking://resources/my-project/",
    "./exports/my-project.ovpack"
)
print(f"Exported to: {path}")

client.close()
```

---

### import_ovpack()

Import a `.ovpack` file.

**Signature**

```python
def import_ovpack(
    self,
    file_path: str,
    parent: str,
    force: bool = False,
    vectorize: bool = True
) -> str
```

**Parameters**

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| file_path | str | Yes | - | Local `.ovpack` file path |
| parent | str | Yes | - | Target parent URI |
| force | bool | No | False | Overwrite existing resources |
| vectorize | bool | No | True | Trigger vectorization after import |

**Returns**

| Type | Description |
|------|-------------|
| str | Root URI of imported resources |

**Example**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# Import a package
uri = client.import_ovpack(
    "./exports/my-project.ovpack",
    "viking://resources/imported/",
    force=True,
    vectorize=True
)
print(f"Imported to: {uri}")

client.wait_processed()
client.close()
```

---

## Managing Resources

### List Resources

```python
# List all resources
entries = client.ls("viking://resources/")

# List with details
for entry in entries:
    type_str = "dir" if entry['isDir'] else "file"
    print(f"{entry['name']} - {type_str}")

# Simple path list
paths = client.ls("viking://resources/", simple=True)
# Returns: ["project-a/", "project-b/", "shared/"]

# Recursive listing
all_entries = client.ls("viking://resources/", recursive=True)
```

### Read Resource Content

```python
# L0: Abstract
abstract = client.abstract("viking://resources/docs/")

# L1: Overview
overview = client.overview("viking://resources/docs/")

# L2: Full content
content = client.read("viking://resources/docs/api.md")
```

### Move Resources

```python
client.mv(
    "viking://resources/old-project/",
    "viking://resources/new-project/"
)
```

### Delete Resources

```python
# Delete single file
client.rm("viking://resources/docs/old.md")

# Delete directory recursively
client.rm("viking://resources/old-project/", recursive=True)
```

### Create Links

```python
# Link related resources
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
```

### Get Relations

```python
relations = client.relations("viking://resources/docs/auth/")
for rel in relations:
    print(f"{rel['uri']}: {rel['reason']}")
```

### Remove Links

```python
client.unlink(
    "viking://resources/docs/auth/",
    "viking://resources/docs/security/"
)
```

## Best Practices

### Organize by Project

```
viking://resources/
├── project-a/
│   ├── docs/
│   ├── specs/
│   └── references/
├── project-b/
│   └── ...
└── shared/
    └── common-docs/
```

## Related Documentation

- [Retrieval](retrieval.md) - Search resources
- [File System](filesystem.md) - File operations
- [Context Types](../concepts/context-types.md) - Resource concept
