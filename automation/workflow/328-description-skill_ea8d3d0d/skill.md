---
description: Workflows and patterns for exploring unknown structured data files with dasel v3 — discover structure, list keys, find nested values, sample arrays, identify data types across JSON, YAML, TOML, XML, CSV, HCL, INI
---

# Data Exploration with Dasel v3

<when_to_use>

Activate this skill when:

- Exploring unfamiliar structured data files (config, API responses, datasets)
- Discovering the schema or shape of a document before modifying it
- Investigating nested config structures (Kubernetes manifests, CI pipelines, package files)
- Sampling large arrays or deeply nested objects to understand content
- Identifying data types before transformation or extraction

</when_to_use>

## Supported Formats

Dasel auto-detects format from file extension. Override with `-i <format>` when reading from stdin or when extension is ambiguous.

Format identifiers: `json`, `yaml`, `toml`, `xml`, `csv`, `hcl`, `ini`

## Universal Exploration Workflow

Follow this sequence when encountering an unknown structured data file. Each step narrows scope.

### Step 1 — Format Detection

Dasel infers format from file extension. For stdin or non-standard extensions, specify explicitly:

```bash
cat mystery_file | dasel -i json 'keys($this)'
```

### Step 2 — Top-Level Keys

```bash
dasel -f config.yaml 'keys($this)'
```

Output: array of top-level key names. This is always the first exploration command.

### Step 3 — Structure Preview

For small files (configs, manifests), dump the full document:

```bash
dasel -f config.yaml
```

For large files, skip to Step 4.

### Step 4 — Nested Key Discovery

Navigate level by level:

```bash
dasel -f config.yaml 'server'
dasel -f config.yaml 'keys(server)'
dasel -f config.yaml 'keys(server.logging)'
```

Recursive key discovery across all depths:

```bash
dasel -f config.yaml '..keys($this)'
```

### Step 5 — Array Sampling

Preview first few elements without loading entire array:

```bash
dasel -f data.json 'items[0:3]'
```

Single element inspection:

```bash
dasel -f data.json 'items[0]'
```

### Step 6 — Type Inspection

Determine the type of any node:

```bash
dasel -f data.json 'typeOf(settings)'
dasel -f data.json 'typeOf(items[0].count)'
```

Return values: `"string"`, `"array"`, `"bool"`, `"null"`, `"int"`, `"float"`

### Step 7 — Value Extraction

Once path is known, extract specific values:

```bash
dasel -f config.yaml 'database.connection.host'
dasel -f data.json 'users[0].email'
```

## Exploration Patterns

### Breadth-First Exploration

Start at root, enumerate keys at each level before going deeper:

```bash
dasel -f file.json 'keys($this)'           # Level 0
dasel -f file.json 'keys(metadata)'         # Level 1
dasel -f file.json 'keys(metadata.labels)'  # Level 2
```

### Search-Based Exploration (Large Files)

When the file is too large for manual traversal, use `search()` with predicates:

```bash
# Find all objects containing a specific key
dasel -f data.json 'search(has("email"))'

# Find all objects with both "id" and "name" keys
dasel -f data.json 'search(has("id") && has("name"))'

# Find nodes where a value matches
dasel -f data.json 'search($this == 42)'
```

### Count Elements

```bash
dasel -f data.json 'len(items)'
dasel -f data.json 'len(keys($this))'
```

### Unique Value Discovery

Extract a field from all array elements, then deduplicate in shell:

```bash
dasel -f data.json 'items.map(category)' | dasel -i json '$this...' | sort -u
```

### Recursive Descent

Find all values for a key name at any depth:

```bash
dasel -f data.json '..name'
```

Get first element of every nested array:

```bash
dasel -f data.json '..[0]'
```

## Format-Specific Recipes

For detailed per-format exploration commands, see [Format-Specific Recipes](./references/format-recipes.md).

## References

- [Dasel v3 Documentation](https://daseldocs.tomwright.me) (fetched 2026-02-19)
- [Dasel Functions Index](https://daseldocs.tomwright.me/functions) (fetched 2026-02-19)
- [Dasel GitHub Repository](https://github.com/TomWright/dasel) (fetched 2026-02-19)
