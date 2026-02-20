# Format-Specific Exploration Recipes

Per-format commands for exploring structured data with dasel v3. Each section covers the format's unique characteristics and provides concrete commands.

## Table of Contents

- [JSON Exploration](#json-exploration)
- [YAML Exploration](#yaml-exploration)
- [TOML Exploration](#toml-exploration)
- [XML Exploration](#xml-exploration)
- [CSV Exploration](#csv-exploration)
- [HCL Exploration](#hcl-exploration)
- [INI Exploration](#ini-exploration)

---

## JSON Exploration

JSON is dasel's native format. All operations work without caveats.

### Discover Structure

```bash
# Top-level keys
dasel -f data.json 'keys($this)'

# Nested object keys
dasel -f data.json 'keys(config.database)'
```

### Inspect Arrays of Objects

```bash
# Preview first element to understand object shape
dasel -f data.json 'users[0]'

# Keys of first element (reveals schema of array items)
dasel -f data.json 'keys(users[0])'

# Count elements
dasel -f data.json 'len(users)'

# Sample first 3 elements
dasel -f data.json 'users[0:3]'
```

### Search Nested Objects

```bash
# Find all objects containing "email" at any depth
dasel -f data.json 'search(has("email"))'

# Find objects where "status" equals "active"
dasel -f data.json 'search(status == "active")'

# Combined predicate
dasel -f data.json 'search(has("id") && has("name"))'
```

### Extract Specific Fields from Arrays

```bash
# Get all names from users array
dasel -f data.json 'users.map(name)'

# Get all unique roles
dasel -f data.json 'users.map(role)' | dasel -i json '$this...' | sort -u
```

### Mixed-Type Detection

```bash
# Check type of ambiguous field
dasel -f data.json 'typeOf(config.port)'
# Returns: "int", "string", "float", etc.

# Check if field is array or object
dasel -f data.json 'typeOf(metadata)'
```

---

## YAML Exploration

YAML supports anchors, multi-document files, and deep nesting common in Kubernetes and CI configs.

### Discover Structure

```bash
# Top-level keys
dasel -f config.yaml 'keys($this)'

# Deeply nested config (common in Kubernetes)
dasel -f deployment.yaml 'keys(spec.template.spec)'
```

### Deep Nesting Navigation

```bash
# Step through Kubernetes manifest
dasel -f deployment.yaml 'keys($this)'
# -> ["apiVersion", "kind", "metadata", "spec"]

dasel -f deployment.yaml 'keys(spec)'
# -> ["replicas", "selector", "template"]

dasel -f deployment.yaml 'spec.template.spec.containers[0].image'
```

### Inspect Nested Arrays

```bash
# List container names in a pod spec
dasel -f deployment.yaml 'spec.template.spec.containers.map(name)'

# Count containers
dasel -f deployment.yaml 'len(spec.template.spec.containers)'

# Get all environment variable names from first container
dasel -f deployment.yaml 'spec.template.spec.containers[0].env.map(name)'
```

### Recursive Key Search

```bash
# Find all keys named "image" at any depth
dasel -f deployment.yaml '..image'

# Find all objects with a "name" field
dasel -f deployment.yaml 'search(has("name"))'
```

---

## TOML Exploration

TOML uses tables (sections) and arrays of tables. Common in Rust (Cargo.toml), Python (pyproject.toml), and Go configs.

### Enumerate Tables

```bash
# Top-level tables and keys
dasel -f pyproject.toml 'keys($this)'

# Nested tables
dasel -f Cargo.toml 'keys(dependencies)'
dasel -f pyproject.toml 'keys(tool.ruff)'
```

### Inspect Arrays of Tables

```bash
# pyproject.toml optional-dependencies or similar array sections
dasel -f Cargo.toml 'keys(bin[0])'
dasel -f Cargo.toml 'len(bin)'
```

### Value Extraction

```bash
# Package metadata
dasel -f pyproject.toml 'project.name'
dasel -f pyproject.toml 'project.version'

# Dependency list
dasel -f pyproject.toml 'keys(project.dependencies)'
```

### Type Checking

```bash
# Verify field types (TOML distinguishes int, float, string, bool, datetime)
dasel -f config.toml 'typeOf(server.port)'
dasel -f config.toml 'typeOf(server.debug)'
```

---

## XML Exploration

XML has elements, attributes, and text content. Dasel represents these as nested objects.

### Element Listing

```bash
# Top-level element (XML has single root)
dasel -f data.xml 'keys($this)'

# Child elements
dasel -f pom.xml 'keys(project)'
dasel -f pom.xml 'keys(project.dependencies)'
```

### Inspect Repeated Elements

```bash
# Count dependency entries in Maven POM
dasel -f pom.xml 'len(project.dependencies.dependency)'

# First dependency
dasel -f pom.xml 'project.dependencies.dependency[0]'
```

### Attribute Access

XML attributes are accessible as properties on the element. The exact representation depends on dasel's XML mapping. Inspect the element first:

```bash
# View full element to see attribute representation
dasel -f data.xml 'root.element[0]'
```

### Recursive Search

```bash
# Find all elements with a specific child
dasel -f data.xml 'search(has("version"))'
```

---

## CSV Exploration

CSV files are represented as arrays of objects (header row becomes keys).

### Header Extraction

```bash
# Get column names from first row
dasel -f data.csv 'keys($this[0])'
```

### Row Count

```bash
dasel -f data.csv 'len($this)'
```

### Sample Rows

```bash
# First row (as object with column keys)
dasel -f data.csv '$this[0]'

# First 3 rows
dasel -f data.csv '$this[0:3]'
```

### Column Sampling

```bash
# Extract single column values
dasel -f data.csv '$this.map(name)'

# Unique values in a column
dasel -f data.csv '$this.map(category)' | dasel -i json '$this...' | sort -u
```

### Type Inspection

CSV values are typically strings. Check with:

```bash
dasel -f data.csv 'typeOf($this[0].age)'
```

---

## HCL Exploration

HCL (HashiCorp Configuration Language) uses blocks with labels. Common in Terraform files.

### Block Discovery

```bash
# Top-level block types
dasel -f main.tf 'keys($this)'
# -> ["resource", "variable", "output", "provider"]
```

### Inspect Block Labels

```bash
# Resource types
dasel -f main.tf 'keys(resource)'

# Specific resource
dasel -f main.tf 'resource.aws_instance'
dasel -f main.tf 'keys(resource.aws_instance)'
```

### Nested Block Inspection

```bash
# Terraform resource attributes
dasel -f main.tf 'resource.aws_instance.web'
dasel -f main.tf 'keys(resource.aws_instance.web)'
```

### Variable Discovery

```bash
# List all variable names
dasel -f variables.tf 'keys(variable)'

# Get variable default value
dasel -f variables.tf 'variable.region.default'
```

---

## INI Exploration

INI files have sections and key-value pairs. Common in systemd, PHP, Git config.

### Section Listing

```bash
# List all sections
dasel -f config.ini 'keys($this)'
```

### Key Listing Within Section

```bash
# Keys in a specific section
dasel -f config.ini 'keys(database)'
dasel -f config.ini 'keys(server)'
```

### Value Extraction

```bash
dasel -f config.ini 'database.host'
dasel -f config.ini 'server.port'
```

### Full Section Dump

```bash
# View all keys and values in a section
dasel -f config.ini 'database'
```

---

## Cross-Format Tips

- **Type ambiguity**: CSV and INI store everything as strings. Use `typeOf()` to confirm before numeric operations.
- **XML root**: XML always has exactly one root element. `keys($this)` returns a single-element array.
- **TOML datetimes**: TOML natively supports datetime types. `typeOf()` reports these distinctly.
- **HCL labels**: HCL block labels become nested object keys. `resource.aws_instance.web` navigates type -> label -> label.
- **Large files**: Use `search()` instead of recursive descent (`..`) when you need predicate-based filtering, not exhaustive traversal.

## References

- [Dasel Query Syntax](https://daseldocs.tomwright.me/syntax/query-syntax.md) (fetched 2026-02-19)
- [Dasel Functions Index](https://daseldocs.tomwright.me/functions) (fetched 2026-02-19)
- [Dasel Arrays and Slices](https://daseldocs.tomwright.me/syntax/arrays-slices.md) (fetched 2026-02-19)
