# Format-Specific Patterns

Practical dasel v3 patterns for each supported data format, plus conversion caveats.

## JSON Patterns

### Nested Object Queries

```bash
# Deep path access
cat package.json | dasel -i json 'dependencies.react'

# Array of objects -- get all names
cat data.json | dasel -i json 'users.map(name)'

# Filter and extract
cat data.json | dasel -i json 'users.filter(role == "admin").map(email)'
```

### Array Manipulation

```bash
# Get array length
cat data.json | dasel -i json 'len(items)'

# Slice first 3
cat data.json | dasel -i json 'items[0:2]'

# Sum numeric field
cat data.json | dasel -i json 'orders.map(total).sum($this)'
```

### Modify JSON

```bash
# Set nested value
cat config.json | dasel -i json --root 'settings.theme = "dark"' > config_new.json

# Add to array
cat data.json | dasel -i json --root 'tags = [tags..., "new-tag"]'

# Remove key by reconstruction
echo '{"keep": 1, "drop": 2, "also_keep": 3}' \
  | dasel -i json --root '{ keep, also_keep }'
```

### Recursive Search in JSON

```bash
# Find all objects with "id" field anywhere in tree
cat complex.json | dasel -i json 'search(has("id"))'

# Find all string values
cat data.json | dasel -i json 'search(typeOf($this) == "string")'
```

## YAML Patterns

### Basic Queries

```bash
# Read nested key
cat config.yaml | dasel -i yaml 'database.host'

# Read list item
cat config.yaml | dasel -i yaml 'servers[0].name'

# Get all server names
cat config.yaml | dasel -i yaml 'servers.map(name)'
```

### Modify YAML

```bash
# Update a value
cat config.yaml | dasel -i yaml --root 'database.port = 5433' > config_new.yaml

# YAML preserves structure when using same input/output format
cat config.yaml | dasel -i yaml -o yaml --root 'app.debug = false'
```

### Multi-Document YAML

dasel processes the first document in a multi-document YAML stream by default. To work with specific documents, split them first or pipe individual documents:

```bash
# Process first document
cat multi.yaml | dasel -i yaml 'metadata.name'
```

### YAML Anchors and Aliases

dasel resolves YAML anchors/aliases before query execution. Queries operate on the resolved values, not the anchor references:

```bash
# Given YAML with anchors:
# defaults: &defaults
#   timeout: 30
# production:
#   <<: *defaults
#   timeout: 60

cat config.yaml | dasel -i yaml 'production.timeout'
# Output: 60 (resolved value, not anchor)
```

## TOML Patterns

### Table Access

```bash
# Read from a table
cat config.toml | dasel -i toml 'database.server'

# Read from nested table
cat config.toml | dasel -i toml 'servers.production.host'
```

### Arrays of Tables

```bash
# TOML [[products]] becomes an array
cat config.toml | dasel -i toml 'products[0].name'

# Filter array of tables
cat config.toml | dasel -i toml 'products.filter(price > 10)'

# Get all product names
cat config.toml | dasel -i toml 'products.map(name)'
```

### Modify TOML

```bash
# Set value in TOML table
cat config.toml | dasel -i toml --root 'database.port = 5432' > config_new.toml

# TOML type preservation
cat config.toml | dasel -i toml --root 'settings.max_connections = 100'
```

### Inline Tables

```bash
# Access inline table values
# Given: point = { x = 1, y = 2 }
cat data.toml | dasel -i toml 'point.x'
```

## XML Patterns

XML has two parsing modes controlled by `--read-flag xml-mode=structured`:

- **Friendly mode** (default): Attributes stored with `-` prefix (e.g., `-id`), text content as `#text`. Source: `parsing/xml/reader.go:98-165`
- **Structured mode**: Explicit keys: `name`, `attrs`, `content`, `children`. Source: `parsing/xml/reader.go:58-96`

### XML Security Limits

Source: `parsing/xml/reader.go:15-21`

- Maximum input size: 10 MB (`maxXMLSize = 10_000_000`)
- Maximum bytes per comment: 10 KB (`maxCommentLength = 10_000`)
- Maximum comments per document: 1,000 (`maxTotalComments = 1_000`)

### Element Access

```bash
# Read element text
echo '<root><name>Alice</name></root>' | dasel -i xml 'root.name'

# Nested elements
cat data.xml | dasel -i xml 'catalog.book[0].title'

# All book titles
cat data.xml | dasel -i xml 'catalog.book.map(title)'
```

### Attributes (Friendly Mode)

XML attributes are accessed via a special `-` prefix convention in friendly mode:

```bash
# Access attribute
# Given: <item id="123">content</item>
echo '<root><item id="123">content</item></root>' \
  | dasel -i xml 'root.item.-id'

# Filter by attribute
cat data.xml | dasel -i xml 'items.item.filter(-type == "active")'
```

### Text Nodes (Friendly Mode)

```bash
# Element text content accessed directly (simple elements)
echo '<root><message>Hello World</message></root>' \
  | dasel -i xml 'root.message'
# Output: Hello World

# Mixed content: text accessible via #text key
# Given: <item id="1">content</item>
echo '<root><item id="1">content</item></root>' \
  | dasel -i xml 'root.item.#text'
```

### Structured Mode

```bash
# Access with structured mode
echo '<root><item id="1">content</item></root>' \
  | dasel -i xml --read-flag xml-mode=structured 'root.children[0].attrs.id'

# Structured mode keys: name, attrs, content, children
echo '<root><item id="1">content</item></root>' \
  | dasel -i xml --read-flag xml-mode=structured 'root.children[0].content'
```

### Modify XML

```bash
# Set element value
cat data.xml | dasel -i xml --root 'config.timeout = 30' > data_new.xml
```

## CSV Patterns

### Column Access

CSV data is treated as an array of arrays (or array of objects if headers present):

```bash
# Access by row and column index
echo 'name,age
Alice,30
Bob,25' | dasel -i csv '[0].name'
# Output: Alice

# Get all values from a column
echo 'name,age
Alice,30
Bob,25' | dasel -i csv 'map(name)'
# Output: ["Alice", "Bob"]
```

### Row Filtering

```bash
# Filter rows by column value
echo 'name,age,active
Alice,30,true
Bob,25,false
Carol,35,true' | dasel -i csv 'filter(active == "true")'

# Note: CSV values are strings; compare as strings
echo 'name,score
Alice,90
Bob,75' | dasel -i csv 'filter(toInt(score) > 80)'
```

### CSV to JSON

```bash
# Convert CSV to JSON array of objects
cat data.csv | dasel -i csv -o json
```

## HCL Patterns

HCL (HashiCorp Configuration Language) support was added in v3.0.0.

### Block Access

```bash
# Read resource block
# Given HCL:
# resource "aws_instance" "web" {
#   ami = "ami-12345"
#   instance_type = "t2.micro"
# }
cat main.tf | dasel -i hcl 'resource.aws_instance.web.ami'
```

### Labels

HCL blocks with labels are accessed as nested maps:

```bash
# Access by block type and label
cat main.tf | dasel -i hcl 'variable.region.default'

# List all resource types
cat main.tf | dasel -i hcl 'keys(resource)'
```

### Nested Blocks

```bash
# Access nested block
cat main.tf | dasel -i hcl 'resource.aws_security_group.web.ingress[0].from_port'
```

## INI Patterns

INI support was added in v3.0.0.

### Section and Key Access

```bash
# Read key from section
# Given INI:
# [database]
# host = localhost
# port = 5432
cat config.ini | dasel -i ini 'database.host'
# Output: localhost

cat config.ini | dasel -i ini 'database.port'
# Output: 5432
```

### Modify INI

```bash
# Update a value
cat config.ini | dasel -i ini --root 'database.port = 3306' > config_new.ini
```

### Global Keys

Keys outside any section are accessed at the top level:

```bash
# Given INI:
# version = 2
# [server]
# host = localhost
cat config.ini | dasel -i ini 'version'
```

## Format Conversion Matrix

| From / To | JSON | YAML | TOML | XML | CSV |
|-----------|------|------|------|-----|-----|
| **JSON** | -- | Clean | Clean | Caveats | Flat only |
| **YAML** | Clean | -- | Clean | Caveats | Flat only |
| **TOML** | Clean | Clean | -- | Caveats | Flat only |
| **XML** | Caveats | Caveats | Caveats | -- | No |
| **CSV** | Clean | Clean | Limited | No | -- |

### Conversion Syntax

```bash
# Generic pattern
cat input_file | dasel -i <input_format> -o <output_format>

# JSON to YAML
cat data.json | dasel -i json -o yaml

# YAML to TOML
cat config.yaml | dasel -i yaml -o toml

# TOML to JSON (compact)
cat config.toml | dasel -i toml -o json --compact

# CSV to JSON
cat data.csv | dasel -i csv -o json
```

### Conversion Caveats

**XML conversions**:

- XML attributes (`-attr`) may not round-trip cleanly to JSON/YAML/TOML
- XML text content mixed with child elements requires careful handling
- XML namespaces are not preserved in conversion

**CSV limitations** (source: `parsing/csv/writer.go:16`, `parsing/csv/csv.go:24`):

- Root value must be a slice/array for CSV writing: `"csv writer expects root output to be a slice/array"`
- CSV to structured formats works only for flat tabular data
- Nested structures in JSON/YAML/TOML cannot represent well as CSV
- All CSV values are strings; numeric/boolean types must be cast after conversion
- Custom delimiter via `--rw-flag csv-delimiter=;` (source: `parsing/csv/csv.go:24`)

**TOML limitations** (source: `parsing/toml/toml_writer.go:28-29`):

- TOML does not support `null` values; null values from JSON/YAML are dropped or converted to empty strings
- TOML requires all keys to be strings
- TOML datetime types may not convert cleanly to other formats
- Does not preserve formatting metadata (multiline strings, etc.) on round-trip

**INI limitations** (source: `parsing/ini/ini_writer.go:26-28, 73-74, 88-89`):

- Root must be a map: `"ini can only represent map values"`
- Cannot represent slice/array values directly: `"ini writer cannot represent slice values directly"`
- Sections can only represent map values
- All values are stored as strings

**HCL limitations** (source: `parsing/hcl/writer.go:49-50`):

- Body must be a map: `"hcl body is expected to be a map"`
- HCL and INI conversion to other formats may lose format-specific semantics (HCL block labels, INI section structure)
- Round-trip fidelity is best when staying within the same format

**YAML notes** (source: `parsing/yaml/yaml_reader.go:24-54`):

- Multi-document YAML is supported: single document returns the value directly; multiple documents return a branched slice
- Empty YAML returns null

## Extended Format Options

Format-specific behavior can be customized via `--read-flag`, `--write-flag`, and `--rw-flag`. Each accepts `name=value` pairs.

Source: `internal/cli/read_write_flag.go`, `internal/cli/query.go:7-9`

| Flag | Format | Effect | Source |
|------|--------|--------|--------|
| `csv-delimiter` | CSV | Custom field delimiter | `parsing/csv/csv.go:24` |
| `xml-mode=structured` | XML | Structured parsing (name/attrs/content/children) | `parsing/xml/reader.go:25` |
| `hcl-block-format=array` | HCL | Always read block labels as slices | `parsing/hcl/reader.go:16` |

```bash
# CSV with semicolon delimiter
echo 'a;b;c\n1;2;3' | dasel -i csv --rw-flag csv-delimiter=';' 'map(a)'

# XML structured mode
cat data.xml | dasel -i xml --read-flag xml-mode=structured 'root.children[0].name'

# HCL array block format
cat main.tf | dasel -i hcl --read-flag hcl-block-format=array 'resource'
```

## Sources

- dasel README: <https://raw.githubusercontent.com/TomWright/dasel/master/README.md> (fetched 2026-02-19)
- dasel docs: <https://daseldocs.tomwright.me> (fetched 2026-02-19)
- Arrays/slices: <https://daseldocs.tomwright.me/syntax/arrays-slices.md> (fetched 2026-02-19)
- Objects/maps: <https://daseldocs.tomwright.me/syntax/objects-maps.md> (fetched 2026-02-19)
- CHANGELOG (format support): <https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md> (fetched 2026-02-19)
- Source code analysis: `parsing/xml/reader.go`, `parsing/csv/writer.go`, `parsing/csv/csv.go`, `parsing/ini/ini_writer.go`, `parsing/hcl/writer.go`, `parsing/toml/toml_writer.go`, `parsing/yaml/yaml_reader.go`, `internal/cli/read_write_flag.go` (analyzed 2026-02-19)
