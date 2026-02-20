# Dasel CLI Tool - Research Findings

## What Dasel Is and What Problem It Solves

Dasel (Data-Select) is a single-binary CLI tool and Go library for querying, modifying, and transforming structured data files. It solves the problem of needing separate tools (`jq` for JSON, `yq` for YAML, `xmllint` for XML) by providing one unified query syntax that works across all supported formats.

Primary use cases:
- Query values from structured config/data files in shell scripts
- Modify structured data files in-place without full-language scripts
- Convert between data formats (JSON to YAML, TOML to JSON, etc.)
- Recursive search through nested documents

Current stable version: **v3.2.2** (released 2026-02-13). v3 is a major rewrite with revised syntax.

SOURCES:
- README: <https://raw.githubusercontent.com/TomWright/dasel/master/README.md> (fetched 2026-02-19)
- Docs: <https://daseldocs.tomwright.me> (fetched 2026-02-19)
- Releases: <https://github.com/TomWright/dasel/releases> (fetched 2026-02-19)
- CHANGELOG: <https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md> (fetched 2026-02-19)

---

## Supported Data Formats

Confirmed from `cmd/dasel/main.go` import list and documentation:

- JSON (`json`)
- YAML (`yaml`)
- TOML (`toml`)
- XML (`xml`)
- CSV (`csv`)
- HCL (`hcl`) — added in v3.0.0
- INI (`ini`) — added in v3.0.0

Format identifiers (used with `-i`/`-o` flags): `json`, `yaml`, `toml`, `xml`, `csv`, `hcl`, `ini`

---

## Binary Release URL Pattern

Base URL pattern:

```text
https://github.com/TomWright/dasel/releases/download/v{VERSION}/dasel_{platform}
```

### Available Platforms (v3.2.2, 18 assets total)

| Asset Name | Platform | Notes |
|---|---|---|
| `dasel_linux_amd64` | Linux x86_64 | Primary Linux target |
| `dasel_linux_amd64.gz` | Linux x86_64 compressed | |
| `dasel_linux_arm64` | Linux ARM64 | |
| `dasel_linux_arm64.gz` | Linux ARM64 compressed | |
| `dasel_linux_arm32` | Linux ARM32 | |
| `dasel_linux_arm32.gz` | Linux ARM32 compressed | |
| `dasel_linux_386` | Linux x86 32-bit | |
| `dasel_linux_386.gz` | Linux x86 32-bit compressed | |
| `dasel_darwin_amd64` | macOS Intel | |
| `dasel_darwin_amd64.gz` | macOS Intel compressed | |
| `dasel_darwin_arm64` | macOS Apple Silicon | |
| `dasel_darwin_arm64.gz` | macOS Apple Silicon compressed | |
| `dasel_windows_amd64.exe` | Windows x64 | |
| `dasel_windows_amd64.exe.gz` | Windows x64 compressed | |
| `dasel_windows_386.exe` | Windows 32-bit | |
| `dasel_windows_386.exe.gz` | Windows 32-bit compressed | |

### Concrete Examples (v3.2.2)

```text
https://github.com/TomWright/dasel/releases/download/v3.2.2/dasel_linux_amd64
https://github.com/TomWright/dasel/releases/download/v3.2.2/dasel_linux_arm64
https://github.com/TomWright/dasel/releases/download/v3.2.2/dasel_windows_amd64.exe
https://github.com/TomWright/dasel/releases/download/v3.2.2/dasel_darwin_arm64
```

---

## Checksum Verification

Checksums are embedded directly in the GitHub release page per-asset. There is **no separate SHA256SUMS file** in releases. Each asset listing shows its sha256 hash inline.

Observed sha256 values for v3.2.2:

```text
dasel_linux_amd64:    313f055bd8c7c59bb8a9cad863f0685769d6f8a93b5f9b9a0986c95c76513832
dasel_linux_arm64:    135ac22a599fca49d48af44bd5dfa83bccdd91cd813c34c018e3574f76e842e6
dasel_linux_arm32:    6c43ff668c8205465e5509f94f39095d1ba86275e2e9d2d85976abe4c0cf1083
dasel_linux_386:      06c4024d39494ea46bff086ff3622cb2e8a722c3f68b604fa0cf6cf03365e887
dasel_darwin_amd64:   cb064c256cc33e90ddc0630935bf6c6d640a34cf07eed5a47fdd365c13514239
dasel_darwin_arm64:   03f90053514080d57ed29b8558829e7e51f5ef8e62f25aa0029f67d2a9050e39
dasel_windows_amd64:  bf5ccd37700525f7b55592064cb74fa7fae3d86b76f8a7b7fc7ef2c237c9c2c0
dasel_windows_386:    299ae7e3c889f13f9e26ee323c50cb6a3395284b7d625c9c2fc43e1695d1c87f
```

To verify after download:

```bash
echo "313f055bd8c7c59bb8a9cad863f0685769d6f8a93b5f9b9a0986c95c76513832  dasel" | sha256sum -c
```

Checksums per release can be fetched dynamically via GitHub API:

```bash
curl -s https://api.github.com/repos/tomwright/dasel/releases/latest \
  | grep -E '"browser_download_url"|"sha256"'
```

---

## Version Detection Command

```bash
dasel --version
```

To get the latest available version from GitHub API:

```bash
curl -s https://api.github.com/repos/tomwright/dasel/releases/latest \
  | grep '"tag_name"' \
  | cut -d'"' -f4
```

---

## Install Instructions

### Linux (64-bit amd64) — Auto-detect latest

```bash
curl -sSLf "$(curl -sSLf https://api.github.com/repos/tomwright/dasel/releases/latest \
  | grep browser_download_url \
  | grep linux_amd64 \
  | grep -v .gz \
  | cut -d'"' -f4)" -L -o dasel && chmod +x dasel
mv ./dasel /usr/local/bin/dasel
```

### Linux (manual, pinned version)

```bash
VERSION=v3.2.2
ARCH=linux_amd64   # or linux_arm64, linux_arm32, linux_386
curl -sSLf "https://github.com/TomWright/dasel/releases/download/${VERSION}/dasel_${ARCH}" \
  -o /usr/local/bin/dasel
chmod +x /usr/local/bin/dasel
```

### Windows (PowerShell)

```powershell
$releases = curl -sSLf https://api.github.com/repos/tomwright/dasel/releases/latest
Invoke-WebRequest -Uri (($releases | ConvertFrom-Json).assets `
                    | Where-Object { $_.name -eq "dasel_windows_amd64.exe" } `
                    | Select-Object -ExpandProperty browser_download_url) `
                    -OutFile dasel.exe
```

### Other install methods

- Homebrew (macOS/Linux): `brew install dasel`
- Docker: `docker run -i --rm ghcr.io/tomwright/dasel:latest`
- Go: `go install github.com/tomwright/dasel/v3/cmd/dasel@master`
- Scoop (Windows): `scoop bucket add extras && scoop install dasel`
- ASDF: `asdf plugin add dasel && asdf install dasel <version>`
- Mise: `mise install dasel@<version>`
- Nix: `nix-env -iA nixpkgs.dasel`

---

## Core CLI Syntax (v3)

### Flags

```text
dasel [flags] [QUERY]

Key flags:
  -f, --file <path>     Input file path (reads stdin if omitted)
  -i, --input <format>  Input parser: json, yaml, toml, xml, csv, hcl, ini
  -o, --output <format> Output parser (defaults to input format)
  --root                Output the full document after modification
  -w, --write <path>    Write result to file (modifies in-place if same as -f)
  --var <name>=<value>  Pass variables into query (repeatable)
  --compact             Compact output (no pretty-printing)
  --colour              Enable colourised output
  --indent <string>     Set indentation string
```

### Select (Query) — Read a Value

```bash
# From piped stdin
echo '{"foo": {"bar": "baz"}}' | dasel -i json 'foo.bar'
# Output: "baz"

# From file
dasel -f data.json 'foo.bar'
dasel -f config.yaml 'database.host'
dasel -f data.toml 'server.port'
```

### Modify in Place — Write a Value

In v3, `put` and `delete` subcommands were removed. Modifications use assignment expressions with `--root` to output the full document:

```bash
# Update a value and output full document
echo '{"foo": {"bar": "baz"}}' | dasel -i json --root 'foo.bar = "newvalue"'

# Modify file in place
dasel -f config.yaml --root 'server.port = 9090' > config_new.yaml

# Numeric assignment
echo '{"count": 1}' | dasel -i json --root 'count = 42'

# Boolean assignment
echo '{"enabled": false}' | dasel -i json --root 'enabled = true'
```

### Delete a Field

In v3, deletion is done by constructing a new object/map without the unwanted key:

```bash
# Remove a key using map spread (keeping all except the deleted key)
# Build new object without 'unwanted_key':
echo '{"keep": "yes", "remove": "no"}' | dasel -i json --root \
  '{ keep }'
```

### Array Operations

```bash
# Access by index (zero-based)
echo '[10, 20, 30]' | dasel -i json '[1]'
# Output: 20

# Access last element
echo '[10, 20, 30]' | dasel -i json '[len($this)-1]'

# Range slice [start:end]
echo '[1,2,3,4,5]' | dasel -i json '[0:2]'
# Output: [1, 2, 3]

# Append to array
echo '[1,2,3]' | dasel -i json --root '[$this..., 4]'
# Output: [1, 2, 3, 4]
```

### Batch Modify Array Elements — `each`

```bash
# Multiply all numbers by 2
echo '[1,2,3,4,5]' | dasel -i json --root 'each($this = $this*2)'
# Output: [2, 4, 6, 8, 10]

# Increment all integers
echo '[1,2,3]' | dasel -i json 'each($this = $this+1)'
```

### Format Conversion

```bash
# JSON to YAML
cat data.json | dasel -i json -o yaml

# YAML to TOML
cat config.yaml | dasel -i yaml -o toml

# JSON to XML
cat data.json | dasel -i json -o xml
```

### Recursive Descent (`..`)

```bash
# Find all values with key "name" at any depth
echo '{"user": {"name": "Alice"}}' | dasel -i json '..name'
# Output: ["Alice"]

# Get first element of every nested array
dasel -f data.json '..[0]'

# Get all values at any depth
dasel -f data.json '..*'
```

### Search Function

More powerful than `..` — supports arbitrary predicates:

```bash
# Find all objects that have a "name" key
dasel -f data.json 'search(has("name"))'

# Find nodes where value equals 42
dasel -f data.json 'search($this == 42)'

# Combined predicate
dasel -f data.json 'search(has("id") && has("name"))'
```

### Filter and Map

```bash
# Filter array where active == true
echo '...' | dasel -i json 'users.filter(active == true)'

# Map: extract single field from array of objects
echo '...' | dasel -i json 'users.map(name)'

# Map with transformation
echo '[1,2,3]' | dasel -i json '[1,2,3].map($this * 2)'
```

### Conditionals

```bash
# if/else expression
dasel -i json -f input.json '.foo.if(bar == "baz") { bong } else { qux }'

# Conditional with literal output
echo '{ "count": 7 }' | dasel -i json 'if(count > 5) { "many" } else { "few" }'
```

### Variable Assignment (v3 feature)

```bash
# Multi-statement query with variables
dasel -f data.json '$activeUsers = users.filter(active == true); $activeUsers.map(name)'
```

### Type Casting Functions

```bash
# typeOf: get type as string
echo '"hello"' | dasel -i json 'typeOf($this)'
# Output: "string"

# toString: convert to string
echo '123' | dasel -i json 'toString($this)'
# Output: "123"

# toInt: convert to integer
echo '"42"' | dasel -i json 'toInt($this)'
# Output: 42

# toFloat: convert to float
echo '"3.14"' | dasel -i json 'toFloat($this)'
```

Type strings returned by `typeOf`: `"string"`, `"array"`, `"bool"`, `"null"`, `"int"`, `"float"`

### Other Useful Functions

```bash
# len: length of array or string
echo '[1,2,3]' | dasel -i json 'len($this)'
# Output: 3

# has: check key/index exists (returns bool)
echo '{"foo": "bar"}' | dasel -i json 'has("foo")'
# Output: true

# join: join array elements into string
echo '["a","b","c"]' | dasel -i json 'join(",")'
# Output: "a,b,c"

# sum: sum numeric array
echo '[1,2,3,4,5]' | dasel -i json 'sum($this)'
# Output: 15

# keys: get map keys as array
echo '{"a": 1, "b": 2}' | dasel -i json 'keys($this)'

# reverse: reverse an array
echo '[1,2,3]' | dasel -i json 'reverse($this)'
```

---

## v3 Breaking Changes from v2 (Important)

From CHANGELOG v3.0.0:

- `put` and `delete` subcommands **removed**. Modifications done inline in query with `--root` flag.
- Query/selector syntax **revamped** — v2 syntax is not compatible with v3.
- Go module path changed to `github.com/tomwright/dasel/v3`.
- CLI framework changed from Cobra to Kong.

If a plugin needs to support v2 as well, the `put` command was used in v2:

```bash
# v2 syntax (NOT v3):
dasel put string -f data.json -o - '.foo.bar' 'newvalue'
dasel delete -f data.json '.foo.bar'
```

---

## Special Syntax Notes

### $root and $this Variables

- `$root` — references the root document
- `$this` — references the current node (used inside `each`, `map`, `filter`, `search`)

### Spread Operator

```bash
# Spread array into arguments
doSomething([1, 2, 3]...)  # equivalent to doSomething(1, 2, 3)

# Merge objects
{ {"firstName": "Tom"}..., "lastName": "Wright" }
# resolves to { "firstName": "Tom", "lastName": "Wright" }

# Output array elements as separate documents
[1, 2, 3]...
# 1
# 2
# 3
```

### Multi-statement Queries (semicolon separator)

```bash
# Statements separated by semicolons; final statement is output
dasel -f data.json '$active = users.filter(active == true); $active.map(name)'
```

### Comments in Queries

v3 supports comments in multi-statement queries (syntax not yet documented in official docs but noted in v3.0.0 release notes).

---

## References

1. README (raw): <https://raw.githubusercontent.com/TomWright/dasel/master/README.md> (fetched 2026-02-19)
2. Installation docs: <https://daseldocs.tomwright.me/getting-started/installation> (fetched 2026-02-19)
3. Releases page: <https://github.com/TomWright/dasel/releases> (fetched 2026-02-19)
4. Release v3.2.2 assets: <https://github.com/TomWright/dasel/releases/tag/v3.2.2> (fetched 2026-02-19)
5. CHANGELOG: <https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md> (fetched 2026-02-19)
6. Query syntax docs: <https://daseldocs.tomwright.me/syntax/query-syntax.md> (fetched 2026-02-19)
7. Types/Literals: <https://daseldocs.tomwright.me/syntax/types-literals.md> (fetched 2026-02-19)
8. Recursive descent: <https://daseldocs.tomwright.me/syntax/recursive-descent.md> (fetched 2026-02-19)
9. Arrays/slices: <https://daseldocs.tomwright.me/syntax/arrays-slices.md> (fetched 2026-02-19)
10. Objects/maps: <https://daseldocs.tomwright.me/syntax/objects-maps.md> (fetched 2026-02-19)
11. Conditionals: <https://daseldocs.tomwright.me/syntax/conditionals.md> (fetched 2026-02-19)
12. Spread: <https://daseldocs.tomwright.me/syntax/spread.md> (fetched 2026-02-19)
13. Functions index: <https://daseldocs.tomwright.me/functions> (fetched 2026-02-19)
14. search(): <https://daseldocs.tomwright.me/functions/search.md> (fetched 2026-02-19)
15. each(): <https://daseldocs.tomwright.me/functions/each.md> (fetched 2026-02-19)
16. filter(): <https://daseldocs.tomwright.me/functions/filter.md> (fetched 2026-02-19)
17. map(): <https://daseldocs.tomwright.me/functions/map.md> (fetched 2026-02-19)
18. typeOf(): <https://daseldocs.tomwright.me/functions/typeof.md> (fetched 2026-02-19)
19. toString(): <https://daseldocs.tomwright.me/functions/tostring.md> (fetched 2026-02-19)
20. toInt(): <https://daseldocs.tomwright.me/functions/toint.md> (fetched 2026-02-19)
21. has(): <https://daseldocs.tomwright.me/functions/has.md> (fetched 2026-02-19)
22. join(): <https://daseldocs.tomwright.me/functions/join.md> (fetched 2026-02-19)
23. main.go (format parsers): <https://github.com/TomWright/dasel/blob/master/cmd/dasel/main.go> (fetched 2026-02-19)
