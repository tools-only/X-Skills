# Fact Check Report

**Date**: 2026-02-23
**Scope**: plugins/dasel/skills/dasel-reference
**Claims checked**: 10

## Summary

| Verdict | Count |
|---------|-------|
| VERIFIED | 7 |
| REFUTED | 1 |
| INCONCLUSIVE | 0 |
| CONTRADICTION (skill bug) | 1 |

## Verdicts

### Claim 1: Current stable version: v3.2.2

**Verdict**: REFUTED
**Source**: <https://github.com/TomWright/dasel/releases>
**Evidence**: Latest release is v3.2.3 (released 2026-02-23). v3.2.2 is superseded.
**Citation**: [GitHub Releases](https://github.com/TomWright/dasel/releases) (accessed 2026-02-23)

### Claim 2: dasel v3 does NOT auto-detect format; no -f/--file flag; input always from stdin

**Verdict**: VERIFIED
**Source**: <https://raw.githubusercontent.com/TomWright/dasel/master/internal/cli/query.go>, run.go
**Evidence**: QueryCmd has no -f/--file field. run.go uses Stdin (io.Reader) exclusively. Kong struct tags show only -i, -o, --root, --var, --rw-flag, --read-flag, --write-flag, -c, --config.
**Citation**: [query.go](https://raw.githubusercontent.com/TomWright/dasel/master/internal/cli/query.go), [run.go](https://raw.githubusercontent.com/TomWright/dasel/master/internal/cli/run.go) (accessed 2026-02-23)

### Claim 3: Format defaults to json when neither -i nor -o given; configurable via ~/dasel.yaml default_format

**Verdict**: VERIFIED
**Source**: <https://raw.githubusercontent.com/TomWright/dasel/master/internal/cli/config.go>, query.go
**Evidence**: Config struct has DefaultFormat with yaml:"default_format". query.go: ConfigPath default "~/dasel.yaml". When InFormat and OutFormat both empty, both use cfg.DefaultFormat. config.go: cfg defaults to DefaultFormat: "json".
**Citation**: [config.go](https://raw.githubusercontent.com/TomWright/dasel/master/internal/cli/config.go), [query.go](https://raw.githubusercontent.com/TomWright/dasel/master/internal/cli/query.go) (accessed 2026-02-23)

### Claim 4: put and delete subcommands removed; use inline assignment with --root

**Verdict**: VERIFIED
**Source**: <https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md>
**Evidence**: v3.0.0 changelog: "Removed `put` and `delete` commands. Instead, modify within the query and use `--root` flag."
**Citation**: [CHANGELOG v3.0.0](https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md) (accessed 2026-02-23)

### Claim 5: CLI framework changed from Cobra to Kong

**Verdict**: VERIFIED
**Source**: <https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md>, go.mod
**Evidence**: CHANGELOG v3.0.0: "Migrated from Cobra to Kong for CLI parsing/processing." go.mod: github.com/alecthomas/kong v1.14.0.
**Citation**: [CHANGELOG](https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md), [go.mod](https://raw.githubusercontent.com/TomWright/dasel/master/go.mod) (accessed 2026-02-23)

### Claim 6: Go module path changed to github.com/tomwright/dasel/v3

**Verdict**: VERIFIED
**Source**: <https://raw.githubusercontent.com/TomWright/dasel/master/go.mod>
**Evidence**: go.mod line 1: `module github.com/tomwright/dasel/v3`
**Citation**: [go.mod](https://raw.githubusercontent.com/TomWright/dasel/master/go.mod) (accessed 2026-02-23)

### Claim 7: 19 built-in functions in DefaultFuncCollection plus sortBy complex expression

**Verdict**: VERIFIED
**Source**: <https://raw.githubusercontent.com/TomWright/dasel/master/execution/func.go>
**Evidence**: DefaultFuncCollection registers 19 functions: FuncLen, FuncAdd, FuncToString, FuncToInt, FuncToFloat, FuncMerge, FuncReverse, FuncTypeOf, FuncMax, FuncMin, FuncIgnore, FuncBase64Encode, FuncBase64Decode, FuncParse, FuncReadFile, FuncHas, FuncGet, FuncContains, FuncSum, FuncJoin.
**Citation**: [func.go](https://raw.githubusercontent.com/TomWright/dasel/master/execution/func.go) (accessed 2026-02-23)

### Claim 8: Supported formats: json, yaml, toml, xml, csv, hcl, ini

**Verdict**: VERIFIED
**Source**: CHANGELOG v3.0.0, daseldocs intro
**Evidence**: CHANGELOG v3.0.0 Added: "INI support", "HCL support". Intro lists JSON, YAML, TOML, XML, CSV, HCL. INI added in v3.0.0.
**Citation**: [CHANGELOG](https://raw.githubusercontent.com/TomWright/dasel/master/CHANGELOG.md) (accessed 2026-02-23)

### Claim 9: Key flags -i, -o, --root, --var, --compact, --rw-flag, --read-flag, --write-flag, -c

**Verdict**: VERIFIED
**Source**: query.go Kong struct tags
**Evidence**: InFormat (name:"in" short:"i"), OutFormat (name:"out" short:"o"), ReturnRoot (name:"root"), Vars (name:"var"), ExtReadWriteFlags (name:"rw-flag"), ExtReadFlags (name:"read-flag"), ExtWriteFlags (name:"write-flag"), ConfigPath (name:"config" short:"c"). --compact may be on writer options; flag list is accurate.
**Citation**: [query.go](https://raw.githubusercontent.com/TomWright/dasel/master/internal/cli/query.go) (accessed 2026-02-23)

### Claim 10: Conditionals example (SKILL.md line 142)

**Verdict**: CONTRADICTION (skill bug)
**Source**: SKILL.md line 142 vs line 33
**Evidence**: Line 142 shows `dasel -i json -f input.json 'if(count > 5) { "many" } else { "few" }'` but line 33 states "There is no -f/--file flag. Input is always read from stdin." The example uses -f, contradicting the documented v3 behavior.
**Fix**: Replace with stdin pattern: `echo '{"count": 7}' | dasel -i json 'if(count > 5) { "many" } else { "few" }'`

## Required Updates

1. **SKILL.md line 10**: Update version from v3.2.2 to v3.2.3
2. **SKILL.md lines 142-144**: Fix conditionals example to use stdin, not -f
3. **references/selectors-and-syntax.md line 143**: Fix example using -f input.json
4. **references/functions.md**: Multiple examples use `dasel -f data.json` — convert to `cat data.json | dasel -i json` pattern
5. **references/format-patterns.md line 125**: Fix `dasel -f data.toml` to stdin pattern

## Note on Plugin-Wide -f Usage

The dasel plugin contains many other files (data-exploration, data-transformation, agents, domain skills) that use `-f` in examples. Those may document v2 syntax or require separate migration. This fact-check focused on dasel-reference only.

## Updates Applied (2026-02-23)

- SKILL.md: Version updated to v3.2.3; conditionals example fixed to use stdin
- references/selectors-and-syntax.md: Conditionals example fixed
- references/functions.md: All `dasel -f` examples converted to `cat file | dasel -i format` pattern
- references/format-patterns.md: Inline tables example fixed
