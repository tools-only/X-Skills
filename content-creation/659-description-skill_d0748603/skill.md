---
description: "Dasel v3 query patterns for InstallAnywhere .iap_xml installer definitions — use when querying action sequences, discovering variables, resolving platform conditions, navigating panels, or comparing installer variants. Files are 2.5+ MB, 65,000+ lines — too large for context reads, requires structural dasel queries."
---

# InstallAnywhere .iap_xml Query Patterns

<when_to_use>

Load this skill when working with any `.iap_xml` file — auditing install actions, listing file copy operations, discovering variables, finding platform-specific branches, or comparing two installer variants (e.g., RA vs NON_RA builds).

</when_to_use>

## Operational Constraints

- Always pass `-i xml` explicitly — auto-detection fails on `.iap_xml` extensions
- Write all intermediate results to `/tmp/` — never to the source tree
- Files are 2.5–2.6 MB, 65,000+ lines — use structural dasel queries, never read the file into context

## Root Structure

InstallAnywhere `.iap_xml` files have a fixed root path: `InstallAnywhere.project`.

```bash
# Top-level element names at the document root
dasel -f installer.iap_xml -i xml 'keys($this)'

# Navigate to the install action sequence
dasel -f installer.iap_xml -i xml 'InstallAnywhere.project.installSequence'

# Children of the project node
dasel -f installer.iap_xml -i xml 'InstallAnywhere.project.keys($this)'
```

Recursive descent — finds elements at any depth without knowing the full path:

```bash
# All elements named "action" at any depth
dasel -f installer.iap_xml -i xml '..action'

# All elements named "panel" at any depth
dasel -f installer.iap_xml -i xml '..panel'

# Any node at any depth where a specific attribute exists
dasel -f installer.iap_xml -i xml 'search(has("-name"))'
```

`..elementName` descends recursively. `search(predicate)` matches nodes at any depth.

## Action Sequence Analysis

```bash
# All action names in execution order
dasel -f installer.iap_xml -i xml '..action.map("-name")'

# All action types
dasel -f installer.iap_xml -i xml '..action.map("-type")'

# Count all action elements
dasel -f installer.iap_xml -i xml 'len(..action)'

# Write for analysis
dasel -f installer.iap_xml -i xml '..action.map("-name")' > /tmp/action_names.txt
dasel -f installer.iap_xml -i xml '..action.map("-type")' > /tmp/action_types.txt
```

## Variable Discovery

```bash
# All variable definitions
dasel -f installer.iap_xml -i xml '..variable.map("-name")' > /tmp/installer_vars.txt
```

## Platform Conditions

```bash
# All platform conditions (Windows vs Linux branches)
dasel -f installer.iap_xml -i xml 'search(has("-platform")).map("-platform")' > /tmp/platforms.txt
```

## Cross-Variant Comparison

When comparing two builds of the same installer (e.g., RA vs NON_RA, Windows vs Linux variant):

```bash
# Compare action sets between two installer files
dasel -f installer_a.iap_xml -i xml '..action.map("-name")' | sort > /tmp/actions_a.txt
dasel -f installer_b.iap_xml -i xml '..action.map("-name")' | sort > /tmp/actions_b.txt
diff /tmp/actions_a.txt /tmp/actions_b.txt > /tmp/installer_diff.txt

# Compare panel sets
dasel -f installer_a.iap_xml -i xml '..panel.map("-name")' | sort > /tmp/panels_a.txt
dasel -f installer_b.iap_xml -i xml '..panel.map("-name")' | sort > /tmp/panels_b.txt
diff /tmp/panels_a.txt /tmp/panels_b.txt
```

## File Copy Operations

```bash
# All file copy sources
dasel -f installer.iap_xml -i xml '..installFile.map("-source")' > /tmp/file_copies.txt
```

## Panel Navigation

```bash
# Panel names and count
dasel -f installer.iap_xml -i xml '..panel.map("-name")'
dasel -f installer.iap_xml -i xml 'len(..panel)'
```

## Attribute Syntax Reference

XML attributes use `-` prefix in dasel friendly mode:

```bash
# Access specific attribute on first match
dasel -f installer.iap_xml -i xml '..action[0].-name'

# Discover attributes and children of a node
dasel -f installer.iap_xml -i xml '..action[0].keys($this)'

# Extract attribute values across all matching nodes
dasel -f installer.iap_xml -i xml '..elementName.map("-attributeName")'

# Count matching elements
dasel -f installer.iap_xml -i xml 'len(..elementName)'

# Search with attribute predicate at any depth
dasel -f installer.iap_xml -i xml 'search(has("-attributeName"))'
```
