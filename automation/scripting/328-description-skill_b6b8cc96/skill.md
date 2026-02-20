---
description: "Dasel v3 query patterns for Hibernate .hbm.xml mapping files — entity-table binding, Java property-to-column extraction, one-to-many set/list/bag relationship tracing, many-to-one foreign key discovery, batch scanning across 60+ HBM files. Use when querying Hibernate ORM class mappings, extracting schema metadata from Java persistence layer, or auditing entity-column relationships in enterprise legacy codebases."
---

# Hibernate HBM Mapping Queries

<when_to_use>

Load this skill when querying Hibernate `.hbm.xml` mapping files — extracting entity-table bindings, property-column mappings, collection relationships (set/list/bag), foreign key discovery, or running batch scans across an enterprise persistence layer with 60+ HBM files.

</when_to_use>

Domain skill for dasel v3 queries against Hibernate `.hbm.xml` mapping files.

**Attribute syntax (required):** XML attributes use `-` prefix in dasel friendly mode (default). `name` → `-name`, `table` → `-table`, `column` → `-column`.

**Parser flag (required):** Always pass `-i xml` explicitly. Do not rely on auto-detection for `.hbm.xml` files.

## Entity Discovery

```bash
# Fully qualified Java entity class name
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.-name'

# Database table the entity maps to
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.-table'
```

## Property Column Mapping

```bash
# All mapped Java property names in the entity
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.property.map(-name)'

# All database column names for mapped properties
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.property.map(-column)'
```

## Collection Relationships

`<set>`, `<list>`, `<bag>` — one-to-many relationships.

```bash
# All set collection names (one-to-many)
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.set.map(-name)'

# All list collection names
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.list.map(-name)'
```

## Foreign Key Discovery

`<many-to-one>` — foreign key columns pointing to other entities.

```bash
# All foreign key column names
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.many-to-one.map(-column)'

# All referenced entity class names (the FK target)
dasel -f User.hbm.xml -i xml 'hibernate-mapping.class.many-to-one.map(-class)'
```

## Batch Entity-Table Mapping

Extract entity→table pairs across all `.hbm.xml` files. Write to `/tmp/` — never to the source tree.

```bash
for f in $(fdfind -e hbm.xml .); do
  entity=$(dasel -f "$f" -i xml 'hibernate-mapping.class.-name' 2>/dev/null)
  table=$(dasel -f "$f" -i xml 'hibernate-mapping.class.-table' 2>/dev/null)
  echo "$entity -> $table"
done > /tmp/entity_table_map.txt
```

Malformed files and missing attributes produce empty values via `2>/dev/null` suppression.
