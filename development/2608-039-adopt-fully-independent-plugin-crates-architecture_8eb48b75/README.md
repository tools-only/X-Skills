# 039 Adopt Fully Independent Plugin Crates Architecture

| Property | Value |
|----------|-------|
| **Name** | 039 Adopt Fully Independent Plugin Crates Architecture |
| **Repository** | [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/039-adopt-fully-independent-plugin-crates-architecture.md) (⭐ 3.3k) |
| **Original Path** | `docs/docs/architecture/adr/039-adopt-fully-independent-plugin-crates-architecture.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-02-14 |
| **Updated** | 2026-02-14 |
| **File Hash** | `8eb48b7559548132...` |

## Description

The current pii_filter plugin is not a separate crate and embeds PyO3 dependencies and macros directly. This couples plugin logic to Python bindings, making it difficult to add new plugins and increasing longterm maintenance costs as we expand support for both Rust and Python implementations.

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [IBM/mcp-context-forge](https://raw.githubusercontent.com/IBM/mcp-context-forge/main/docs/docs/architecture/adr/039-adopt-fully-independent-plugin-crates-architecture.md)*
