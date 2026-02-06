# Memory Architecture

| Property | Value |
|----------|-------|
| **Name** | Memory Architecture |
| **Repository** | [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/docs/memory_architecture.md) (⭐ 24) |
| **Original Path** | `docs/memory_architecture.md` |
| **Category** | data-analysis |
| **Subcategory** | visualization |
| **Tags** | data analysis |
| **Created** | 2026-01-31 |
| **Updated** | 2026-02-05 |
| **File Hash** | `095a8125806cc74a...` |

## Description

| 存储 | 用途 | 更新频率 | 注入方式 |
|||||
| ChromaDB | 向量索引，语义搜索 | 实时（有新记忆就索引） | 按需搜索相关记忆 |
| memories.json | 完整记忆库 | 实时 | 作为 ChromaDB 的元数据 |
| MEMORY.md | 精华摘要（最重要的 1015 条） | 每日凌晨刷新 | 每次系统提示都带上 |

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [openakita/openakita](https://raw.githubusercontent.com/openakita/openakita/main/docs/memory_architecture.md)*
