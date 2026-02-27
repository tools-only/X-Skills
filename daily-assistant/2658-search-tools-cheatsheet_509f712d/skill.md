---
title: "Search Tools Cheatsheet"
description: "Quick reference for choosing between rg, grepai, Serena and ast-grep"
tags: [cheatsheet, search, reference]
---

# Search Tools Cheatsheet (1-Page Reference)

> **Quick reference for choosing between rg, grepai, Serena & ast-grep**

---

## ⚡ Quick Decision (5 Seconds)

```
Know exact text?     → rg
Know exact name?     → rg or Serena
Know concept?        → grepai
Know structure?      → ast-grep
Need dependencies?   → grepai trace
```

---

## 📊 Speed vs Intelligence

```
Fast ←─────────────────────────────────→ Smart

rg        Serena       ast-grep      grepai
~20ms     ~100ms       ~200ms        ~500ms

Exact     Symbols      Structure     Meaning
```

---

## 🎯 Use Cases (Choose the Right Tool)

| Task | Tool | Command |
|------|------|---------|
| "Find TODO comments" | `rg` | `rg "TODO"` |
| "Find login function" | `rg`/`Serena` | `rg "login"` or `serena find_symbol` |
| "Find auth code" | `grepai` | `grepai search "authentication"` |
| "Who calls login?" | `grepai` | `grepai trace callers "login"` |
| "Get file structure" | `Serena` | `serena get_symbols_overview` |
| "Async without try/catch" | `ast-grep` | `ast-grep "async function $F"` |
| "React class → hooks" | `ast-grep` | Migration pattern |
| "Remember decision" | `Serena` | `serena write_memory` |

---

## 🔍 Tool Capabilities Matrix

| Feature | rg | grepai | Serena | ast-grep |
|---------|:--:|:------:|:------:|:--------:|
| **Exact match** | ✅ | ❌ | ✅ | ✅ |
| **Semantic** | ❌ | ✅ | ❌ | ❌ |
| **Call graph** | ❌ | ✅ | ❌ | ❌ |
| **Symbols** | ❌ | ❌ | ✅ | ❌ |
| **AST patterns** | ❌ | ❌ | ❌ | ✅ |
| **Memory** | ❌ | ❌ | ✅ | ❌ |
| **Setup** | ✅ | ⚠️ | ⚠️ | ⚠️ |

---

## 🚀 Combined Workflow (5 Steps)

**Example: Refactor authentication**

```bash
# 1. DISCOVER (grepai - semantic)
grepai search "authentication flow"

# 2. STRUCTURE (Serena - symbols)
serena get_symbols_overview --file auth.service.ts

# 3. DEPENDENCIES (grepai - call graph)
grepai trace callers "login"

# 4. PATTERNS (ast-grep - structure)
ast-grep "async function login"

# 5. VERIFY (rg - exact)
rg "validateSession" --type ts
```

---

## ⚠️ Common Mistakes

| ❌ Wrong | ✅ Right |
|---------|---------|
| `grepai search "createSession"` | `rg "createSession"` |
| `rg "auth.*login.*session"` | `grepai search "auth flow"` |
| `rg + sed` for refactoring | `Serena find_symbol` |
| Refactor without checking callers | `grepai trace callers` first |

---

## 📈 When to Escalate

```
Start: rg (90% of searches)
  ↓
  Need meaning? → grepai
  Need symbols? → Serena
  Need AST? → ast-grep
```

---

## 🎓 Full Guide

📖 [Search Tools Mastery](./workflows/search-tools-mastery.md) - Complete workflows, scenarios, benchmarks

---

**Print this page for quick reference** | Last updated: January 2026
