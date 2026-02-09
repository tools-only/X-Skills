# Search Tool Preference

Use Exa MCP tools for all web search:
- `web_search_exa` — general web search
- `get_code_context_exa` — code, GitHub, docs, Stack Overflow
- `company_research_exa` — company/vendor info

Do not use the built-in WebSearch tool. If Exa tools are not visible, discover them with `ToolSearch(query: 'exa')`.

# Documentation Search (QMD) — Always Use First

**QMD is the primary way to find documentation.** Do NOT manually read `docs/index.md`, `TECHNICAL_OVERVIEW.md`, or subdirectory CLAUDE.md files for doc lookup — use QMD search instead.

```bash
# Search for docs relevant to your task
qmd_search "authentication flow"

# Get a specific document by path
qmd_get "motium/cortex/docs/architecture/authentication.md"
```

**Why QMD over manual reads:**
- Searches across all indexed documentation in one query
- No need to read index.md first — search finds the right doc directly
- Token-efficient — returns excerpts, not full files
- `.gitignore` respected — no node_modules pollution

**Only use Read tool for:** `CLAUDE.md` (root) and `.claude/MEMORIES.md` at session start.
For everything else, search first with QMD.

If QMD tools are not available, fall back to reading `docs/index.md` manually.
