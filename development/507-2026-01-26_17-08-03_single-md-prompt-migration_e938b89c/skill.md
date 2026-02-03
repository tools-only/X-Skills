# Research - Single MD Prompt Migration

**Date:** 2026-01-26
**Owner:** agent
**Phase:** Research

## Goal

Research migrating from the complex multi-section XML prompt system back to a single MD file for the agent system prompt, while keeping tool prompts separate (they are already independent).

**IMPORTANT:** Local mode is being deleted in this change, so only ONE prompt file is needed (`system_prompt.md`), not two.

---

## Additional Search

```bash
grep -ri "section" /home/tuna/tunacode/.claude/ 2>/dev/null | head -20
```

---

## Findings

### Current Architecture Summary

Tunacode uses a **two-layer prompt composition system**:

1. **Agent System Prompt**: Built from 11 modular XML sections composed via templates
2. **Tool Prompts**: Separate XML files loaded independently via `@base_tool` decorator

**Tool prompts are already separate** from agent system prompts and require NO changes.

### Current Prompt Flow

```
load_system_prompt() [agent_config.py:182]
│
├──> SectionLoader.__init__(sections_dir) [loader.py:17]
│
├──> loader.load_section(each SystemPromptSection) [loader.py:21]
│    └──> Check {section_name}.xml/.md/.txt in prompts/sections/
│         └──> _read_with_cache(path) with mtime caching
│
├──> compose_prompt(template, sections) [prompting_engine.py:82]
│    └──> Replace {{SECTION_NAME}} with section content
│
└──> resolve_prompt(prompt) [prompting_engine.py:77]
     └──> Replace {{CWD}}, {{OS}}, {{DATE}} with runtime values
```

### Files to DELETE (6 files + local_prompt.md)

| File | Lines | Purpose |
|------|-------|---------|
| `src/tunacode/core/prompting/loader.py` | 67 | `SectionLoader` class - loads section files |
| `src/tunacode/core/prompting/sections.py` | 51 | `SystemPromptSection` enum - defines 12 sections |
| `src/tunacode/core/prompting/templates.py` | 86 | Template definitions (MAIN, LOCAL, RESEARCH) |
| `src/tunacode/prompts/sections/` | dir | 11 XML section files (~26,710 chars total) |
| `src/tunacode/core/prompting/local_prompt.md` | 72 | Local mode prompt (local mode being deleted) |
| `src/tunacode/prompts/sections/agent_role.xml` | 9 | Agent identity |
| `src/tunacode/prompts/sections/critical_rules.xml` | 36 | Behavior rules |
| `src/tunacode/prompts/sections/tool_use.xml` | 81 | Tool access rules |
| `src/tunacode/prompts/sections/search_pattern.xml` | 101 | GLOB->GREP->READ funnel |
| `src/tunacode/prompts/sections/completion.xml` | 10 | Submit tool signaling |
| `src/tunacode/prompts/sections/parallel_exec.xml` | 104 | Parallel execution rules |
| `src/tunacode/prompts/sections/output_style.xml` | 94 | Communication style |
| `src/tunacode/prompts/sections/examples.xml` | 186 | Few-shot examples |
| `src/tunacode/prompts/sections/advanced_patterns.xml` | 217 | Advanced workflows |
| `src/tunacode/prompts/sections/system_info.xml` | 7 | Dynamic placeholders |
| `src/tunacode/prompts/sections/user_instructions.xml` | 4 | User context placeholder |

### Files to MODIFY (Core Layer - 3 files)

| File | Lines | Changes |
|------|-------|---------|
| `src/tunacode/core/prompting/__init__.py` | 30 | Remove exports for deleted modules |
| `src/tunacode/core/prompting/prompting_engine.py` | 99 | Keep `resolve_prompt()` for {{CWD}}, {{OS}}, {{DATE}} |
| `src/tunacode/core/agents/agent_components/agent_config.py` | 475 | Replace `load_system_prompt()` with simple file read |

### Files to CREATE (1 file)

| File | Purpose |
|------|---------|
| `src/tunacode/prompts/system_prompt.md` | Single MD file with all agent prompt content |

### Files to MODIFY (Documentation - 13 files)

| File | Lines | Impact |
|------|-------|--------|
| `docs/codebase-map/modules/core-prompting.md` | 112 | Rewrite - section composition removed |
| `docs/codebase-map/modules/prompts.md` | 145 | Rewrite - modular sections removed |
| `docs/codebase-map/modules/templates.md` | 113 | Delete or rewrite - templates deprecated |
| `docs/codebase-map/architecture/architecture.md` | 426+ | Update - modular prompts references |
| `docs/codebase-map/MAP.md` | 756 | Update - SectionLoader references |
| `docs/codebase-map/modules/core-agents.md` | 154 | Update - template selection logic |
| `docs/codebase-map/state/state.md` | 80 | Update - SectionLoader caching note |
| `docs/codebase-map/modules/utils-limits.md` | 24 | DELETE or update - local mode being removed |
| `docs/codebase-map/modules/INDEX.md` | 107 | Update - module descriptions |
| `docs/codebase-map/modules/00-overview.md` | 38 | Update - prompts/ directory reference |
| `docs/codebase-map/structure/tree-structure.txt` | 163 | Update - directory structure |
| `docs/codebase-map/modules/tools-overview.md` | 230 | Update - XML prompt loading |
| `.claude/JOURNAL.md` | 161 | Update - LOCAL_TEMPLATE references |

### Tests (Mostly Stable - 1 file)

| File | Status |
|------|--------|
| `tests/unit/core/test_prompting_engine.py` | **KEEP** - tests {{placeholder}} resolution which is still needed |

---

## Tool Prompts - NO CHANGES NEEDED

Tool prompts are **already separate** from agent system prompts:

### Tool Prompt Flow (Independent - No Changes)

```
Tool function decorated with @base_tool [decorators.py:89]
│
├──> load_prompt_from_xml(func.__name__) [decorators.py:119]
│    └──> Load tools/prompts/{func_name}_prompt.xml
│         └──> Parse <description> tag, return text
│
└──> wrapper.__doc__ = xml_prompt [decorators.py:121]
     └──> pydantic-ai reads __doc__ when Tool() wraps function
```

### Tool Prompt Files (Keep Unchanged)

All files in `src/tunacode/tools/prompts/` remain unchanged:
- `bash_prompt.xml`
- `glob_prompt.xml`
- `grep_prompt.xml`
- `list_dir_prompt.xml`
- `read_file_prompt.xml`
- `write_file_prompt.xml`
- `update_file_prompt.xml`
- `web_fetch_prompt.xml`
- `submit_prompt.xml`

### Related Files (No Changes)

- `src/tunacode/tools/xml_helper.py` - `load_prompt_from_xml()` function
- `src/tunacode/tools/decorators.py` - `@base_tool` decorator (line 119-121)

---

## Key Patterns / Solutions Found

### Reference Implementation: `local_prompt.md` (BEING DELETED)

The file `src/tunacode/core/prompting/local_prompt.md` (72 lines) was the **reference implementation** for the single-file approach:

```python
# Current local mode loading (agent_config.py:228-232) - BEING REMOVED
prompting_dir = Path(__file__).parent.parent.parent / "prompting"
tunacode_path = prompting_dir / "local_prompt.md"
tunacode_content = tunacode_path.read_text(encoding="utf-8")
```

This demonstrates the pattern to follow:
- Single MD file with all content inline
- Uses `{{USER_INSTRUCTIONS}}` placeholder for context injection
- Direct file loading, no SectionLoader, no compose_prompt
- Simple `Path.read_text()` approach

**NOTE:** `local_prompt.md` itself will be deleted as part of local mode removal, but it serves as the structural reference for `system_prompt.md`.

### Proposed New `load_system_prompt()`

```python
def load_system_prompt(base_path: Path, model: str | None = None) -> str:
    """Load the system prompt from a single MD file.

    Local mode is being deleted - only one prompt file needed.

    Args:
        base_path: Base path to the tunacode package
        model: Optional model name (for future per-model prompts)

    Returns:
        System prompt with dynamic placeholders resolved
    """
    prompt_file = base_path / "prompts" / "system_prompt.md"

    if not prompt_file.exists():
        raise FileNotFoundError(
            f"Required prompt file not found: {prompt_file}"
        )

    # Read single file and resolve dynamic placeholders
    content = prompt_file.read_text(encoding="utf-8")
    return resolve_prompt(content)  # Handles {{CWD}}, {{OS}}, {{DATE}}
```

**NOTE:** Local mode is being deleted - no conditional logic needed. Just ONE prompt file.

### Placeholder Resolution - KEEP

The `{{CWD}}`, `{{OS}}`, `{{DATE}}` placeholders are **still needed** and should be preserved:

```python
# prompting_engine.py - KEEP THIS
def resolve_prompt(template: str) -> str:
    """Resolve {{CWD}}, {{OS}}, {{DATE}} placeholders."""
    # ... existing implementation ...
```

### What Goes Into `system_prompt.md`

Content from the 11 section files, in this order (from `MAIN_TEMPLATE`):

1. `## AGENT_ROLE` - Content from `agent_role.xml`
2. `## SEARCH_PATTERN` - Content from `search_pattern.xml`
3. `## CRITICAL_RULES` - Content from `critical_rules.xml`
4. `## TOOL_USE` - Content from `tool_use.xml`
5. `## COMPLETION` - Content from `completion.xml`
6. `## PARALLEL_EXEC` - Content from `parallel_exec.xml`
7. `## OUTPUT_STYLE` - Content from `output_style.xml`
8. `## EXAMPLES` - Content from `examples.xml`
9. `## ADVANCED_PATTERNS` - Content from `advanced_patterns.xml`
10. `## SYSTEM_INFO` - Content from `system_info.xml`
11. `{{USER_INSTRUCTIONS}}` - Placeholder for AGENTS.md context

---

## Knowledge Gaps

1. **RESEARCH_TEMPLATE deletion**: Current research shows `RESEARCH_TEMPLATE` is being removed (ticket tun-5a52). Need to verify this is complete before single MD migration.

2. **Template override system**: `TEMPLATE_OVERRIDES` dict allows model-specific templates. Need to decide if this is kept with single MD approach (would require separate files per model).

3. **Section dependencies**: Verify no external code depends on specific section names (e.g., tests that mock sections).

4. **Migration verification**: Need to confirm the composed single MD produces identical output to current section-based composition.

---

## Deletion Strategy (4 Phases)

### Phase 1: Create Single MD File
- Create `src/tunacode/prompts/system_prompt.md` with all 11 sections
- Verify content matches current composed output
- Keep section files for now (A/B testing)

### Phase 2: Modify Core Loading
- Update `load_system_prompt()` to read single file
- Remove `SectionLoader` import and usage
- Keep `resolve_prompt()` for dynamic placeholders

### Phase 3: Delete Unused Files
- Delete `src/tunacode/core/prompting/loader.py`
- Delete `src/tunacode/core/prompting/sections.py`
- Delete `src/tunacode/core/prompting/templates.py`
- Delete `src/tunacode/prompts/sections/` directory
- Update `src/tunacode/core/prompting/__init__.py` exports

### Phase 4: Update Documentation
- Rewrite `core-prompting.md` and `prompts.md`
- Update all 13 documentation files
- Update `tree-structure.txt`

---

## References

### Core Implementation Files
- `src/tunacode/core/prompting/__init__.py:1-30` - Public API exports
- `src/tunacode/core/prompting/loader.py:1-67` - SectionLoader class
- `src/tunacode/core/prompting/sections.py:1-51` - SystemPromptSection enum
- `src/tunacode/core/prompting/templates.py:1-86` - Template definitions
- `src/tunacode/core/prompting/prompting_engine.py:82-98` - compose_prompt function
- `src/tunacode/core/prompting/prompting_engine.py:77-79` - resolve_prompt function
- `src/tunacode/core/agents/agent_components/agent_config.py:182-217` - load_system_prompt function
- `src/tunacode/core/prompting/local_prompt.md:1-72` - Reference implementation (BEING DELETED with local mode)

### Tool Prompt Files (Separate - No Changes)
- `src/tunacode/tools/xml_helper.py:10-34` - load_prompt_from_xml function
- `src/tunacode/tools/decorators.py:119-121` - Tool prompt attachment

### Test Files
- `tests/unit/core/test_prompting_engine.py:1-90` - Placeholder resolution tests

### Documentation Files
- `docs/codebase-map/modules/core-prompting.md` - Core module documentation
- `docs/codebase-map/modules/prompts.md` - Prompt system overview
- `docs/codebase-map/modules/templates.md` - Template documentation

### GitHub Permalinks
- [agent_config.py:182-217](https://github.com/alchemiststudiosDOTai/tunacode/blob/b6eb3990926b44c95b5314f866ca856bc737a12f/src/tunacode/core/agents/agent_components/agent_config.py#L182-L217)
- [loader.py:1-67](https://github.com/alchemiststudiosDOTai/tunacode/blob/b6eb3990926b44c95b5314f866ca856bc737a12f/src/tunacode/core/prompting/loader.py)
- [sections.py:1-51](https://github.com/alchemiststudiosDOTai/tunacode/blob/b6eb3990926b44c95b5314f866ca856bc737a12f/src/tunacode/core/prompting/sections.py)
- [templates.py:1-86](https://github.com/alchemiststudiosDOTai/tunacode/blob/b6eb3990926b44c95b5314f866ca856bc737a12f/src/tunacode/core/prompting/templates.py)
