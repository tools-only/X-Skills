---
title: "Single MD Prompt Migration – Plan"
phase: Plan
date: "2026-01-26T17:11:28"
owner: "agent"
parent_research: "memory-bank/research/2026-01-26_17-08-03_single-md-prompt-migration.md"
git_commit_at_plan: "b6eb3990"
tags: [plan, prompt-migration, coding]
---

## Goal

Migrate from the complex multi-section XML prompt composition system to a single `system_prompt.md` file for agent system prompts. Tool prompts remain separate and unchanged.

**Singular outcome:** One unified `prompts/system_prompt.md` file replaces 11 section files + 3 template modules + section loader logic.

**Non-goals:**
- No tool prompt changes (already independent)
- No new placeholder types
- No documentation rewrites beyond updating references to deleted files
- No per-model prompt variants (TEMPLATE_OVERRIDES deleted)

---

## Scope & Assumptions

### In Scope
- Deleting `src/tunacode/core/prompting/loader.py`
- Deleting `src/tunacode/core/prompting/sections.py`
- Deleting `src/tunacode/core/prompting/templates.py`
- Deleting `src/tunacode/prompts/sections/` directory (11 XML files)
- Deleting `src/tunacode/core/prompting/local_prompt.md` (local mode removal)
- Creating `src/tunacode/prompts/system_prompt.md`
- Modifying `agent_config.py` to load single file
- Updating `__init__.py` exports

### Out of Scope
- Tool prompts (already independent)
- New dynamic placeholder types
- Per-model prompt customization (TEMPLATE_OVERRIDES removed)
- Full documentation rewrites (minimal updates only)

### Assumptions
- Local mode is already being deleted (per ticket tun-5a52)
- RESEARCH_TEMPLATE already removed
- `resolve_prompt()` function remains for {{CWD}}, {{OS}}, {{DATE}} placeholders
- Python 3.11+ Path.read_text() available
- Existing test `test_prompting_engine.py` remains valid (tests placeholder resolution)

---

## Deliverables

### Source Code
1. `src/tunacode/prompts/system_prompt.md` - Single unified agent prompt
2. Modified `src/tunacode/core/agents/agent_components/agent_config.py` - Simplified `load_system_prompt()`
3. Modified `src/tunacode/core/prompting/__init__.py` - Removed exports

### Deletions
- `src/tunacode/core/prompting/loader.py`
- `src/tunacode/core/prompting/sections.py`
- `src/tunacode/core/prompting/templates.py`
- `src/tunacode/prompts/sections/` (entire directory)
- `src/tunacode/core/prompting/local_prompt.md`

### Documentation (minimal updates)
- Update `docs/codebase-map/modules/core-prompting.md` - remove section references
- Update `docs/codebase-map/modules/prompts.md` - single file reference
- Delete `docs/codebase-map/modules/templates.md` - no longer applicable

---

## Readiness

### Preconditions
- Current git state: `b6eb3990` (uncommitted changes exist in working directory)
- All 11 section files exist in `prompts/sections/`
- `agent_config.py` has existing `load_system_prompt()` function
- `resolve_prompt()` in `prompting_engine.py` is tested and working

### Dependencies
- None (pure Python stdlib)

### Data Schema
- Input: 11 XML section files in `prompts/sections/`
- Output: 1 MD file `prompts/system_prompt.md`

---

## Milestones

### M1: Single MD File Creation
Create `system_prompt.md` by composing content from all 11 section files in MAIN_TEMPLATE order.

### M2: Core Loading Logic Update
Modify `agent_config.py` to read single file instead of composing sections.

### M3: Cleanup Unused Code
Delete loader.py, sections.py, templates.py, prompts/sections/, local_prompt.md.

### M4: Module API Updates
Update `__init__.py` exports and minimal documentation fixes.

---

## Work Breakdown (Tasks)

### Task T1: Create system_prompt.md from Section Files
**Owner:** agent
**Estimate:** Small
**Dependencies:** None
**Milestone:** M1

**Description:**
Create `src/tunacode/prompts/system_prompt.md` by concatenating content from all 11 section files in MAIN_TEMPLATE order with `====` separators. Remove XML wrapper tags (`###Instruction###`, `<instructions>`, etc.) and keep only the inner content.

**Section order (from MAIN_TEMPLATE):**
1. `agent_role.xml` - Agent identity
2. `search_pattern.xml` - GLOB->GREP->READ funnel
3. `critical_rules.xml` - Behavior rules
4. `tool_use.xml` - Tool access rules
5. `completion.xml` - Submit tool signaling
6. `parallel_exec.xml` - Parallel execution rules
7. `output_style.xml` - Communication style
8. `examples.xml` - Few-shot examples
9. `advanced_patterns.xml` - Advanced workflows
10. `system_info.xml` - Dynamic placeholders
11. `{{USER_INSTRUCTIONS}}` - Placeholder (keep as-is)

**Acceptance Test:**
- File `src/tunacode/prompts/system_prompt.md` exists
- File contains all 11 sections in correct order with `====` separators
- `{{USER_INSTRUCTIONS}}` placeholder present
- No XML tags (`###Instruction###`, `<instructions>`, etc.) in output

**Files touched:**
- `src/tunacode/prompts/system_prompt.md` (CREATE)

---

### Task T2: Simplify load_system_prompt() in agent_config.py
**Owner:** agent
**Estimate:** Small
**Dependencies:** T1
**Milestone:** M2

**Description:**
Replace `load_system_prompt()` function with simple file read. Remove SectionLoader, compose_prompt, template selection, and local mode logic. Keep `resolve_prompt()` call for dynamic placeholders.

**New implementation:**
```python
def load_system_prompt(base_path: Path, model: str | None = None) -> str:
    """Load the system prompt from a single MD file.

    Local mode is being deleted - only one prompt file needed.

    Args:
        base_path: Base path to the tunacode package
        model: Optional model name (reserved for future use)

    Returns:
        System prompt with dynamic placeholders resolved
    """
    prompt_file = base_path / "prompts" / "system_prompt.md"

    if not prompt_file.exists():
        raise FileNotFoundError(
            f"Required prompt file not found: {prompt_file}"
        )

    content = prompt_file.read_text(encoding="utf-8")
    return resolve_prompt(content)
```

**Acceptance Test:**
- Function exists in `agent_config.py`
- No imports of SectionLoader, compose_prompt, templates, is_local_mode
- Function reads `prompts/system_prompt.md` directly
- Function calls `resolve_prompt()` on content
- Function raises FileNotFoundError if file missing

**Files touched:**
- `src/tunacode/core/agents/agent_components/agent_config.py` (MODIFY)

---

### Task T3: Update __init__.py Exports
**Owner:** agent
**Estimate:** Trivial
**Dependencies:** None
**Milestone:** M2

**Description:**
Remove exports for deleted modules from `src/tunacode/core/prompting/__init__.py`. Keep `resolve_prompt` export.

**Remove exports for:**
- `SectionLoader`
- `SystemPromptSection`
- `compose_prompt`
- `MAIN_TEMPLATE`, `LOCAL_TEMPLATE`, `RESEARCH_TEMPLATE`
- `TEMPLATE_OVERRIDES`

**Acceptance Test:**
- `__init__.py` exports only `resolve_prompt`
- No references to deleted modules

**Files touched:**
- `src/tunacode/core/prompting/__init__.py` (MODIFY)

---

### Task T4: Delete Unused Module Files
**Owner:** agent
**Estimate:** Trivial
**Dependencies:** T2, T3
**Milestone:** M3

**Description:**
Delete three module files: loader.py, sections.py, templates.py from `src/tunacode/core/prompting/`.

**Acceptance Test:**
- `src/tunacode/core/prompting/loader.py` deleted
- `src/tunacode/core/prompting/sections.py` deleted
- `src/tunacode/core/prompting/templates.py` deleted

**Files touched:**
- `src/tunacode/core/prompting/loader.py` (DELETE)
- `src/tunacode/core/prompting/sections.py` (DELETE)
- `src/tunacode/core/prompting/templates.py` (DELETE)

---

### Task T5: Delete prompts/sections/ Directory
**Owner:** agent
**Estimate:** Trivial
**Dependencies:** T1 (content migrated to system_prompt.md)
**Milestone:** M3

**Description:**
Delete entire `src/tunacode/prompts/sections/` directory containing 11 XML section files.

**Files to delete:**
- `agent_role.xml`
- `search_pattern.xml`
- `critical_rules.xml`
- `tool_use.xml`
- `completion.xml`
- `parallel_exec.xml`
- `output_style.xml`
- `examples.xml`
- `advanced_patterns.xml`
- `system_info.xml`
- `user_instructions.xml`

**Acceptance Test:**
- `src/tunacode/prompts/sections/` directory deleted
- All 11 XML files deleted

**Files touched:**
- `src/tunacode/prompts/sections/` (DELETE directory)

---

### Task T6: Delete local_prompt.md
**Owner:** agent
**Estimate:** Trivial
**Dependencies:** None
**Milestone:** M3

**Description:**
Delete `src/tunacode/core/prompting/local_prompt.md` as part of local mode removal.

**Acceptance Test:**
- `src/tunacode/core/prompting/local_prompt.md` deleted

**Files touched:**
- `src/tunacode/core/prompting/local_prompt.md` (DELETE)

---

### Task T7: Update core-prompting.md Documentation
**Owner:** agent
**Estimate:** Small
**Dependencies:** T4, T5, T6 (deletions complete)
**Milestone:** M4

**Description:**
Update `docs/codebase-map/modules/core-prompting.md` to reflect single-file architecture. Remove references to SectionLoader, sections, templates. Describe new simple file loading approach.

**Acceptance Test:**
- No references to SectionLoader, SystemPromptSection, compose_prompt
- Documents `load_system_prompt()` reading single file
- Notes that `resolve_prompt()` remains for placeholder resolution

**Files touched:**
- `docs/codebase-map/modules/core-prompting.md` (MODIFY)

---

### Task T8: Update prompts.md Documentation
**Owner:** agent
**Estimate:** Small
**Dependencies:** T4, T5, T6 (deletions complete)
**Milestone:** M4

**Description:**
Update `docs/codebase-map/modules/prompts.md` to reference single `system_prompt.md` file instead of 11 section files. Note that tool prompts remain unchanged.

**Acceptance Test:**
- References `prompts/system_prompt.md` as single agent prompt file
- No references to prompts/sections/ directory
- Clarifies tool prompts are separate and unchanged

**Files touched:**
- `docs/codebase-map/modules/prompts.md` (MODIFY)

---

### Task T9: Delete templates.md Documentation
**Owner:** agent
**Estimate:** Trivial
**Dependencies:** T4 (templates.py deleted)
**Milestone:** M4

**Description:**
Delete `docs/codebase-map/modules/templates.md` as template system is being removed.

**Acceptance Test:**
- `docs/codebase-map/modules/templates.md` deleted

**Files touched:**
- `docs/codebase-map/modules/templates.md` (DELETE)

---

## Risks & Mitigations

### Risk 1: Section Content Loss During Migration
**Mitigation:** Task T1 explicitly requires reading from existing section files. Do not delete `prompts/sections/` until after `system_prompt.md` is created and verified.

### Risk 2: Dynamic Placeholders Stop Working
**Mitigation:** Task T2 keeps `resolve_prompt()` call. Existing test `test_prompting_engine.py` validates placeholder resolution.

### Risk 3: Hidden Dependencies on Section Names
**Mitigation:** Search codebase for references to `SystemPromptSection`, `SectionLoader`, `compose_prompt` before starting T4.

### Risk 4: Template Overrides in Active Use
**Mitigation:** Current `TEMPLATE_OVERRIDES` dict is empty. If non-empty at execution time, document which models need custom prompts before deletion.

---

## Test Strategy

### Validation Test (for Task T2)
After completing T2, run:
```bash
# Verify load_system_prompt works
uv run python -c "
from pathlib import Path
from tunacode.core.agents.agent_components.agent_config import load_system_prompt
prompt = load_system_prompt(Path('src/tunacode'))
assert '{{USER_INSTRUCTIONS}}' in prompt
assert '{{CWD}}' not in prompt  # Should be resolved
print('OK')
"
```

### Existing Tests
- `tests/unit/core/test_prompting_engine.py` - KEEP, tests `resolve_prompt()` which is unchanged
- No new tests needed for simple file read operation

---

## References

### Research Document
- `memory-bank/research/2026-01-26_17-08-03_single-md-prompt-migration.md`

### Key Code Files
- `src/tunacode/core/agents/agent_components/agent_config.py:182-217` - Current `load_system_prompt()`
- `src/tunacode/core/prompting/templates.py:1-86` - MAIN_TEMPLATE definition
- `src/tunacode/core/prompting/sections.py:1-51` - SystemPromptSection enum
- `src/tunacode/core/prompting/prompting_engine.py:77-79` - `resolve_prompt()` (KEEP)

### Section Files (11 files to migrate)
- `src/tunacode/prompts/sections/agent_role.xml`
- `src/tunacode/prompts/sections/search_pattern.xml`
- `src/tunacode/prompts/sections/critical_rules.xml`
- `src/tunacode/prompts/sections/tool_use.xml`
- `src/tunacode/prompts/sections/completion.xml`
- `src/tunacode/prompts/sections/parallel_exec.xml`
- `src/tunacode/prompts/sections/output_style.xml`
- `src/tunacode/prompts/sections/examples.xml`
- `src/tunacode/prompts/sections/advanced_patterns.xml`
- `src/tunacode/prompts/sections/system_info.xml`
- `src/tunacode/prompts/sections/user_instructions.xml`

---

## Final Gate

**Output Summary:**
- Plan file: `memory-bank/plan/2026-01-26_17-11-28_single-md-prompt-migration.md`
- Milestones: 4 (M1: File Creation, M2: Loading Logic, M3: Cleanup, M4: Docs)
- Tasks: 9 (T1-T9)
- Files to delete: 16 (3 modules + 11 sections + 1 local prompt + 1 doc)
- Files to create: 1 (system_prompt.md)
- Files to modify: 4 (agent_config.py, __init__.py, 2 docs)

**Next Command:**
```bash
/context-engineer:execute "memory-bank/plan/2026-01-26_17-11-28_single-md-prompt-migration.md"
```

**Ready for Coding:** Yes - all tasks have explicit acceptance tests and file targets.
