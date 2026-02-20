---
description: Process all research files to add/update Integration Opportunities sections. Orchestrates the research-context-agent to cross-reference each file with existing skills, agents, and plugins.
---

# Process Research Integration Opportunities

Batch process all research markdown files to discover and document integration opportunities with existing repository capabilities.

---

## What This Does

1. Finds all markdown files in `research/` (excluding README.md files)
2. For each research file, spawns the **research-context-agent**
3. Agent reads the file, searches for connections, appends Integration Opportunities section
4. Reports summary of all processed files

---

## Usage

```bash
/process-research-integration
```

**Options**:

```bash
# Process a specific category
/process-research-integration --category developer-tools

# Process a single file
/process-research-integration --file research/developer-tools/loguru.md

# Force reprocess (replace existing Integration Opportunities)
/process-research-integration --force

# Dry run (show what would be processed without making changes)
/process-research-integration --dry-run
```

---

## Implementation

<implementation>

This command uses the batch processing script at `./scripts/process-research-integration.py`.

**Script Features**:
- Cross-platform Python script using `uv` with inline script metadata
- CLI built with Typer for command-line options
- Rich console output for progress tracking and results visualization
- Supports category filtering, single file processing, force reprocessing, and dry-run mode

**For manual invocation**:
```bash
# Using uv (recommended - handles dependencies automatically)
uv run ./scripts/process-research-integration.py --help

# Or make executable and run directly
chmod +x ./scripts/process-research-integration.py
./scripts/process-research-integration.py --help
```

**Processing workflow**:
1. Find research files matching criteria (category/file/all)
2. For each file, spawn research-context-agent via Task tool
3. Agent reads, analyzes, validates with WebSearch, appends Integration Opportunities
4. Track and display results in formatted tables

**Agent orchestration** (to be implemented):
- Use Task tool to spawn research-context-agent instances
- Pass file path as parameter
- Sequential processing (one file at a time)
- Collect results and generate summary

</implementation>

---

## Expected Output

```
🔍 Processing Research Integration Opportunities

Found 58 research files across 24 categories

Processing: research/developer-tools/loguru.md
✅ Completed (2 enhancements, 0 skills, 1 cross-ref)

Processing: research/developer-tools/traycer.md
✅ Completed (1 enhancement, 0 skills, 0 cross-refs)

...

📊 Summary
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Files processed:        58
Enhancements found:     87
New skill candidates:   12
New MCP candidates:     8
Cross-references:       23
Failures:               0
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

✅ All research files processed
```

---

## Notes

- **Processing time**: With 58 files, expect ~30-60 minutes for full batch
- **Incremental updates**: Can re-run safely; existing sections will be replaced
- **Review required**: Auto-generated opportunities should be reviewed before implementation
- **Cross-reference sync**: When updating one file, related files may need updates too

---

## Quality Assurance

The research-context-agent applies these quality gates:

✅ Concrete, actionable descriptions (not vague suggestions)
✅ Specific connection between research and target
✅ Empty sections omitted (not filled with "None")
✅ Existing content preserved (append-only)
✅ Idempotent (can re-run without duplication)

---

## Related

- Agent: `.claude/agents/research-context-agent.md` — The agent that does the actual processing
- Script: `./scripts/process-research-integration.py` — Python batch processing implementation
- Skill: `.claude/skills/research-curator/` — Related research management workflow
- Directory: `research/` — All research entries
