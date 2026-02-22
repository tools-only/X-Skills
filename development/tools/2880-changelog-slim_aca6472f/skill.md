# Memory Copilot Slim Initiative Changes

## Summary

Refactored Memory Copilot to store minimal initiative data (pointers only), with detailed task data expected to live in Task Copilot. This reduces context bloat when `/continue` loads initiatives.

## Changes Made

### 1. Database Schema (v2 Migration)

**New columns added to `initiatives` and `initiatives_archive`:**
- `task_copilot_linked` (INTEGER) - Whether initiative is linked to Task Copilot
- `active_prd_ids` (TEXT) - JSON array of active PRD IDs
- `current_focus` (TEXT) - Current work focus (max 100 chars)
- `next_action` (TEXT) - Next action to take (max 100 chars)

**Existing columns retained for backward compatibility:**
- `completed`, `in_progress`, `blocked`, `resume_instructions` (marked DEPRECATED)

**Migration:** Automatic on first run. Version bumped from 1 to 2.

### 2. Type Definitions

**Updated `Initiative` interface** (`src/types.ts`):
```typescript
interface Initiative {
  // NEW: Task Copilot integration
  taskCopilotLinked?: boolean;
  activePrdIds?: string[];

  // NEW: Slim resume context
  currentFocus?: string;   // Max 100 chars
  nextAction?: string;     // Max 100 chars

  // DEPRECATED (but kept)
  completed: string[];
  inProgress: string[];
  blocked: string[];
  resumeInstructions?: string;
}
```

**Added new input/output types:**
- `InitiativeSlimInput` - Input for slim operation
- `InitiativeSlimOutput` - Output with savings metrics

### 3. New Tool: `initiative_slim`

**Purpose:** Slim down bloated initiatives by removing task lists

**Input:**
```json
{
  "archiveDetails": true  // Archive removed data to file (default: true)
}
```

**Output:**
```json
{
  "initiativeId": "abc-123",
  "archived": true,
  "archivePath": "~/.claude/memory/archives/abc-123_2025-01-15_pre_slim.json",
  "removedFields": ["completed", "inProgress", "blocked", "resumeInstructions"],
  "beforeSize": 2400,
  "afterSize": 600,
  "savings": "75% reduction"
}
```

**Behavior:**
- Archives removed data to `~/.claude/memory/archives/{id}_pre_slim.json`
- Removes: `completed`, `inProgress`, `blocked`, `resumeInstructions`
- Keeps: `decisions`, `lessons`, `keyFiles` (permanent knowledge)
- Resets Task Copilot fields to defaults

### 4. Enhanced `initiative_update`

**New fields accepted:**
- `taskCopilotLinked` (boolean)
- `activePrdIds` (string[])
- `currentFocus` (string, auto-truncated to 100 chars)
- `nextAction` (string, auto-truncated to 100 chars)

**Backward compatibility:**
- Still accepts deprecated fields (`completed`, `inProgress`, etc.)
- Warns to console if `completed` array exceeds 10 items

### 5. Enhanced `initiative_get`

**New behavior:**
- Adds hint to output if initiative is bloated (>10 tasks or >200 char resume)
- Hint text: "HINT: This initiative has bloated task lists. Consider using initiative_slim to reduce context usage."

### 6. Documentation Updates

**README.md:**
- Added "Initiative Structure" section explaining slim vs legacy modes
- Added usage examples for slim workflow
- Updated tool descriptions

## Files Modified

1. `src/types.ts` - Added new fields and types
2. `src/db/schema.ts` - Added new columns, v2 migration
3. `src/db/client.ts` - Updated upsert/archive to handle new columns, migration logic
4. `src/tools/initiative-tools.ts` - Added `initiativeSlim`, updated conversion functions
5. `src/tools/index.ts` - Exported new function
6. `src/index.ts` - Added tool definition and handler
7. `README.md` - Documentation updates

## Backward Compatibility

- Existing initiatives continue to work unchanged
- Old fields (`completed`, etc.) are retained in database
- Migration is automatic and non-destructive
- `initiative_slim` provides opt-in migration path

## Token Savings

**Before (bloated initiative):**
- 50+ completed tasks
- 5-10 in-progress tasks
- Verbose resume instructions
- Total: ~2000-3000 tokens

**After (slim initiative):**
- Task Copilot pointers only
- 100 char focus/action strings
- Decisions/lessons retained
- Total: ~500-800 tokens

**Typical savings: 60-75% reduction in context usage**

## Testing

To test the changes:

1. Start initiative: `initiative_start`
2. Add some tasks using legacy fields: `initiative_update` with `completed`, `inProgress`
3. Check for hint: `initiative_get` (should show bloat warning if >10 items)
4. Slim it down: `initiative_slim`
5. Verify reduction: Check `beforeSize` vs `afterSize` in output
6. Confirm archive: Check file at `archivePath`
7. Test new fields: `initiative_update` with `currentFocus`, `nextAction`

## Migration Guide

### For Existing Initiatives

```bash
# 1. Get current initiative to see size
claude tool initiative_get

# 2. Slim it down (archives to file first)
claude tool initiative_slim '{"archiveDetails": true}'

# 3. Update with new slim fields
claude tool initiative_update '{
  "taskCopilotLinked": true,
  "activePrdIds": ["PRD-001"],
  "currentFocus": "Implementing Phase 2",
  "nextAction": "Run task TASK-123"
}'
```

### For New Initiatives

```bash
# Use slim mode from the start
claude tool initiative_start '{"name": "New Feature", "goal": "Build it"}'
claude tool initiative_update '{
  "currentFocus": "Starting Phase 1",
  "nextAction": "Create PRD"
}'
```

## Next Steps

1. Update `/continue` command to use slim fields when resuming
2. Integrate with Task Copilot for seamless task management
3. Consider deprecation timeline for legacy fields (future v3)
