# One-Task-Per-File System Implementation

**Date**: 2026-02-05
**Status**: Complete

## Overview

Implemented directory-based task organization for SAM (Stateless Agent Methodology) task management, allowing one task per file instead of multi-task files.

## Implementation Summary

### 1. Parser Updates

**File**: `plugins/python3-development/skills/implementation-manager/scripts/implementation_manager.py`

**Changes**:
- Added `parse_task_directory()` function to discover and parse all `.md` files in a directory
- Updated `parse_task_file()` to detect directories and delegate to `parse_task_directory()`
- Updated `find_task_file_by_slug()` to support directory paths in addition to file paths
- Added natural numeric sorting for task IDs (T1, T2, T10 instead of T1, T10, T2)

**Key Features**:
- Automatic format detection (single file vs directory)
- Backward compatibility with existing single-file format
- Each file must have YAML frontmatter
- Files without frontmatter are skipped with warning
- Parse errors are logged but don't stop processing of other files

### 2. Split Script

**File**: `plugins/python3-development/scripts/split_task_file.py`

**Purpose**: Convert multi-task files into directory structure

**Features**:
- Reads existing task files (supports both YAML frontmatter and legacy markdown formats)
- Creates one file per task with naming convention `{task-id}-{slug}.md`
- Preserves all frontmatter metadata
- Generates basic body structure (Context, Objective, Requirements sections)
- Prompts before overwriting existing files (unless `--force` used)
- Provides clear success/failure feedback with Rich output

**Usage**:
```bash
# Split into tasks/ subdirectory
split_task_file.py plugin-validator-tasks.md

# Split into custom directory
split_task_file.py tasks.md ./my-tasks/

# Force overwrite
split_task_file.py --force tasks.md
```

### 3. Documentation Updates

**File**: `.claude/docs/TASK_FILE_FORMAT.md`

**Added Sections**:
- File Organization Options (single file vs directory)
- Directory-Based Task Organization
- File Naming Convention
- Parser Behavior
- Splitting Multi-Task Files
- Updated migration strategy to reflect completed phases

## File Naming Convention

Format: `{task-id}-{slug}.md`

Components:
- `{task-id}`: Exact task ID from frontmatter (T1, 1.1, T15, etc.)
- `{slug}`: URL-friendly task title (lowercase, hyphens, max 50 chars)

Examples:
- `T1-data-models.md` (task T1: "Data Models and Error Codes")
- `1.1-prepare-host.md` (task 1.1: "Prepare Host Environment")
- `T15-cli-tests.md` (task T15: "CLI Integration Tests")

## Testing

### Test Case 1: Split Existing File

```bash
split_task_file.py \
  /home/ubuntulinuxqa2/repos/claude_skills/plugins/plugin-creator/planning/plugin-validator-tasks.md \
  /home/ubuntulinuxqa2/repos/claude_skills/plugins/plugin-creator/planning/tasks \
  --force
```

**Result**: ✅ Successfully split 23 tasks into individual files

### Test Case 2: Directory Discovery

```bash
implementation_manager.py status \
  /home/ubuntulinuxqa2/repos/claude_skills/plugins/plugin-creator/planning \
  tasks
```

**Result**: ✅ Discovered and parsed all 23 tasks correctly

### Test Case 3: Ready Tasks Query

```bash
implementation_manager.py ready-tasks \
  /home/ubuntulinuxqa2/repos/claude_skills/plugins/plugin-creator/planning \
  tasks
```

**Result**: ✅ Correctly identified T23 as the only ready task

### Test Case 4: Backward Compatibility

```bash
implementation_manager.py status \
  /home/ubuntulinuxqa2/repos/claude_skills/plugins/plugin-creator/planning \
  plugin-validator-tasks.md
```

**Result**: ✅ Original multi-task file still works correctly

## Parser Behavior Details

### Automatic Format Detection

The parser automatically detects the format:

1. **Path is a directory**: Use `parse_task_directory()`
   - Discovers all `*.md` files
   - Parses each as a single task
   - Requires YAML frontmatter
   - Sorts tasks numerically by ID

2. **Path is a file**: Use `parse_task_file()` → `parse_task_content()`
   - Checks for YAML frontmatter format
   - Falls back to legacy markdown format if no frontmatter
   - Can contain multiple tasks

### Task Sorting

Tasks are sorted numerically by extracting integer components from task IDs:

- T1 → (1,) → sorts before T10 → (10,)
- 1.1 → (1, 1) → sorts before 1.2 → (1, 2)
- Natural ordering for human readability

### Error Handling

- **Directory not found**: Raises `FileNotFoundError`
- **Not a directory**: Raises `NotADirectoryError`
- **No .md files**: Returns empty list with warning
- **Invalid frontmatter**: Logs warning, skips file, continues parsing
- **Parse error**: Logs warning, skips file, continues parsing

## Benefits

### For Multi-Agent Workflows

1. **Parallel Execution**: Each agent can work on their own task file without git conflicts
2. **Granular History**: Git blame/log shows per-task changes instead of bulk file updates
3. **Clear Ownership**: File path directly indicates task being worked on
4. **Easy Reorganization**: Move/rename individual tasks without affecting others

### For Human Developers

1. **Better Navigation**: Jump directly to specific task file by ID
2. **Clearer Structure**: Each file contains only relevant task information
3. **Easier Review**: Review one task at a time in PRs
4. **Flexible Organization**: Group tasks into subdirectories if needed

### For Task Management

1. **Scalability**: Handles 100+ tasks without performance issues
2. **Concurrent Updates**: Multiple agents can update different tasks simultaneously
3. **Atomic Operations**: Task status updates don't require parsing entire file
4. **Simple Discovery**: List directory to see all tasks

## Migration Path

### For Existing Projects

1. **Run split script** on existing multi-task file:
   ```bash
   split_task_file.py tasks-feature-name.md
   ```

2. **Review generated files** in `tasks/` directory

3. **Fill in task details** (Context, Objective, Requirements) for each file

4. **Update references** to point to directory instead of file:
   ```bash
   # Old
   implementation_manager.py status /project tasks-feature-name.md

   # New
   implementation_manager.py status /project tasks
   ```

5. **Optional**: Remove original multi-task file after verifying migration

### Backward Compatibility

The parser maintains full backward compatibility:
- Existing single-file task files continue to work
- No changes required to existing workflows
- Mixed usage is supported (some features use files, others use directories)

## Files Modified

1. `plugins/python3-development/skills/implementation-manager/scripts/implementation_manager.py`
   - Added 95 lines (parse_task_directory function)
   - Modified 2 functions (parse_task_file, find_task_file_by_slug)
   - All tests pass
   - All linting passes

2. `plugins/python3-development/scripts/split_task_file.py`
   - New file, 238 lines
   - PEP 723 compliant
   - Executable with proper shebang
   - All linting passes

3. `.claude/docs/TASK_FILE_FORMAT.md`
   - Added directory organization section
   - Updated migration guide
   - Added examples and usage patterns

## Verification

- ✅ All linting passes (ruff, mypy, basedpyright)
- ✅ Parser discovers 23 tasks from directory
- ✅ Parser correctly identifies ready tasks
- ✅ Backward compatibility with single-file format maintained
- ✅ Natural numeric sorting works correctly (T1, T2, T10)
- ✅ Documentation updated with examples and usage

## Next Steps

1. Migrate existing multi-task files to directory structure as needed
2. Update task creation workflows to use one-task-per-file format
3. Consider adding directory support to `find_task_files()` for auto-discovery
4. Add validation for file naming convention compliance
