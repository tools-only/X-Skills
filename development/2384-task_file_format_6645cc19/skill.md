# Task File Format Specification v2.0

**Status**: Design Document
**Created**: 2026-02-04
**Purpose**: Define YAML frontmatter-based task file format to replace regex-based parsing

---

## Problem Statement

### Current Issues

The existing task file format uses markdown with bold field markers that require complex regex parsing:

```markdown
## Task 1.1: Create Skill Directory Structure

**Status**: NOT STARTED
**Dependencies**: None
**Priority**: 1
**Complexity**: Low
```

**Problems**:

1. **Fragile parsing**: Regex patterns break with formatting variations
2. **Multiple heading levels**: Parser fails when task sections use different heading depths
3. **No schema validation**: Fields can have any value without type checking
4. **Ambiguous structure**: Hard to distinguish task metadata from content
5. **Error-prone editing**: Easy to break format accidentally

### Requirements

1. **YAML frontmatter**: Store structured metadata in validated YAML block
2. **Standard parsing**: Use YAML parser library instead of regex
3. **Type safety**: Validate field types and values with schema
4. **Human-readable**: Maintain easy editing in text editors
5. **Backwards compatibility**: Support migration from old format
6. **SAM methodology**: Continue supporting Stateless Agent Methodology workflows

---

## Format Specification

### File Organization Options

**Option 1: Single File (Multi-Task)**

All tasks in one file, each with YAML frontmatter separated by `---` delimiters.
Use for small projects or when tasks have simple dependencies.

**Option 2: Directory (One-Task-Per-File)**

Each task is a separate `.md` file in a directory.
Recommended for larger projects with many tasks or complex parallel workflows.

File naming convention: `{task-id}-{slug}.md`
- Examples: `T1-data-models.md`, `1.1-prepare-host.md`, `T15-cli-tests.md`
- Slug is derived from task title (lowercase, hyphens, max 50 chars)

The parser automatically detects and supports both formats.

### Task File Structure

```markdown
---
task: T1
title: Data Models and Error Codes
status: complete
agent: python-cli-architect
dependencies: []
priority: 1
complexity: medium
created: 2026-02-02T15:00:00Z
started: 2026-02-02T15:15:00Z
completed: 2026-02-02T15:30:00Z
---

## Context

Implement core data structures for validation results, issues, and error codes.

## Objective

Create type-safe data models for ValidationResult, ValidationIssue, ComplexityMetrics.

## Requirements

1. Create ValidationResult dataclass with passed/errors/warnings/info
2. Create ValidationIssue dataclass with field/severity/message/code
3. All dataclasses must be frozen or use __post_init__ validation

## Constraints

- Use Python 3.11+ syntax (str | None, not Optional[str])
- Error codes must remain stable (no code reuse)

## Expected Outputs

- File created: plugins/plugin-creator/scripts/plugin-validator.py
- Models: ValidationResult, ValidationIssue, ComplexityMetrics

## Acceptance Criteria

1. All dataclasses type-check with mypy strict mode
2. ValidationIssue.format() produces expected output format

## Verification Steps

```bash
uv run mypy --strict plugins/plugin-creator/scripts/plugin-validator.py
uv run pytest tests/test_data_models.py -v
```

## Can Parallelize With

T2 (after data models complete)

## Handoff

Report:
- Data model file path
- All error codes implemented (count 23)
- mypy strict mode status
```

### Field Definitions

#### Required Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `task` | string | Unique task identifier | `"T1"`, `"1.1"`, `"T15"` |
| `title` | string | Brief task description | `"Create Data Models"` |
| `status` | enum | Task state | `"not-started"`, `"in-progress"`, `"complete"`, `"blocked"` |

#### Optional Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `agent` | string | Agent responsible for task | `"python-cli-architect"` |
| `dependencies` | array | Task IDs that must complete first | `["T1", "T2"]` |
| `priority` | integer | Priority level (1-5, 1=highest) | `1` |
| `complexity` | enum | Complexity estimate | `"low"`, `"medium"`, `"high"` |
| `created` | datetime | ISO 8601 timestamp when created | `"2026-02-02T15:00:00Z"` |
| `started` | datetime | ISO 8601 timestamp when started | `"2026-02-02T15:15:00Z"` |
| `completed` | datetime | ISO 8601 timestamp when completed | `"2026-02-02T15:30:00Z"` |
| `blocked-by` | array | External blockers (not task IDs) | `["API access", "Design approval"]` |
| `parallelize-with` | array | Tasks that can run concurrently | `["T2", "T3"]` |

#### Status Values

| Status | Description | Usage |
|--------|-------------|-------|
| `not-started` | Task not yet begun | Default for new tasks |
| `in-progress` | Task actively being worked on | Set when agent claims task |
| `complete` | Task finished successfully | Set after verification passes |
| `blocked` | Task cannot proceed | Set when waiting on external dependency |

#### Complexity Values

| Complexity | Description | Typical Duration |
|------------|-------------|------------------|
| `low` | Simple, well-defined task | <1 hour |
| `medium` | Moderate complexity | 1-4 hours |
| `high` | Complex, requires significant work | >4 hours |

---

## Markdown Body Sections

### Recommended Sections

The markdown body SHOULD include these sections for SAM methodology compliance:

1. **Context**: Why this task exists, background information
2. **Objective**: Clear statement of what this task achieves
3. **Requirements**: Numbered list of specific requirements
4. **Constraints**: Technical or policy constraints that must be respected
5. **Expected Outputs**: Files, artifacts, or deliverables produced
6. **Acceptance Criteria**: Testable conditions for task completion
7. **Verification Steps**: Commands or procedures to verify completion
8. **Can Parallelize With**: Tasks that can run concurrently (optional)
9. **Handoff**: What information to report to orchestrator (optional)

### Section Guidelines

**Context** - Provides historical and situational awareness:

```markdown
## Context

Implement core data structures for validation results, issues, and error codes. These models are used by all validators and reporters.
```

**Objective** - Single sentence describing the goal:

```markdown
## Objective

Create type-safe data models for ValidationResult, ValidationIssue, and ComplexityMetrics with complete error code catalog.
```

**Requirements** - Numbered list of specific deliverables:

```markdown
## Requirements

1. Create ValidationResult dataclass with passed/errors/warnings/info
2. Create ValidationIssue dataclass with field/severity/message/code/line
3. Create ComplexityMetrics dataclass with token counts and thresholds
```

**Constraints** - Hard limits and non-negotiable requirements:

```markdown
## Constraints

- Use Python 3.11+ syntax (str | None, not Optional[str])
- All dataclasses must be frozen or use __post_init__ validation
- Error codes must remain stable (no code reuse)
```

**Expected Outputs** - Concrete artifacts:

```markdown
## Expected Outputs

- File created: plugins/plugin-creator/scripts/plugin-validator.py
- Models: ValidationResult, ValidationIssue, ComplexityMetrics
- Constants: ERROR_CODE_BASE_URL, token thresholds
```

**Acceptance Criteria** - Testable pass/fail conditions:

```markdown
## Acceptance Criteria

1. All dataclasses type-check with mypy strict mode
2. Error code constants match architecture catalog exactly
3. ValidationIssue.format() produces expected output format
```

**Verification Steps** - Executable commands:

````markdown
## Verification Steps

```bash
# Type checking
uv run mypy --strict plugins/plugin-creator/scripts/plugin-validator.py

# Unit test data models
uv run pytest tests/test_data_models.py -v
```
````

---

## YAML Schema

### JSON Schema Definition

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "object",
  "required": ["task", "title", "status"],
  "properties": {
    "task": {
      "type": "string",
      "pattern": "^[A-Z]?\\d+(\\.\\d+)?$",
      "description": "Unique task identifier (e.g., T1, 1.1, T15)"
    },
    "title": {
      "type": "string",
      "minLength": 5,
      "maxLength": 100,
      "description": "Brief task description"
    },
    "status": {
      "type": "string",
      "enum": ["not-started", "in-progress", "complete", "blocked"],
      "description": "Current task state"
    },
    "agent": {
      "type": "string",
      "description": "Agent responsible for executing this task"
    },
    "dependencies": {
      "type": "array",
      "items": {
        "type": "string",
        "pattern": "^[A-Z]?\\d+(\\.\\d+)?$"
      },
      "default": [],
      "description": "Task IDs that must complete before this task"
    },
    "priority": {
      "type": "integer",
      "minimum": 1,
      "maximum": 5,
      "default": 3,
      "description": "Priority level (1=highest, 5=lowest)"
    },
    "complexity": {
      "type": "string",
      "enum": ["low", "medium", "high"],
      "default": "medium",
      "description": "Complexity estimate"
    },
    "created": {
      "type": "string",
      "format": "date-time",
      "description": "ISO 8601 timestamp when task created"
    },
    "started": {
      "type": "string",
      "format": "date-time",
      "description": "ISO 8601 timestamp when work began"
    },
    "completed": {
      "type": "string",
      "format": "date-time",
      "description": "ISO 8601 timestamp when task finished"
    },
    "blocked-by": {
      "type": "array",
      "items": {
        "type": "string"
      },
      "default": [],
      "description": "External blockers preventing progress"
    },
    "parallelize-with": {
      "type": "array",
      "items": {
        "type": "string",
        "pattern": "^[A-Z]?\\d+(\\.\\d+)?$"
      },
      "default": [],
      "description": "Tasks that can run concurrently with this one"
    }
  }
}
```

### Validation Rules

1. **task**: Must match pattern `^[A-Z]?\d+(\.\d+)?$` (e.g., `"T1"`, `"1.1"`, `"15"`)
2. **status**: Must be one of: `"not-started"`, `"in-progress"`, `"complete"`, `"blocked"`
3. **priority**: Integer 1-5 (1=highest, 5=lowest)
4. **complexity**: Must be one of: `"low"`, `"medium"`, `"high"`
5. **dependencies**: Array of task IDs matching task pattern
6. **timestamps**: ISO 8601 format with timezone (e.g., `"2026-02-02T15:00:00Z"`)

---

## Directory-Based Task Organization

### When to Use Directory Structure

Use one-task-per-file organization when:
- Project has >10 tasks
- Tasks can run in parallel
- Multiple agents work on different tasks concurrently
- Need granular git history per task
- Tasks may be reorganized or reprioritized frequently

### File Naming Convention

Format: `{task-id}-{slug}.md`

Components:
- `{task-id}`: Matches the `task:` field in frontmatter (e.g., T1, 1.1, T15)
- `{slug}`: URL-friendly version of task title (lowercase, hyphens, max 50 chars)

Examples:
- `T1-data-models.md` (task T1: "Data Models and Error Codes")
- `1.1-prepare-host.md` (task 1.1: "Prepare Host Environment")
- `T15-cli-tests.md` (task T15: "CLI Integration Tests")

### Parser Behavior

The `implementation_manager.py` parser automatically detects format:

**Single File**: Parses all tasks from one file
```bash
implementation_manager.py status /path/to/project plugin-validator-tasks.md
```

**Directory**: Discovers all .md files, parses each as one task
```bash
implementation_manager.py status /path/to/project tasks/
```

**Task Ordering**: Tasks are sorted numerically by ID (T1, T2, T10, not T1, T10, T2)

### Splitting Multi-Task Files

Convert existing multi-task files to directory structure:

```bash
# Split into tasks/ subdirectory
split-task-file.py plugin-validator-tasks.md

# Split into custom directory
split-task-file.py tasks.md ./my-tasks/

# Force overwrite existing files
split-task-file.py --force tasks.md
```

The script:
1. Parses tasks from input file (supports both YAML and markdown formats)
2. Creates output directory if needed
3. Generates one file per task with naming convention
4. Preserves all frontmatter metadata
5. Creates basic body structure (Context, Objective, Requirements sections)

After splitting, update task content in individual files.

## Migration Guide

### Migration Strategy

**Phase 1: Parser Update** (✅ COMPLETE)

1. ✅ Updated `implementation_manager.py` to parse YAML frontmatter
2. ✅ Maintained backwards compatibility with markdown format
3. ✅ Added directory discovery for one-task-per-file organization
4. ✅ Return same Task dataclass regardless of format

**Phase 2: Template Creation** (✅ COMPLETE)

1. ✅ Created task template at `plugins/python3-development/templates/sam-task-template.md`
2. ✅ Documented frontmatter fields and body sections
3. ✅ Provided examples for common task types

**Phase 3: Migration Script** (✅ COMPLETE)

1. ✅ Created `plugins/python3-development/scripts/split-task-file.py`
2. ✅ Converts existing markdown tasks to YAML frontmatter
3. ✅ Preserves all metadata and content
4. ✅ Supports directory-based organization

**Phase 4: Adoption** (in progress)

1. ✅ Update task creation workflows to use new format
2. ⏳ Migrate existing task files to YAML format or directories
3. ⏳ Mark old format as deprecated in documentation
4. ⏳ Remove markdown parsing support after migration complete

### Conversion Example

**Old Format** (markdown with bold fields):

```markdown
## Task T1: Data Models and Error Codes

**Status**: COMPLETE
**Started**: 2026-02-02T15:15:00Z
**Completed**: 2026-02-02T15:30:00Z
**Agent**: python-cli-architect
**Dependencies**: None
**Priority**: 1 (Foundational)
**Complexity**: Medium

Implement core data structures...
```

**New Format** (YAML frontmatter):

```markdown
---
task: T1
title: Data Models and Error Codes
status: complete
started: 2026-02-02T15:15:00Z
completed: 2026-02-02T15:30:00Z
agent: python-cli-architect
dependencies: []
priority: 1
complexity: medium
---

## Context

Implement core data structures...
```

### Field Mapping

| Old Format | New Format | Transformation |
|------------|------------|----------------|
| `**Status**: COMPLETE` | `status: complete` | Lowercase, hyphenated |
| `**Dependencies**: None` | `dependencies: []` | Empty array instead of "None" |
| `**Priority**: 1 (Foundational)` | `priority: 1` | Integer only, drop description |
| `**Complexity**: Medium` | `complexity: medium` | Lowercase |
| `**Agent**: python-cli-architect` | `agent: python-cli-architect` | No change |
| `**Started**: <timestamp>` | `started: <timestamp>` | No change |
| `**Completed**: <timestamp>` | `completed: <timestamp>` | No change |

### Status Value Mapping

| Old Status | New Status |
|------------|------------|
| `NOT STARTED` | `not-started` |
| `IN PROGRESS` | `in-progress` |
| `COMPLETE` | `complete` |
| `:x:` (emoji) | `not-started` |
| `:white_check_mark:` (emoji) | `complete` |
| `:arrows_counterclockwise:` (emoji) | `in-progress` |

---

## Implementation Changes

### Parser Updates Required

**File**: `plugins/python3-development/skills/implementation-manager/scripts/implementation_manager.py`

**Changes**:

1. Add YAML parser import: `import yaml`
2. Add frontmatter detection function
3. Create YAML parsing path
4. Maintain backwards compatibility with markdown parsing
5. Update validation to use schema

**Functions to Add**:

```python
def has_yaml_frontmatter(content: str) -> bool:
    """Detect if file uses YAML frontmatter format."""
    return content.strip().startswith('---\n')

def parse_yaml_frontmatter(content: str) -> tuple[dict[str, object], str]:
    """Extract YAML frontmatter and markdown body.

    Returns:
        Tuple of (frontmatter_dict, body_content)
    """
    # Split on frontmatter delimiters
    parts = content.split('---\n', 2)
    if len(parts) < 3:
        raise ValueError("Invalid frontmatter format")

    frontmatter = yaml.safe_load(parts[1])
    body = parts[2].strip()

    return frontmatter, body

def parse_task_with_frontmatter(task_section: str) -> Task:
    """Parse task from YAML frontmatter format."""
    frontmatter, body = parse_yaml_frontmatter(task_section)

    # Validate against schema
    validate_task_frontmatter(frontmatter)

    # Convert to Task dataclass
    return Task(
        id=frontmatter['task'],
        name=frontmatter['title'],
        status=TaskStatus(frontmatter['status'].upper().replace('-', '_')),
        dependencies=frontmatter.get('dependencies', []),
        agent=frontmatter.get('agent'),
        priority=TaskPriority(frontmatter.get('priority', 3)),
        complexity=frontmatter.get('complexity', 'medium').capitalize(),
        started=frontmatter.get('started'),
        completed=frontmatter.get('completed'),
    )
```

**Schema Validation**:

```python
import jsonschema

TASK_SCHEMA = {
    "type": "object",
    "required": ["task", "title", "status"],
    "properties": {
        "task": {"type": "string", "pattern": "^[A-Z]?\\d+(\\.\\d+)?$"},
        "title": {"type": "string", "minLength": 5, "maxLength": 100},
        "status": {
            "type": "string",
            "enum": ["not-started", "in-progress", "complete", "blocked"]
        },
        "priority": {"type": "integer", "minimum": 1, "maximum": 5},
        "complexity": {"type": "string", "enum": ["low", "medium", "high"]},
        "dependencies": {"type": "array", "items": {"type": "string"}},
        # ... rest of schema
    }
}

def validate_task_frontmatter(frontmatter: dict[str, object]) -> None:
    """Validate frontmatter against JSON schema."""
    try:
        jsonschema.validate(instance=frontmatter, schema=TASK_SCHEMA)
    except jsonschema.ValidationError as e:
        raise ValueError(f"Invalid task frontmatter: {e.message}") from e
```

---

## Template File

**Location**: `plugins/python3-development/templates/sam-task-template.md`

**Purpose**: Provide reusable template for creating new SAM-compliant task files

**Template Contents**:

```markdown
---
task: ""  # REQUIRED: Unique identifier (e.g., T1, 1.1, T15)
title: ""  # REQUIRED: Brief task description (5-100 chars)
status: not-started  # REQUIRED: not-started | in-progress | complete | blocked
agent: ""  # OPTIONAL: Agent responsible for execution
dependencies: []  # OPTIONAL: Task IDs that must complete first
priority: 3  # OPTIONAL: 1-5 (1=highest, 5=lowest)
complexity: medium  # OPTIONAL: low | medium | high
created: ""  # OPTIONAL: ISO 8601 timestamp
started: ""  # OPTIONAL: ISO 8601 timestamp
completed: ""  # OPTIONAL: ISO 8601 timestamp
blocked-by: []  # OPTIONAL: External blockers
parallelize-with: []  # OPTIONAL: Tasks that can run concurrently
---

## Context

Provide background and rationale for this task. Explain why it exists and what problem it solves.

## Objective

Single clear sentence describing what this task achieves.

## Requirements

1. First specific requirement
2. Second specific requirement
3. Third specific requirement

## Constraints

- Technical constraint or limitation
- Policy or architectural constraint
- Performance or quality constraint

## Expected Outputs

- File created: path/to/file.py
- Artifact generated: description
- Configuration updated: description

## Acceptance Criteria

1. First testable pass/fail condition
2. Second testable pass/fail condition
3. Third testable pass/fail condition

## Verification Steps

```bash
# Command to verify first criterion
command1 --verify

# Command to verify second criterion
command2 --test
```

## Can Parallelize With

List task IDs that can run concurrently: T2, T3, T4

**Reason**: Explanation of why parallelization is safe

## Handoff

Report to orchestrator:
- Key information needed for next phase
- Artifacts produced and their locations
- Any blockers or issues encountered
```

---

## Migration Script Specification

**Script**: `plugins/python3-development/scripts/migrate-task-format.py`

**Purpose**: Convert existing markdown-format task files to YAML frontmatter format

**Requirements**:

1. Read task file
2. Parse using existing regex-based parser
3. Extract all Task objects
4. Convert to YAML frontmatter format
5. Preserve markdown body content
6. Validate output against schema
7. Write updated file
8. Report conversion statistics

**CLI Interface**:

```bash
# Migrate single file
uv run plugins/python3-development/scripts/migrate-task-format.py tasks-refactor-plugin.md

# Migrate all task files in directory
uv run plugins/python3-development/scripts/migrate-task-format.py .claude/plan/

# Dry run (report changes without writing)
uv run plugins/python3-development/scripts/migrate-task-format.py --dry-run tasks.md

# Validate migrated files
uv run plugins/python3-development/scripts/migrate-task-format.py --validate tasks.md
```

**Output**:

```text
Migrating: tasks-refactor-plugin.md
  ✅ Task T1: Data Models and Error Codes
  ✅ Task T2: Validator Protocol Definition
  ✅ Task T3: Port FrontmatterValidator

✅ Migrated 3 tasks successfully
✅ All tasks validated against schema
✅ File written: tasks-refactor-plugin.md
```

---

## Benefits Summary

### Developer Benefits

1. **Easier Parsing**: Use `yaml.safe_load()` instead of complex regex
2. **Type Safety**: Validate fields with JSON schema
3. **Better Errors**: Schema validation provides clear error messages
4. **Extensibility**: Add new fields without breaking parser
5. **Standard Tools**: YAML editors provide syntax highlighting and validation

### Agent Benefits

1. **Reliable Parsing**: No regex edge cases or formatting variations
2. **Schema Validation**: Catch invalid data before processing
3. **Clear Structure**: Frontmatter vs body separation is explicit
4. **Maintainability**: Easier to extend with new fields
5. **Debugging**: YAML syntax errors are clear and specific

### Orchestrator Benefits

1. **Consistency**: All tasks use identical metadata format
2. **Validation**: Schema ensures all required fields present
3. **Automation**: Easier to generate tasks programmatically
4. **Querying**: Parse once, query many times efficiently
5. **Migration**: Clear migration path from old format

---

## Validation Examples

### Valid Task File

```markdown
---
task: T1
title: Create Data Models
status: in-progress
agent: python-cli-architect
dependencies: []
priority: 1
complexity: medium
created: 2026-02-02T15:00:00Z
started: 2026-02-02T15:15:00Z
---

## Context

This task creates foundational data models.

## Requirements

1. Create ValidationResult dataclass
2. Create ValidationIssue dataclass
```

✅ Valid - all required fields present, valid values

### Invalid Task Files

**Missing required field**:

```markdown
---
task: T1
status: in-progress
---
```

❌ Invalid - missing required `title` field

**Invalid status value**:

```markdown
---
task: T1
title: Create Models
status: STARTED
---
```

❌ Invalid - status must be one of: not-started, in-progress, complete, blocked

**Invalid priority**:

```markdown
---
task: T1
title: Create Models
status: not-started
priority: 10
---
```

❌ Invalid - priority must be 1-5

**Invalid task ID pattern**:

```markdown
---
task: Task-1
title: Create Models
status: not-started
---
```

❌ Invalid - task ID must match pattern `^[A-Z]?\d+(\.\d+)?$`

**Invalid timestamp format**:

```markdown
---
task: T1
title: Create Models
status: complete
completed: 2026-02-02
---
```

❌ Invalid - timestamp must be ISO 8601 with timezone (e.g., 2026-02-02T15:30:00Z)

---

## Next Steps

1. **Parser Implementation**: Update `implementation_manager.py` to parse YAML frontmatter
2. **Template Creation**: Create `sam-task-template.md` with documentation
3. **Migration Script**: Implement `migrate-task-format.py` converter
4. **Testing**: Validate migration on existing task files
5. **Documentation**: Update SAM methodology docs to reference new format

---

## References

- **SAM Methodology**: `plugins/python3-development/skills/implementation-manager/SKILL.md`
- **Current Parser**: `plugins/python3-development/skills/implementation-manager/scripts/implementation_manager.py`
- **Existing Task Files**: `.claude/plan/tasks-*.md`, `plugins/*/planning/*-tasks.md`
- **YAML Specification**: <https://yaml.org/spec/1.2.2/>
- **JSON Schema**: <https://json-schema.org/draft/2020-12/schema>

---

## Appendix: Field Evolution

### Possible Future Fields

These fields may be added in future versions without breaking existing parsers:

- `estimate-hours`: Estimated effort in hours
- `actual-hours`: Actual time spent
- `assignee`: Person or team assigned (vs agent)
- `labels`: Tags for categorization
- `epic`: Parent epic or feature ID
- `sprint`: Sprint or iteration ID
- `verification-agent`: Agent responsible for validation
- `retry-count`: Number of execution attempts
- `last-error`: Most recent error message

### Version Compatibility

All parsers MUST ignore unknown fields to maintain forward compatibility. New fields added in future versions will be optional and will not break existing task files.
