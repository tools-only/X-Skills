# WM-004: Update state schema with metrics definitions

## Overview
Extend the state schema to include metrics fields for validation.

## Files to Modify
- `.zerg/schemas/state.schema.json`

## Requirements

### Updated Schema

Replace the entire schema file with:

```json
{
  "$schema": "http://json-schema.org/draft-07/schema#",
  "type": "object",
  "required": ["feature", "started_at", "current_level", "tasks"],
  "properties": {
    "feature": {"type": "string"},
    "started_at": {"type": "string", "format": "date-time"},
    "current_level": {"type": "integer", "minimum": 0},
    "paused": {"type": "boolean"},
    "error": {"type": ["string", "null"]},
    "tasks": {
      "type": "object",
      "additionalProperties": {"$ref": "#/definitions/task"}
    },
    "workers": {
      "type": "object",
      "additionalProperties": {"$ref": "#/definitions/worker"}
    },
    "levels": {
      "type": "object",
      "additionalProperties": {"$ref": "#/definitions/level"}
    },
    "execution_log": {
      "type": "array",
      "items": {"$ref": "#/definitions/event"}
    },
    "checkpoints": {
      "type": "array",
      "items": {"type": "string"}
    },
    "metrics": {
      "oneOf": [
        {"type": "null"},
        {"$ref": "#/definitions/metrics"}
      ]
    }
  },
  "definitions": {
    "task": {
      "type": "object",
      "properties": {
        "id": {"type": "string"},
        "status": {"enum": ["pending", "todo", "claimed", "in_progress", "complete", "failed", "blocked"]},
        "worker_id": {"type": ["integer", "null"]},
        "created_at": {"type": ["string", "null"]},
        "claimed_at": {"type": ["string", "null"]},
        "started_at": {"type": ["string", "null"]},
        "completed_at": {"type": ["string", "null"]},
        "updated_at": {"type": ["string", "null"]},
        "duration_ms": {"type": ["integer", "null"]},
        "retry_count": {"type": "integer", "minimum": 0},
        "error": {"type": ["string", "null"]}
      }
    },
    "worker": {
      "type": "object",
      "properties": {
        "worker_id": {"type": "integer"},
        "status": {"type": "string"},
        "current_task": {"type": ["string", "null"]},
        "port": {"type": ["integer", "null"]},
        "container_id": {"type": ["string", "null"]},
        "worktree_path": {"type": ["string", "null"]},
        "branch": {"type": ["string", "null"]},
        "started_at": {"type": ["string", "null"]},
        "ready_at": {"type": ["string", "null"]},
        "health_check_at": {"type": ["string", "null"]},
        "last_task_completed_at": {"type": ["string", "null"]},
        "tasks_completed": {"type": "integer", "minimum": 0},
        "context_usage": {"type": "number", "minimum": 0, "maximum": 1}
      }
    },
    "level": {
      "type": "object",
      "properties": {
        "name": {"type": "string"},
        "status": {"type": "string"},
        "started_at": {"type": ["string", "null"]},
        "completed_at": {"type": ["string", "null"]},
        "updated_at": {"type": ["string", "null"]},
        "merge_status": {"type": ["string", "null"]},
        "merge_commit": {"type": ["string", "null"]},
        "merge_completed_at": {"type": ["string", "null"]}
      }
    },
    "event": {
      "type": "object",
      "properties": {
        "timestamp": {"type": "string"},
        "event": {"type": "string"},
        "data": {"type": "object"}
      }
    },
    "metrics": {
      "type": "object",
      "properties": {
        "computed_at": {"type": "string", "format": "date-time"},
        "total_duration_ms": {"type": ["integer", "null"]},
        "workers_used": {"type": "integer", "minimum": 0},
        "tasks_total": {"type": "integer", "minimum": 0},
        "tasks_completed": {"type": "integer", "minimum": 0},
        "tasks_failed": {"type": "integer", "minimum": 0},
        "levels_completed": {"type": "integer", "minimum": 0},
        "worker_metrics": {
          "type": "array",
          "items": {"$ref": "#/definitions/worker_metrics"}
        },
        "level_metrics": {
          "type": "array",
          "items": {"$ref": "#/definitions/level_metrics"}
        }
      }
    },
    "worker_metrics": {
      "type": "object",
      "properties": {
        "worker_id": {"type": "integer"},
        "initialization_ms": {"type": ["integer", "null"]},
        "uptime_ms": {"type": "integer", "minimum": 0},
        "tasks_completed": {"type": "integer", "minimum": 0},
        "tasks_failed": {"type": "integer", "minimum": 0},
        "total_task_duration_ms": {"type": "integer", "minimum": 0},
        "avg_task_duration_ms": {"type": "number", "minimum": 0}
      }
    },
    "level_metrics": {
      "type": "object",
      "properties": {
        "level": {"type": "integer"},
        "duration_ms": {"type": ["integer", "null"]},
        "task_count": {"type": "integer", "minimum": 0},
        "completed_count": {"type": "integer", "minimum": 0},
        "failed_count": {"type": "integer", "minimum": 0},
        "avg_task_duration_ms": {"type": "number", "minimum": 0},
        "p50_duration_ms": {"type": "integer", "minimum": 0},
        "p95_duration_ms": {"type": "integer", "minimum": 0}
      }
    }
  }
}
```

## Verification

```bash
python -c "import json; d = json.load(open('.zerg/schemas/state.schema.json')); assert 'metrics' in str(d); print('Schema OK')"
```

## Acceptance Criteria
- [ ] task definition includes claimed_at property
- [ ] task definition includes duration_ms property
- [ ] worker definition includes ready_at property
- [ ] worker definition includes last_task_completed_at property
- [ ] metrics definition exists with all required fields
- [ ] worker_metrics definition exists
- [ ] level_metrics definition exists
- [ ] Schema is valid JSON
