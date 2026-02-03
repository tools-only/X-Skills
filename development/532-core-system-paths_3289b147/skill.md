---
title: Core System Paths
path: src/tunacode/core/system_paths.py
type: file
depth: 1
description: Core facade for system path helpers used by the UI
exports: [get_project_id, delete_session_file, check_for_updates]
seams: [M]
---

# Core System Paths

## Purpose
Provide core-layer access to session and update helpers without UI importing utils directly.

## Key Functions

### get_project_id
Returns a stable project identifier based on the current repository or working directory.

### delete_session_file
Removes a persisted session file for a project/session pair.

### check_for_updates
Checks whether a newer TunaCode CLI version is available.

## Integration Points

- **ui/app.py** - session metadata initialization
- **ui/main.py** - background update checks
- **ui/commands/__init__.py** - update and resume commands
