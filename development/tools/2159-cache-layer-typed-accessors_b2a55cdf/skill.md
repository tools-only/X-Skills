---
title: Cache layer: CacheManager + typed accessors (mtime/version)
link: cache-layer-typed-accessors
type: doc
path: src/tunacode/infrastructure/cache/
depth: 0
seams: [A, D, M]
ontological_relations:
  - relates_to: [[cache]]
  - affects: [[tunacode.infrastructure.cache]]
  - affects: [[tunacode.infrastructure.cache.caches]]
tags:
  - cache
  - invalidation
  - testing
created_at: 2026-02-06T19:17:21.647084+00:00
updated_at: 2026-02-06T19:17:21.647084+00:00
uuid: 3dbefb03-11aa-47bb-a20a-7c546dfde6cb
---

## Summary

Tunacode uses a strict cache registry (`CacheManager`) plus *typed cache accessors* as the only supported public boundary for cache reads/writes. This removes ad-hoc module-level dict caches, enforces fail-fast behavior for unregistered caches, and makes invalidation strategies (manual vs mtime-driven) testable and deterministic.

## Context

Implemented as part of epic `tun-gqgx` (cache unification). The initial migration targets:
- agent reuse + versioning (`agent_config.py`)
- `AGENTS.md` project context loading
- `.gitignore`-driven ignore manager caching

## Changes

- Added `tunacode.infrastructure.cache`:
  - `CacheManager` singleton with explicit `register_cache()` and strict `get_cache()` (raises if missing)
  - per-cache key metadata (`set_metadata`/`get_metadata`)
  - strategies:
    - `ManualStrategy`: never auto-invalidates
    - `MtimeStrategy`: invalidates when `os.stat(metadata.path).st_mtime_ns != metadata.mtime_ns` (missing file treated as `mtime_ns=0`)
  - deterministic `clear_all()`/`clear_cache()` helpers for tests

- Added typed accessor modules under `tunacode.infrastructure.cache.caches/`:
  - `agents.py`: version-aware get/set/invalidate/clear
  - `tunacode_context.py`: mtime-aware `get_context(path: Path) -> str`

- Added typed accessor modules under `tunacode.tools.cache_accessors/`:
  - `ignore_manager_cache.py`: mtime-aware `get_ignore_manager(root: Path) -> IgnoreManager`
  - `ripgrep_cache.py`: manual cache for platform identifier + binary path
  - `xml_prompts_cache.py`: mtime-aware XML prompt caching (including cached None)

## Behavioral Impact

- No module should expose cache dicts as globals or import them in tests.
- Call sites should import accessors (e.g. `tunacode.infrastructure.cache.caches.agents`) and avoid touching `CacheManager` directly.
- Tests should clear caches via `tunacode.infrastructure.cache.clear_all()` or per-accessor `clear_*()` helpers (no singleton reset).

## Related Cards

- [[cache]]
