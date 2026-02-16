# Requirements: Performance Analysis

**Feature**: performance-analysis
**Status**: APPROVED
**Created**: 2026-01-30

## Context

Backlog item #12: Add `--performance` option to `zerg analyze` that runs a comprehensive performance audit covering 140 factors across 16 categories. Uses optional external static analysis tools with graceful degradation. Dynamically selects semgrep registry rules based on detected project stack.

## Functional Requirements

### FR-1: Performance audit via `zerg analyze --performance`
- Add `--performance` boolean flag to CLI
- Add `"performance"` to `--check` choices
- EXCLUDE performance from `--check all` (heavyweight)
- Run comprehensive audit using available static analysis tools

### FR-2: Factor catalog
- Load 140 factors from JSON catalog (package data)
- Filter to static-analysis-only factors (skip runtime-only)
- Map each tool to the factors it covers

### FR-3: Stack detection
- Detect project languages from file extensions
- Detect frameworks from package files (package.json, requirements.txt, go.mod, etc.)
- Detect infrastructure (Dockerfile, docker-compose, k8s manifests)

### FR-4: Tool registry with graceful degradation
- All 11 tools are optional externals (no pip deps added)
- Check availability via `shutil.which()` + version command
- Skip unavailable tools gracefully, report what was skipped
- Parallel availability checking via ThreadPoolExecutor

### FR-5: Tool adapters
- 11 adapters: semgrep, radon, lizard, vulture, jscpd, deptry, pipdeptree, dive, hadolint, trivy, cloc
- Semgrep: dynamically pull rules from registry based on detected stack
- Each adapter maps tool output to PerformanceFinding objects
- Adapters run in parallel via ThreadPoolExecutor

### FR-6: Scoring and aggregation
- Per-category score: 100 - weighted_penalties (CRITICAL=25, HIGH=10, MEDIUM=5, LOW=2)
- Overall score: weighted average of category scores by factor count
- Categories with no tools available get score=None (not 100)

### FR-7: Output formats (4 formats, parity with existing analyze)
- Rich CLI: category-grouped tables, color-coded scores, tool availability panel
- JSON: full PerformanceReport serialization
- SARIF 2.1.0: standard tool/rule/result mapping
- Markdown: category headers, finding tables, summary section

## Non-Functional Requirements

### NFR-1: No new pip dependencies
### NFR-2: Testable without external tools (mock all tool outputs)
### NFR-3: Follow existing BaseChecker pattern for integration
### NFR-4: Follow diagnostics/ submodule pattern for code organization
