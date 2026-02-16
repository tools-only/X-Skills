# Technical Design: Performance Analysis

## Metadata
- **Feature**: performance-analysis
- **Status**: APPROVED
- **Created**: 2026-01-30

## 1. Overview

### 1.1 Summary
Add `zerg analyze --performance` as a comprehensive static performance audit. A new `zerg/performance/` submodule (following the `zerg/diagnostics/` pattern) houses 11 tool adapters, a factor catalog, stack detection, scoring aggregation, and 4 output formatters. Integration into the existing `analyze.py` is minimal: one new CheckType enum value, one bridge checker class, and CLI flag additions.

### 1.2 Goals
- Cover 140 performance factors across 16 categories via static analysis
- Graceful degradation when tools are missing
- Dynamic semgrep rule selection based on project stack
- Rich multi-format output (text, JSON, SARIF, markdown)

### 1.3 Non-Goals
- Runtime profiling or dynamic analysis
- Advisory checklists for non-static factors
- Adding new pip dependencies

## 2. Architecture

### 2.1 High-Level Design

```
CLI (analyze.py)
  └─ PerformanceChecker (bridge)
       └─ PerformanceAuditor (aggregator.py)
            ├─ FactorCatalog (catalog.py) ─── factors.json
            ├─ StackDetector (stack_detector.py)
            ├─ ToolRegistry (tool_registry.py)
            └─ Adapters[] (adapters/*.py)
                 ├─ SemgrepAdapter
                 ├─ RadonAdapter
                 ├─ LizardAdapter
                 ├─ VultureAdapter
                 ├─ JscpdAdapter
                 ├─ DeptryAdapter
                 ├─ PipdeptreeAdapter
                 ├─ DiveAdapter
                 ├─ HadolintAdapter
                 ├─ TrivyAdapter
                 └─ ClocAdapter
```

### 2.2 Component Breakdown

| Component | Responsibility | Files |
|-----------|---------------|-------|
| types | Data models for findings, reports, scores | `zerg/performance/types.py` |
| catalog | Factor loading, filtering, tool-to-factor mapping | `zerg/performance/catalog.py` |
| stack_detector | Detect languages, frameworks, infra | `zerg/performance/stack_detector.py` |
| tool_registry | Check tool availability, advisory output | `zerg/performance/tool_registry.py` |
| adapters | 11 tool adapters + base ABC | `zerg/performance/adapters/` |
| aggregator | Orchestrate adapters, score categories | `zerg/performance/aggregator.py` |
| formatters | Rich, JSON, SARIF, Markdown output | `zerg/performance/formatters.py` |

### 2.3 Data Flow

1. CLI invokes `PerformanceChecker.check(files)`
2. Bridge creates `PerformanceAuditor` and calls `run(files)`
3. Auditor loads factor catalog, detects stack, checks tool availability
4. Auditor selects applicable adapters (based on stack + tool availability)
5. Adapters run in parallel via ThreadPoolExecutor
6. Each adapter returns `list[PerformanceFinding]`
7. Auditor aggregates findings by category, computes scores
8. Returns `PerformanceReport` to bridge
9. Bridge converts to `AnalysisResult` for standard flow
10. CLI uses perf-specific formatters when `check == performance`

## 3. Detailed Design

### 3.1 Data Models (types.py)

```python
class Severity(Enum):
    CRITICAL = "critical"   # Weight: 25 points
    HIGH = "high"           # Weight: 10 points
    MEDIUM = "medium"       # Weight: 5 points
    LOW = "low"             # Weight: 2 points
    INFO = "info"           # Weight: 0 points

@dataclass
class PerformanceFactor:
    id: int
    category: str
    factor: str
    description: str
    cli_tools: list[str]
    security_note: str | None = None

@dataclass
class PerformanceFinding:
    factor_id: int
    factor_name: str
    category: str
    severity: Severity
    message: str
    file: str = ""
    line: int = 0
    tool: str = ""
    rule_id: str = ""
    suggestion: str = ""

@dataclass
class ToolStatus:
    name: str
    available: bool
    version: str = ""
    factors_covered: int = 0

@dataclass
class CategoryScore:
    category: str
    score: float | None  # None if no tools available
    findings: list[PerformanceFinding]
    factors_checked: int
    factors_total: int

@dataclass
class PerformanceReport:
    overall_score: float | None
    categories: list[CategoryScore]
    tool_statuses: list[ToolStatus]
    findings: list[PerformanceFinding]
    factors_checked: int
    factors_total: int
    detected_stack: dict[str, list[str]]

    def to_dict(self) -> dict: ...
    def top_issues(self, limit: int = 20) -> list[str]: ...
```

### 3.2 Adapter Interface (adapters/base.py)

```python
class BaseToolAdapter(ABC):
    name: str
    tool_name: str
    factors_covered: list[int]

    @abstractmethod
    def run(self, files: list[str], project_path: str, stack: DetectedStack) -> list[PerformanceFinding]: ...

    def is_applicable(self, stack: DetectedStack) -> bool:
        return True  # Override for language/infra-specific adapters
```

### 3.3 Semgrep Adapter Registry Selection

```python
# Map detected stack to semgrep config sets
SEMGREP_CONFIGS = {
    "python": ["p/python", "p/python-best-practices"],
    "javascript": ["p/javascript", "p/nodejs"],
    "typescript": ["p/typescript"],
    "go": ["p/golang"],
    "rust": ["p/rust"],
    "docker": ["p/dockerfile"],
}
# Always include: p/performance (if available)
```

### 3.4 Scoring Algorithm

Per-category:
```
penalty = sum(SEVERITY_WEIGHTS[f.severity] for f in category_findings)
max_penalty = factors_checked * 25  # Assume worst case
score = max(0, 100 - (penalty / max(1, factors_checked)) * 10)
```

Overall:
```
weighted_sum = sum(cat.score * cat.factors_checked for cat in categories if cat.score is not None)
total_weight = sum(cat.factors_checked for cat in categories if cat.score is not None)
overall = weighted_sum / max(1, total_weight)
```

## 4. Key Decisions

### 4.1 Submodule vs Inline
**Decision**: New `zerg/performance/` submodule
**Rationale**: 11 adapters + orchestrator + formatters would bloat analyze.py past 2000 lines. Follows established `zerg/diagnostics/` pattern.

### 4.2 Tool Dependencies
**Decision**: All optional externals, no pip deps
**Rationale**: Users install only what they need. Graceful degradation preserves usability.

### 4.3 Semgrep Rules
**Decision**: Dynamic registry pull based on detected stack
**Rationale**: Project-specific rules provide better coverage than shipping static rules.

### 4.4 Non-Static Factors
**Decision**: Skip entirely (no advisory checklist)
**Rationale**: User preference. Only report what tools can measure.

### 4.5 Performance Excluded from --check all
**Decision**: Require explicit --performance or --check performance
**Rationale**: Performance audit is heavyweight; don't slow down standard quick checks.

## 5. Implementation Plan

### 5.1 Phase Summary

| Phase | Level | Tasks | Parallel | Description |
|-------|-------|-------|----------|-------------|
| Foundation | 1 | 4 | Yes | Types, catalog, stack detection, tool registry |
| Adapters | 2 | 4 | Yes | 11 tool adapters grouped by domain |
| Aggregation | 3 | 2 | Yes | Orchestrator + formatters |
| Integration | 4 | 1 | No | CLI modifications to analyze.py |
| Testing | 5 | 2 | Yes | Unit tests with mocked tools |

### 5.2 Recommended Workers
- Optimal: **4 workers** (matches max parallelization at L1 and L2)

## 6. Testing Strategy

### 6.1 Unit Tests
- types.py: Serialization, to_dict(), severity ordering
- catalog.py: Load, filter static-only, tool mapping
- stack_detector.py: Detect Python, JS, Go, Rust, Docker, K8s projects
- tool_registry.py: Available/missing tool detection (mock shutil.which)
- adapters: Each adapter with mocked subprocess output → findings
- aggregator: Scoring math with synthetic findings
- formatters: JSON schema, SARIF structure, Rich output, Markdown headings

### 6.2 Integration
- CLI `--performance` flag invocation (mocked tools)
- All 4 output formats produce valid output
- Graceful degradation with zero tools available

## 7. Risk Assessment

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Semgrep registry rules change | Medium | Low | Map by category, not exact rule ID |
| Tool output format changes | Low | Medium | Version-pin expected output format per adapter |
| Too many findings overwhelm output | Medium | Low | top_issues() limit, collapsible Rich panels |
