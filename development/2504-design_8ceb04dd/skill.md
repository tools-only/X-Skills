# Technical Design: troubleshoot-enhancement

## Metadata
- **Feature**: troubleshoot-enhancement
- **Status**: DRAFT
- **Created**: 2026-01-30

---

## 1. Overview

### 1.1 Summary
Enhance the ZERG troubleshoot system from basic regex error parsing to a world-class diagnostic engine with multi-language error intelligence, cross-worker log correlation, Bayesian hypothesis ranking, code-aware recovery planning, deep environment diagnostics, and interactive resolution workflows.

### 1.2 Goals
- Make `/zerg:troubleshoot` the most comprehensive code debugging tool available
- Provide actionable, evidence-backed root cause analysis
- Automate the entire diagnostic workflow from symptom to resolution

### 1.3 Non-Goals
- AI/LLM-powered analysis (this is deterministic code analysis)
- Real-time monitoring (this is post-hoc diagnosis)
- Support for non-text-based errors (binary/image)

---

## 2. Architecture

### 2.1 High-Level Design

```
┌─────────────────────────────────────────────────────────────────┐
│                    TroubleshootCommand (enhanced)                │
│  ┌───────────┐ ┌──────────────┐ ┌────────────┐ ┌────────────┐  │
│  │ ErrorIntel │ │ LogCorrelator│ │ Hypothesis │ │ Recovery   │  │
│  │ Engine     │ │              │ │ Engine     │ │ Planner    │  │
│  └─────┬─────┘ └──────┬───────┘ └─────┬──────┘ └─────┬──────┘  │
│        │               │               │              │         │
│  ┌─────┴─────┐ ┌──────┴───────┐ ┌─────┴──────┐ ┌────┴───────┐ │
│  │ MultiLang │ │  Timeline    │ │  Bayesian  │ │ CodeAware  │ │
│  │ Parser    │ │  Builder     │ │  Scorer    │ │ Fixer      │ │
│  └─────┬─────┘ └──────┬───────┘ └─────┬──────┘ └─────┬──────┘ │
│        │               │               │              │         │
│  ┌─────┴─────┐ ┌──────┴───────┐ ┌─────┴──────┐ ┌────┴───────┐ │
│  │ ErrorChain│ │  CrossWorker │ │  Knowledge │ │ DependGraph│ │
│  │ Analyzer  │ │  Correlator  │ │  Base      │ │ Analyzer   │ │
│  └───────────┘ └──────────────┘ └────────────┘ └────────────┘  │
│                                                                  │
│  ┌─────────────┐ ┌──────────────┐ ┌──────────────────────────┐  │
│  │ EnvDiag     │ │ Interactive  │ │ ReportGenerator          │  │
│  │ (enhanced)  │ │ Wizard       │ │ (md/json/html)           │  │
│  └─────────────┘ └──────────────┘ └──────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

### 2.2 Component Breakdown

| Component | Responsibility | Files |
|-----------|---------------|-------|
| ErrorIntelEngine | Multi-language error parsing, fingerprinting, chain analysis | `zerg/diagnostics/error_intel.py` |
| LogCorrelator | Timeline reconstruction, cross-worker correlation, temporal clustering | `zerg/diagnostics/log_correlator.py` |
| HypothesisEngine | Bayesian scoring, evidence tracking, auto-testing, knowledge base | `zerg/diagnostics/hypothesis_engine.py` |
| CodeAwareFixer | Dependency graph, git-aware recovery, context-aware fix suggestions | `zerg/diagnostics/code_fixer.py` |
| EnvDiagnostics | Python env, Docker, network, resources, config validation | `zerg/diagnostics/env_diagnostics.py` |
| InteractiveWizard | Guided troubleshooting, progressive disclosure, session persistence | `zerg/diagnostics/interactive.py` |
| ReportGenerator | Markdown/JSON/HTML report generation | `zerg/diagnostics/report_generator.py` |
| DiagnosticTypes | Shared types, enums, protocols | `zerg/diagnostics/types.py` |
| TroubleshootCommand | Enhanced orchestrator integrating all engines | `zerg/commands/troubleshoot.py` (modify) |

### 2.3 Data Flow

1. User invokes `zerg troubleshoot` with error/feature/flags
2. ErrorIntelEngine parses and classifies the error across languages
3. LogCorrelator builds timeline and finds cross-worker patterns
4. HypothesisEngine generates ranked hypotheses with evidence
5. HypothesisEngine auto-tests hypotheses and updates confidence
6. CodeAwareFixer generates recovery plan with dependency awareness
7. EnvDiagnostics checks environment health
8. ReportGenerator produces structured output
9. InteractiveWizard guides user through resolution (if `--interactive`)

---

## 3. Detailed Design

### 3.1 Data Models (types.py)

```python
class ErrorSeverity(Enum):
    CRITICAL = "critical"   # System down, data loss
    ERROR = "error"         # Feature broken
    WARNING = "warning"     # Degraded functionality
    INFO = "info"           # Informational

class ErrorCategory(Enum):
    WORKER_FAILURE = "worker_failure"
    TASK_FAILURE = "task_failure"
    STATE_CORRUPTION = "state_corruption"
    INFRASTRUCTURE = "infrastructure"
    CODE_ERROR = "code_error"
    DEPENDENCY = "dependency"
    MERGE_CONFLICT = "merge_conflict"
    ENVIRONMENT = "environment"
    CONFIGURATION = "configuration"
    UNKNOWN = "unknown"

@dataclass
class ErrorFingerprint:
    hash: str               # Dedupe key
    language: str            # python, javascript, go, rust, java, cpp
    error_type: str
    message_template: str    # Parameterized message
    file: str
    line: int
    column: int
    function: str
    module: str
    chain: list[ErrorFingerprint]  # caused-by chain

@dataclass
class TimelineEvent:
    timestamp: str
    worker_id: int
    event_type: str          # error, warning, info, state_change
    message: str
    source_file: str
    line_number: int
    correlation_id: str      # For grouping related events

@dataclass
class Evidence:
    description: str
    source: str              # log, state, git, system, code
    confidence: float        # 0.0 - 1.0
    data: dict[str, Any]

@dataclass
class ScoredHypothesis:
    description: str
    category: ErrorCategory
    prior_probability: float
    evidence_for: list[Evidence]
    evidence_against: list[Evidence]
    posterior_probability: float  # After evidence
    test_command: str
    test_result: str | None      # None = untested
    suggested_fix: str
```

---

## 4. Key Decisions

### Decision: Bayesian Hypothesis Scoring

**Context**: Current system uses simple high/medium/low likelihood strings.

**Options**:
1. Simple weighted scoring
2. Bayesian probability updates
3. ML-based classification

**Decision**: Bayesian probability updates

**Rationale**: Provides principled evidence integration, interpretable confidence scores, and doesn't require training data. Simple enough to implement in pure Python.

### Decision: New modules vs. refactoring existing

**Context**: Existing `diagnostics/` has 4 modules. Enhancement is substantial.

**Options**:
1. Refactor existing modules in-place
2. Create new modules alongside existing
3. Replace existing entirely

**Decision**: Create new modules alongside existing, then update TroubleshootCommand to use them. Existing modules continue working, new modules add capabilities.

**Rationale**: Zero-risk migration path. Old tests keep passing. New functionality is additive.

---

## 5. Implementation Plan

### Phase Summary

| Phase | Level | Tasks | Parallel |
|-------|-------|-------|----------|
| Foundation | 1 | 2 | Yes |
| Core Engines | 2 | 5 | Yes |
| Integration | 3 | 2 | Yes |
| Testing | 4 | 5 | Yes |
| Quality | 5 | 1 | No |

### File Ownership

| File | Task ID | Operation |
|------|---------|-----------|
| `zerg/diagnostics/types.py` | TE-L1-001 | create |
| `zerg/diagnostics/knowledge_base.py` | TE-L1-002 | create |
| `zerg/diagnostics/error_intel.py` | TE-L2-001 | create |
| `zerg/diagnostics/log_correlator.py` | TE-L2-002 | create |
| `zerg/diagnostics/hypothesis_engine.py` | TE-L2-003 | create |
| `zerg/diagnostics/code_fixer.py` | TE-L2-004 | create |
| `zerg/diagnostics/env_diagnostics.py` | TE-L2-005 | create |
| `zerg/commands/troubleshoot.py` | TE-L3-001 | modify |
| `zerg/diagnostics/__init__.py` | TE-L3-002 | modify |
| `tests/unit/test_diagnostics/test_types.py` | TE-L4-001 | create |
| `tests/unit/test_diagnostics/test_error_intel.py` | TE-L4-002 | create |
| `tests/unit/test_diagnostics/test_log_correlator.py` | TE-L4-003 | create |
| `tests/unit/test_diagnostics/test_hypothesis_engine.py` | TE-L4-004 | create |
| `tests/unit/test_diagnostics/test_code_fixer.py` | TE-L4-005 | create |
| `zerg/data/commands/zerg:troubleshoot.md` | TE-L5-001 | modify |

### Dependency Graph

```
L1: TE-L1-001 (types) ──┐
    TE-L1-002 (kb)    ──┤
                         │
L2: TE-L2-001 (error)  ─┤─ depends on L1
    TE-L2-002 (logs)   ─┤
    TE-L2-003 (hypo)   ─┤
    TE-L2-004 (fixer)  ─┤
    TE-L2-005 (env)    ─┘
                         │
L3: TE-L3-001 (cmd)   ──┤─ depends on L2
    TE-L3-002 (init)   ─┘
                         │
L4: TE-L4-001..005     ──┤─ depends on L3
                         │
L5: TE-L5-001 (docs)  ──┘─ depends on L4
```

---

## 6. Risk Assessment

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Backward compat break | Low | High | New modules additive, old API preserved |
| Test coverage gap | Medium | Medium | Each L4 task covers one L2 module |
| Subprocess commands fail in CI | Medium | Low | Mock all subprocess calls in tests |

---

## 7. Parallel Execution Notes

### Recommended Workers
- Minimum: 2 workers
- Optimal: 5 workers (widest level is L2 with 5 tasks)
- Maximum: 5 workers

### Estimated Duration
- Single worker: ~15 tasks sequential
- With 5 workers: 5 levels of parallelized work

---

## 8. Approval

| Role | Name | Date | Signature |
|------|------|------|-----------|
| Architecture | | | PENDING |
