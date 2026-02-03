# sf-ai-agentforce-observability Validation Tracking

## Overview

This document tracks validation runs for the `sf-ai-agentforce-observability` skill. Validation ensures the skill's CLI commands and API integrations work correctly against live Salesforce Data 360 environments.

**Validation Cycle:** Every 30 days or after significant changes

---

## Current Status

| Metric | Value |
|--------|-------|
| Last Validation | _Not yet validated_ |
| Overall Score | _Pending_ |
| Status | 🟡 PENDING |
| Next Validation Due | _After initial run_ |

---

## Scoring Rubric (100 Points)

| Tier | Category | Points | Description |
|------|----------|--------|-------------|
| T1 | Auth & Connectivity | 25 | JWT auth, consumer key resolution, API access |
| T2 | Extraction Commands | 30 | extract, extract-tree, extract-incremental, count |
| T3 | Analysis Commands | 20 | analyze, debug-session, topics, quality-report |
| T4 | Data Model/Schema | 15 | DMO names, field casing, Parquet structure |
| T5 | Negative/Error Cases | 10 | Auth failures, invalid arguments |

**Pass Threshold:** ≥80 points

---

## Validation History

### Run #0 - Initial (Pending)

| Tier | Score | Status | Notes |
|------|-------|--------|-------|
| T1 | -/25 | ⏳ | Pending |
| T2 | -/30 | ⏳ | Pending |
| T3 | -/20 | ⏳ | Pending |
| T4 | -/15 | ⏳ | Pending |
| T5 | -/10 | ⏳ | Pending |
| **Total** | **-/100** | ⏳ | **Pending** |

---

## Test Mapping to SKILL.md Claims

### T1: Auth & Connectivity (25 pts)

| Test | SKILL.md Section | Points |
|------|------------------|--------|
| Key path resolution order | "Key Path Resolution Order" | 5 |
| Consumer key resolution | "Auth Setup" | 5 |
| v65.0 API connectivity | "Prerequisites Checklist" | 10 |
| Error message quality | "Common Issues & Fixes" | 5 |

### T2: Extraction Commands (30 pts)

| Test | SKILL.md Section | Points |
|------|------------------|--------|
| extract creates 4 directories | "Output Directory Structure" | 10 |
| Date filtering works | "CLI Quick Reference" | 5 |
| extract-tree for session IDs | "Extract Session Tree" | 5 |
| Incremental watermark behavior | "Incremental Extraction" | 5 |
| count command | "Count Records" | 5 |

### T3: Analysis Commands (20 pts)

| Test | SKILL.md Section | Points |
|------|------------------|--------|
| analyze generates stats | "Analysis Examples" | 5 |
| debug-session timeline | "Debug Session Timeline" | 5 |
| topics routing analysis | "Topic Analysis" | 5 |
| quality-report metrics | "Quality Report" | 5 |

### T4: Data Model/Schema (15 pts)

| Test | SKILL.md Section | Points |
|------|------------------|--------|
| DMO names correct | "Session Tracing Data Model" | 5 |
| AiAgent field casing | "Key Schema Notes" | 5 |
| Parquet directory structure | "Output Directory Structure" | 5 |

### T5: Negative/Error Cases (10 pts)

| Test | SKILL.md Section | Points |
|------|------------------|--------|
| Auth failure handling | "Common Issues & Fixes" | 5 |
| Invalid args rejected | "CLI Quick Reference" | 5 |

---

## Known Issues

_None recorded yet._

---

## Remediation Log

| Date | Issue | Resolution | Commit |
|------|-------|------------|--------|
| _None yet_ | | | |

---

## Running Validation

```bash
# Quick offline test (T3, T4, T5 only - uses fixtures)
cd sf-ai-agentforce-observability
python3 validation/scripts/run_validation.py --offline

# Full validation against live org
python3 validation/scripts/run_validation.py --org Vivint-DevInt

# Generate report and update this file
python3 validation/scripts/run_validation.py --org Vivint-DevInt --report

# Run specific tier
python3 validation/scripts/run_validation.py --org Vivint-DevInt --tier T1
```

---

## Integration with test_v65.py

The existing `scripts/test_v65.py` provides a quick smoke test with 7 progressive levels. The validation framework wraps these and adds:

- **Scoring** - Quantified pass/fail with points
- **Isolation** - Offline tests using fixtures
- **Coverage** - Tests CLI commands, not just API
- **Tracking** - Historical validation runs

| test_v65.py Level | Validation Tier |
|-------------------|-----------------|
| Level 1-2 (Connectivity) | T1 |
| Level 3-5 (Query API) | T1 |
| Level 6-7 (Extraction) | T2 |

---

## Next Validation Checklist

- [ ] Ensure test org has STDM data (Agentforce sessions)
- [ ] Verify JWT key exists at `~/.sf/jwt/{org}.key`
- [ ] Verify consumer key exists at `~/.sf/jwt/{org}.consumer-key`
- [ ] Run: `python3 validation/scripts/run_validation.py --org <alias> --report`
- [ ] Update this file with results
- [ ] Commit changes
