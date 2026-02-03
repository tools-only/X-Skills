# Incident Response & Root Cause Analysis Skill

**🎯 Zweck**: Systematische Fehleranalyse und Dokumentation zur Verhinderung von Regressionen.

**⚠️ KRITISCH**: Basiert auf 3 Production Incidents (2025-11-12, 2025-11-18, 2025-11-21).

## Wann wird dieser Skill aktiviert?

- Production Outage oder kritischer Bug aufgetreten
- Du debuggst komplexe, mehrstufige Fehler
- Du erstellst Incident Reports oder Post-Mortems
- Du implementierst Prevention Measures nach einem Incident
- Keywords: incident, outage, root cause, post-mortem, lessons learned

## Root Cause Analysis Framework

### 5 Whys Methode

**Pattern**: Frage 5× "Warum?" um zur Root Cause zu gelangen.

**Beispiel - OAuth 404 Error (2025-11-21):**

1. **Warum** bekommen User 404 nach Login?
   → Weil `GOOGLE_REDIRECT_URI` nicht definiert ist

2. **Warum** ist `GOOGLE_REDIRECT_URI` nicht definiert?
   → Weil Vercel Env-Var heißt `GOOGLE_REDIRECT_URL` (falsch)

3. **Warum** wurde falscher Name verwendet?
   → Weil kein Standard/RFC befolgt (OAuth 2.0 spezifiziert "URI")

4. **Warum** fiel das nicht früher auf?
   → Weil hardcoded fallback den Fehler maskierte: `process.env.GOOGLE_REDIRECT_URI || 'https://...'`

5. **Warum** gab es hardcoded fallback?
   → Weil keine Fail-Fast Validation beim Startup

**Root Cause**: Fehlende Environment Variable Validation + Silent Fallbacks
**Prevention**: Zod-basierte Startup Validation + Automated Tests

---

### Ishikawa (Fishbone) Diagram

**Pattern**: Kategorisiere Ursachen in: People, Process, Tools, Environment

**Beispiel - CSP Violation (2025-11-18):**

```
                        Website Styling Broken
                                 │
                  ┌──────────────┼──────────────┐
                  │              │              │
              People          Process         Tools
                  │              │              │
            Dev added      No PR Review    Vite didn't
            CDN link      checklist for    prevent CDN
            manually      CSP violations      injection
                  │              │              │
                  └──────────────┼──────────────┘
                                 │
                          Environment
                                 │
                         CSP blocked
                         external CDN
```

**Root Causes Identified**:
1. **People**: Manual CDN link in `index.html`
2. **Process**: No CSP validation in PR review
3. **Tools**: Vite config allowed external URLs
4. **Environment**: Production CSP stricter than dev

**Prevention**: Pre-commit hook + GitHub Actions CSP check

---

### Timeline of Events Template

**Pattern**: Chronologische Auflistung mit Zeitstempeln.

```markdown
| Time | Event | Impact |
|------|-------|--------|
| 13:00 | PR #20 merged (OEE Dashboard) | None |
| 13:30 | Deployment to Vercel completed | None |
| 14:00 | User reports: "Website hat kein Design" | 🔴 HIGH |
| 14:15 | GitHub Actions test workflow fails | 🟡 MEDIUM |
| 15:00 | Investigation started | - |
| 15:10 | Root Cause #1 identified (CSP) | - |
| 15:15 | Root Cause #2 identified (Schema) | - |
| 15:20 | Fixes deployed | - |
| 15:31 | Production verified (fra1 region) | ✅ RESOLVED |
```

**Regel**: Timestamp ALLES (detection, investigation, fix, verification).

---

## Incident Severity Classification

### Severity Matrix

| Severity | Impact | Response Time | Examples |
|----------|--------|---------------|----------|
| 🔴 **CRITICAL** | Complete outage, data loss | Immediate (<30min) | Website unstyled, Login broken, Data leak |
| 🟡 **HIGH** | Major feature broken, workaround exists | <2 hours | API timeout, Payment failing |
| 🟢 **MEDIUM** | Minor feature broken, no user impact | <24 hours | Typo in UI, Slow query |
| 🔵 **LOW** | Cosmetic issue, enhancement | <1 week | Color inconsistency, Missing tooltip |

**Regel**: Überschätze severity (besser 🔴 als 🟢 zu niedrig einschätzen).

---

## Incident Report Template

```markdown
# Incident Report: [Title] ([Date])

## Incident Summary

**Date**: YYYY-MM-DD
**Duration**: X hours
**Severity**: [CRITICAL/HIGH/MEDIUM/LOW]
**Status**: [IN PROGRESS/RESOLVED/MONITORING]

## Symptoms

1. [Observable user-facing issue]
2. [Error messages or logs]
3. [Affected systems or features]

---

## Root Cause Analysis

### Root Cause #1: [Title] 🚨 SEVERITY

**Problem Location**:
- File: `path/to/file.ts:123`
- Component: [API/Frontend/Database/etc.]

**Technical Details**:
```typescript
// ❌ BROKEN CODE
const value = process.env.VAR || 'fallback';
```

**Why This Happened**:
1. [Immediate cause]
2. [Underlying cause]
3. [Systemic issue]

**Impact**:
- [User impact]
- [Business impact]
- [Technical debt created]

**Fix Applied**:
```typescript
// ✅ FIXED CODE
if (!process.env.VAR) {
  throw new Error('VAR not configured');
}
const value = process.env.VAR;
```

---

## Timeline of Events

[See template above]

---

## Resolution

### Fix #1: [Title] (Commit [hash])

**Changes**:
1. [Specific change made]
2. [Files modified]
3. [Configuration updated]

**Result**:
- [Verification step]
- [Test results]
- [Deployment confirmed]

---

## Prevention Measures

### Immediate Actions (Completed ✅)

- [x] [Action taken immediately]
- [x] [Emergency fix deployed]

### Process Improvements (To Implement 📋)

#### 1. [Improvement Title]

**Implementation**:
```yaml
# Example: GitHub Actions check
- name: Check for CSP violations
  run: |
    if grep -r "cdn.tailwindcss.com" *.html; then
      echo "ERROR: External CDN detected"
      exit 1
    fi
```

**Timeline**: [Date]
**Owner**: [Team/Person]

---

## Key Learnings

### What Went Wrong

1. [Mistake #1]
2. [Mistake #2]

### What Went Right

1. [Good practice #1]
2. [Fast detection/response]

---

## Technical Debt Identified

1. [Debt item #1] - Priority: [HIGH/MEDIUM/LOW]
2. [Debt item #2] - Priority: [HIGH/MEDIUM/LOW]

---

## Follow-Up Actions

| Priority | Action | Owner | Deadline | Status |
|----------|--------|-------|----------|--------|
| 🔴 HIGH | [Action] | [Owner] | [Date] | [TODO/IN PROGRESS/DONE] |

---

## Stakeholder Communication

**Internal**: [How team was notified]
**External**: [Customer communication, if any]

---

## Conclusion

[1-2 paragraph summary of incident, fixes, and prevention]

---

**Report Author**: [Name]
**Review Date**: [Date]
**Status**: ✅ Resolved & Documented
```

---

## Automated Prevention Measures

### 1. Pre-Commit Hook - Schema Change Detection

**Pattern**: Husky hook warnt bei Schema-Änderungen ohne Test-Updates.

**File**: `.husky/pre-commit`
```bash
#!/bin/bash

# Schema Change Detection
if git diff --cached --name-only | grep -q "utils/validation.ts"; then
  echo "⚠️  SCHEMA CHANGE DETECTED"
  echo ""
  echo "   File changed: utils/validation.ts"
  echo "   Action required: Update utils/validation.test.ts"
  echo ""

  if ! git diff --cached --name-only | grep -q "utils/validation.test.ts"; then
    echo "   ❌ ERROR: utils/validation.test.ts NOT staged"
    echo ""
    printf "   Did you update the tests? (y/n): "
    read REPLY
    echo ""

    case "$REPLY" in
      [Yy]) ;;
      *)
        echo "   ❌ COMMIT BLOCKED: Update tests first"
        exit 1
        ;;
    esac
  else
    echo "   ✅ utils/validation.test.ts is staged"
  fi
fi
```

**Lesson Learned**: OEE Schema Migration (2025-11-18) broke tests weil Updates fehlten.

---

### 2. GitHub Actions - CSP Validation

**Pattern**: CI/CD prüft static HTML files auf external CDN links.

**File**: `.github/workflows/test.yml`
```yaml
jobs:
  csp-validation:
    name: CSP Compliance Check
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Check for CSP violations
        run: |
          echo "🔍 Checking static HTML files for external CDN..."

          PROHIBITED_CDNS="cdn.tailwindcss.com cdn.jsdelivr.net unpkg.com cdnjs.cloudflare.com"
          VIOLATIONS=0

          for FILE in index.html public/*.html; do
            for CDN in $PROHIBITED_CDNS; do
              if grep -q "$CDN" "$FILE" 2>/dev/null; then
                echo "❌ CSP VIOLATION in $FILE: $CDN"
                VIOLATIONS=1
              fi
            done
          done

          if [ "$VIOLATIONS" -eq 1 ]; then
            echo ""
            echo "❌ COMMIT BLOCKED: External CDN usage detected"
            echo ""
            echo "Fix: Remove external CDN links, use inline CSS or Vite bundling"
            exit 1
          fi

          echo "✅ No CSP violations detected"
```

**Lesson Learned**: Tailwind CDN (2025-11-18) broke production website styling.

---

### 3. Smoke Tests - Post-Deployment Verification

**Pattern**: Automatische Health Checks nach jedem Deployment.

**File**: `scripts/smoke-test.sh`
```bash
#!/bin/bash
set -e

DEPLOYMENT_URL="${1:-https://manufacturing-insights.vercel.app}"

echo "🧪 Running smoke tests for $DEPLOYMENT_URL"
echo ""

# Test 1: Homepage returns 200
echo "Test 1: Homepage accessibility..."
STATUS=$(curl -s -o /dev/null -w "%{http_code}" "$DEPLOYMENT_URL")
if [ "$STATUS" -ne 200 ]; then
  echo "❌ FAILED: Homepage returned $STATUS"
  exit 1
fi
echo "✅ Homepage: 200 OK"

# Test 2: API health endpoint
echo "Test 2: API health check..."
STATUS=$(curl -s -o /dev/null -w "%{http_code}" "$DEPLOYMENT_URL/api/health")
if [ "$STATUS" -ne 200 ]; then
  echo "❌ FAILED: API health returned $STATUS"
  exit 1
fi
echo "✅ API health: 200 OK"

# Test 3: OAuth init endpoint
echo "Test 3: OAuth initialization..."
STATUS=$(curl -s -o /dev/null -w "%{http_code}" "$DEPLOYMENT_URL/api/auth?action=google")
if [ "$STATUS" -ne 302 ] && [ "$STATUS" -ne 200 ]; then
  echo "❌ FAILED: OAuth init returned $STATUS"
  exit 1
fi
echo "✅ OAuth init: $STATUS (redirect)"

# Test 4: DSGVO region verification
echo "Test 4: DSGVO region check..."
if command -v vercel &> /dev/null; then
  REGION=$(npx vercel inspect "$DEPLOYMENT_URL" --wait 2>&1 | grep -o "\[fra1\]")
  if [ -z "$REGION" ]; then
    echo "❌ CRITICAL: Deployment NOT in Frankfurt (fra1)!"
    exit 1
  fi
  echo "✅ DSGVO compliant: Frankfurt (fra1)"
else
  echo "⚠️  Vercel CLI not available, skipping region check"
fi

echo ""
echo "✅ All smoke tests passed"
```

**Usage**:
```bash
# After deployment
bash scripts/smoke-test.sh https://manufacturing-insights.vercel.app
```

---

## Pull Request Template with Incident Prevention

**Pattern**: Checklist zwingt Developers zu systematischem Review.

**File**: `.github/PULL_REQUEST_TEMPLATE.md`
```markdown
## Description

[Brief description of changes]

## Type of Change

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Configuration/Infrastructure change

## Pre-Merge Checklist

### Testing ✅
- [ ] All tests passing locally (`npm run test`)
- [ ] Build successful (`npm run build`)
- [ ] No TypeScript errors (`npx tsc --noEmit`)

### Security & Compliance 🔒
- [ ] **NO hardcoded fallbacks** for environment variables
- [ ] **NO console.logs** with PII (emails, user data)
- [ ] **NO secrets** in code (use .env or Vercel env vars)
- [ ] OAuth redirects use "URI" not "URL" (RFC 6749)

### Performance & UX 🚀
- [ ] **Max 2 nested layout divs** (no conflicting styles)
- [ ] **Static HTML uses inline CSS** (no external CDN)
- [ ] **Processing time < 10s** for API endpoints
- [ ] **Mobile responsive** checked

### Schema Changes 📊
- [ ] If `utils/validation.ts` changed → `utils/validation.test.ts` updated
- [ ] Backwards compatibility maintained OR migration guide provided
- [ ] API contracts not broken

### DSGVO Compliance 🇪🇺
- [ ] `vercel.json` has `regions: ["fra1"]` (if infrastructure change)
- [ ] No data stored server-side (except session with TTL)
- [ ] Presidio anonymization for user inputs (if applicable)

### Deployment Verification 🚢
- [ ] Environment variables documented in `.env.example`
- [ ] Smoke tests will pass (homepage, API health, OAuth)
- [ ] No breaking changes to production

## Incident Prevention

**Does this PR address a previous incident?**
- [ ] Yes (link incident report): [INCIDENT_REPORT_YYYY-MM-DD.md](...)
- [ ] No

**Prevention measures implemented:**
- [x] [Specific prevention measure]
- [x] [Automated test added]

## Verification Steps

1. [Step to reproduce/verify]
2. [Expected behavior]

## Screenshots (if applicable)

[Add screenshots]

---

**⚠️ CRITICAL**: If ANY checkbox is unchecked, explain why in comments!
```

---

## Monitoring & Alerting Setup

### Sentry Integration (Error Tracking)

**Pattern**: Capture errors + CSP violations automatisch.

```typescript
// utils/sentry.ts
import * as Sentry from '@sentry/nextjs';

Sentry.init({
  dsn: process.env.SENTRY_DSN,
  environment: process.env.VERCEL_ENV || 'development',

  // Capture 100% of errors in production
  tracesSampleRate: 1.0,

  // Capture CSP violations
  beforeSend(event) {
    // Log CSP violations
    if (event.exception?.values?.[0]?.type === 'CSPViolation') {
      console.error('[CSP VIOLATION]', event);
    }
    return event;
  },

  // Don't send PII
  beforeBreadcrumb(breadcrumb) {
    if (breadcrumb.data?.email || breadcrumb.data?.user) {
      delete breadcrumb.data.email;
      delete breadcrumb.data.user;
    }
    return breadcrumb;
  },
});
```

---

## Lessons Learned Documentation Format

**Pattern**: Sofort nach Fix dokumentieren (nicht später!).

**File**: `LESSONS_LEARNED_[TOPIC]_[DATE].md`

```markdown
# Lessons Learned: [Topic] ([Date])

**Status**: ✅ RESOLVED

---

## 🔍 Root Causes (X Critical Issues)

### 1. [Root Cause Title]

**Symptom**: [Observable behavior]

**Root Cause**: [Technical explanation]

**What Happened**:
```typescript
// ❌ BROKEN
[broken code]

// ✅ FIXED
[fixed code]
```

**Files Fixed**:
- `path/to/file.ts:123`
- `path/to/other.ts:456`

**Key Insight**: [1-sentence takeaway]

---

## 📋 Complete Solution

[Detailed implementation with code examples]

---

## 🚫 Common Mistakes

### Mistake 1: [Title]
❌ [What not to do]
✅ [What to do instead]

---

## ✅ Testing Checklist

- [ ] [Verification step 1]
- [ ] [Verification step 2]

---

## 💡 Key Takeaways

1. [Lesson 1]
2. [Lesson 2]
3. [Lesson 3]

---

**Final Commit**: `[hash]` - "[commit message]"
**Status**: ✅ [Feature/Bug] now works consistently
```

---

## Post-Incident Review Meeting

**Pattern**: Team Review innerhalb 48h nach Incident Resolution.

**Agenda**:
1. **Timeline Review** (10 min)
   - Was geschah wann?
   - Wie lange bis Detection?
   - Wie lange bis Resolution?

2. **Root Cause Analysis** (15 min)
   - 5 Whys walkthrough
   - Contributing factors
   - Systemic issues identified

3. **What Went Well** (5 min)
   - Fast detection
   - Good communication
   - Effective debugging

4. **What Can Improve** (10 min)
   - Process gaps
   - Tool limitations
   - Knowledge gaps

5. **Action Items** (10 min)
   - Assign owners
   - Set deadlines
   - Prioritize (HIGH/MEDIUM/LOW)

6. **Documentation** (5 min)
   - Incident report completed?
   - Lessons learned documented?
   - Runbooks updated?

**Output**: Action Items mit Deadlines

---

## Incident Metrics to Track

```typescript
// Tracking metrics
interface IncidentMetrics {
  // Detection
  timeToDetect: number;        // How long until we knew?
  detectionMethod: string;      // User report, monitoring, logs?

  // Response
  timeToAcknowledge: number;   // How long until someone started?
  timeToIdentify: number;      // How long to find root cause?
  timeToResolve: number;       // How long to deploy fix?

  // Impact
  usersAffected: number;       // How many users impacted?
  severity: 'CRITICAL' | 'HIGH' | 'MEDIUM' | 'LOW';

  // Prevention
  preventable: boolean;        // Could automation have caught this?
  recurringIssue: boolean;     // Have we seen this before?
}
```

**Goal**: Reduce MTTR (Mean Time To Resolution) über Zeit.

---

## Ressourcen

- [Google SRE Book - Postmortem Culture](https://sre.google/sre-book/postmortem-culture/)
- [Atlassian Incident Handbook](https://www.atlassian.com/incident-management/handbook)
- [PagerDuty Incident Response](https://response.pagerduty.com/)
- [Etsy Debriefing Facilitation Guide](https://extfiles.etsy.com/DebriefingFacilitationGuide.pdf)

---

## Projektspezifische Anwendung

### ManufacturingInsideAnalyzer

**Documented Incidents:**
1. OAuth Cookie Loop (2025-11-12) → `LESSONS_LEARNED_OAUTH.md`
2. CSP + Schema Migration (2025-11-18) → `INCIDENT_REPORT_2025-11-18.md`
3. OAuth Domain Mismatch (2025-11-21) → Integrated in prevention measures

**Prevention Measures Implemented:**
- ✅ Pre-commit hooks (schema detection, CSP checks)
- ✅ GitHub Actions (env-var consistency, CSP validation)
- ✅ Smoke tests (post-deployment verification)
- ✅ PR template (security checklist)
- ✅ Zod validation (environment variables)

**MTTR Reduction:**
- 2025-11-12: 4 hours (3 root causes)
- 2025-11-18: 2 hours (2 root causes)
- 2025-11-21: 1 hour (1 root cause + immediate prevention)

---

**🎯 Ziel**: Learn from EVERY incident, prevent ALL regressions!
