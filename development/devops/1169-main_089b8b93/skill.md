# Environment Validation Skill

**🎯 Zweck**: Verhindert Silent Failures durch fehlende oder falsch benannte Environment Variables.

**⚠️ KRITISCH**: Dieser Skill basiert auf Lessons Learned vom 2025-11-21 (OAuth Domain Mismatch Bug).

## Wann wird dieser Skill aktiviert?

- Du arbeitest mit `process.env.*` in TypeScript/Node.js
- Du erstellst oder änderst `.env.example` Dateien
- Du erwähnst Keywords: environment variables, config, secrets, API keys
- Du debuggst OAuth, Authentication oder API-Configuration Probleme

## Goldene Regeln

### 1. Fail-Fast Pattern (KEINE Silent Fallbacks!)

**Problem**: Silent fallbacks verschleiern Konfigurationsfehler.

**SCHLECHT:**
```typescript
// ❌ Bug versteckt: Wenn GOOGLE_REDIRECT_URI fehlt, verwendet Code falsche Domain
const redirectUri = process.env.GOOGLE_REDIRECT_URI || 'https://wrong-domain.com/callback';

// ❌ Sicherheitsrisiko: Admin-Liste mit Fallback
const adminEmails = (process.env.ADMIN_EMAILS || 'default@admin.com').split(',');
```

**GUT:**
```typescript
// ✅ Expliziter Fehler bei fehlender Config
const redirectUri = process.env.GOOGLE_REDIRECT_URI;
if (!redirectUri) {
  throw new Error('FATAL: GOOGLE_REDIRECT_URI not configured in environment variables');
}

// ✅ Fail-Fast mit aussagekräftigem Error
if (!process.env.ADMIN_EMAILS) {
  console.error('CRITICAL: ADMIN_EMAILS missing. Server cannot start.');
  process.exit(1);
}
const adminEmails = process.env.ADMIN_EMAILS.split(',');
```

**Ausnahme - Operational Defaults erlaubt:**
```typescript
// ✅ Safe default für DSGVO-konforme Region
const geminiRegion = process.env.GEMINI_API_REGION || 'europe-west1';

// ✅ Safe default für Timeout
const timeoutMs = parseInt(process.env.VERTEX_TIMEOUT_MS || '55000', 10);

// ✅ Safe default für Model Name
const model = process.env.GEMINI_MODEL_NAME || 'gemini-2.5-flash-lite';
```

**Regel**: Nur DSGVO-sichere, operationelle Defaults erlauben. Niemals für:
- Secrets (API Keys, Passwords)
- Email-Adressen (Admin-Listen, FROM-Adressen)
- URLs/URIs (Redirect-URLs, Webhooks)

---

### 2. RFC-Compliant Naming (OAuth = URI, nicht URL)

**Hintergrund**: OAuth 2.0 RFC 6749 spezifiziert `redirect_uri`, nicht `redirect_url`.

**SCHLECHT:**
```bash
# .env
GOOGLE_REDIRECT_URL=https://...  # ❌ Falsch nach RFC 6749
OAUTH_CALLBACK_URL=https://...   # ❌ Inkonsistent
```

**GUT:**
```bash
# .env
GOOGLE_REDIRECT_URI=https://manufacturing-insights.vercel.app/api/auth/google/callback  # ✅ RFC-konform
MAGIC_LINK_REDIRECT_URI=https://...  # ✅ Konsistent
```

**Regel**: Folge RFCs und Standards:
- OAuth: `redirect_uri`, `client_id`, `client_secret`
- Webhooks: `webhook_url` (kein Standard, aber üblich)
- API Endpoints: `*_url` oder `*_endpoint`

---

### 3. Service-First Naming Convention

**Pattern**: `{SERVICE}_{TYPE}_{DESCRIPTOR}`

**SCHLECHT:**
```bash
# ❌ Unklar welcher Service
API_KEY=xxx
SECRET=yyy
URL=zzz
```

**GUT:**
```bash
# ✅ Service zuerst, dann Typ
GEMINI_API_KEY=xxx
GOOGLE_CLIENT_SECRET=yyy
RESEND_API_KEY=zzz

# ✅ Hierarchisch gruppiert
CHAT_MODEL_NAME=gemini-2.5-flash-lite
CHAT_TIMEOUT_MS=10000
CHAT_USER_AUTO_BLOCK_THRESHOLD=20

# ✅ R2-Variablen mit Ziffern
R2_ACCOUNT_ID=abc123
R2_ACCESS_KEY_ID=xyz789
R2_BUCKET_NAME=security-telemetry
```

**Regel**: Alphabetische Sortierung nach Service ermöglicht Gruppierung in `.env.example`.

---

### 4. Zod-basierte Validation (Startup Check)

**Pattern**: Validiere ALLE env vars beim App-Start, nicht zur Runtime.

**Beispiel - `utils/envValidation.ts`:**
```typescript
import { z } from 'zod';

const envSchema = z.object({
  // Required Secrets
  GOOGLE_CLIENT_ID: z.string().min(1, 'GOOGLE_CLIENT_ID is required'),
  GOOGLE_CLIENT_SECRET: z.string().min(1, 'GOOGLE_CLIENT_SECRET is required'),

  // OAuth URI Validation
  GOOGLE_REDIRECT_URI: z.string().url('Must be valid URL').refine(
    (url) => url.includes('/api/auth/google/callback'),
    'GOOGLE_REDIRECT_URI must end with /api/auth/google/callback'
  ),

  // DSGVO Region Validation
  GEMINI_API_REGION: z.enum(['europe-west1', 'europe-west4', 'europe-north1'])
    .default('europe-west1'),

  // JSON Credentials Validation
  GOOGLE_APPLICATION_CREDENTIALS: z.string().min(100).refine(
    (str) => {
      try {
        const parsed = JSON.parse(str);
        return parsed.type === 'service_account' && parsed.project_id;
      } catch { return false; }
    },
    'GOOGLE_APPLICATION_CREDENTIALS must be valid service account JSON'
  ),
});

export type ValidatedEnv = z.infer<typeof envSchema>;

export function validateEnv(): ValidatedEnv {
  const result = envSchema.safeParse(process.env);

  if (!result.success) {
    const errors = result.error.errors.map(e =>
      `  - ${e.path.join('.')}: ${e.message}`
    ).join('\n');

    throw new Error(
      '🚨 ENVIRONMENT VALIDATION FAILED:\n\n' + errors +
      '\n\nFix .env.local or Vercel Environment Variables.'
    );
  }

  return result.data;
}
```

**Startup Integration - `utils/startup.ts`:**
```typescript
import { validateEnv } from './envValidation';

export function runStartupChecks(): ValidatedEnv {
  console.log('🚀 Running application startup checks...');

  // 1. Environment Validation
  const env = validateEnv();
  console.log('✅ Environment variables validated');

  // 2. DSGVO Compliance Check
  if (!['europe-west1', 'europe-west4', 'europe-north1'].includes(env.GEMINI_API_REGION)) {
    throw new Error('❌ DSGVO VIOLATION: Region must be in EU!');
  }
  console.log(`✅ DSGVO compliant region: ${env.GEMINI_API_REGION}`);

  // 3. Critical Services Reachability (optional)
  // await checkRedisConnection();
  // await checkDatabaseConnection();

  console.log('✅ All startup checks passed');
  return env;
}
```

**API Entry Point - `api/analyze.ts`:**
```typescript
import { runStartupChecks } from '../utils/startup';

// Run checks ONCE at module load (serverless cold start)
const ENV = runStartupChecks();

export default async function handler(req: VercelRequest, res: VercelResponse) {
  // Use validated ENV instead of process.env
  const apiKey = ENV.GEMINI_API_KEY;
  const region = ENV.GEMINI_API_REGION;

  // ...rest of handler
}
```

**Regel**: Fail at startup, not at runtime!

---

### 5. Automated Consistency Checks

**Pattern**: CI/CD prüft automatisch Code vs `.env.example` Konsistenz.

**Script - `scripts/check-env-vars.sh`:**
```bash
#!/bin/bash
set -e

echo "🔍 Checking environment variable consistency..."

# Extract all process.env.VAR_NAME from code (INCLUDE DIGITS for R2_*)
CODE_VARS=$(grep -rh "process\.env\.[A-Z0-9_]\+" api/ utils/ \
  --include="*.ts" | grep -o "process\.env\.[A-Z0-9_]\+" | \
  sed 's/process\.env\.//' | sort -u)

# Extract all VAR_NAME= from .env.example (INCLUDE DIGITS)
EXAMPLE_VARS=$(grep "^[A-Z0-9_]\+=" .env.example | cut -d'=' -f1 | sort -u)

# Find vars in code but not in .env.example
MISSING=""
for var in $CODE_VARS; do
  if ! echo "$EXAMPLE_VARS" | grep -q "^$var$"; then
    MISSING="$MISSING\n  - $var"
  fi
done

if [ -n "$MISSING" ]; then
  echo "❌ Variables used in code but missing from .env.example:"
  echo -e "$MISSING"
  echo ""
  echo "Action: Add these variables to .env.example"
  exit 1
fi

echo "✅ All code variables are documented in .env.example"
```

**Lesson Learned**: Regex MUSS Ziffern matchen (`[A-Z0-9_]+`), sonst wird `R2_ACCOUNT_ID` als `R` erkannt!

**GitHub Actions - `.github/workflows/test.yml`:**
```yaml
jobs:
  env-var-check:
    name: Environment Variable Consistency
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Check env var consistency
        run: bash scripts/check-env-vars.sh
```

---

### 6. .env.example Best Practices

**Struktur:**
```bash
# Manufacturing Insights Analyzer - Environment Variables
# Copy to .env.local: cp .env.example .env.local

# =============================================================================
# GEMINI AI (REQUIRED)
# =============================================================================
GEMINI_API_KEY=your_gemini_api_key_here
GEMINI_API_REGION=europe-west1  # EU Region (DSGVO)
GEMINI_MODEL_NAME=gemini-2.5-flash-lite

# =============================================================================
# GOOGLE OAUTH (REQUIRED)
# =============================================================================
GOOGLE_CLIENT_ID=your_google_client_id.apps.googleusercontent.com
GOOGLE_CLIENT_SECRET=your_google_client_secret
GOOGLE_REDIRECT_URI=http://localhost:3000/api/auth/google/callback

# =============================================================================
# VERCEL AUTO-PROVIDED (DO NOT SET MANUALLY)
# =============================================================================
# Auto-provided by Vercel in production:
VERCEL_ENV=development
VERCEL_URL=localhost:3000
VERCEL_GIT_COMMIT_SHA=abc123def456

# =============================================================================
# OPTIONAL FEATURES
# =============================================================================
# Anonymizer (Presidio Backend)
PRESIDIO_ANALYZE_URL=https://your-presidio-instance.com/analyze
PRESIDIO_API_KEY=your_presidio_api_key

# Cloudflare R2 (Security Telemetry Export)
R2_ACCOUNT_ID=your_cloudflare_account_id
R2_ACCESS_KEY_ID=your_r2_access_key
R2_SECRET_ACCESS_KEY=your_r2_secret_key
```

**Regeln:**
1. Gruppiere logisch mit Kommentaren (`# ===...`)
2. Markiere REQUIRED vs OPTIONAL
3. Gib Beispielwerte (niemals echte Secrets!)
4. Dokumentiere Auto-Provided Variables (VERCEL_*, NODE_ENV)
5. Erkläre Defaults in Kommentaren

---

## Vitest Tests - Hardcoded Fallback Detection

**Datei - `tests/envVarConsistency.test.ts`:**
```typescript
import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';

describe('Environment Variable Consistency', () => {
  it('should NOT have hardcoded fallbacks for environment variables', () => {
    const files = ['api/analyze.ts', 'api/auth.ts', 'api/chat.ts'];
    const fallbackRegex = /process\.env\.([A-Z_]+)\s*\|\|\s*['"`]([^'"`]+)['"`]/g;

    // Allowed operational defaults (DSGVO-safe)
    const allowedDefaults = new Set([
      'GEMINI_API_REGION',
      'GEMINI_MODEL_NAME',
      'VERTEX_TIMEOUT_MS',
      'CHAT_TIMEOUT_MS',
    ]);

    const criticalFallbacks = [];

    for (const file of files) {
      const content = readFileSync(file, 'utf-8');
      let match;

      while ((match = fallbackRegex.exec(content)) !== null) {
        const varName = match[1];
        if (!allowedDefaults.has(varName)) {
          criticalFallbacks.push({ file, varName, line: match[0] });
        }
      }
    }

    if (criticalFallbacks.length > 0) {
      const message = criticalFallbacks.map(f =>
        `  ❌ ${f.file}\n     Variable: ${f.varName}\n     ${f.line}`
      ).join('\n\n');

      throw new Error(
        `🚨 CRITICAL: Hardcoded fallbacks detected!\n\n` +
        `Lesson Learned: Silent fallbacks caused OAuth domain mismatch bug.\n` +
        `Solution: Use explicit error handling instead of || operator.\n\n` +
        `Found ${criticalFallbacks.length} critical fallback(s):\n\n${message}`
      );
    }

    expect(criticalFallbacks).toEqual([]);
  });
});
```

---

## Pre-Commit Hook Integration

**Datei - `.husky/pre-commit`:**
```bash
#!/bin/sh

# Environment Variable Consistency Check
echo "🔍 Checking environment variables..."

if ! bash scripts/check-env-vars.sh; then
  echo ""
  echo "❌ COMMIT BLOCKED: Environment variable issues detected"
  echo ""
  echo "Fix: Add missing variables to .env.example"
  exit 1
fi

echo "✅ Environment variables consistent"
```

---

## Debugging Checklist

Wenn env-var Probleme auftreten:

- [ ] **Vercel Dashboard**: Richtig geschrieben? Typos in Variable-Namen?
- [ ] **Case-Sensitive**: `GOOGLE_REDIRECT_URI` ≠ `google_redirect_uri`
- [ ] **Trailing Spaces**: `"value "` ≠ `"value"`
- [ ] **Quotes**: Vercel braucht KEINE Quotes (nicht `"value"`, nur `value`)
- [ ] **Deployment**: Nach Env-Var-Änderung → Re-Deploy triggern
- [ ] **Logs**: `npx vercel logs --prod` zeigt missing env var errors
- [ ] **Local vs Production**: `.env.local` vs Vercel Dashboard unterschiedlich?

---

## Common Pitfalls

### Pitfall 1: GOOGLE_REDIRECT_URL statt GOOGLE_REDIRECT_URI
```typescript
// ❌ Variable falsch benannt
process.env.GOOGLE_REDIRECT_URI  // undefined!

// Vercel Dashboard hat:
GOOGLE_REDIRECT_URL=https://...  // Falsch!

// ✅ Fix: Rename in Vercel
GOOGLE_REDIRECT_URI=https://...
```

### Pitfall 2: Regex matcht keine Ziffern
```bash
# ❌ Erkennt R2_ACCOUNT_ID nur als "R"
CODE_VARS=$(grep -o "process\.env\.[A-Z_]\+")

# ✅ Inkludiere Ziffern
CODE_VARS=$(grep -o "process\.env\.[A-Z0-9_]\+")
```

### Pitfall 3: Hardcoded Fallback maskiert Missing Var
```typescript
// ❌ Wenn RESEND_FROM_EMAIL fehlt → sendet von falscher Adresse
from: process.env.RESEND_FROM_EMAIL || 'onboarding@resend.dev'

// ✅ Fail explizit
if (!process.env.RESEND_FROM_EMAIL) {
  throw new Error('RESEND_FROM_EMAIL not configured');
}
```

---

## Ressourcen

- [OAuth 2.0 RFC 6749](https://datatracker.ietf.org/doc/html/rfc6749) - redirect_uri Spezifikation
- [12-Factor App](https://12factor.net/config) - Environment-basierte Config
- [Zod Documentation](https://zod.dev/) - TypeScript Schema Validation
- [Vercel Environment Variables](https://vercel.com/docs/environment-variables) - Best Practices

---

## Projektspezifische Anwendung

### ManufacturingInsideAnalyzer
- ✅ Zod Validation in `utils/envValidation.ts`
- ✅ Startup Checks in `utils/startup.ts`
- ✅ Automated CI/CD Check in `.github/workflows/test.yml`
- ✅ Pre-Commit Hook in `.husky/pre-commit`
- ✅ 68 Variablen dokumentiert in `.env.example`

### Lesson Learned Origin
**Date**: 2025-11-21
**Incident**: OAuth 404 after login
**Root Cause**: `GOOGLE_REDIRECT_URL` (wrong) vs `GOOGLE_REDIRECT_URI` (correct)
**Impact**: Login completely broken, masked by hardcoded fallback
**Prevention**: This skill + automated tests

---

**🎯 Ziel**: Zero Silent Failures durch Environment Variables!
