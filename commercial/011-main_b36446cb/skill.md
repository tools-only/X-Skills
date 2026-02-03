# DSGVO-Compliance Skill

**⚠️ KRITISCH**: Dieser Skill muss bei allen DSGVO-relevanten Änderungen aktiv sein!

## Goldene Regeln

### 1. Frankfurt europe-west3 Region PFLICHT
```json
// Cloud Run
{
  "regions": ["europe-west3"]  // ✅ IMMER Frankfurt
}
```

### 2. Kein PII-Logging
```ts
// ❌ NIEMALS
console.log('User email:', user.email);

// ✅ Nur Metadaten
console.log('Analysis completed:', { duration, status });
```

### 3. Anonymisierung VOR Verarbeitung
```ts
// ✅ Richtige Reihenfolge
const anonymized = await anonymizeData(rawData);
const result = await analyzeData(anonymized);
```

### 4. Klare Datenretenzen
- Session-Daten: Firestore mit TTL (24h max)
- Analysen: Nicht gespeichert
- Logs: Nur Metadaten, 7 Tage Retention

## ManufacturingInsideAnalyzer Spezifika

### Presidio Backend
```python
# api/presidio-analyze.py
from presidio_analyzer import AnalyzerEngine

analyzer = AnalyzerEngine()
results = analyzer.analyze(
    text=input_text,
    language='de',  # Deutsche PII-Patterns
    entities=['PERSON', 'EMAIL', 'PHONE_NUMBER']
)
```

### Rate Limiting
```ts
// 5 Analysen pro E-Mail (Free Tier)
const remaining = await firestore.get(`rate:${email}`);
if (remaining <= 0) throw new Error('Rate limit exceeded');
```

### Audit Logging (nur Metadaten)
```ts
await logAuditEvent({
  action: 'analysis_completed',
  userId: hashedUserId,  // ✅ Gehasht
  timestamp: Date.now(),
  status: 'success'
  // ❌ KEIN fileContent, KEINE userData
});
```

## Checkliste vor Deployment
- [ ] `Cloud Run` hat `regions: ["europe-west3"]`
- [ ] Keine console.logs mit PII
- [ ] Presidio aktiviert für alle User-Inputs
- [ ] Rate Limits getestet
- [ ] Audit Logs nur Metadaten

## Deployment Verification (KRITISCH - Lesson Learned 2025-11-18)

### Problem: Region Drift nach Re-Deployment

**Symptom**: Deployment erfolgreich, aber Funktionen laufen plötzlich in US-Region (iad1, sfo1).

**Root Cause**: Cloud Run kann Region aendern bei:
- Git-basierten Deployments ohne explicit regions
- Rollbacks zu älteren Deployments
- Manual changes via Console

**PFLICHT-Checks nach JEDEM Deployment**:

```bash
# 1. Deploy to production
gcloud run deploy fabrikiq-backend --source backend/ --region europe-west3 --platform managed

# 2. IMMEDIATELY verify region (DO NOT SKIP!)
gcloud run services describe fabrikiq-backend --region europe-west3

# Expected output:
# Service Info:
# ├── λ api/analyze (64.13KB) [europe-west3] ✅
# ├── λ api/auth (12.75KB) [europe-west3] ✅
# ├── λ api/chat (48.92KB) [europe-west3] ✅
# └── λ api/admin (23.45KB) [europe-west3] ✅

# ❌ If you see [iad1], [sfo1], [sin1] → DSGVO VIOLATION!
```

**Automated Verification Script**:
```bash
# scripts/verify-dsgvo-region.sh
#!/bin/bash
set -e

echo "🔍 Verifying DSGVO compliance (Frankfurt region)..."

DEPLOYMENT_URL="https://app.your-domain.com"
INSPECT_OUTPUT=$(gcloud run services describe fabrikiq-backend --region europe-west3 "$DEPLOYMENT_URL" --wait 2>&1)

# Check if ALL lambda functions are in europe-west3
if echo "$INSPECT_OUTPUT" | grep -q "\[europe-west3\]"; then
  echo "✅ All functions in Frankfurt (europe-west3) region"
else
  echo "❌ CRITICAL: Functions NOT in Frankfurt!"
  echo "$INSPECT_OUTPUT"
  exit 1
fi

# Check for non-EU regions
if echo "$INSPECT_OUTPUT" | grep -E "\[(iad1|sfo1|sin1|hnd1|bom1)\]"; then
  echo "❌ DSGVO VIOLATION: Non-EU region detected!"
  echo "$INSPECT_OUTPUT"
  exit 1
fi

echo "✅ DSGVO compliance verified"
```

**GitHub Actions Integration**:
```yaml
# .github/workflows/dsgvo-verification.yml
name: DSGVO Region Verification

on:
  deployment_status:

jobs:
  verify-region:
    if: github.event.deployment_status.state == 'success'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Verify Frankfurt region
        run: bash scripts/verify-dsgvo-region.sh
        env:
          GCP_SA_KEY: ${{ secrets.GCP_SA_KEY }}

      - name: Notify on violation
        if: failure()
        uses: actions/github-script@v7
        with:
          script: |
            github.rest.issues.create({
              owner: context.repo.owner,
              repo: context.repo.repo,
              title: '🚨 DSGVO VIOLATION: Non-EU Region Detected',
              body: 'Deployment verification failed. Functions NOT in Frankfurt region.',
              labels: ['critical', 'dsgvo', 'security']
            })
```

**Regel**:
- **NIEMALS** Deployment als fertig betrachten ohne Region-Verification
- **Automated alerts** bei Region-Drift (GitHub Issues + Email)
- **Rollback** bei Non-EU Regions sofort durchführen

---

## CSP (Content Security Policy) Compliance

### Problem: External CDNs blockieren DSGVO-konforme Hosting

**Hintergrund**: DSGVO Art. 44 verbietet Datenübertragung in Drittländer ohne Angemessenheitsbeschluss.

**Problematische CDNs**:
- ❌ `https://cdn.tailwindcss.com` (Cloudflare, global)
- ❌ `https://unpkg.com` (Cloudflare, global)
- ❌ `https://cdn.jsdelivr.net` (Cloudflare, global)
- ❌ `https://fonts.googleapis.com` (Google, USA)
- ❌ `https://ajax.googleapis.com` (Google, USA)

**Erlaubte Quellen**:
- ✅ `'self'` (eigene Domain)
- ✅ `run.app` (EU-Region konfiguriert)
- ✅ Lokale npm packages

**CSP Header Configuration**:

> ⚠️ **SECURITY**: Verwende NIEMALS `unsafe-inline` oder `unsafe-eval` in Production!
> Nutze stattdessen Nonces oder Hashes für Inline-Scripts.

```json
// Cloud Run - PRODUCTION (sicher)
{
  "regions": ["europe-west3"],
  "headers": [
    {
      "source": "/(.*)",
      "headers": [
        {
          "key": "Content-Security-Policy",
          "value": "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data: https:; font-src 'self' data:; connect-src 'self' https://*.run.app https://*.google.com; frame-ancestors 'none'; base-uri 'self'; form-action 'self'"
        },
        {
          "key": "X-Content-Type-Options",
          "value": "nosniff"
        },
        {
          "key": "X-Frame-Options",
          "value": "DENY"
        },
        {
          "key": "Referrer-Policy",
          "value": "strict-origin-when-cross-origin"
        }
      ]
    }
  ]
}
```

**Explanation**:
- `default-src 'self'` → Nur eigene Domain (kein CDN!)
- `script-src 'self'` → Keine unsafe-* in Production! Für Dev: Vite Config anpassen
- `connect-src ... https://*.google.com` → Google OAuth + Gemini API (EU-hosted)
- `frame-ancestors 'none'` → Prevent clickjacking
- `X-Frame-Options: DENY` → Double protection gegen iframes

**Für Vite Dev-Mode** (nur lokal, NICHT in Production):
```ts
// vite.config.ts - nur für lokale Entwicklung
export default defineConfig({
  server: {
    headers: {
      'Content-Security-Policy': "script-src 'self' 'unsafe-inline' 'unsafe-eval'"
    }
  }
})
```

**Verification Commands**:
```bash
# 1. Check CSP headers in production
curl -I https://app.your-domain.com | grep -i content-security-policy

# Expected output:
# content-security-policy: default-src 'self'; ...

# 2. Verify no external CDN requests
curl -s https://app.your-domain.com | grep -i "cdn\."

# Expected: NO output (keine CDN-Links)

# 3. Test with browser DevTools
# Open DevTools → Network Tab → Reload
# Filter: "cdn" → Should be EMPTY
```

**Pre-Commit Hook für CDN-Detection**:
```bash
# .husky/pre-commit
#!/bin/sh

echo "🔍 Checking for external CDN links..."

# Search for CDN links in HTML/TSX files
if grep -r "cdn\." src/ index.html public/ 2>/dev/null | grep -v node_modules; then
  echo ""
  echo "❌ DSGVO VIOLATION: External CDN detected!"
  echo ""
  echo "Found external CDN links in the following files:"
  grep -rn "cdn\." src/ index.html public/ 2>/dev/null | grep -v node_modules
  echo ""
  echo "Action: Remove CDN links and use local npm packages instead."
  echo "Example: npm install -D tailwindcss (instead of <script src='https://cdn.tailwindcss.com'>)"
  exit 1
fi

echo "✅ No external CDN links found"
```

**Regel**:
- **NIEMALS** externe CDNs in Production
- **CSP Headers** nach jedem Config-Change prüfen
- **Pre-commit hooks** für automatische Detection
- **Browser DevTools** Network Tab nach Deployment checken

---

## DSGVO-konforme AI Model Configuration

### Google Gemini API (Vertex AI) - EU Region

**PFLICHT-Konfiguration**:
```typescript
// api/analyze.ts
import { VertexAI } from '@google-cloud/vertexai';

// ✅ IMMER Europe-West1 Region (Belgien)
const vertexAI = new VertexAI({
  project: process.env.GOOGLE_PROJECT_ID,
  location: 'europe-west1',  // CRITICAL: EU Region!
});

// Get Gemini model
const model = vertexAI.getGenerativeModel({
  model: 'gemini-2.0-flash',
  generationConfig: {
    temperature: 0.7,
    maxOutputTokens: 2048,
  }
});
```

**Erlaubte EU Regions** (Google Cloud):
- ✅ `europe-west1` (Belgien) - Standard
- ✅ `europe-west4` (Niederlande)
- ✅ `europe-north1` (Finnland)
- ❌ `us-central1` (Iowa, USA) - VERBOTEN
- ❌ `asia-southeast1` (Singapur) - VERBOTEN

**Environment Variable Validation**:
```typescript
// utils/envValidation.ts
import { z } from 'zod';

const envSchema = z.object({
  GEMINI_API_REGION: z.enum(['europe-west1', 'europe-west4', 'europe-north1'])
    .default('europe-west1'),
});

export function validateEnv() {
  const result = envSchema.safeParse(process.env);

  if (!result.success) {
    throw new Error(
      '🚨 DSGVO VIOLATION: Invalid Gemini API region!\n' +
      'Allowed regions: europe-west1, europe-west4, europe-north1'
    );
  }

  return result.data;
}
```

### Anthropic Claude (Vertex AI) - EU Region

**PFLICHT-Konfiguration**:
```typescript
// api/chat.ts
import Anthropic from '@anthropic-ai/sdk';

// ✅ Claude via Google Vertex AI (EU Region)
const anthropic = new Anthropic({
  baseURL: `https://europe-west1-aiplatform.googleapis.com/v1/projects/${process.env.GOOGLE_PROJECT_ID}/locations/europe-west1/publishers/anthropic/models`,
  apiKey: process.env.GOOGLE_APPLICATION_CREDENTIALS,
});

const message = await anthropic.messages.create({
  model: 'claude-3-5-haiku@20241022',
  max_tokens: 500,
  messages: [{ role: 'user', content: userQuestion }],
});
```

**Wichtig**:
- **Nie** direkt `https://api.anthropic.com` verwenden (USA!)
- **Immer** über Google Vertex AI (EU Region)
- **baseURL** MUSS `europe-west1` enthalten

**Startup Check**:
```typescript
// utils/startup.ts
export function runStartupChecks() {
  console.log('🚀 Running DSGVO compliance checks...');

  const env = validateEnv();

  // 1. Check Gemini Region
  if (!['europe-west1', 'europe-west4', 'europe-north1'].includes(env.GEMINI_API_REGION)) {
    throw new Error('❌ DSGVO VIOLATION: Gemini API not in EU region!');
  }
  console.log(`✅ Gemini API: ${env.GEMINI_API_REGION}`);

  // 2. Check Claude Vertex Region (from baseURL)
  if (process.env.CLAUDE_VERTEX_REGION !== 'europe-west1') {
    throw new Error('❌ DSGVO VIOLATION: Claude Vertex not in EU region!');
  }
  console.log(`✅ Claude Vertex: ${process.env.CLAUDE_VERTEX_REGION}`);

  console.log('✅ All DSGVO compliance checks passed');
}
```

**Regel**:
- **Google Gemini**: Nur `europe-west1`, `europe-west4`, `europe-north1`
- **Anthropic Claude**: Nur via Vertex AI mit EU Region
- **OpenAI GPT**: NICHT erlaubt (keine EU-Region verfügbar!)
- **Startup Checks**: Region-Validation BEFORE erste API-Calls

---

## Data Processing Transparency (DSGVO Art. 13/14)

### Informationspflichten beim Datenerhebung

**Minimale Angaben auf UI** (datenschutz.html):
```html
<!-- public/datenschutz.html -->
<h2>1. Datenverarbeitung bei der Analyse</h2>

<h3>Verarbeitete Daten</h3>
<ul>
  <li><strong>Hochgeladene Dateien</strong>: CSV, Excel, JSON (maximal 5 MB)</li>
  <li><strong>Verarbeitungszweck</strong>: KI-gestützte Qualitätsanalyse, Fehlererkennung, Predictive Maintenance</li>
  <li><strong>Rechtsgrundlage</strong>: Art. 6 Abs. 1 lit. b DSGVO (Vertragserfüllung)</li>
  <li><strong>Speicherdauer</strong>: Keine serverseitige Speicherung. Ergebnisse nur in Ihrem Browser (LocalStorage).</li>
</ul>

<h3>KI-Modell und Region</h3>
<ul>
  <li><strong>Anbieter</strong>: Google Gemini 2.0 Flash (via Vertex AI)</li>
  <li><strong>Region</strong>: europe-west1 (Belgien, EU)</li>
  <li><strong>Datenübertragung</strong>: Keine Übertragung außerhalb der EU</li>
  <li><strong>Anonymisierung</strong>: Optional verfügbar (Microsoft Presidio)</li>
</ul>

<h3>Session-Daten</h3>
<ul>
  <li><strong>Gespeichert in</strong>: Firestore (Redis, Frankfurt Region)</li>
  <li><strong>Speicherdauer</strong>: 24 Stunden (automatische Löschung)</li>
  <li><strong>Zweck</strong>: Authentifizierung, Rate Limiting (5 Analysen pro E-Mail)</li>
</ul>

<h3>Ihre Rechte (Art. 15-22 DSGVO)</h3>
<ul>
  <li>Auskunft über gespeicherte Daten (Art. 15)</li>
  <li>Berichtigung falscher Daten (Art. 16)</li>
  <li>Löschung ("Recht auf Vergessenwerden", Art. 17)</li>
  <li>Einschränkung der Verarbeitung (Art. 18)</li>
  <li>Datenübertragbarkeit (Art. 20)</li>
  <li>Widerspruch gegen Verarbeitung (Art. 21)</li>
</ul>

<p><strong>Kontakt</strong>: <a href="mailto:privacy@your-domain.com">privacy@your-domain.com</a></p>
```

**User-Consent Banner** (NICHT Cookie-Banner!):
```tsx
// components/AnalysisConsent.tsx
import { useState } from 'react';

export function AnalysisConsent({ onAccept }: { onAccept: () => void }) {
  const [accepted, setAccepted] = useState(() => {
    return localStorage.getItem('analysis-consent') === 'true';
  });

  if (accepted) return null;

  const handleAccept = () => {
    localStorage.setItem('analysis-consent', 'true');
    setAccepted(true);
    onAccept();
  };

  return (
    <div className="fixed bottom-4 left-4 right-4 bg-slate-800 border border-slate-600 rounded-lg p-4 shadow-lg max-w-2xl mx-auto">
      <h3 className="text-lg font-semibold text-white mb-2">Datenverarbeitung</h3>
      <p className="text-sm text-slate-300 mb-4">
        Ihre hochgeladenen Dateien werden mit <strong>Google Gemini AI (EU-Region)</strong> analysiert.
        Keine serverseitige Speicherung. Ergebnisse nur in Ihrem Browser.
      </p>
      <div className="flex gap-2">
        <button
          onClick={handleAccept}
          className="px-4 py-2 bg-blue-600 text-white rounded hover:bg-blue-700"
        >
          Verstanden
        </button>
        <a
          href="/datenschutz.html"
          target="_blank"
          className="px-4 py-2 border border-slate-600 text-slate-200 rounded hover:bg-slate-700"
        >
          Mehr erfahren
        </a>
      </div>
    </div>
  );
}
```

**Regel**:
- **Transparenz**: Klare Angabe zu Zweck, Rechtsgrundlage, Speicherdauer
- **EU-Region**: Explizit erwähnen (Vertrauen schaffen)
- **User Rights**: Alle DSGVO-Rechte (Art. 15-22) klar benennen
- **Kontakt**: Datenschutzbeauftragter oder Privacy-Email angeben

---

## Audit Logging Best Practices

### Was DARF geloggt werden (Metadaten)

**Erlaubte Log-Daten**:
```typescript
// utils/auditLog.ts
interface AuditEvent {
  timestamp: number;              // ✅ Zeitstempel
  action: string;                 // ✅ z.B. 'analysis_completed'
  userId: string;                 // ✅ Gehashter User-ID (SHA-256)
  duration: number;               // ✅ Verarbeitungszeit in ms
  status: 'success' | 'error';    // ✅ Ergebnisstatus
  errorCode?: string;             // ✅ z.B. '504_TIMEOUT'
  fileSize?: number;              // ✅ Dateigröße in Bytes
  rowCount?: number;              // ✅ Anzahl Zeilen
  columnCount?: number;           // ✅ Anzahl Spalten
  analysisType?: string;          // ✅ z.B. 'quality_analysis'
}
```

**VERBOTENE Log-Daten**:
```typescript
// ❌ NIEMALS loggen:
interface ForbiddenData {
  email: string;           // ❌ E-Mail-Adresse (PII)
  name: string;            // ❌ Benutzername (PII)
  fileContent: string;     // ❌ Dateiinhalt
  analysisResult: string;  // ❌ KI-Ergebnisse
  ipAddress: string;       // ❌ IP-Adresse (personenbezogen)
  sessionToken: string;    // ❌ Session-Cookie
}
```

**Implementation mit Hash**:

> ⚠️ **SECURITY**: Mindestens 32 Zeichen (128 bit) für ausreichende Kollisionsresistenz!

```typescript
// utils/auditLog.ts
import crypto from 'crypto';

function hashUserId(email: string): string {
  return crypto
    .createHash('sha256')
    .update(email.toLowerCase() + process.env.HASH_SALT)
    .digest('hex')
    .slice(0, 32); // 32 chars = 128 bit (Kollisionsresistent)
}

export async function logAnalysis(email: string, metadata: Omit<AuditEvent, 'userId' | 'timestamp'>) {
  const auditEvent: AuditEvent = {
    timestamp: Date.now(),
    userId: hashUserId(email), // ✅ Hashed
    ...metadata,
  };

  // Log to Firestore with 7-day TTL
  await redis.setex(
    `audit:${auditEvent.timestamp}:${auditEvent.userId}`,
    7 * 24 * 60 * 60, // 7 days
    JSON.stringify(auditEvent)
  );

  console.log('[AUDIT]', JSON.stringify(auditEvent));
}
```

**Usage Example**:
```typescript
// api/analyze.ts
try {
  const startTime = Date.now();
  const result = await geminiAPI.generateContent(prompt);
  const duration = Date.now() - startTime;

  await logAnalysis(userEmail, {
    action: 'analysis_completed',
    duration,
    status: 'success',
    fileSize: file.size,
    rowCount: data.length,
    columnCount: Object.keys(data[0]).length,
    analysisType: 'quality',
  });

} catch (error) {
  await logAnalysis(userEmail, {
    action: 'analysis_failed',
    duration: Date.now() - startTime,
    status: 'error',
    errorCode: error.code,
  });
}
```

**Regel**:
- **Metadaten only**: Keine PII (E-Mail, Name, IP)
- **Hashed User IDs**: SHA-256 mit Salt
- **TTL 7 days**: Automatische Löschung nach 1 Woche
- **Purpose Limitation**: Nur für Fehleranalyse, nicht für Marketing

---

## Post-Deployment Smoke Tests

### Automated DSGVO Smoke Test Suite

```bash
# scripts/dsgvo-smoke-test.sh
#!/bin/bash
set -e

PROD_URL="https://app.your-domain.com"

echo "🧪 Running DSGVO Smoke Tests..."
echo ""

# Test 1: Region Verification
echo "Test 1: Verifying Frankfurt region..."
INSPECT_OUTPUT=$(gcloud run services describe fabrikiq-backend --region europe-west3 "$PROD_URL" --wait 2>&1)
if echo "$INSPECT_OUTPUT" | grep -q "\[europe-west3\]"; then
  echo "✅ PASS: All functions in Frankfurt"
else
  echo "❌ FAIL: Functions NOT in Frankfurt!"
  exit 1
fi

# Test 2: CSP Headers
echo "Test 2: Checking CSP headers..."
CSP_HEADER=$(curl -sI "$PROD_URL" | grep -i content-security-policy)
if echo "$CSP_HEADER" | grep -q "default-src 'self'"; then
  echo "✅ PASS: CSP header present"
else
  echo "❌ FAIL: CSP header missing or incorrect!"
  exit 1
fi

# Test 3: No External CDNs
echo "Test 3: Checking for external CDNs..."
HTML_CONTENT=$(curl -s "$PROD_URL")
if echo "$HTML_CONTENT" | grep -qi "cdn\."; then
  echo "❌ FAIL: External CDN detected!"
  echo "$HTML_CONTENT" | grep -i "cdn\."
  exit 1
else
  echo "✅ PASS: No external CDNs"
fi

# Test 4: Privacy Policy Accessible
echo "Test 4: Checking datenschutz.html..."
PRIVACY_STATUS=$(curl -s -o /dev/null -w "%{http_code}" "$PROD_URL/datenschutz.html")
if [ "$PRIVACY_STATUS" -eq 200 ]; then
  echo "✅ PASS: Privacy policy accessible"
else
  echo "❌ FAIL: Privacy policy returns $PRIVACY_STATUS"
  exit 1
fi

# Test 5: Impressum Accessible
echo "Test 5: Checking impressum.html..."
IMPRESSUM_STATUS=$(curl -s -o /dev/null -w "%{http_code}" "$PROD_URL/impressum.html")
if [ "$IMPRESSUM_STATUS" -eq 200 ]; then
  echo "✅ PASS: Impressum accessible"
else
  echo "❌ FAIL: Impressum returns $IMPRESSUM_STATUS"
  exit 1
fi

echo ""
echo "✅ All DSGVO smoke tests passed!"
```

**GitHub Actions Integration**:
```yaml
# .github/workflows/dsgvo-smoke-tests.yml
name: DSGVO Smoke Tests

on:
  deployment_status:

jobs:
  smoke-tests:
    if: github.event.deployment_status.state == 'success'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Run DSGVO smoke tests
        run: bash scripts/dsgvo-smoke-test.sh
        env:
          GCP_SA_KEY: ${{ secrets.GCP_SA_KEY }}

      - name: Notify on failure
        if: failure()
        run: |
          echo "🚨 DSGVO smoke tests failed!"
          exit 1
```

**Regel**:
- **Run after EVERY deployment** (automatisch via GitHub Actions)
- **Block releases** bei Smoke-Test-Failures
- **Manual verification** als Fallback (wenn CI/CD fehlt)
- **Document results** in deployment notes

---

## Bei Unsicherheit: DSGVO-Scan

```bash
# Nutze den gdpr-compliance-scanner
/gdpr-compliance-scanner:scan-gdpr
```

---

## Lesson Learned Origins

**Deployment Verification**: Lesson 2025-11-18 (Region drift prevention)
**CSP Compliance**: Incident 2025-11-18 (Tailwind CDN blocked)
**AI Model Configuration**: Feature 2025-11-12 (Gemini EU region setup)
**Audit Logging**: Best practice from DSGVO Art. 30 (Verarbeitungsverzeichnis)
**Cloud Run Migration**: Q1 2026 (Vercel -> Google Cloud Run)

---

## Cloud Run Specific Commands

### Deployment
```bash
# Deploy to Frankfurt region
gcloud run deploy fabrikiq-backend   --source backend/   --region europe-west3   --platform managed   --allow-unauthenticated

# Deploy with environment variables from Secret Manager
gcloud run deploy fabrikiq-backend   --source backend/   --region europe-west3   --set-secrets=GEMINI_API_KEY=gemini-api-key:latest,JWT_SECRET=jwt-secret:latest
```

### Verification
```bash
# Check service region
gcloud run services describe fabrikiq-backend --region europe-west3

# View logs
gcloud run logs read fabrikiq-backend --region europe-west3 --limit 50

# Check traffic routing
gcloud run services describe fabrikiq-backend   --region europe-west3   --format="table(status.traffic)"
```

### Rollback
```bash
# List revisions
gcloud run revisions list --service fabrikiq-backend --region europe-west3

# Rollback to previous revision
gcloud run services update-traffic fabrikiq-backend   --region europe-west3   --to-revisions REVISION_NAME=100
```




