# Security & Input Validation Skill

Dieser Skill hilft bei der Absicherung deiner Anwendungen gegen häufige Sicherheitsl

ücken: XSS, SQL Injection, CSRF, und mehr.

## Wann wird dieser Skill aktiviert?

- Du erwähnst: security, validation, Zod, sanitize, XSS, CSRF, injection
- Du arbeitest mit User-Input, Forms, API-Endpoints
- Du fragst nach Secret-Management oder Rate-Limiting

## Goldene Sicherheitsregeln

### 1. NIEMALS Secrets im Code
```ts
// ❌ NIEMALS!
const API_KEY = 'sk-1234567890';

// ✅ Immer via Environment Variables
const API_KEY = import.meta.env.VITE_API_KEY;

// ✅ Backend: process.env
const SECRET = process.env.GEMINI_API_KEY;
```

### 2. IMMER Input validieren
```ts
// ❌ Direkt nutzen
app.post('/user', (req, res) => {
  db.insert(req.body); // GEFAHR!
});

// ✅ Mit Zod validieren
app.post('/user', async (req, res) => {
  const result = userSchema.safeParse(req.body);
  if (!result.success) {
    return res.status(400).json({ error: result.error });
  }
  db.insert(result.data);
});
```

### 3. IMMER Output escapen
```tsx
// ❌ XSS-Gefahr
<div dangerouslySetInnerHTML={{ __html: userInput }} />

// ✅ React escaped automatisch
<div>{userInput}</div>

// ✅ Für HTML: DOMPurify nutzen
import DOMPurify from 'isomorphic-dompurify';
<div dangerouslySetInnerHTML={{ __html: DOMPurify.sanitize(html) }} />
```

---

## 1. Input Validation mit Zod

### Setup

```bash
npm install zod
```

### Basis-Patterns

```ts
import { z } from 'zod';

// Email & Password
const loginSchema = z.object({
  email: z.string().email('Ungültige E-Mail'),
  password: z.string().min(8, 'Mindestens 8 Zeichen')
});

// Optional Fields
const userSchema = z.object({
  name: z.string(),
  age: z.number().int().positive().optional()
});

// Enums
const roleSchema = z.enum(['admin', 'user', 'guest']);

// Custom Validation
const passwordSchema = z.string().refine(
  (val) => /^(?=.*[A-Z])(?=.*[0-9])/.test(val),
  { message: 'Password muss Großbuchstaben und Zahlen enthalten' }
);
```

### API Endpoint Protection

```ts
// Vercel Serverless Function
import { z } from 'zod';

const analyzeRequestSchema = z.object({
  file: z.string(),
  options: z.object({
    anonymize: z.boolean().default(true),
    language: z.enum(['de', 'en']).default('de')
  })
});

export default async function handler(req, res) {
  // Parse & Validate
  const result = analyzeRequestSchema.safeParse(req.body);

  if (!result.success) {
    return res.status(400).json({
      error: 'Invalid request',
      details: result.error.flatten()
    });
  }

  const { file, options } = result.data;
  // Sicher: data ist typisiert & validiert!
}
```

### Form Validation (React)

```tsx
import { z } from 'zod';
import { useState } from 'react';

const contactSchema = z.object({
  name: z.string().min(2, 'Name zu kurz'),
  email: z.string().email('Ungültige E-Mail'),
  message: z.string().min(10, 'Nachricht zu kurz')
});

export function ContactForm() {
  const [errors, setErrors] = useState<z.ZodFormattedError<any>>();

  const handleSubmit = (e: React.FormEvent<HTMLFormElement>) => {
    e.preventDefault();
    const formData = new FormData(e.currentTarget);
    const data = Object.fromEntries(formData);

    const result = contactSchema.safeParse(data);

    if (!result.success) {
      setErrors(result.error.format());
      return;
    }

    // Sicher: Submit validated data
    submitForm(result.data);
  };

  return (
    <form onSubmit={handleSubmit}>
      <input name="name" />
      {errors?.name && <span>{errors.name._errors[0]}</span>}

      <input name="email" type="email" />
      {errors?.email && <span>{errors.email._errors[0]}</span>}

      <textarea name="message" />
      {errors?.message && <span>{errors.message._errors[0]}</span>}

      <button>Submit</button>
    </form>
  );
}
```

---

## 2. XSS Prevention

### Automatisches Escaping (React)

React escaped automatisch - **außer** bei `dangerouslySetInnerHTML`:

```tsx
// ✅ Sicher - React escaped
<div>{userInput}</div>

// ⚠️ GEFAHR!
<div dangerouslySetInnerHTML={{ __html: userComment }} />
```

### HTML Sanitization

```bash
npm install isomorphic-dompurify
```

```ts
import DOMPurify from 'isomorphic-dompurify';

// Server-Side (Node.js) oder Client-Side
const sanitized = DOMPurify.sanitize(userHtml, {
  ALLOWED_TAGS: ['b', 'i', 'em', 'strong', 'a'],
  ALLOWED_ATTR: ['href']
});

// React
<div dangerouslySetInnerHTML={{ __html: sanitized }} />
```

### URL Validation

```ts
const urlSchema = z.string().url().refine((url) => {
  // Nur HTTPS erlauben
  return url.startsWith('https://');
}, { message: 'Nur HTTPS-URLs erlaubt' });

// Oder: Whitelist
const safeUrl = (url: string) => {
  const allowed = ['https://example.com', 'https://api.example.com'];
  return allowed.some(domain => url.startsWith(domain));
};
```

---

## 3. SQL Injection Prevention

### Drizzle ORM (Prepared Statements)

```ts
import { drizzle } from 'drizzle-orm/vercel-postgres';
import { eq } from 'drizzle-orm';

// ✅ Sicher - Parameterized Query
const user = await db.select()
  .from(users)
  .where(eq(users.email, userInput)); // Automatisch escaped!

// ❌ NIEMALS!
const query = `SELECT * FROM users WHERE email = '${userInput}'`;
```

### Validation vor DB-Insert

```ts
const userInsertSchema = z.object({
  email: z.string().email(),
  name: z.string().max(255),
  age: z.number().int().min(0).max(150)
});

const result = userInsertSchema.safeParse(req.body);
if (!result.success) throw new Error('Invalid data');

await db.insert(users).values(result.data);
```

---

## 4. Rate Limiting

### Vercel KV (Redis) Pattern

> ⚠️ **DSGVO**: Niemals Klartext-Email als Key! Immer hashen.

```ts
import { kv } from '@vercel/kv';
import crypto from 'crypto';

// Hash email for privacy (DSGVO-konform)
function hashEmail(email: string): string {
  return crypto
    .createHash('sha256')
    .update(email.toLowerCase() + process.env.HASH_SALT)
    .digest('hex')
    .slice(0, 32); // 32 chars = 128 bit entropy
}

async function checkRateLimit(email: string): Promise<boolean> {
  const key = `rate:${hashEmail(email)}`; // ✅ Hashed, nicht Klartext!
  const limit = 5; // Max 5 requests
  const window = 3600; // Per hour

  const current = await kv.get<number>(key) || 0;

  if (current >= limit) {
    return false; // Rate limit exceeded
  }

  await kv.set(key, current + 1, { ex: window });
  return true;
}

// API Endpoint
export default async function handler(req, res) {
  const email = req.body.email;

  if (!await checkRateLimit(email)) {
    return res.status(429).json({ error: 'Rate limit exceeded' });
  }

  // Process request
}
```

### IP-Based Rate Limiting

```ts
function getClientIP(req: Request): string {
  return req.headers.get('x-forwarded-for')?.split(',')[0] ||
         req.headers.get('x-real-ip') ||
         'unknown';
}

async function rateLimitByIP(req: Request): Promise<boolean> {
  const ip = getClientIP(req);
  const key = `ip-rate:${ip}`;

  const requests = await kv.incr(key);
  if (requests === 1) {
    await kv.expire(key, 60); // 60 seconds window
  }

  return requests <= 100; // Max 100 requests/minute
}
```

---

## 5. CSRF Protection

### Token-Based (Next.js/Vercel)

```ts
import { randomBytes } from 'crypto';

// Generate CSRF Token
export function generateCSRFToken(): string {
  return randomBytes(32).toString('hex');
}

// Verify Token
export function verifyCSRFToken(token: string, expected: string): boolean {
  return token === expected;
}

// API Route
export default async function handler(req, res) {
  const csrfToken = req.headers['x-csrf-token'];
  const sessionToken = req.cookies.csrfToken;

  if (!verifyCSRFToken(csrfToken, sessionToken)) {
    return res.status(403).json({ error: 'Invalid CSRF token' });
  }

  // Process request
}
```

### SameSite Cookies

```ts
res.setHeader('Set-Cookie', [
  `session=${sessionId}; HttpOnly; Secure; SameSite=Strict; Path=/`
]);
```

---

## 6. Environment Variables & Secrets

### Struktur

```
.env.local          # Lokal, gitignored
.env.production     # Production (in Vercel/CI)
```

### Nutzung

```ts
// Vite Frontend
const GEMINI_KEY = import.meta.env.VITE_GEMINI_API_KEY;

// Node.js Backend
const SECRET = process.env.GEMINI_API_KEY;

// Validation beim Start
if (!process.env.DATABASE_URL) {
  throw new Error('DATABASE_URL not set');
}
```

### Vercel Deployment

```bash
# CLI
vercel env add GEMINI_API_KEY

# Dashboard: Settings → Environment Variables
```

**NIEMALS committen:**
```
# .gitignore
.env
.env.local
.env*.local
```

---

## 7. Headers & Security

### Security Headers (Vercel)

```json
// vercel.json
{
  "headers": [
    {
      "source": "/(.*)",
      "headers": [
        {
          "key": "X-Content-Type-Options",
          "value": "nosniff"
        },
        {
          "key": "X-Frame-Options",
          "value": "DENY"
        },
        {
          "key": "X-XSS-Protection",
          "value": "1; mode=block"
        },
        {
          "key": "Strict-Transport-Security",
          "value": "max-age=31536000; includeSubDomains"
        },
        {
          "key": "Content-Security-Policy",
          "value": "default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline';"
        }
      ]
    }
  ]
}
```

---

## 8. Projekt-spezifische Patterns

### ManufacturingInsideAnalyzer (DSGVO + PII)

```ts
// Kein PII-Logging!
// ❌ NIEMALS
console.log('User email:', user.email);

// ✅ Nur Metadata
console.log('Analysis completed', {
  duration: elapsedMs,
  status: 'success',
  userId: hashUserId(user.id) // Gehasht!
});

// Rate Limit pro Email
const canAnalyze = await checkRateLimit(email);
if (!canAnalyze) {
  throw new Error('Rate limit: 5 analyses per 24h');
}
```

### digitalTwin (Camera Access)

```ts
// Validate before accessing hardware
const cameraPermissionSchema = z.object({
  requestCamera: z.boolean(),
  resolution: z.enum(['720p', '1080p']).default('720p')
});

async function requestCameraAccess(options: unknown) {
  const validated = cameraPermissionSchema.parse(options);

  if (!validated.requestCamera) return null;

  try {
    const stream = await navigator.mediaDevices.getUserMedia({
      video: { width: validated.resolution === '1080p' ? 1920 : 1280 }
    });
    return stream;
  } catch (err) {
    // Fallback: Mock camera
    return getMockStream();
  }
}
```

---

## 9. Audit & Monitoring

### Security Checklist

- [ ] Alle Secrets in `.env` (nicht im Code)
- [ ] Input-Validation mit Zod
- [ ] Output-Sanitization (DOMPurify)
- [ ] Rate-Limiting aktiv
- [ ] CSRF-Protection (wenn Forms)
- [ ] Security Headers gesetzt
- [ ] HTTPS erzwungen
- [ ] Kein PII-Logging

### Tools

```bash
# npm Audit
npm audit

# Dependency Check
npm install -D @cyclonedx/cdxgen
npx cdxgen -o bom.json

# Snyk (Security Scanner)
npx snyk test
```

---

## Ressourcen

- [OWASP Top 10](https://owasp.org/www-project-top-ten/)
- [Zod Docs](https://zod.dev/)
- [DOMPurify](https://github.com/cure53/DOMPurify)
- [DSGVO Compliance Guide](../dsgvo-compliance/main.md)
