# OAuth Cookie Patterns Skill

**🎯 Zweck**: Verhindert OAuth Login Loops durch Cookie-Handling-Fehler.

**⚠️ KRITISCH**: Basiert auf Lessons Learned vom 2025-11-12 (3 kritische Root Causes).

## Wann wird dieser Skill aktiviert?

- Du implementierst OAuth 2.0 (Google, GitHub, etc.)
- Du arbeitest mit Session Cookies (`Set-Cookie`, `credentials`)
- Du debuggst Login Loops oder "Session not found" Fehler
- Du erwähnst Keywords: OAuth, authentication, cookies, redirect_uri

## Die 3 Kritischen Root Causes (Lessons Learned)

### Root Cause #1: Missing `credentials: 'include'` in Frontend 🚨

**Problem**: Browser sendet Cookies NICHT automatisch in `fetch()` Requests.

**Symptom**:
- Session Cookie existiert in DevTools → Application → Cookies ✅
- Session Cookie NICHT in Network → Request Headers ❌
- Endless Login Loop: Login → Redirect → Check Session → Not Found → Login

**SCHLECHT:**
```typescript
// ❌ Cookie exists but NEVER sent to server
async function checkSession() {
  const response = await fetch('/api/auth/session');
  const data = await response.json();
  // data.session is ALWAYS null, even though cookie exists!
}
```

**GUT:**
```typescript
// ✅ Cookie is included in request
async function checkSession() {
  const response = await fetch('/api/auth/session', {
    credentials: 'include'  // 🎯 CRITICAL: Required for same-origin cookies
  });
  const data = await response.json();
  // data.session now correctly populated
}
```

**Regel**: JEDER `fetch()` Call zu auth endpoints braucht `credentials: 'include'`.

**Affected Files (ManufacturingInsideAnalyzer):**
```typescript
// App.tsx - 4 locations fixed
const checkSession = async () => {
  const res = await fetch('/api/auth/session', { credentials: 'include' });
};

const handleLoginSuccess = async () => {
  const res = await fetch('/api/auth/session', { credentials: 'include' });
};

const persistPreferences = async () => {
  await fetch('/api/preferences', { credentials: 'include' });
};

const handleLogout = async () => {
  await fetch('/api/auth/logout', { credentials: 'include' });
};
```

---

### Root Cause #2: `Domain=.vercel.app` Blocked by Browser 🚨

**Problem**: Setting `Domain=.vercel.app` causes browser to REJECT cookie entirely.

**Symptom**:
- Server logs show `Set-Cookie` header sent ✅
- Browser DevTools → Application → Cookies shows NOTHING ❌
- Silent failure, no error message

**SCHLECHT:**
```typescript
// ❌ Browser rejects this cookie (security restriction)
function buildSessionCookie(token: string): string {
  return [
    `session=${token}`,
    'Domain=.vercel.app',  // 🚨 CAUSES REJECTION
    'HttpOnly',
    'Secure',
    'SameSite=Lax',
    'Path=/',
    `Max-Age=${30 * 24 * 60 * 60}`
  ].join('; ');
}
```

**GUT:**
```typescript
// ✅ NO Domain attribute → browser scopes automatically
function buildSessionCookie(token: string): string {
  const attributes = [
    `session=${token}`,
    'HttpOnly',           // Prevent XSS attacks
    'SameSite=Lax',       // Allow OAuth redirects, prevent CSRF
    'Path=/',             // Available on all routes
    `Max-Age=${30 * 24 * 60 * 60}` // 30 days
  ];

  if (process.env.NODE_ENV === 'production') {
    attributes.push('Secure'); // HTTPS only
    // NO Domain attribute!
  }

  return attributes.join('; ');
}
```

**Warum funktioniert das?**
- Ohne Domain → Cookie ist scoped zu `manufacturing-insights.vercel.app`
- Mit `Domain=.vercel.app` → Browser lehnt ab (security policy)
- Vercel Production Aliases funktionieren perfekt ohne explizites Domain

**Regel**: NIEMALS `Domain` setzen auf Vercel. Lass Browser automatisch scopen.

---

### Root Cause #3: Invalid Google OAuth Credentials 🚨

**Problem**: Even 1 character wrong in `GOOGLE_CLIENT_SECRET` → total failure.

**Symptom**:
- Redirect to `/?error=token_exchange_failed`
- Vercel Logs zeigen: `{ "error": "invalid_client", "error_description": "Unauthorized" }`

**Debugging Steps:**
```bash
# 1. Check Vercel Environment Variables
npx vercel env ls

# 2. Verify exact match with Google Cloud Console
# Go to: https://console.cloud.google.com/apis/credentials
# - Click your OAuth 2.0 Client
# - Copy EXACT Client ID and Secret
# - Compare character-by-character with Vercel

# 3. Update Vercel
GOOGLE_CLIENT_ID=123456789-abc.apps.googleusercontent.com
GOOGLE_CLIENT_SECRET=GOCSPX-xyz123abc456  # EXACT match required!
GOOGLE_REDIRECT_URI=https://manufacturing-insights.vercel.app/api/auth/google/callback

# 4. Trigger redeploy
git commit --allow-empty -m "chore: trigger redeploy for env vars"
git push
```

**Regel**: Copy-Paste Credentials, never type manually!

---

## Cookie Attribute Best Practices

### Übersicht

| Attribute | Wert | Zweck | Pflicht? |
|-----------|------|-------|----------|
| `HttpOnly` | (flag) | Verhindert JavaScript-Zugriff (XSS Protection) | ✅ JA |
| `Secure` | (flag) | Nur HTTPS (production) | ✅ JA (prod) |
| `SameSite` | `Lax` | OAuth-Redirects erlauben, CSRF verhindern | ✅ JA |
| `Path` | `/` | Cookie auf allen Routes verfügbar | ✅ JA |
| `Max-Age` | `2592000` | 30 Tage Lebensdauer | ✅ JA |
| `Domain` | (omit) | Browser scoped automatisch | ❌ NEIN |

### SameSite Attribute Explained

**`SameSite=Strict`** ❌ Blockiert OAuth!
```
Browser: "Cookie nur bei same-site requests senden"
Problem: OAuth redirect von accounts.google.com → Blocked!
```

**`SameSite=Lax`** ✅ OAuth funktioniert!
```
Browser: "Cookie bei top-level navigation (GET) senden"
OAuth: accounts.google.com → GET redirect → Cookie sent ✅
CSRF: POST from evil.com → Cookie NOT sent ✅
```

**`SameSite=None`** ⚠️ Nur mit 3rd-party cookies!
```
Requires: SameSite=None; Secure
Use Case: iframes, cross-site embeds
OAuth: Nicht nötig!
```

**Regel**: `SameSite=Lax` für OAuth flows!

---

## Browser DevTools Debugging

### Schritt 1: Check if Cookie is Set

**Chrome DevTools → Application → Cookies → `https://your-domain.com`**

✅ **Expected after login:**
```
Name: session
Value: eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9... (JWT token)
Domain: manufacturing-insights.vercel.app  (NO leading dot!)
Path: /
Secure: ✓  (production only)
HttpOnly: ✓
SameSite: Lax
Max-Age: 2592000  (30 days)
```

❌ **If NO cookie** → Server isn't setting it OR browser is blocking it.

**Debugging:**
```bash
# Check server logs
npx vercel logs --prod --since 5m

# Look for:
[AUTH] Google callback: Creating session for user@example.com
[AUTH] Setting cookie: session=eyJ...
```

If logs show cookie set, but browser doesn't receive it:
→ **Domain attribute problem** (Remove Domain!)

---

### Schritt 2: Check if Cookie is Sent

**Chrome DevTools → Network → Select request → Headers → Request Headers**

✅ **Expected:**
```
Cookie: session=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...
```

❌ **If cookie exists in Application tab but NOT in Network tab:**
→ **Missing `credentials: 'include'`** in fetch!

**Quick Fix:**
```typescript
// Before
fetch('/api/auth/session')

// After
fetch('/api/auth/session', { credentials: 'include' })
```

---

### Schritt 3: Check Server-Side Cookie Parsing

**Server Logs:**
```typescript
// api/auth.ts - Session check endpoint
export default async function handler(req: VercelRequest, res: VercelResponse) {
  console.log('[AUTH SESSION] Received cookies:', req.cookies);

  const sessionToken = req.cookies.session;
  if (!sessionToken) {
    console.log('[AUTH SESSION] No session cookie found');
    return res.json({ session: null });
  }

  // Verify JWT...
}
```

**Expected Logs:**
```
[AUTH SESSION] Received cookies: { session: 'eyJ...' }
[AUTH SESSION] Session valid for user@example.com
```

**If cookies empty:**
```
[AUTH SESSION] Received cookies: {}
```
→ Frontend not sending cookies (`credentials: 'include'` missing)

---

## Complete Working Implementation

### Server Side - `api/auth.ts`

```typescript
import { sign, verify } from 'jsonwebtoken';
import type { VercelRequest, VercelResponse } from '@vercel/node';

const JWT_SECRET = process.env.JWT_SECRET!;

// Build session cookie with correct attributes
function buildSessionCookie(token: string): string {
  const attributes = [
    `session=${token}`,
    'HttpOnly',
    'SameSite=Lax',
    'Path=/',
    `Max-Age=${30 * 24 * 60 * 60}` // 30 days
  ];

  if (process.env.NODE_ENV === 'production') {
    attributes.push('Secure');
    // NO Domain attribute!
  }

  return attributes.join('; ');
}

// Google OAuth Callback
export default async function handler(req: VercelRequest, res: VercelResponse) {
  const { code } = req.query;

  // 1. Exchange code for tokens
  const tokenResponse = await fetch('https://oauth2.googleapis.com/token', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      code,
      client_id: process.env.GOOGLE_CLIENT_ID,
      client_secret: process.env.GOOGLE_CLIENT_SECRET,
      redirect_uri: process.env.GOOGLE_REDIRECT_URI,
      grant_type: 'authorization_code'
    })
  });

  if (!tokenResponse.ok) {
    const error = await tokenResponse.json();
    console.error('[OAUTH] Token exchange failed:', error);
    return res.redirect('/?error=token_exchange_failed');
  }

  const tokens = await tokenResponse.json();

  // 2. Get user info
  const userResponse = await fetch('https://www.googleapis.com/oauth2/v2/userinfo', {
    headers: { Authorization: `Bearer ${tokens.access_token}` }
  });

  const user = await userResponse.json();

  // 3. Create session token (JWT)
  const sessionToken = sign(
    {
      email: user.email,
      name: user.name,
      picture: user.picture,
      exp: Math.floor(Date.now() / 1000) + (30 * 24 * 60 * 60) // 30 days
    },
    JWT_SECRET
  );

  // 4. Set cookie
  res.setHeader('Set-Cookie', buildSessionCookie(sessionToken));

  console.log('[AUTH] Google callback: Creating session for', user.email);

  // 5. Redirect to success
  return res.redirect('/?login=success');
}
```

### Client Side - `App.tsx`

```typescript
import { useState, useEffect } from 'react';

interface Session {
  email: string;
  name: string;
  picture: string;
}

export function App() {
  const [session, setSession] = useState<Session | null>(null);
  const [loading, setLoading] = useState(true);

  // Check session on app load
  useEffect(() => {
    checkSession();
  }, []);

  // Check session with credentials: 'include'
  async function checkSession() {
    try {
      const response = await fetch('/api/auth/session', {
        credentials: 'include'  // 🎯 CRITICAL!
      });

      const data = await response.json();
      setSession(data.session);
    } catch (error) {
      console.error('Session check failed:', error);
    } finally {
      setLoading(false);
    }
  }

  // Handle OAuth redirect
  useEffect(() => {
    const params = new URLSearchParams(window.location.search);
    if (params.get('login') === 'success') {
      // Remove query param
      window.history.replaceState({}, '', '/');
      // Refresh session
      checkSession();
    }
  }, []);

  // Logout handler
  async function handleLogout() {
    await fetch('/api/auth/logout', {
      method: 'POST',
      credentials: 'include'  // 🎯 CRITICAL!
    });

    setSession(null);
  }

  if (loading) return <div>Loading...</div>;

  if (!session) {
    return (
      <div>
        <h1>Sign In</h1>
        <a href="/api/auth?action=google">
          <button>Sign in with Google</button>
        </a>
      </div>
    );
  }

  return (
    <div>
      <h1>Welcome {session.name}</h1>
      <img src={session.picture} alt={session.name} />
      <button onClick={handleLogout}>Logout</button>
    </div>
  );
}
```

---

## Testing Checklist

✅ **Test in Incognito Window** (fresh state, no old cookies)

**Steps:**
1. Open incognito window
2. Open DevTools → Application → Cookies
3. Navigate to your app
4. Click "Sign in with Google"
5. After redirect, check:
   - [ ] `session` cookie exists
   - [ ] Domain: `your-domain.com` (no leading dot)
   - [ ] Secure: ✓ (production)
   - [ ] HttpOnly: ✓
   - [ ] SameSite: `Lax`
6. Open DevTools → Network
7. Refresh page
8. Click on request to `/api/auth/session`
9. Verify **Cookie header is present** in Request Headers
10. User should stay logged in ✓

---

## Common Mistakes & Fixes

### Mistake 1: `SameSite=Strict` for OAuth
```typescript
// ❌ Blocks OAuth redirects from accounts.google.com
attributes.push('SameSite=Strict');

// ✅ Use Lax - allows top-level GET navigation
attributes.push('SameSite=Lax');
```

### Mistake 2: Setting `Domain=.vercel.app`
```typescript
// ❌ Browser rejects the cookie
attributes.push('Domain=.vercel.app');

// ✅ Omit Domain - browser scopes to current host
// (no Domain attribute)
```

### Mistake 3: Forgetting `credentials: 'include'`
```typescript
// ❌ Cookie never sent to server
fetch('/api/auth/session')

// ✅ Add to EVERY fetch that needs cookies
fetch('/api/auth/session', { credentials: 'include' })
```

### Mistake 4: Not testing in Incognito
```typescript
// ❌ Old cookies with wrong attributes stay in browser
// Test in normal browser window

// ✅ Always test OAuth in fresh incognito window
// Cmd+Shift+N (Chrome) / Cmd+Shift+P (Firefox)
```

---

## Debugging Flowchart

```
Session not found?
│
├─ Cookie exists in DevTools → Application → Cookies?
│  │
│  ├─ YES → Cookie sent in Network → Request Headers?
│  │  │
│  │  ├─ YES → Server-side issue (check JWT validation)
│  │  │
│  │  └─ NO → Missing credentials: 'include' in fetch
│  │
│  └─ NO → Cookie not set by server
│     │
│     ├─ Check server logs: Set-Cookie header sent?
│     │  │
│     │  ├─ YES → Browser blocking (remove Domain attribute)
│     │  │
│     │  └─ NO → OAuth callback not creating session
│     │
│     └─ Verify Google credentials match exactly
```

---

## Ressourcen

- [MDN: Using Fetch - Credentials](https://developer.mozilla.org/en-US/docs/Web/API/Fetch_API/Using_Fetch#sending_a_request_with_credentials_included)
- [MDN: Set-Cookie Header](https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Set-Cookie)
- [OWASP: SameSite Cookie Attribute](https://owasp.org/www-community/SameSite)
- [OAuth 2.0 RFC 6749](https://datatracker.ietf.org/doc/html/rfc6749)
- [Vercel: Understanding Cookies](https://vercel.com/blog/understanding-cookies)

---

## Lesson Learned Origin

**Date**: 2025-11-12
**Incident**: Google OAuth login loop
**Duration**: ~4 hours debugging
**Root Causes**: 3 separate issues compounded
**Impact**: Login completely broken
**Prevention**: This skill + automated testing

**Files Fixed:**
- `App.tsx`: Added `credentials: 'include'` to 4 fetch calls
- `api/auth.ts`: Removed `Domain=.vercel.app` from cookie attributes
- Vercel Dashboard: Fixed `GOOGLE_CLIENT_SECRET` typo

**Final Commit**: `34a0e7a` - "fix: remove Domain=.vercel.app from cookies"

---

**🎯 Ziel**: Zero OAuth Login Loops!
