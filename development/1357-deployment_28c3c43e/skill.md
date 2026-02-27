# Deployment Reference

Production deployment checklist and common mistakes for TanStack Start on Cloudflare Workers.

## Pre-Deploy Checklist

- [ ] `wrangler.jsonc` has correct `account_id`
- [ ] D1 database created and `database_id` in wrangler.jsonc
- [ ] `main` field is `"@tanstack/react-start/server-entry"`
- [ ] `nodejs_compat` in `compatibility_flags`
- [ ] `.dev.vars` is in `.gitignore`
- [ ] No hardcoded secrets in source code

## Deploy Sequence

### 1. Set Production Secrets

```bash
# Generate a DIFFERENT secret than dev
openssl rand -hex 32 | npx wrangler secret put BETTER_AUTH_SECRET

# Use your actual Worker URL (update after first deploy if needed)
echo "https://PROJECT_NAME.SUBDOMAIN.workers.dev" | npx wrangler secret put BETTER_AUTH_URL

# Include all valid origins
echo "http://localhost:3000,https://PROJECT_NAME.SUBDOMAIN.workers.dev" | npx wrangler secret put TRUSTED_ORIGINS
```

### 2. Google OAuth (if using)

```bash
echo "your-client-id" | npx wrangler secret put GOOGLE_CLIENT_ID
echo "your-client-secret" | npx wrangler secret put GOOGLE_CLIENT_SECRET
```

Add production redirect URI in Google Cloud Console:
```
https://PROJECT_NAME.SUBDOMAIN.workers.dev/api/auth/callback/google
```

### 3. Migrate Remote Database

```bash
pnpm db:migrate:remote
```

**Always run migrations on BOTH local AND remote** before testing.

### 4. Build and Deploy

```bash
pnpm build && npx wrangler deploy
```

### 5. Post-Deploy Verification

- [ ] App loads at production URL
- [ ] Auth login/register works
- [ ] Database operations work (create, read, update, delete)
- [ ] Theme persists across page loads
- [ ] No console errors

## Custom Domain (Optional)

After initial deploy:

1. Add custom domain in Cloudflare Dashboard → Workers → your worker → Triggers → Custom Domains
2. Update secrets:
   ```bash
   echo "https://yourdomain.com" | npx wrangler secret put BETTER_AUTH_URL
   echo "http://localhost:3000,https://yourdomain.com" | npx wrangler secret put TRUSTED_ORIGINS
   ```
3. Update Google OAuth redirect URI to use custom domain
4. Redeploy: `npx wrangler deploy`

## Common Mistakes

### 1. TRUSTED_ORIGINS Not Set

**Symptom**: User signs in but gets silently redirected to homepage. Auth cookie not accepted.

**Fix**: Set the TRUSTED_ORIGINS secret with ALL valid URLs (comma-separated, no spaces):
```bash
echo "http://localhost:3000,https://your-app.workers.dev" | npx wrangler secret put TRUSTED_ORIGINS
```

### 2. BETTER_AUTH_URL Mismatch

**Symptom**: Auth works locally but fails in production. Redirect loops or "invalid callback" errors.

**Fix**: The URL must match the actual production URL exactly — including `https://`, no trailing slash. After first deploy, check the actual Worker URL and update:
```bash
echo "https://actual-url.workers.dev" | npx wrangler secret put BETTER_AUTH_URL
npx wrangler deploy  # Must redeploy after setting secret
```

### 3. Secrets Set But Not Deployed

**Symptom**: Changed a secret but the app still uses the old value.

**Fix**: `wrangler secret put` does NOT automatically redeploy. Always run `npx wrangler deploy` after setting secrets.

### 4. Database Migrations Not Applied to Remote

**Symptom**: API returns 500 errors. "Table not found" or "no such table" in logs.

**Fix**: Run migrations on remote:
```bash
pnpm db:migrate:remote
```

### 5. Google OAuth Redirect URI Missing

**Symptom**: "Error 400: redirect_uri_mismatch" when signing in with Google.

**Fix**: Add your production URL to Google Cloud Console → APIs & Services → Credentials → your OAuth client → Authorised redirect URIs:
```
https://your-app.workers.dev/api/auth/callback/google
```

### 6. Wrong Vite Plugin Order

**Symptom**: Build fails or SSR doesn't work. Cryptic Vite errors about "environment not found".

**Fix**: Plugin order in `vite.config.ts` must be: `cloudflare()` → `tailwindcss()` → `tanstackStart()` → `viteReact()`. Cloudflare MUST be first.

### 7. Module-Level Env Access

**Symptom**: `env` is undefined, or "Cannot access env outside of request context".

**Fix**: Never import and use `env` from `cloudflare:workers` at module level. Always access inside server function handlers, route loaders, or API route handlers.

```typescript
// ❌ Module level — fails
import { env } from "cloudflare:workers";
const db = drizzle(env.DB);

// ✅ Inside handler — works
export const getItems = createServerFn({ method: "GET" }).handler(async () => {
  const db = drizzle(env.DB); // env accessed during request
});
```

### 8. Auth on Server Function Instead of API Route

**Symptom**: better-auth endpoints return 404 or auth flow doesn't work.

**Fix**: better-auth needs a catch-all API route (`createAPIFileRoute`), not a server function. It handles its own routing internally:

```typescript
// src/routes/api/auth/$.ts
import { createAPIFileRoute } from "@tanstack/react-start/api";

export const APIRoute = createAPIFileRoute("/api/auth/$")({
  GET: ({ request }) => getAuth().handler(request),
  POST: ({ request }) => getAuth().handler(request),
});
```
