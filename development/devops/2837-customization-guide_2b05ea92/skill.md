# Customisation Guide

Complete rebranding and customisation reference for vite-flare-starter forks.

## Security Fingerprints

If you don't change these, attackers can identify your site uses this starter:

| Location | Default Value | How to Change |
|----------|---------------|---------------|
| Page title | "Vite Flare Starter" | `index.html` |
| App name in UI | "Vite Flare Starter" | `VITE_APP_NAME` env var |
| localStorage keys | `vite-flare-starter-theme` | `VITE_APP_ID` env var |
| API tokens | `vfs_` prefix | `VITE_TOKEN_PREFIX` env var |
| Sentry release | `vite-flare-starter@x.x.x` | `VITE_APP_ID` env var |
| GitHub links | jezweb repo | `VITE_GITHUB_URL` (set empty to hide) |
| Worker name | `vite-flare-starter` | `wrangler.jsonc` |
| Database name | `vite-flare-starter-db` | `wrangler.jsonc` |
| R2 buckets | `vite-flare-starter-*` | `wrangler.jsonc` |

## Environment Variables

### Branding Variables (VITE_ prefix = available in frontend)

| Variable | Purpose | Example |
|----------|---------|---------|
| `VITE_APP_NAME` | Display name in UI | "My Cool App" |
| `VITE_APP_ID` | localStorage prefix, Sentry | "mycoolapp" |
| `VITE_TOKEN_PREFIX` | API token prefix | "mca_" |
| `VITE_GITHUB_URL` | GitHub link (empty = hidden) | "" |
| `VITE_FOOTER_TEXT` | Footer copyright text | "2026 My Company" |
| `APP_NAME` | Server-side app name | "My Cool App" |

### Auth Variables

| Variable | Purpose | Notes |
|----------|---------|-------|
| `BETTER_AUTH_SECRET` | Session encryption | Generate with `openssl rand -hex 32` |
| `BETTER_AUTH_URL` | Auth base URL | Must match actual URL exactly |
| `TRUSTED_ORIGINS` | Allowed origins | Comma-separated, include localhost + prod |
| `GOOGLE_CLIENT_ID` | Google OAuth | From Google Cloud Console |
| `GOOGLE_CLIENT_SECRET` | Google OAuth | From Google Cloud Console |
| `ENABLE_EMAIL_LOGIN` | Enable email/password | "true" to enable |
| `ENABLE_EMAIL_SIGNUP` | Enable email signup | "true" to enable (requires ENABLE_EMAIL_LOGIN) |

### Email Variables (Optional)

| Variable | Purpose | Notes |
|----------|---------|-------|
| `EMAIL_FROM` | Sender address | For verification/password reset |
| `EMAIL_API_KEY` | Email service API key | Resend recommended |

## Common Customisations

### Adding a New Database Table

1. Add schema in `src/server/db/schema.ts`
2. Generate migration: `pnpm db:generate`
3. Apply locally: `pnpm db:migrate:local`
4. Apply to production: `pnpm db:migrate:remote`

### Adding a New API Route

1. Create route file in `src/server/routes/`
2. Register in `src/server/index.ts`
3. Add TanStack Query hook in `src/client/hooks/`

### Changing Auth Providers

better-auth supports multiple providers. Edit `src/server/auth.ts`:
- Add provider to `socialProviders`
- Add credentials to `.dev.vars` and production secrets
- Update client-side login buttons

### Custom Theme

Edit Tailwind theme in `src/client/index.css`:
- Update CSS variables in `:root` and `.dark`
- Use `tailwind-theme-builder` skill for guided setup

### Feature Flags

Control features via environment variables:
- `VITE_FEATURE_STYLE_GUIDE=true` — enable style guide page
- `VITE_FEATURE_COMPONENTS=true` — enable component showcase
- Add your own in `src/client/lib/features.ts`

## Production Deployment Checklist

- [ ] `BETTER_AUTH_SECRET` set (different from dev!)
- [ ] `BETTER_AUTH_URL` matches actual Worker URL
- [ ] `TRUSTED_ORIGINS` includes all valid URLs
- [ ] Google OAuth redirect URI includes production URL
- [ ] Remote database migrated (`pnpm db:migrate:remote`)
- [ ] No `vite-flare-starter` references in config files
- [ ] Favicon replaced
- [ ] CLAUDE.md updated
- [ ] `.dev.vars` is NOT committed (check `.gitignore`)

## Common Mistakes

### 1. TRUSTED_ORIGINS Not Set

**Symptom**: User signs in but gets redirected to homepage (auth cookie not accepted).

**Fix**:
```bash
echo "http://localhost:5173,https://your-app.workers.dev" | npx wrangler secret put TRUSTED_ORIGINS
```

### 2. BETTER_AUTH_URL Mismatch

**Symptom**: Auth fails in production but works locally.

**Fix**: The URL must match exactly — including protocol (`https://`) and no trailing slash.

### 3. Database Migrations Not Applied

**Symptom**: API returns 500 errors on database operations.

**Fix**: Run both local and remote migrations:
```bash
pnpm db:migrate:local
pnpm db:migrate:remote
```

### 4. Secrets Set But Not Deployed

**Symptom**: Changes to secrets don't take effect.

**Fix**: Setting a secret via `wrangler secret put` does NOT automatically redeploy. Run `pnpm deploy` after.

### 5. Google OAuth Redirect URI Missing

**Symptom**: "Error 400: redirect_uri_mismatch" when signing in with Google.

**Fix**: Add your production URL to Google Cloud Console > OAuth Client > Authorised redirect URIs:
```
https://your-app.workers.dev/api/auth/callback/google
```
