<CODE_GUIDE_ADDON_HACKATHON>
> ⚠️ **READ BASELINE FIRST**
> This add-on modifies `CODE_GUIDE.md` for fast prototypes and demos.
> Rules here override the Baseline where they conflict.

## Interpretation
- **Relax**: Baseline rule exists, but simpler path is acceptable.
- **Skip**: Ignore that baseline guidance for hackathon scope.
- **Override**: Use this guidance instead of baseline.

## Security Floor (never relax)
- Never commit secrets
- Hash passwords if handling real users
- Verify webhook signatures if present
- Don't expose stack traces to clients

---

## Frontend

### Structure
- **Override**: Flat layout (`/components`, `/routes`, `/lib`) is fine; refactor later if the project continues.

### Data Fetching
- **Relax**: Client-side fetching is acceptable; use server fetch only when equally easy.

### Components
- **Relax**: Larger components are fine temporarily; extract only for obvious reuse.

### Styling
- **Override**: Use the fastest path (e.g., Tailwind + simple components); skip tokens and theming.

### Assets
- **Relax**: Default `next/image` settings; skip font optimization.

### Testing
- **Relax**: One E2E happy-path test is enough; manual testing acceptable for demos.

### Error Handling
- **Relax**: `console.error` is fine; add error tracking only if trivial.

---

## Backend

### Architecture
- **Override**: Single service; avoid layering and microservices.

### Config
- **Relax**: Single `.env`; validate only required values; document in README.

### Database
- **Override**: SQLite or hosted dev DB; schema push acceptable; provide a seed script.

### Auth & Integrations
- **Override**: Use hosted providers for auth, payments, file storage rather than custom builds.

### API
- **Relax**: Minimal endpoints; return simple, consistent error shape; stub unimplemented features with 501.

### Observability
- **Relax**: Console logs with request path/status; add `/health`; skip structured logging.

### Background Work
- **Skip**: Avoid external queues/schedulers; use in-process timers only if essential.

### Webhooks
- **Relax**: Shared secret; retries optional; log deliveries.

### Testing
- **Relax**: One integration test for the main flow if time allows.

---

## Graduation Checklist

When moving beyond hackathon mode, address:
- [ ] Design tokens and consistent styling
- [ ] Proper error handling and logging
- [ ] Test coverage for critical paths
- [ ] Production security headers
- [ ] Caching strategy
</CODE_GUIDE_ADDON_HACKATHON>
