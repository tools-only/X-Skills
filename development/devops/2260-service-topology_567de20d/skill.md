# Service Topology

This file documents the services in your application stack for the `/appfix` workflow.

**IMPORTANT**: Fill this out before using `/appfix`. The workflow needs these URLs to check health and verify fixes.

## Services

### Frontend

| Property | Value |
|----------|-------|
| Name | `[YOUR_FRONTEND_NAME]` |
| Staging URL | `https://staging.example.com` |
| Production URL | `https://example.com` |
| Health Endpoint | `/api/health` |
| Technology | Next.js / React / Vue / etc. |

### Backend API

| Property | Value |
|----------|-------|
| Name | `[YOUR_API_NAME]` |
| Staging URL | `https://api-staging.example.com` |
| Production URL | `https://api.example.com` |
| Health Endpoint | `/health` |
| Technology | FastAPI / Express / etc. |

### Workers (Optional)

| Property | Value |
|----------|-------|
| Name | `[YOUR_WORKER_NAME]` |
| Staging URL | N/A (background worker) |
| Health Endpoint | N/A |
| Technology | Celery / Bull / etc. |

## Deployment

| Environment | Trigger Command | Watch Command |
|-------------|-----------------|---------------|
| Staging | `gh workflow run deploy.yml -f environment=staging` | `gh run watch --exit-status` |
| Production | `gh workflow run deploy.yml -f environment=production` | `gh run watch --exit-status` |

## Health Check Commands

```bash
# Frontend health
curl -sf https://staging.example.com/api/health

# Backend health
curl -sf https://api-staging.example.com/health

# All services quick check
for url in "https://staging.example.com/api/health" "https://api-staging.example.com/health"; do
  curl -sf "$url" && echo " OK: $url" || echo " FAIL: $url"
done
```

## Log Sources

| Source | Command |
|--------|---------|
| Azure Container Logs | `az containerapp logs show --name [app] --resource-group [rg] --type console --tail 100` |
| LogFire (if configured) | `curl -H "Authorization: Bearer $LOGFIRE_READ_TOKEN" "https://logfire-api.pydantic.dev/v1/query?level=error&since=1h"` |
| Browser Console | Chrome MCP `read_console_messages` |

## Critical Paths to Verify

List the user flows that must work after any fix:

1. **Login Flow**: Navigate to `/login`, enter credentials, verify redirect to dashboard
2. **Data Display**: Check that tables/lists load actual data (not spinners)
3. **API Calls**: Verify no 500 errors in Network tab

## Infrastructure Repository

**IMPORTANT**: If you make direct infrastructure changes using `az cli` commands, you MUST create a PR to this repo to sync the IaC files with the actual state.

| Property | Value |
|----------|-------|
| Repo URL | `https://github.com/[ORG]/[INFRA-REPO]` |
| Branch | `main` |
| IaC Type | Terraform / Bicep / ARM |
| Path | `environments/staging/` |

### Why This Matters

When you run `az containerapp update` or similar commands:
1. The change is applied immediately to Azure
2. But the IaC files (Terraform/Bicep) still have the OLD values
3. Next time someone runs `terraform apply`, YOUR fix gets overwritten
4. Infrastructure drift causes confusion and outages

### Required Steps After az CLI Changes

1. Document changes in `.claude/infra-changes.md`
2. Clone this infra repo
3. Update the relevant IaC files to match your az CLI changes
4. Create a PR with description of what changed and why

## Notes

Add any project-specific notes here:

-
