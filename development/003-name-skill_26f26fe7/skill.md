---
name: lokalise-security-basics
description: |
  Apply Lokalise security best practices for API tokens and access control.
  Use when securing API tokens, implementing least privilege access,
  or auditing Lokalise security configuration.
  Trigger with phrases like "lokalise security", "lokalise secrets",
  "secure lokalise", "lokalise API token security".
allowed-tools: Read, Write, Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Lokalise Security Basics

## Overview
Security best practices for Lokalise API tokens, webhook secrets, and access control.

## Prerequisites
- Lokalise SDK installed
- Understanding of environment variables
- Access to Lokalise profile settings

## Instructions

### Step 1: Configure Environment Variables
```bash
# .env (NEVER commit to git)
LOKALISE_API_TOKEN=abc123...
LOKALISE_WEBHOOK_SECRET=whsec_***
LOKALISE_PROJECT_ID=123456789.abcdef

# .gitignore (MUST include)
.env
.env.local
.env.*.local
*.env
```

### Step 2: Token Types and Permissions
| Token Type | Use Case | Permissions |
|------------|----------|-------------|
| Read-only | CI builds, reporting | Read projects, keys, translations |
| Read-write | Full integration | All CRUD operations |
| OAuth2 | User-facing apps | Scoped per user consent |

```typescript
// Use appropriate token type
const readOnlyClient = new LokaliseApi({
  apiKey: process.env.LOKALISE_READ_TOKEN!,
});

const fullAccessClient = new LokaliseApi({
  apiKey: process.env.LOKALISE_WRITE_TOKEN!,
});
```

### Step 3: Implement Token Rotation
```bash
# Token rotation procedure:
# 1. Generate new token in Profile Settings -> API tokens
# 2. Update environment variable in all environments
# 3. Test new token
# 4. Revoke old token

# Verify new token works
curl -X GET "https://api.lokalise.com/api2/projects" \
  -H "X-Api-Token: $LOKALISE_API_TOKEN_NEW" \
  | jq '.projects | length'

# If successful, update and revoke old token
export LOKALISE_API_TOKEN="$LOKALISE_API_TOKEN_NEW"
```

### Step 4: Webhook Secret Verification
```typescript
import crypto from "crypto";

function verifyWebhookSignature(
  payload: string,
  receivedSecret: string,
  expectedSecret: string
): boolean {
  // Lokalise sends secret in X-Secret header
  return crypto.timingSafeEqual(
    Buffer.from(receivedSecret),
    Buffer.from(expectedSecret)
  );
}

// Express middleware
app.post("/webhooks/lokalise", (req, res) => {
  const receivedSecret = req.headers["x-secret"] as string;
  const expectedSecret = process.env.LOKALISE_WEBHOOK_SECRET!;

  if (!verifyWebhookSignature(req.body, receivedSecret, expectedSecret)) {
    console.error("Invalid webhook signature");
    return res.status(401).json({ error: "Invalid signature" });
  }

  // Process webhook
  handleWebhook(req.body);
  res.status(200).json({ received: true });
});
```

### Step 5: Apply Least Privilege
```typescript
// Project-level access control
// Use Team settings to limit user access to specific projects

// Environment-specific tokens
const tokenConfig = {
  development: {
    token: process.env.LOKALISE_DEV_TOKEN,
    projectId: process.env.LOKALISE_DEV_PROJECT,
  },
  staging: {
    token: process.env.LOKALISE_STAGING_TOKEN,
    projectId: process.env.LOKALISE_STAGING_PROJECT,
  },
  production: {
    token: process.env.LOKALISE_PROD_TOKEN,
    projectId: process.env.LOKALISE_PROD_PROJECT,
  },
};

function getLokaliseConfig() {
  const env = process.env.NODE_ENV || "development";
  return tokenConfig[env as keyof typeof tokenConfig];
}
```

## Output
- Secure API token storage
- Environment-specific access controls
- Webhook signature verification
- Token rotation procedure documented

## Error Handling
| Security Issue | Detection | Mitigation |
|----------------|-----------|------------|
| Exposed API token | Git scanning, logs | Rotate immediately |
| Missing .gitignore | Code review | Add .env to ignore |
| Webhook without verification | Security audit | Add signature check |
| Single token for all envs | Config review | Separate per environment |

## Examples

### Service Account Pattern
```typescript
// Use different clients for different operations
const clients = {
  reader: new LokaliseApi({
    apiKey: process.env.LOKALISE_READ_TOKEN!,
  }),
  writer: new LokaliseApi({
    apiKey: process.env.LOKALISE_WRITE_TOKEN!,
  }),
};

// Read operations use read-only token
async function getTranslations(projectId: string) {
  return clients.reader.translations().list({ project_id: projectId });
}

// Write operations use write token
async function updateTranslation(projectId: string, translationId: number, text: string) {
  return clients.writer.translations().update(translationId, {
    project_id: projectId,
    translation: text,
  });
}
```

### Secret Scanning Prevention
```typescript
// Pre-commit hook to prevent token leaks
// .husky/pre-commit

const dangerousPatterns = [
  /LOKALISE_API_TOKEN\s*=\s*["']?[a-zA-Z0-9]{40,}/,
  /X-Api-Token:\s*[a-zA-Z0-9]{40,}/,
];

function scanForSecrets(content: string): boolean {
  return dangerousPatterns.some(pattern => pattern.test(content));
}
```

### Audit Logging
```typescript
interface AuditEntry {
  timestamp: Date;
  action: string;
  userId?: string;
  projectId: string;
  resource: string;
  result: "success" | "failure";
  metadata?: Record<string, any>;
}

async function auditLog(entry: Omit<AuditEntry, "timestamp">): Promise<void> {
  const log: AuditEntry = { ...entry, timestamp: new Date() };

  // Log locally for compliance
  console.log("[AUDIT]", JSON.stringify(log));

  // Optional: Send to logging service
  await sendToLoggingService(log);
}

// Usage
await auditLog({
  action: "lokalise.keys.create",
  projectId: "123456.abc",
  resource: "/keys",
  result: "success",
  metadata: { keyCount: 5 },
});
```

### Security Checklist
```markdown
## Lokalise Security Checklist

### Token Management
- [ ] API tokens stored in environment variables
- [ ] .env files in .gitignore
- [ ] Different tokens for dev/staging/prod
- [ ] Read-only tokens where possible
- [ ] Token rotation schedule (quarterly)

### Webhook Security
- [ ] HTTPS-only webhook endpoints
- [ ] X-Secret header validated
- [ ] Timing-safe comparison for secrets

### Access Control
- [ ] Minimal team member permissions
- [ ] Project-level access restrictions
- [ ] Regular access audits

### Monitoring
- [ ] Audit logging enabled
- [ ] Failed auth attempts monitored
- [ ] Unusual API usage alerts
```

### CI/CD Secret Management
```yaml
# GitHub Actions
jobs:
  deploy:
    steps:
      - name: Download translations
        env:
          LOKALISE_API_TOKEN: ${{ secrets.LOKALISE_API_TOKEN }}
        run: npm run i18n:pull

# GitLab CI
deploy:
  variables:
    LOKALISE_API_TOKEN: $LOKALISE_API_TOKEN  # From CI/CD settings
  script:
    - npm run i18n:pull
```

## Resources
- [Lokalise API Authentication](https://developers.lokalise.com/reference/api-authentication)
- [Lokalise Team Management](https://docs.lokalise.com/en/articles/1400472-team-management)
- [OAuth2 Guide](https://developers.lokalise.com/reference/oauth2-authentication)

## Next Steps
For production deployment, see `lokalise-prod-checklist`.
