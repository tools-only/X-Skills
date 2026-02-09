---
name: langfuse-security-basics
description: |
  Implement Langfuse security best practices for API keys and data privacy.
  Use when securing Langfuse integration, protecting API keys,
  or implementing data privacy controls for LLM observability.
  Trigger with phrases like "langfuse security", "langfuse API key security",
  "langfuse data privacy", "secure langfuse", "langfuse PII".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Security Basics

## Overview
Security best practices for Langfuse API keys, data privacy, and compliance.

## Prerequisites
- Langfuse SDK installed
- Understanding of API key management
- Knowledge of data privacy requirements (GDPR, etc.)

## Instructions

### Step 1: Secure API Key Management

```bash
# NEVER commit API keys to git
echo '.env*' >> .gitignore
echo '!.env.example' >> .gitignore

# Create .env.example with placeholders
cat > .env.example << 'EOF'
# Langfuse Configuration
LANGFUSE_PUBLIC_KEY=pk-lf-your-public-key
LANGFUSE_SECRET_KEY=sk-lf-your-secret-key
LANGFUSE_HOST=https://cloud.langfuse.com
EOF

# Verify no secrets in git history
git log --all --full-history -- '*/.env*' '*.env'
```

### Step 2: Use Environment-Specific Keys

```typescript
// config/langfuse.ts
interface LangfuseConfig {
  publicKey: string;
  secretKey: string;
  host: string;
  enabled: boolean;
}

function getLangfuseConfig(): LangfuseConfig {
  const env = process.env.NODE_ENV || "development";

  // Different keys per environment
  const configs: Record<string, LangfuseConfig> = {
    development: {
      publicKey: process.env.LANGFUSE_PUBLIC_KEY_DEV!,
      secretKey: process.env.LANGFUSE_SECRET_KEY_DEV!,
      host: process.env.LANGFUSE_HOST_DEV || "https://cloud.langfuse.com",
      enabled: true,
    },
    staging: {
      publicKey: process.env.LANGFUSE_PUBLIC_KEY_STAGING!,
      secretKey: process.env.LANGFUSE_SECRET_KEY_STAGING!,
      host: process.env.LANGFUSE_HOST_STAGING || "https://cloud.langfuse.com",
      enabled: true,
    },
    production: {
      publicKey: process.env.LANGFUSE_PUBLIC_KEY_PROD!,
      secretKey: process.env.LANGFUSE_SECRET_KEY_PROD!,
      host: process.env.LANGFUSE_HOST || "https://cloud.langfuse.com",
      enabled: true,
    },
  };

  return configs[env] || configs.development;
}
```

### Step 3: Implement PII Scrubbing

```typescript
interface ScrubConfig {
  patterns: RegExp[];
  replacement: string;
}

const DEFAULT_SCRUB_CONFIG: ScrubConfig = {
  patterns: [
    // Email addresses
    /\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b/g,
    // Phone numbers
    /\b\d{3}[-.]?\d{3}[-.]?\d{4}\b/g,
    // SSN
    /\b\d{3}-\d{2}-\d{4}\b/g,
    // Credit card numbers
    /\b\d{4}[-\s]?\d{4}[-\s]?\d{4}[-\s]?\d{4}\b/g,
    // API keys
    /\b(sk|pk|api|key|token|secret)[-_]?[a-zA-Z0-9]{20,}\b/gi,
    // IP addresses
    /\b\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3}\b/g,
  ],
  replacement: "[REDACTED]",
};

function scrubPII(data: any, config = DEFAULT_SCRUB_CONFIG): any {
  if (typeof data === "string") {
    let scrubbed = data;
    for (const pattern of config.patterns) {
      scrubbed = scrubbed.replace(pattern, config.replacement);
    }
    return scrubbed;
  }

  if (Array.isArray(data)) {
    return data.map((item) => scrubPII(item, config));
  }

  if (typeof data === "object" && data !== null) {
    const scrubbed: Record<string, any> = {};
    for (const [key, value] of Object.entries(data)) {
      // Scrub sensitive field names
      if (/password|secret|token|key|auth/i.test(key)) {
        scrubbed[key] = "[REDACTED]";
      } else {
        scrubbed[key] = scrubPII(value, config);
      }
    }
    return scrubbed;
  }

  return data;
}

// Usage
const trace = langfuse.trace({
  name: "user-query",
  input: scrubPII(userInput),
  metadata: scrubPII(metadata),
});
```

### Step 4: Create Secure Trace Wrapper

```typescript
interface SecureTraceOptions {
  scrubInputs?: boolean;
  scrubOutputs?: boolean;
  excludeFields?: string[];
  allowedMetadataFields?: string[];
}

function createSecureLangfuse(
  config: ConstructorParameters<typeof Langfuse>[0],
  securityOptions: SecureTraceOptions = {}
) {
  const langfuse = new Langfuse(config);

  const {
    scrubInputs = true,
    scrubOutputs = true,
    excludeFields = ["password", "token", "apiKey", "secret"],
    allowedMetadataFields = ["userId", "sessionId", "environment"],
  } = securityOptions;

  function sanitize(data: any, shouldScrub: boolean): any {
    if (!shouldScrub) return data;

    // Remove excluded fields
    if (typeof data === "object" && data !== null) {
      const sanitized = { ...data };
      for (const field of excludeFields) {
        delete sanitized[field];
      }
      return scrubPII(sanitized);
    }

    return scrubPII(data);
  }

  function filterMetadata(
    metadata: Record<string, any>
  ): Record<string, any> {
    const filtered: Record<string, any> = {};
    for (const field of allowedMetadataFields) {
      if (metadata[field] !== undefined) {
        filtered[field] = metadata[field];
      }
    }
    return filtered;
  }

  return {
    trace(params: Parameters<typeof langfuse.trace>[0]) {
      return langfuse.trace({
        ...params,
        input: sanitize(params.input, scrubInputs),
        output: params.output ? sanitize(params.output, scrubOutputs) : undefined,
        metadata: params.metadata
          ? filterMetadata(params.metadata)
          : undefined,
      });
    },

    flush: () => langfuse.flushAsync(),
    shutdown: () => langfuse.shutdownAsync(),
  };
}
```

### Step 5: Implement Key Rotation

```typescript
// scripts/rotate-langfuse-keys.ts
import { Langfuse } from "langfuse";

async function rotateKeys() {
  // 1. Generate new keys in Langfuse dashboard
  console.log("1. Create new API keys in Langfuse dashboard");
  console.log("   https://cloud.langfuse.com/settings/api-keys\n");

  // 2. Update secrets in your secret manager
  console.log("2. Update secrets in your environment:");
  console.log("   - AWS Secrets Manager");
  console.log("   - GCP Secret Manager");
  console.log("   - HashiCorp Vault");
  console.log("   - GitHub Secrets\n");

  // 3. Deploy with new keys
  console.log("3. Deploy application with new keys\n");

  // 4. Verify new keys work
  const langfuse = new Langfuse({
    publicKey: process.env.LANGFUSE_PUBLIC_KEY_NEW!,
    secretKey: process.env.LANGFUSE_SECRET_KEY_NEW!,
  });

  const trace = langfuse.trace({
    name: "key-rotation-test",
    metadata: { test: true },
  });

  await langfuse.flushAsync();
  console.log("4. Verified new keys work\n");

  // 5. Revoke old keys
  console.log("5. Revoke old keys in Langfuse dashboard");
  console.log("   https://cloud.langfuse.com/settings/api-keys\n");

  console.log("Key rotation complete!");
}

rotateKeys().catch(console.error);
```

## Output
- Secure API key management
- Environment-specific key separation
- PII scrubbing for inputs/outputs
- Secure trace wrapper
- Key rotation procedure

## Security Checklist

| Item | Status | Notes |
|------|--------|-------|
| API keys in .gitignore | Required | Never commit keys |
| Environment separation | Required | Dev/staging/prod keys |
| PII scrubbing enabled | Recommended | Before sending traces |
| Secrets in secret manager | Required | Not in code/config |
| Key rotation schedule | Recommended | Quarterly minimum |
| Audit logging enabled | Recommended | Track API key usage |

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Key exposed in logs | Logging full requests | Scrub before logging |
| PII in traces | No scrubbing | Enable PII scrubbing |
| Key in git history | Committed .env | Rotate keys immediately |
| Cross-env data leak | Shared keys | Use separate keys per env |

## Examples

### Secrets Manager Integration
```typescript
import { SecretsManager } from "@aws-sdk/client-secrets-manager";

const secretsManager = new SecretsManager({ region: "us-east-1" });

async function getLangfuseSecrets() {
  const response = await secretsManager.getSecretValue({
    SecretId: "langfuse/production",
  });

  const secrets = JSON.parse(response.SecretString!);

  return {
    publicKey: secrets.LANGFUSE_PUBLIC_KEY,
    secretKey: secrets.LANGFUSE_SECRET_KEY,
  };
}

// Initialize Langfuse with secrets from manager
const secrets = await getLangfuseSecrets();
const langfuse = new Langfuse(secrets);
```

### Audit Key Usage
```typescript
// Track which keys are used where
const trace = langfuse.trace({
  name: "operation",
  metadata: {
    keyPrefix: process.env.LANGFUSE_PUBLIC_KEY?.slice(0, 10),
    environment: process.env.NODE_ENV,
    service: process.env.SERVICE_NAME,
  },
});
```

## Resources
- [Langfuse Security](https://langfuse.com/docs/security)
- [GDPR Compliance](https://langfuse.com/docs/security/gdpr)
- [Langfuse Self-Hosting](https://langfuse.com/docs/deployment/self-host)

## Next Steps
For production readiness, see `langfuse-prod-checklist`.
