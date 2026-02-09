---
name: lokalise-multi-env-setup
description: |
  Configure Lokalise across development, staging, and production environments.
  Use when setting up multi-environment deployments, configuring per-environment secrets,
  or implementing environment-specific Lokalise configurations.
  Trigger with phrases like "lokalise environments", "lokalise staging",
  "lokalise dev prod", "lokalise environment setup", "lokalise config by env".
allowed-tools: Read, Write, Edit, Bash(aws:*), Bash(gcloud:*), Bash(vault:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Lokalise Multi-Environment Setup

## Overview
Configure Lokalise across development, staging, and production environments with proper isolation.

## Prerequisites
- Separate Lokalise projects or branches per environment
- Secret management solution (Vault, AWS Secrets Manager, etc.)
- CI/CD pipeline with environment variables
- Environment detection in application

## Environment Strategy

| Environment | Lokalise Setup | Use Case |
|-------------|---------------|----------|
| Development | Separate project or `dev` branch | Local testing, rapid iteration |
| Staging | Separate project or `staging` branch | Pre-prod validation, QA |
| Production | Main project/branch | Live users |

## Instructions

### Step 1: Create Configuration Structure
```
config/
├── lokalise/
│   ├── base.json           # Shared config
│   ├── development.json    # Dev overrides
│   ├── staging.json        # Staging overrides
│   └── production.json     # Prod overrides
```

### base.json
```json
{
  "enableCompression": true,
  "cacheEnabled": true,
  "cacheTtlSeconds": 300,
  "syncOnBuild": true,
  "fallbackLocale": "en",
  "supportedLocales": ["en", "es", "fr", "de"]
}
```

### development.json
```json
{
  "projectId": "${LOKALISE_PROJECT_ID_DEV}",
  "branch": "development",
  "cacheEnabled": false,
  "syncOnBuild": false,
  "debug": true
}
```

### staging.json
```json
{
  "projectId": "${LOKALISE_PROJECT_ID}",
  "branch": "staging",
  "cacheEnabled": true,
  "cacheTtlSeconds": 60,
  "syncOnBuild": true
}
```

### production.json
```json
{
  "projectId": "${LOKALISE_PROJECT_ID}",
  "branch": null,
  "cacheEnabled": true,
  "cacheTtlSeconds": 3600,
  "syncOnBuild": true,
  "webhookSecret": "${LOKALISE_WEBHOOK_SECRET}"
}
```

### Step 2: Environment Detection and Config Loading
```typescript
// src/config/lokalise.ts
import baseConfig from "../../config/lokalise/base.json";

type Environment = "development" | "staging" | "production";

interface LokaliseConfig {
  projectId: string;
  branch: string | null;
  apiToken: string;
  enableCompression: boolean;
  cacheEnabled: boolean;
  cacheTtlSeconds: number;
  syncOnBuild: boolean;
  fallbackLocale: string;
  supportedLocales: string[];
  debug?: boolean;
  webhookSecret?: string;
}

function detectEnvironment(): Environment {
  const env = process.env.NODE_ENV || "development";
  if (["development", "staging", "production"].includes(env)) {
    return env as Environment;
  }
  return "development";
}

export function getLokaliseConfig(): LokaliseConfig {
  const env = detectEnvironment();
  const envConfig = require(`../../config/lokalise/${env}.json`);

  const config = {
    ...baseConfig,
    ...envConfig,
    apiToken: process.env.LOKALISE_API_TOKEN!,
  };

  // Resolve environment variable placeholders
  config.projectId = config.projectId.replace(
    /\${(\w+)}/g,
    (_, key) => process.env[key] || ""
  );

  if (config.webhookSecret) {
    config.webhookSecret = process.env.LOKALISE_WEBHOOK_SECRET;
  }

  return config;
}

// Get project ID with branch suffix
export function getProjectIdWithBranch(): string {
  const config = getLokaliseConfig();
  if (config.branch) {
    return `${config.projectId}:${config.branch}`;
  }
  return config.projectId;
}
```

### Step 3: Environment-Specific Clients
```typescript
// src/services/lokalise/client.ts
import { LokaliseApi } from "@lokalise/node-api";
import { getLokaliseConfig, getProjectIdWithBranch } from "../../config/lokalise";

let client: LokaliseApi | null = null;

export function getLokaliseClient(): LokaliseApi {
  if (!client) {
    const config = getLokaliseConfig();
    client = new LokaliseApi({
      apiKey: config.apiToken,
      enableCompression: config.enableCompression,
    });
  }
  return client;
}

export function getProjectId(): string {
  return getProjectIdWithBranch();
}

// For debugging - log current environment
export function logEnvironmentInfo(): void {
  const config = getLokaliseConfig();
  console.log({
    environment: process.env.NODE_ENV,
    projectId: config.projectId,
    branch: config.branch,
    cacheEnabled: config.cacheEnabled,
    debug: config.debug,
  });
}
```

### Step 4: Secret Management by Environment

#### Local Development
```bash
# .env.local (git-ignored)
LOKALISE_API_TOKEN=dev_token_here
LOKALISE_PROJECT_ID=123456789.abcdef
LOKALISE_PROJECT_ID_DEV=987654321.fedcba
```

#### CI/CD (GitHub Actions)
```yaml
jobs:
  deploy:
    environment: ${{ github.ref == 'refs/heads/main' && 'production' || 'staging' }}
    env:
      LOKALISE_API_TOKEN: ${{ secrets.LOKALISE_API_TOKEN }}
      LOKALISE_PROJECT_ID: ${{ secrets.LOKALISE_PROJECT_ID }}
```

#### Production (Secret Manager)
```bash
# AWS Secrets Manager
aws secretsmanager get-secret-value \
  --secret-id lokalise/production \
  --query SecretString --output text

# GCP Secret Manager
gcloud secrets versions access latest --secret=lokalise-api-token

# HashiCorp Vault
vault kv get -field=api_token secret/lokalise/production
```

### Step 5: Environment Isolation Guards
```typescript
// Prevent dangerous operations in wrong environment
function guardProductionOperation(operation: string): void {
  const config = getLokaliseConfig();
  const env = process.env.NODE_ENV;

  if (env !== "production" && config.projectId.includes("prod")) {
    throw new Error(
      `SAFETY: ${operation} targeting production project from ${env} environment`
    );
  }
}

// Example usage
async function deleteAllKeys(projectId: string) {
  guardProductionOperation("deleteAllKeys");
  // ... dangerous operation
}

// Prevent accidental production writes from dev
function isProductionSafe(): boolean {
  const env = process.env.NODE_ENV;
  return env === "production" || process.env.ALLOW_PROD_WRITE === "true";
}
```

## Output
- Multi-environment config structure
- Environment detection logic
- Secure secret management
- Production safeguards enabled

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Wrong environment | Missing NODE_ENV | Set environment variable |
| Secret not found | Wrong secret path | Verify secret manager config |
| Branch not found | Typo in branch name | Check Lokalise project branches |
| Production guard | Safety mechanism | Use correct environment |

## Examples

### Branch-Based Environment Workflow
```typescript
// Use Lokalise branches instead of separate projects
const envToBranch: Record<string, string | null> = {
  development: "dev",
  staging: "staging",
  production: null,  // Main branch
};

function getProjectIdForEnv(env: string): string {
  const projectId = process.env.LOKALISE_PROJECT_ID!;
  const branch = envToBranch[env];
  return branch ? `${projectId}:${branch}` : projectId;
}
```

### CI Pipeline with Environment Matrix
```yaml
# .github/workflows/sync-translations.yml
jobs:
  sync:
    strategy:
      matrix:
        environment: [staging, production]
    environment: ${{ matrix.environment }}
    steps:
      - name: Sync translations
        env:
          LOKALISE_API_TOKEN: ${{ secrets.LOKALISE_API_TOKEN }}
          NODE_ENV: ${{ matrix.environment }}
        run: npm run i18n:pull
```

### Quick Environment Check
```typescript
const config = getLokaliseConfig();
console.log(`
Environment: ${process.env.NODE_ENV}
Project ID: ${config.projectId}
Branch: ${config.branch || 'main'}
Cache: ${config.cacheEnabled ? `enabled (${config.cacheTtlSeconds}s)` : 'disabled'}
`);
```

## Resources
- [Lokalise Branching](https://docs.lokalise.com/en/articles/3183265-project-branching)
- [12-Factor App Config](https://12factor.net/config)
- [AWS Secrets Manager](https://aws.amazon.com/secrets-manager/)

## Next Steps
For observability setup, see `lokalise-observability`.
