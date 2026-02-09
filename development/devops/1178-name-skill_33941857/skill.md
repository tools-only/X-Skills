---
name: openevidence-multi-env-setup
description: |
  Configure OpenEvidence across development, staging, and production environments.
  Use when setting up multiple environments, managing environment-specific configurations,
  or implementing environment promotion strategies for clinical AI applications.
  Trigger with phrases like "openevidence environments", "openevidence staging",
  "openevidence dev setup", "multi-environment openevidence".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Multi-Environment Setup

## Overview
Configure and manage OpenEvidence integrations across development, staging, and production environments with proper isolation and promotion strategies.

## Prerequisites
- OpenEvidence accounts for each environment
- Separate API keys per environment
- Infrastructure for each environment (cloud or on-premise)
- CI/CD pipeline configured

## Environment Strategy

| Environment | API Endpoint | Purpose | Data |
|-------------|--------------|---------|------|
| Development | sandbox.openevidence.com | Local development | Synthetic |
| Staging | sandbox.openevidence.com | Integration testing | Synthetic |
| Production | api.openevidence.com | Live clinical use | Real (PHI) |

## Instructions

### Step 1: Environment Configuration Files
```typescript
// config/environments/development.ts
export const developmentConfig = {
  env: 'development',
  openevidence: {
    baseUrl: 'https://api.sandbox.openevidence.com',
    timeout: 60000, // Longer timeout for debugging
    retries: 1,
    rateLimit: {
      enabled: false, // No rate limiting in dev
    },
  },
  cache: {
    enabled: false, // Disable cache for fresh responses
    ttlSeconds: 0,
  },
  logging: {
    level: 'debug',
    prettyPrint: true,
  },
  features: {
    deepConsult: true,
    webhooks: false, // No webhooks in local dev
    auditLogging: false,
  },
};

// config/environments/staging.ts
export const stagingConfig = {
  env: 'staging',
  openevidence: {
    baseUrl: 'https://api.sandbox.openevidence.com',
    timeout: 45000,
    retries: 2,
    rateLimit: {
      enabled: true,
      requestsPerMinute: 60,
    },
  },
  cache: {
    enabled: true,
    ttlSeconds: 1800, // 30 min cache
  },
  logging: {
    level: 'info',
    prettyPrint: false,
  },
  features: {
    deepConsult: true,
    webhooks: true,
    auditLogging: true,
  },
};

// config/environments/production.ts
export const productionConfig = {
  env: 'production',
  openevidence: {
    baseUrl: 'https://api.openevidence.com',
    timeout: 30000,
    retries: 3,
    rateLimit: {
      enabled: true,
      requestsPerMinute: 300, // Higher limit in prod
    },
  },
  cache: {
    enabled: true,
    ttlSeconds: 3600, // 1 hour cache
  },
  logging: {
    level: 'warn',
    prettyPrint: false,
  },
  features: {
    deepConsult: true,
    webhooks: true,
    auditLogging: true,
  },
};
```

### Step 2: Configuration Loader
```typescript
// src/config/index.ts
import { developmentConfig } from './environments/development';
import { stagingConfig } from './environments/staging';
import { productionConfig } from './environments/production';

type Environment = 'development' | 'staging' | 'production';

const configs = {
  development: developmentConfig,
  staging: stagingConfig,
  production: productionConfig,
};

export function loadConfig(env?: Environment) {
  const environment = env || (process.env.NODE_ENV as Environment) || 'development';

  if (!configs[environment]) {
    throw new Error(`Unknown environment: ${environment}`);
  }

  const baseConfig = configs[environment];

  // Override with environment variables
  return {
    ...baseConfig,
    openevidence: {
      ...baseConfig.openevidence,
      apiKey: process.env.OPENEVIDENCE_API_KEY,
      orgId: process.env.OPENEVIDENCE_ORG_ID,
      ...(process.env.OPENEVIDENCE_BASE_URL && {
        baseUrl: process.env.OPENEVIDENCE_BASE_URL,
      }),
      ...(process.env.OPENEVIDENCE_TIMEOUT && {
        timeout: parseInt(process.env.OPENEVIDENCE_TIMEOUT),
      }),
    },
  };
}

// Singleton config instance
let configInstance: ReturnType<typeof loadConfig> | null = null;

export function getConfig() {
  if (!configInstance) {
    configInstance = loadConfig();
  }
  return configInstance;
}

// For testing - reset config
export function resetConfig() {
  configInstance = null;
}
```

### Step 3: Environment-Aware Client Factory
```typescript
// src/openevidence/client-factory.ts
import { OpenEvidenceClient } from '@openevidence/sdk';
import { getConfig } from '../config';
import { ClinicalQueryCache } from './cache';
import { MetricsCollector } from '../monitoring/metrics';

export class OpenEvidenceClientFactory {
  private static instances: Map<string, OpenEvidenceClient> = new Map();

  static getClient(environment?: string): OpenEvidenceClient {
    const config = getConfig();
    const env = environment || config.env;

    if (!this.instances.has(env)) {
      const client = new OpenEvidenceClient({
        apiKey: config.openevidence.apiKey,
        orgId: config.openevidence.orgId,
        baseUrl: config.openevidence.baseUrl,
        timeout: config.openevidence.timeout,
      });

      this.instances.set(env, client);
    }

    return this.instances.get(env)!;
  }

  // For testing - inject mock client
  static setClient(env: string, client: OpenEvidenceClient): void {
    this.instances.set(env, client);
  }

  static clearClients(): void {
    this.instances.clear();
  }
}
```

### Step 4: Secret Management Per Environment
```typescript
// src/secrets/manager.ts
import { SecretManagerServiceClient } from '@google-cloud/secret-manager';

interface SecretPaths {
  apiKey: string;
  orgId: string;
  webhookSecret?: string;
}

const SECRET_PATHS: Record<string, SecretPaths> = {
  development: {
    apiKey: 'local', // Use .env in development
    orgId: 'local',
  },
  staging: {
    apiKey: 'projects/staging-project/secrets/openevidence-api-key/versions/latest',
    orgId: 'projects/staging-project/secrets/openevidence-org-id/versions/latest',
    webhookSecret: 'projects/staging-project/secrets/openevidence-webhook-secret/versions/latest',
  },
  production: {
    apiKey: 'projects/prod-project/secrets/openevidence-api-key/versions/latest',
    orgId: 'projects/prod-project/secrets/openevidence-org-id/versions/latest',
    webhookSecret: 'projects/prod-project/secrets/openevidence-webhook-secret/versions/latest',
  },
};

export class SecretManager {
  private client: SecretManagerServiceClient;

  constructor() {
    this.client = new SecretManagerServiceClient();
  }

  async getSecrets(env: string): Promise<{
    apiKey: string;
    orgId: string;
    webhookSecret?: string;
  }> {
    const paths = SECRET_PATHS[env];

    if (!paths) {
      throw new Error(`Unknown environment: ${env}`);
    }

    // Development uses environment variables
    if (paths.apiKey === 'local') {
      return {
        apiKey: process.env.OPENEVIDENCE_API_KEY!,
        orgId: process.env.OPENEVIDENCE_ORG_ID!,
        webhookSecret: process.env.OPENEVIDENCE_WEBHOOK_SECRET,
      };
    }

    // Staging/Production use Secret Manager
    const [apiKeyVersion] = await this.client.accessSecretVersion({ name: paths.apiKey });
    const [orgIdVersion] = await this.client.accessSecretVersion({ name: paths.orgId });

    const secrets: any = {
      apiKey: apiKeyVersion.payload?.data?.toString(),
      orgId: orgIdVersion.payload?.data?.toString(),
    };

    if (paths.webhookSecret) {
      const [webhookVersion] = await this.client.accessSecretVersion({ name: paths.webhookSecret });
      secrets.webhookSecret = webhookVersion.payload?.data?.toString();
    }

    return secrets;
  }
}
```

### Step 5: Environment Promotion Workflow
```yaml
# .github/workflows/promotion.yml
name: Environment Promotion

on:
  workflow_dispatch:
    inputs:
      source_env:
        description: 'Source environment'
        required: true
        type: choice
        options:
          - staging
      target_env:
        description: 'Target environment'
        required: true
        type: choice
        options:
          - production
      version:
        description: 'Version to promote (e.g., v1.2.3)'
        required: true

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Validate version exists in source
        run: |
          # Check that the version exists and is deployed in source
          SOURCE_VERSION=$(curl -sf "${{ secrets.STAGING_URL }}/health" | jq -r '.version')
          if [ "$SOURCE_VERSION" != "${{ inputs.version }}" ]; then
            echo "Version mismatch: staging has $SOURCE_VERSION, expected ${{ inputs.version }}"
            exit 1
          fi

      - name: Validate staging tests passed
        run: |
          # Check that all tests passed for this version
          gh run list --workflow=openevidence-ci.yml \
            --branch=main --status=success \
            --json conclusion,headBranch | jq -e '.[0]'

  promote:
    needs: validate
    runs-on: ubuntu-latest
    environment: production
    steps:
      - uses: actions/checkout@v4
        with:
          ref: ${{ inputs.version }}

      - name: Setup Cloud SDK
        uses: google-github-actions/setup-gcloud@v2
        with:
          project_id: ${{ secrets.PROD_PROJECT_ID }}

      - name: Deploy to Production
        run: |
          gcloud run deploy clinical-evidence-api \
            --image gcr.io/${{ secrets.STAGING_PROJECT_ID }}/clinical-evidence-api:${{ inputs.version }} \
            --region us-central1 \
            --set-env-vars NODE_ENV=production \
            --no-traffic

      - name: Run smoke tests
        run: |
          REVISION=$(gcloud run revisions list --service clinical-evidence-api \
            --region us-central1 --format 'value(metadata.name)' --limit 1)
          # Test the new revision directly
          ./scripts/smoke-test.sh "$REVISION"

      - name: Gradual traffic shift
        run: |
          # 10% -> 50% -> 100%
          for percent in 10 50 100; do
            gcloud run services update-traffic clinical-evidence-api \
              --region us-central1 \
              --to-latest=$percent
            sleep 300 # Wait 5 minutes between shifts
          done

      - name: Create production release
        uses: softprops/action-gh-release@v1
        with:
          tag_name: ${{ inputs.version }}-production
          name: Production Release ${{ inputs.version }}
          body: |
            Promoted to production from staging.
            Original version: ${{ inputs.version }}
```

### Step 6: Environment Health Checks
```typescript
// src/health/environment-health.ts
import { OpenEvidenceClientFactory } from '../openevidence/client-factory';
import { getConfig } from '../config';

interface EnvironmentHealth {
  environment: string;
  status: 'healthy' | 'degraded' | 'unhealthy';
  checks: {
    openevidence: { status: string; latencyMs?: number };
    cache: { status: string; hitRate?: number };
    database: { status: string };
  };
  version: string;
  timestamp: string;
}

export async function checkEnvironmentHealth(): Promise<EnvironmentHealth> {
  const config = getConfig();
  const checks: EnvironmentHealth['checks'] = {
    openevidence: { status: 'unknown' },
    cache: { status: 'unknown' },
    database: { status: 'unknown' },
  };

  // Check OpenEvidence connectivity
  try {
    const client = OpenEvidenceClientFactory.getClient();
    const start = Date.now();
    await client.health.check();
    checks.openevidence = {
      status: 'healthy',
      latencyMs: Date.now() - start,
    };
  } catch (error: any) {
    checks.openevidence = { status: 'unhealthy' };
  }

  // Check cache
  try {
    if (config.cache.enabled) {
      const stats = await cache.getStats();
      checks.cache = {
        status: 'healthy',
        hitRate: stats.hitRate,
      };
    } else {
      checks.cache = { status: 'disabled' };
    }
  } catch {
    checks.cache = { status: 'unhealthy' };
  }

  // Check database
  try {
    await db.$queryRaw`SELECT 1`;
    checks.database = { status: 'healthy' };
  } catch {
    checks.database = { status: 'unhealthy' };
  }

  // Determine overall status
  const unhealthyCount = Object.values(checks).filter(
    c => c.status === 'unhealthy'
  ).length;

  let status: EnvironmentHealth['status'] = 'healthy';
  if (unhealthyCount > 0) status = 'degraded';
  if (checks.openevidence.status === 'unhealthy') status = 'unhealthy';

  return {
    environment: config.env,
    status,
    checks,
    version: process.env.npm_package_version || 'unknown',
    timestamp: new Date().toISOString(),
  };
}
```

## Output
- Environment-specific configuration files
- Configuration loader with env var overrides
- Secret management per environment
- Promotion workflow
- Health checks per environment

## Environment Checklist

### Development
- [ ] Sandbox API key configured
- [ ] Local .env file created
- [ ] Cache disabled for fresh responses
- [ ] Debug logging enabled

### Staging
- [ ] Separate sandbox API key
- [ ] Secrets in Secret Manager
- [ ] Webhooks configured
- [ ] HIPAA audit logging enabled
- [ ] Performance testing completed

### Production
- [ ] Production API key from OpenEvidence
- [ ] BAA signed
- [ ] All secrets in Secret Manager
- [ ] Monitoring and alerting configured
- [ ] Runbook documented
- [ ] On-call rotation established

## Error Handling
| Environment Issue | Detection | Resolution |
|-------------------|-----------|------------|
| Wrong API endpoint | Health check fails | Verify baseUrl in config |
| Secret not found | Startup failure | Check Secret Manager permissions |
| Config mismatch | Unexpected behavior | Validate config loading |
| Promotion failure | CI/CD error | Check version tagging |

## Resources
- [OpenEvidence Sandbox](https://sandbox.openevidence.com/)
- [Google Cloud Secret Manager](https://cloud.google.com/secret-manager)

## Next Steps
For observability setup, see `openevidence-observability`.
