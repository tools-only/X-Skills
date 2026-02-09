---
name: langfuse-multi-env-setup
description: |
  Configure Langfuse across development, staging, and production environments.
  Use when setting up multi-environment deployments, configuring per-environment keys,
  or implementing environment-specific Langfuse configurations.
  Trigger with phrases like "langfuse environments", "langfuse staging",
  "langfuse dev prod", "langfuse environment setup", "langfuse config by env".
allowed-tools: Read, Write, Edit, Bash(aws:*), Bash(gcloud:*), Bash(vault:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Multi-Environment Setup

## Overview
Configure Langfuse across development, staging, and production environments.

## Prerequisites
- Separate Langfuse projects or API keys per environment
- Secret management solution (Vault, AWS Secrets Manager, etc.)
- CI/CD pipeline with environment variables
- Environment detection in application

## Environment Strategy

| Environment | Purpose | API Keys | Sampling |
|-------------|---------|----------|----------|
| Development | Local dev | Dev project | 100% |
| Staging | Pre-prod validation | Staging project | 50% |
| Production | Live traffic | Production project | 10-100% |

## Instructions

### Step 1: Create Configuration Structure

```
config/
├── langfuse/
│   ├── base.ts           # Shared config
│   ├── development.ts    # Dev overrides
│   ├── staging.ts        # Staging overrides
│   └── production.ts     # Prod overrides
```

### Step 2: Implement Base Configuration

```typescript
// config/langfuse/base.ts
export interface LangfuseConfig {
  publicKey: string;
  secretKey: string;
  host: string;
  enabled: boolean;
  flushAt: number;
  flushInterval: number;
  requestTimeout: number;
  sampling: {
    rate: number;
    alwaysSampleErrors: boolean;
    alwaysSampleSlowTraces: boolean;
    slowTraceThresholdMs: number;
  };
  piiScrubbing: {
    enabled: boolean;
    fields: string[];
  };
}

export const baseConfig: Partial<LangfuseConfig> = {
  host: "https://cloud.langfuse.com",
  requestTimeout: 15000,
  sampling: {
    rate: 1.0,
    alwaysSampleErrors: true,
    alwaysSampleSlowTraces: true,
    slowTraceThresholdMs: 5000,
  },
  piiScrubbing: {
    enabled: true,
    fields: ["email", "password", "ssn", "creditCard", "phone"],
  },
};
```

### Step 3: Create Environment-Specific Configs

```typescript
// config/langfuse/development.ts
import { LangfuseConfig, baseConfig } from "./base";

export const developmentConfig: LangfuseConfig = {
  ...baseConfig,
  publicKey: process.env.LANGFUSE_PUBLIC_KEY_DEV!,
  secretKey: process.env.LANGFUSE_SECRET_KEY_DEV!,
  host: process.env.LANGFUSE_HOST_DEV || "http://localhost:3000",
  enabled: true,
  flushAt: 1,           // Immediate flush for debugging
  flushInterval: 1000,  // 1 second
  sampling: {
    ...baseConfig.sampling!,
    rate: 1.0,          // Sample everything in dev
  },
  piiScrubbing: {
    ...baseConfig.piiScrubbing!,
    enabled: false,     // Disable for debugging
  },
} as LangfuseConfig;
```

```typescript
// config/langfuse/staging.ts
import { LangfuseConfig, baseConfig } from "./base";

export const stagingConfig: LangfuseConfig = {
  ...baseConfig,
  publicKey: process.env.LANGFUSE_PUBLIC_KEY_STAGING!,
  secretKey: process.env.LANGFUSE_SECRET_KEY_STAGING!,
  host: process.env.LANGFUSE_HOST_STAGING || "https://cloud.langfuse.com",
  enabled: true,
  flushAt: 15,
  flushInterval: 5000,
  sampling: {
    ...baseConfig.sampling!,
    rate: 0.5,          // 50% sampling in staging
  },
} as LangfuseConfig;
```

```typescript
// config/langfuse/production.ts
import { LangfuseConfig, baseConfig } from "./base";

export const productionConfig: LangfuseConfig = {
  ...baseConfig,
  publicKey: process.env.LANGFUSE_PUBLIC_KEY_PROD!,
  secretKey: process.env.LANGFUSE_SECRET_KEY_PROD!,
  host: process.env.LANGFUSE_HOST || "https://cloud.langfuse.com",
  enabled: true,
  flushAt: 25,
  flushInterval: 5000,
  sampling: {
    ...baseConfig.sampling!,
    rate: parseFloat(process.env.LANGFUSE_SAMPLE_RATE || "0.1"),
  },
} as LangfuseConfig;
```

### Step 4: Environment Detection and Config Loading

```typescript
// config/langfuse/index.ts
import { LangfuseConfig } from "./base";
import { developmentConfig } from "./development";
import { stagingConfig } from "./staging";
import { productionConfig } from "./production";

type Environment = "development" | "staging" | "production";

const configs: Record<Environment, LangfuseConfig> = {
  development: developmentConfig,
  staging: stagingConfig,
  production: productionConfig,
};

function detectEnvironment(): Environment {
  // Check explicit environment variable
  const explicitEnv = process.env.LANGFUSE_ENV as Environment;
  if (explicitEnv && configs[explicitEnv]) {
    return explicitEnv;
  }

  // Fall back to NODE_ENV
  const nodeEnv = process.env.NODE_ENV;

  if (nodeEnv === "production") {
    // Distinguish staging from production
    if (
      process.env.VERCEL_ENV === "preview" ||
      process.env.DEPLOYMENT_ENV === "staging"
    ) {
      return "staging";
    }
    return "production";
  }

  return "development";
}

export function getLangfuseConfig(): LangfuseConfig & { environment: Environment } {
  const environment = detectEnvironment();
  const config = configs[environment];

  // Validate required fields
  if (!config.publicKey || !config.secretKey) {
    console.warn(
      `Langfuse: Missing credentials for ${environment} environment`
    );
  }

  return {
    ...config,
    environment,
  };
}

// Singleton export
let cachedConfig: ReturnType<typeof getLangfuseConfig> | null = null;

export function getConfig(): ReturnType<typeof getLangfuseConfig> {
  if (!cachedConfig) {
    cachedConfig = getLangfuseConfig();
    console.log(`Langfuse: Loaded ${cachedConfig.environment} configuration`);
  }
  return cachedConfig;
}
```

### Step 5: Create Environment-Aware Langfuse Client

```typescript
// lib/langfuse.ts
import { Langfuse } from "langfuse";
import { getConfig, LangfuseConfig } from "../config/langfuse";

class LangfuseFactory {
  private static instance: Langfuse | null = null;
  private static config: LangfuseConfig | null = null;

  static getInstance(): Langfuse {
    if (!LangfuseFactory.instance) {
      const config = getConfig();
      LangfuseFactory.config = config;

      if (!config.enabled) {
        console.log("Langfuse: Disabled by configuration");
        return LangfuseFactory.createNoOpClient();
      }

      LangfuseFactory.instance = new Langfuse({
        publicKey: config.publicKey,
        secretKey: config.secretKey,
        baseUrl: config.host,
        flushAt: config.flushAt,
        flushInterval: config.flushInterval,
        requestTimeout: config.requestTimeout,
      });

      LangfuseFactory.registerShutdown();
    }

    return LangfuseFactory.instance!;
  }

  static getConfig(): LangfuseConfig | null {
    return LangfuseFactory.config;
  }

  private static createNoOpClient(): Langfuse {
    // Return a mock that does nothing
    return {
      trace: () => LangfuseFactory.createNoOpTrace(),
      flushAsync: async () => {},
      shutdownAsync: async () => {},
    } as unknown as Langfuse;
  }

  private static createNoOpTrace() {
    return {
      id: "noop",
      span: () => LangfuseFactory.createNoOpTrace(),
      generation: () => ({ id: "noop", end: () => {} }),
      update: () => {},
      getTraceUrl: () => "",
    };
  }

  private static registerShutdown() {
    const shutdown = async () => {
      if (LangfuseFactory.instance) {
        console.log("Langfuse: Shutting down...");
        await LangfuseFactory.instance.shutdownAsync();
        LangfuseFactory.instance = null;
      }
    };

    process.on("SIGTERM", shutdown);
    process.on("SIGINT", shutdown);
  }
}

export const langfuse = LangfuseFactory.getInstance();
export const langfuseConfig = LangfuseFactory.getConfig();
```

### Step 6: Secret Management Integration

```typescript
// lib/secrets/langfuse.ts
import { SecretsManagerClient, GetSecretValueCommand } from "@aws-sdk/client-secrets-manager";
import { SecretManagerServiceClient } from "@google-cloud/secret-manager";

type SecretProvider = "aws" | "gcp" | "vault" | "env";

interface LangfuseSecrets {
  publicKey: string;
  secretKey: string;
}

async function getSecretsFromAWS(secretName: string): Promise<LangfuseSecrets> {
  const client = new SecretsManagerClient({});
  const response = await client.send(
    new GetSecretValueCommand({ SecretId: secretName })
  );
  return JSON.parse(response.SecretString!);
}

async function getSecretsFromGCP(secretPath: string): Promise<LangfuseSecrets> {
  const client = new SecretManagerServiceClient();
  const [publicKeyVersion] = await client.accessSecretVersion({
    name: `${secretPath}/langfuse-public-key/versions/latest`,
  });
  const [secretKeyVersion] = await client.accessSecretVersion({
    name: `${secretPath}/langfuse-secret-key/versions/latest`,
  });

  return {
    publicKey: publicKeyVersion.payload?.data?.toString() || "",
    secretKey: secretKeyVersion.payload?.data?.toString() || "",
  };
}

export async function getLangfuseSecrets(
  provider: SecretProvider = "env",
  options?: { secretName?: string; secretPath?: string }
): Promise<LangfuseSecrets> {
  switch (provider) {
    case "aws":
      return getSecretsFromAWS(
        options?.secretName || `langfuse/${process.env.NODE_ENV}`
      );

    case "gcp":
      return getSecretsFromGCP(
        options?.secretPath || `projects/${process.env.GCP_PROJECT}/secrets`
      );

    case "env":
    default:
      return {
        publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
        secretKey: process.env.LANGFUSE_SECRET_KEY!,
      };
  }
}
```

## Output
- Multi-environment configuration structure
- Environment detection logic
- Environment-aware Langfuse client
- Secret management integration
- Graceful degradation when disabled

## Environment Variables Reference

```bash
# Development
LANGFUSE_PUBLIC_KEY_DEV=pk-lf-dev-...
LANGFUSE_SECRET_KEY_DEV=sk-lf-dev-...
LANGFUSE_HOST_DEV=http://localhost:3000

# Staging
LANGFUSE_PUBLIC_KEY_STAGING=pk-lf-staging-...
LANGFUSE_SECRET_KEY_STAGING=sk-lf-staging-...
LANGFUSE_HOST_STAGING=https://cloud.langfuse.com

# Production
LANGFUSE_PUBLIC_KEY_PROD=pk-lf-prod-...
LANGFUSE_SECRET_KEY_PROD=sk-lf-prod-...
LANGFUSE_HOST=https://cloud.langfuse.com
LANGFUSE_SAMPLE_RATE=0.1

# Override environment detection
LANGFUSE_ENV=staging
```

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Wrong environment | Missing NODE_ENV | Set environment variable |
| Secret not found | Wrong secret path | Verify secret manager config |
| Config merge fails | Invalid values | Validate config on load |
| Cross-env data leak | Shared keys | Use separate projects per env |

## Resources
- [Langfuse Projects](https://langfuse.com/docs)
- [12-Factor App Config](https://12factor.net/config)
- [AWS Secrets Manager](https://aws.amazon.com/secrets-manager/)
- [GCP Secret Manager](https://cloud.google.com/secret-manager)

## Next Steps
For observability setup, see `langfuse-observability`.
