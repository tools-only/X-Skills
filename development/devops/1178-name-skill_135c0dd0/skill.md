---
name: twinmind-multi-env-setup
description: |
  Set up TwinMind across multiple environments (dev, staging, production).
  Use when configuring environment-specific settings, managing multiple API keys,
  or implementing environment promotion workflows.
  Trigger with phrases like "twinmind environments", "twinmind multi-env",
  "twinmind staging setup", "twinmind dev vs prod".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Multi-Environment Setup

## Overview
Configure TwinMind for development, staging, and production environments with proper isolation and promotion workflows.

## Prerequisites
- Multiple TwinMind accounts or API keys per environment
- Understanding of environment management
- CI/CD pipeline for deployments

## Instructions

### Step 1: Environment Configuration Structure

```typescript
// src/config/environments.ts
export type Environment = 'development' | 'staging' | 'production';

export interface EnvironmentConfig {
  name: Environment;
  twinmind: {
    apiKey: string;
    baseUrl: string;
    webhookSecret: string;
    model: 'ear-2' | 'ear-3' | 'ear-3-custom';
    features: {
      diarization: boolean;
      autoSummary: boolean;
      actionItems: boolean;
      debug: boolean;
    };
  };
  limits: {
    rateLimit: number;
    concurrentJobs: number;
    dailyQuota?: number;
  };
  integrations: {
    slack?: { webhook: string; channel: string };
    calendar?: { enabled: boolean };
    email?: { enabled: boolean };
  };
}

// Environment-specific configurations
const configs: Record<Environment, Partial<EnvironmentConfig>> = {
  development: {
    name: 'development',
    twinmind: {
      apiKey: process.env.TWINMIND_API_KEY_DEV!,
      baseUrl: 'https://api.twinmind.com/v1',
      webhookSecret: process.env.TWINMIND_WEBHOOK_SECRET_DEV!,
      model: 'ear-2',  // Use cheaper model in dev
      features: {
        diarization: false,  // Faster without diarization
        autoSummary: true,
        actionItems: true,
        debug: true,
      },
    },
    limits: {
      rateLimit: 30,
      concurrentJobs: 1,
      dailyQuota: 10,  // Limit dev usage
    },
    integrations: {
      slack: { webhook: process.env.SLACK_WEBHOOK_DEV!, channel: '#dev-alerts' },
      calendar: { enabled: false },
      email: { enabled: false },
    },
  },

  staging: {
    name: 'staging',
    twinmind: {
      apiKey: process.env.TWINMIND_API_KEY_STAGING!,
      baseUrl: 'https://api.twinmind.com/v1',
      webhookSecret: process.env.TWINMIND_WEBHOOK_SECRET_STAGING!,
      model: 'ear-3',
      features: {
        diarization: true,
        autoSummary: true,
        actionItems: true,
        debug: true,
      },
    },
    limits: {
      rateLimit: 60,
      concurrentJobs: 2,
      dailyQuota: 50,
    },
    integrations: {
      slack: { webhook: process.env.SLACK_WEBHOOK_STAGING!, channel: '#staging-alerts' },
      calendar: { enabled: true },
      email: { enabled: false },  // No real emails in staging
    },
  },

  production: {
    name: 'production',
    twinmind: {
      apiKey: process.env.TWINMIND_API_KEY_PROD!,
      baseUrl: 'https://api.twinmind.com/v1',
      webhookSecret: process.env.TWINMIND_WEBHOOK_SECRET_PROD!,
      model: 'ear-3',
      features: {
        diarization: true,
        autoSummary: true,
        actionItems: true,
        debug: false,  // No debug in production
      },
    },
    limits: {
      rateLimit: 300,  // Enterprise limits
      concurrentJobs: 10,
    },
    integrations: {
      slack: { webhook: process.env.SLACK_WEBHOOK_PROD!, channel: '#production-alerts' },
      calendar: { enabled: true },
      email: { enabled: true },
    },
  },
};

export function getEnvironmentConfig(): EnvironmentConfig {
  const env = (process.env.NODE_ENV || 'development') as Environment;
  const config = configs[env];

  if (!config) {
    throw new Error(`Unknown environment: ${env}`);
  }

  return config as EnvironmentConfig;
}
```

### Step 2: Environment Variable Files

```bash
# .env.development
NODE_ENV=development
TWINMIND_API_KEY_DEV=tm_sk_dev_xxx
TWINMIND_WEBHOOK_SECRET_DEV=whsec_dev_xxx
SLACK_WEBHOOK_DEV=https://hooks.slack.com/services/xxx

# .env.staging
NODE_ENV=staging
TWINMIND_API_KEY_STAGING=tm_sk_staging_xxx
TWINMIND_WEBHOOK_SECRET_STAGING=whsec_staging_xxx
SLACK_WEBHOOK_STAGING=https://hooks.slack.com/services/xxx

# .env.production
NODE_ENV=production
TWINMIND_API_KEY_PROD=tm_sk_prod_xxx
TWINMIND_WEBHOOK_SECRET_PROD=whsec_prod_xxx
SLACK_WEBHOOK_PROD=https://hooks.slack.com/services/xxx

# .env.example (committed to repo)
NODE_ENV=development
TWINMIND_API_KEY_DEV=
TWINMIND_WEBHOOK_SECRET_DEV=
```

### Step 3: Environment-Aware Client Factory

```typescript
// src/twinmind/factory.ts
import { TwinMindClient } from './client';
import { getEnvironmentConfig, Environment } from '../config/environments';

// Client instances per environment
const clients = new Map<Environment, TwinMindClient>();

export function getTwinMindClient(env?: Environment): TwinMindClient {
  const config = getEnvironmentConfig();
  const targetEnv = env || config.name;

  if (!clients.has(targetEnv)) {
    const envConfig = configs[targetEnv];

    clients.set(targetEnv, new TwinMindClient({
      apiKey: envConfig.twinmind.apiKey,
      baseUrl: envConfig.twinmind.baseUrl,
      timeout: 30000,
      debug: envConfig.twinmind.features.debug,
    }));
  }

  return clients.get(targetEnv)!;
}

// Reset client (useful for testing)
export function resetClient(env?: Environment): void {
  if (env) {
    clients.delete(env);
  } else {
    clients.clear();
  }
}
```

### Step 4: Environment Validation

```typescript
// scripts/validate-env.ts
import { Environment } from '../src/config/environments';

interface ValidationResult {
  environment: string;
  valid: boolean;
  issues: string[];
}

async function validateEnvironment(env: Environment): Promise<ValidationResult> {
  const issues: string[] = [];

  // Check required environment variables
  const requiredVars = [
    `TWINMIND_API_KEY_${env.toUpperCase()}`,
    `TWINMIND_WEBHOOK_SECRET_${env.toUpperCase()}`,
  ];

  for (const varName of requiredVars) {
    if (!process.env[varName]) {
      issues.push(`Missing required variable: ${varName}`);
    }
  }

  // Validate API key format
  const apiKey = process.env[`TWINMIND_API_KEY_${env.toUpperCase()}`];
  if (apiKey && !apiKey.startsWith('tm_sk_')) {
    issues.push('API key should start with "tm_sk_"');
  }

  // Test API connectivity
  if (apiKey) {
    try {
      const response = await fetch('https://api.twinmind.com/v1/health', {
        headers: { 'Authorization': `Bearer ${apiKey}` },
      });
      if (!response.ok) {
        issues.push(`API health check failed: ${response.status}`);
      }
    } catch (error: any) {
      issues.push(`API connectivity failed: ${error.message}`);
    }
  }

  return {
    environment: env,
    valid: issues.length === 0,
    issues,
  };
}

async function validateAllEnvironments(): Promise<void> {
  const environments: Environment[] = ['development', 'staging', 'production'];

  console.log('Validating TwinMind environments...\n');

  for (const env of environments) {
    const result = await validateEnvironment(env);

    console.log(`${env.toUpperCase()}: ${result.valid ? 'VALID' : 'INVALID'}`);

    if (result.issues.length > 0) {
      result.issues.forEach(issue => console.log(`  - ${issue}`));
    }
    console.log();
  }
}

validateAllEnvironments();
```

### Step 5: Promotion Workflow

```typescript
// scripts/promote-config.ts
import * as fs from 'fs';
import * as path from 'path';

interface PromotionResult {
  success: boolean;
  fromEnv: string;
  toEnv: string;
  changes: string[];
}

async function promoteConfig(
  fromEnv: 'development' | 'staging',
  toEnv: 'staging' | 'production'
): Promise<PromotionResult> {
  const changes: string[] = [];

  // Read source config
  const sourceConfigPath = path.join(__dirname, `../config/twinmind.${fromEnv}.json`);
  const sourceConfig = JSON.parse(fs.readFileSync(sourceConfigPath, 'utf-8'));

  // Read target config
  const targetConfigPath = path.join(__dirname, `../config/twinmind.${toEnv}.json`);
  const targetConfig = JSON.parse(fs.readFileSync(targetConfigPath, 'utf-8'));

  // Identify changes to promote
  const promotableKeys = ['features', 'webhookEvents', 'integrationSettings'];

  for (const key of promotableKeys) {
    if (JSON.stringify(sourceConfig[key]) !== JSON.stringify(targetConfig[key])) {
      changes.push(`${key}: ${JSON.stringify(sourceConfig[key])}`);
      targetConfig[key] = sourceConfig[key];
    }
  }

  // Write updated target config
  if (changes.length > 0) {
    fs.writeFileSync(targetConfigPath, JSON.stringify(targetConfig, null, 2));
  }

  return {
    success: true,
    fromEnv,
    toEnv,
    changes,
  };
}

// CLI interface
async function main(): Promise<void> {
  const args = process.argv.slice(2);
  const fromEnv = args[0] as 'development' | 'staging';
  const toEnv = args[1] as 'staging' | 'production';

  if (!fromEnv || !toEnv) {
    console.log('Usage: ts-node promote-config.ts <from-env> <to-env>');
    console.log('Example: ts-node promote-config.ts staging production');
    process.exit(1);
  }

  const result = await promoteConfig(fromEnv, toEnv);

  console.log(`\nPromotion: ${fromEnv} -> ${toEnv}`);
  console.log(`Changes: ${result.changes.length}`);

  if (result.changes.length > 0) {
    console.log('\nPromoted settings:');
    result.changes.forEach(c => console.log(`  - ${c}`));
  } else {
    console.log('\nNo changes to promote.');
  }
}

main();
```

### Step 6: Feature Flags per Environment

```typescript
// src/features/flags.ts
import { getEnvironmentConfig, Environment } from '../config/environments';

export interface FeatureFlags {
  enableDiarization: boolean;
  enableAutoSummary: boolean;
  enableActionItems: boolean;
  enableRealTimeTranscription: boolean;
  enableEmailFollowUps: boolean;
  maxAudioDurationMinutes: number;
  enableBetaFeatures: boolean;
}

const defaultFlags: FeatureFlags = {
  enableDiarization: false,
  enableAutoSummary: true,
  enableActionItems: true,
  enableRealTimeTranscription: false,
  enableEmailFollowUps: false,
  maxAudioDurationMinutes: 60,
  enableBetaFeatures: false,
};

const environmentFlags: Record<Environment, Partial<FeatureFlags>> = {
  development: {
    enableBetaFeatures: true,
    maxAudioDurationMinutes: 10,
  },
  staging: {
    enableDiarization: true,
    enableBetaFeatures: true,
    maxAudioDurationMinutes: 30,
  },
  production: {
    enableDiarization: true,
    enableRealTimeTranscription: true,
    enableEmailFollowUps: true,
    maxAudioDurationMinutes: 120,
    enableBetaFeatures: false,
  },
};

export function getFeatureFlags(): FeatureFlags {
  const config = getEnvironmentConfig();
  return {
    ...defaultFlags,
    ...environmentFlags[config.name],
  };
}

export function isFeatureEnabled(feature: keyof FeatureFlags): boolean {
  const flags = getFeatureFlags();
  return Boolean(flags[feature]);
}
```

## Output
- Environment configuration structure
- Environment variable templates
- Client factory with environment awareness
- Environment validation script
- Config promotion workflow
- Feature flags per environment

## Environment Matrix

| Setting | Development | Staging | Production |
|---------|-------------|---------|------------|
| Model | ear-2 | ear-3 | ear-3 |
| Diarization | Off | On | On |
| Debug mode | On | On | Off |
| Rate limit | 30/min | 60/min | 300/min |
| Email sending | Off | Off | On |
| Daily quota | 10 hrs | 50 hrs | Unlimited |

## Error Handling

| Issue | Cause | Solution |
|-------|-------|----------|
| Wrong API key used | Env mismatch | Validate before deploy |
| Feature not available | Wrong environment | Check feature flags |
| Config not loading | Missing file | Verify config paths |
| Promotion failed | Schema mismatch | Validate configs first |

## Resources
- [12-Factor App Config](https://12factor.net/config)
- [Environment Best Practices](https://twinmind.com/docs/environments)

## Next Steps
For monitoring setup, see `twinmind-observability`.
