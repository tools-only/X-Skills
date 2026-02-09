---
name: langfuse-prod-checklist
description: |
  Langfuse production readiness checklist and verification.
  Use when preparing to deploy Langfuse to production,
  validating production configuration, or auditing existing setup.
  Trigger with phrases like "langfuse production", "langfuse prod ready",
  "deploy langfuse", "langfuse checklist", "langfuse go live".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Production Checklist

## Overview
Comprehensive checklist for deploying Langfuse observability to production.

## Prerequisites
- Completed development and staging testing
- Access to production secrets management
- Production Langfuse project created

## Production Readiness Checklist

### 1. Authentication & Security
| Item | Status | Command/Action |
|------|--------|----------------|
| Production API keys | [ ] | Create in Langfuse dashboard |
| Keys in secrets manager | [ ] | AWS/GCP/Vault |
| Keys NOT in code/git | [ ] | `git log --all -- '*key*'` |
| Environment isolation | [ ] | Separate dev/staging/prod keys |
| PII scrubbing enabled | [ ] | Verify scrub functions |

### 2. SDK Configuration
| Item | Status | Recommended Value |
|------|--------|-------------------|
| `flushAt` | [ ] | 15-50 (batch size) |
| `flushInterval` | [ ] | 5000-10000ms |
| `requestTimeout` | [ ] | 10000-30000ms |
| `enabled` | [ ] | `true` for prod |
| Singleton pattern | [ ] | Single client instance |

### 3. Error Handling
| Item | Status | Implementation |
|------|--------|----------------|
| Try/catch around traces | [ ] | All trace operations |
| Span `.end()` in finally | [ ] | Prevent hanging spans |
| Error level traces | [ ] | Log errors to Langfuse |
| Graceful shutdown | [ ] | `shutdownAsync()` on exit |

### 4. Performance
| Item | Status | Verification |
|------|--------|--------------|
| Async tracing (non-blocking) | [ ] | Don't await traces in hot path |
| Batching configured | [ ] | `flushAt` > 1 |
| No memory leaks | [ ] | Monitor heap usage |
| Connection pooling | [ ] | Singleton client |

### 5. Observability
| Item | Status | Setup |
|------|--------|-------|
| Trace URL logging | [ ] | Log trace IDs for debugging |
| SDK error monitoring | [ ] | Catch/log SDK errors |
| Dashboard alerts | [ ] | High error rate, latency |
| Cost tracking | [ ] | Token usage monitoring |

## Instructions

### Step 1: Verify Production Configuration

```typescript
// config/langfuse.prod.ts
import { Langfuse } from "langfuse";

function validateConfig(config: Record<string, any>): void {
  const required = ["publicKey", "secretKey"];
  const missing = required.filter((key) => !config[key]);

  if (missing.length > 0) {
    throw new Error(`Missing Langfuse config: ${missing.join(", ")}`);
  }

  // Verify production keys (not dev/test)
  if (config.publicKey.includes("test") || config.publicKey.includes("dev")) {
    console.warn("WARNING: Using non-production API keys!");
  }
}

export function createProductionLangfuse() {
  const config = {
    publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
    secretKey: process.env.LANGFUSE_SECRET_KEY!,
    baseUrl: process.env.LANGFUSE_HOST || "https://cloud.langfuse.com",
    // Production settings
    flushAt: 25,
    flushInterval: 5000,
    requestTimeout: 15000,
    enabled: process.env.NODE_ENV === "production",
  };

  validateConfig(config);

  return new Langfuse(config);
}
```

### Step 2: Implement Production Error Handling

```typescript
// lib/langfuse-production.ts
import { Langfuse } from "langfuse";

let langfuseInstance: Langfuse | null = null;

export function getLangfuse(): Langfuse {
  if (!langfuseInstance) {
    langfuseInstance = createProductionLangfuse();

    // Register graceful shutdown
    const shutdown = async () => {
      console.log("Flushing Langfuse traces...");
      await langfuseInstance?.shutdownAsync();
      console.log("Langfuse shutdown complete");
    };

    process.on("SIGTERM", shutdown);
    process.on("SIGINT", shutdown);
    process.on("beforeExit", shutdown);
  }

  return langfuseInstance;
}

// Production-safe trace wrapper
export async function safeTrace<T>(
  name: string,
  operation: (trace: ReturnType<typeof langfuse.trace>) => Promise<T>,
  metadata?: Record<string, any>
): Promise<T> {
  const langfuse = getLangfuse();

  let trace;
  try {
    trace = langfuse.trace({
      name,
      metadata: {
        ...metadata,
        environment: "production",
      },
    });
  } catch (error) {
    console.error("Failed to create Langfuse trace:", error);
    // Continue without tracing - don't fail the operation
    return operation(null as any);
  }

  try {
    const result = await operation(trace);
    trace.update({ output: { success: true } });
    return result;
  } catch (error) {
    trace.update({
      output: { error: String(error) },
      level: "ERROR",
    });
    throw error;
  }
}
```

### Step 3: Run Pre-Deployment Verification

```typescript
// scripts/verify-langfuse-prod.ts
import { Langfuse } from "langfuse";

async function verifyProduction() {
  console.log("=== Langfuse Production Verification ===\n");

  const checks: Array<{ name: string; check: () => Promise<boolean> }> = [
    {
      name: "Environment variables set",
      check: async () => {
        const required = [
          "LANGFUSE_PUBLIC_KEY",
          "LANGFUSE_SECRET_KEY",
        ];
        const missing = required.filter((key) => !process.env[key]);
        if (missing.length > 0) {
          console.log(`  Missing: ${missing.join(", ")}`);
          return false;
        }
        return true;
      },
    },
    {
      name: "API keys are production keys",
      check: async () => {
        const key = process.env.LANGFUSE_PUBLIC_KEY || "";
        const isProduction =
          !key.includes("test") &&
          !key.includes("dev") &&
          key.startsWith("pk-lf-");
        if (!isProduction) {
          console.log("  Key appears to be non-production");
        }
        return isProduction;
      },
    },
    {
      name: "Can create trace",
      check: async () => {
        const langfuse = new Langfuse();
        const trace = langfuse.trace({
          name: "production-verification",
          metadata: { test: true },
        });
        await langfuse.flushAsync();
        return true;
      },
    },
    {
      name: "Graceful shutdown works",
      check: async () => {
        const langfuse = new Langfuse();
        langfuse.trace({ name: "shutdown-test" });
        await langfuse.shutdownAsync();
        return true;
      },
    },
  ];

  let passed = 0;
  let failed = 0;

  for (const { name, check } of checks) {
    process.stdout.write(`[ ] ${name}...`);
    try {
      const result = await check();
      if (result) {
        console.log("\r[✓]");
        passed++;
      } else {
        console.log("\r[✗]");
        failed++;
      }
    } catch (error) {
      console.log(`\r[✗] Error: ${error}`);
      failed++;
    }
  }

  console.log(`\n=== Results: ${passed} passed, ${failed} failed ===`);

  if (failed > 0) {
    process.exit(1);
  }
}

verifyProduction();
```

### Step 4: Configure Production Monitoring

```typescript
// Monitor Langfuse health in production
interface LangfuseMetrics {
  tracesCreated: number;
  tracesFailed: number;
  flushLatencyMs: number[];
  lastFlushTime: Date | null;
}

const metrics: LangfuseMetrics = {
  tracesCreated: 0,
  tracesFailed: 0,
  flushLatencyMs: [],
  lastFlushTime: null,
};

// Expose metrics endpoint
app.get("/metrics/langfuse", (req, res) => {
  const avgFlushLatency =
    metrics.flushLatencyMs.length > 0
      ? metrics.flushLatencyMs.reduce((a, b) => a + b) /
        metrics.flushLatencyMs.length
      : 0;

  res.json({
    tracesCreated: metrics.tracesCreated,
    tracesFailed: metrics.tracesFailed,
    avgFlushLatencyMs: avgFlushLatency.toFixed(2),
    lastFlushTime: metrics.lastFlushTime?.toISOString(),
    errorRate: (
      (metrics.tracesFailed / (metrics.tracesCreated || 1)) *
      100
    ).toFixed(2),
  });
});
```

## Output
- Verified production configuration
- Error handling implemented
- Graceful shutdown configured
- Monitoring in place
- All checklist items passed

## Production Configuration Example

```typescript
// Final production Langfuse configuration
export const productionLangfuseConfig = {
  // Authentication
  publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
  secretKey: process.env.LANGFUSE_SECRET_KEY!,
  baseUrl: process.env.LANGFUSE_HOST || "https://cloud.langfuse.com",

  // Batching (optimized for production)
  flushAt: 25,          // Batch 25 events
  flushInterval: 5000,  // Flush every 5 seconds

  // Timeouts
  requestTimeout: 15000, // 15 second timeout

  // Enable/disable
  enabled: process.env.NODE_ENV === "production",

  // Debug (disabled in production)
  debug: false,
};
```

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Missing traces in prod | Flush not called | Verify shutdown handler |
| High latency | Large batches | Reduce `flushAt` |
| Memory growth | Client recreation | Use singleton pattern |
| Lost traces on deploy | No graceful shutdown | Add SIGTERM handler |

## Resources
- [Langfuse Production Guide](https://langfuse.com/docs/deployment)
- [Langfuse Self-Hosting](https://langfuse.com/docs/deployment/self-host)
- [Langfuse SDK Reference](https://langfuse.com/docs/sdk)

## Next Steps
For SDK upgrades, see `langfuse-upgrade-migration`.
