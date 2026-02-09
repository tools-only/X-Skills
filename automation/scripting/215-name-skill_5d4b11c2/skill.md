---
name: langfuse-reference-architecture
description: |
  Production-grade Langfuse architecture patterns and best practices.
  Use when designing LLM observability infrastructure, planning Langfuse deployment,
  or implementing enterprise-grade tracing architecture.
  Trigger with phrases like "langfuse architecture", "langfuse design",
  "langfuse infrastructure", "langfuse enterprise", "langfuse at scale".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Reference Architecture

## Overview
Production-grade architecture patterns for Langfuse LLM observability at scale.

## Prerequisites
- Understanding of distributed systems
- Knowledge of cloud infrastructure
- Familiarity with observability patterns

## Architecture Patterns

### Pattern 1: Basic Cloud Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     Application Layer                        │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐              │
│  │ API      │    │ Worker   │    │ Cron     │              │
│  │ Service  │    │ Service  │    │ Jobs     │              │
│  └────┬─────┘    └────┬─────┘    └────┬─────┘              │
│       │               │               │                     │
│       └───────────────┴───────────────┘                     │
│                       │                                      │
│              ┌────────┴────────┐                            │
│              │ Langfuse SDK    │                            │
│              │ (Singleton)     │                            │
│              └────────┬────────┘                            │
│                       │                                      │
└───────────────────────┼──────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────┐
│                  Langfuse Cloud                              │
│  ┌─────────────────────────────────────────────────────┐   │
│  │ Ingestion API → Processing → PostgreSQL → Dashboard  │   │
│  └─────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Pattern 2: Self-Hosted Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                        VPC                                   │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌─────────────────────────────────────────────────────┐   │
│  │                 Application Cluster                  │   │
│  │  ┌───────┐  ┌───────┐  ┌───────┐  ┌───────┐       │   │
│  │  │ Pod 1 │  │ Pod 2 │  │ Pod 3 │  │ Pod N │       │   │
│  │  └───┬───┘  └───┬───┘  └───┬───┘  └───┬───┘       │   │
│  │      └──────────┴─────┬────┴──────────┘            │   │
│  └─────────────────────────┼──────────────────────────┘   │
│                            │                               │
│                   ┌────────┴────────┐                     │
│                   │ Internal LB     │                     │
│                   └────────┬────────┘                     │
│                            │                               │
│  ┌─────────────────────────┼──────────────────────────┐   │
│  │         Langfuse Self-Hosted Cluster                │   │
│  │  ┌───────────┐    ┌───────────┐    ┌───────────┐  │   │
│  │  │ Langfuse  │    │ Langfuse  │    │ Langfuse  │  │   │
│  │  │ Instance 1│    │ Instance 2│    │ Instance 3│  │   │
│  │  └─────┬─────┘    └─────┬─────┘    └─────┬─────┘  │   │
│  │        └────────────────┼────────────────┘         │   │
│  │                         │                          │   │
│  │                ┌────────┴────────┐                 │   │
│  │                │ PostgreSQL RDS  │                 │   │
│  │                │ (Multi-AZ)      │                 │   │
│  │                └─────────────────┘                 │   │
│  └────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### Pattern 3: High-Scale Architecture with Buffer

```
┌─────────────────────────────────────────────────────────────┐
│                    Application Layer                         │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌────────────────────────────────────────────────────┐    │
│  │              Regional Application Clusters          │    │
│  │  ┌─────────┐  ┌─────────┐  ┌─────────┐            │    │
│  │  │ US-East │  │ EU-West │  │ AP-South│            │    │
│  │  └────┬────┘  └────┬────┘  └────┬────┘            │    │
│  │       └────────────┼────────────┘                  │    │
│  └────────────────────┼───────────────────────────────┘    │
│                       │                                     │
│              ┌────────┴────────┐                           │
│              │ Langfuse SDK    │                           │
│              │ (Batched)       │                           │
│              └────────┬────────┘                           │
│                       │                                     │
│              ┌────────┴────────┐                           │
│              │  Message Queue  │  ← Buffer for high volume │
│              │  (SQS/Kafka)    │                           │
│              └────────┬────────┘                           │
│                       │                                     │
│              ┌────────┴────────┐                           │
│              │ Ingestion       │  ← Async workers          │
│              │ Workers         │                           │
│              └────────┬────────┘                           │
│                       │                                     │
└───────────────────────┼─────────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────────┐
│                  Langfuse (Cloud/Self-Hosted)                │
└─────────────────────────────────────────────────────────────┘
```

## Instructions

### Step 1: Implement Singleton SDK Pattern

```typescript
// lib/langfuse/client.ts
import { Langfuse } from "langfuse";

class LangfuseClient {
  private static instance: Langfuse | null = null;
  private static shutdownRegistered = false;

  static getInstance(): Langfuse {
    if (!LangfuseClient.instance) {
      LangfuseClient.instance = new Langfuse({
        publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
        secretKey: process.env.LANGFUSE_SECRET_KEY!,
        baseUrl: process.env.LANGFUSE_HOST,
        // Production settings
        flushAt: parseInt(process.env.LANGFUSE_FLUSH_AT || "25"),
        flushInterval: parseInt(process.env.LANGFUSE_FLUSH_INTERVAL || "5000"),
        requestTimeout: 15000,
      });

      if (!LangfuseClient.shutdownRegistered) {
        LangfuseClient.registerShutdown();
      }
    }

    return LangfuseClient.instance;
  }

  private static registerShutdown() {
    const shutdown = async (signal: string) => {
      console.log(`${signal} received. Flushing Langfuse...`);
      if (LangfuseClient.instance) {
        await LangfuseClient.instance.shutdownAsync();
        LangfuseClient.instance = null;
      }
    };

    process.on("SIGTERM", () => shutdown("SIGTERM"));
    process.on("SIGINT", () => shutdown("SIGINT"));
    process.on("beforeExit", () => shutdown("beforeExit"));

    LangfuseClient.shutdownRegistered = true;
  }
}

export const langfuse = LangfuseClient.getInstance();
```

### Step 2: Implement Trace Context Propagation

```typescript
// lib/langfuse/context.ts
import { AsyncLocalStorage } from "async_hooks";

interface TraceContext {
  traceId: string;
  parentSpanId?: string;
  userId?: string;
  sessionId?: string;
}

const traceStorage = new AsyncLocalStorage<TraceContext>();

export function withTraceContext<T>(
  context: TraceContext,
  fn: () => T
): T {
  return traceStorage.run(context, fn);
}

export function getTraceContext(): TraceContext | undefined {
  return traceStorage.getStore();
}

// Middleware for Express
export function langfuseMiddleware() {
  return (req: Request, res: Response, next: NextFunction) => {
    const trace = langfuse.trace({
      name: `${req.method} ${req.path}`,
      userId: req.user?.id,
      sessionId: req.session?.id,
      metadata: {
        method: req.method,
        path: req.path,
        userAgent: req.headers["user-agent"],
      },
    });

    const context: TraceContext = {
      traceId: trace.id,
      userId: req.user?.id,
      sessionId: req.session?.id,
    };

    withTraceContext(context, () => {
      // Attach trace to request for easy access
      req.langfuseTrace = trace;

      // Finish trace on response
      res.on("finish", () => {
        trace.update({
          output: { statusCode: res.statusCode },
          level: res.statusCode >= 400 ? "ERROR" : undefined,
        });
      });

      next();
    });
  };
}
```

### Step 3: Implement Queue-Based Ingestion

```typescript
// lib/langfuse/queue.ts
import { SQSClient, SendMessageCommand } from "@aws-sdk/client-sqs";
import { Langfuse } from "langfuse";

interface QueuedTrace {
  name: string;
  input?: any;
  output?: any;
  metadata?: Record<string, any>;
  userId?: string;
  sessionId?: string;
  timestamp: string;
}

// Producer: Application sends to queue
class QueuedLangfuseProducer {
  private sqs: SQSClient;
  private queueUrl: string;

  constructor() {
    this.sqs = new SQSClient({});
    this.queueUrl = process.env.LANGFUSE_QUEUE_URL!;
  }

  async trace(params: Omit<QueuedTrace, "timestamp">) {
    const message: QueuedTrace = {
      ...params,
      timestamp: new Date().toISOString(),
    };

    await this.sqs.send(
      new SendMessageCommand({
        QueueUrl: this.queueUrl,
        MessageBody: JSON.stringify(message),
        MessageGroupId: params.sessionId || "default",
      })
    );
  }
}

// Consumer: Worker processes queue
class QueuedLangfuseConsumer {
  private langfuse: Langfuse;

  constructor() {
    this.langfuse = new Langfuse();
  }

  async processMessage(message: QueuedTrace) {
    const trace = this.langfuse.trace({
      name: message.name,
      input: message.input,
      output: message.output,
      metadata: message.metadata,
      userId: message.userId,
      sessionId: message.sessionId,
      timestamp: new Date(message.timestamp),
    });

    return trace.id;
  }

  async processBatch(messages: QueuedTrace[]) {
    for (const message of messages) {
      await this.processMessage(message);
    }
    await this.langfuse.flushAsync();
  }
}
```

### Step 4: Multi-Environment Configuration

```typescript
// config/langfuse.ts
type Environment = "development" | "staging" | "production";

interface LangfuseEnvironmentConfig {
  publicKey: string;
  secretKey: string;
  host: string;
  flushAt: number;
  flushInterval: number;
  enabled: boolean;
  sampling: {
    rate: number;
    alwaysSampleErrors: boolean;
  };
}

const ENVIRONMENT_CONFIGS: Record<Environment, LangfuseEnvironmentConfig> = {
  development: {
    publicKey: process.env.LANGFUSE_PUBLIC_KEY_DEV!,
    secretKey: process.env.LANGFUSE_SECRET_KEY_DEV!,
    host: process.env.LANGFUSE_HOST_DEV || "http://localhost:3000",
    flushAt: 1,
    flushInterval: 1000,
    enabled: true,
    sampling: { rate: 1.0, alwaysSampleErrors: true },
  },
  staging: {
    publicKey: process.env.LANGFUSE_PUBLIC_KEY_STAGING!,
    secretKey: process.env.LANGFUSE_SECRET_KEY_STAGING!,
    host: process.env.LANGFUSE_HOST_STAGING || "https://cloud.langfuse.com",
    flushAt: 15,
    flushInterval: 5000,
    enabled: true,
    sampling: { rate: 0.5, alwaysSampleErrors: true },
  },
  production: {
    publicKey: process.env.LANGFUSE_PUBLIC_KEY_PROD!,
    secretKey: process.env.LANGFUSE_SECRET_KEY_PROD!,
    host: process.env.LANGFUSE_HOST || "https://cloud.langfuse.com",
    flushAt: 25,
    flushInterval: 5000,
    enabled: true,
    sampling: { rate: 0.1, alwaysSampleErrors: true },
  },
};

export function getLangfuseConfig(): LangfuseEnvironmentConfig {
  const env = (process.env.NODE_ENV || "development") as Environment;
  return ENVIRONMENT_CONFIGS[env] || ENVIRONMENT_CONFIGS.development;
}
```

### Step 5: Implement Service Mesh Tracing

```typescript
// For microservices: propagate trace context across services
interface TraceHeaders {
  "x-langfuse-trace-id": string;
  "x-langfuse-parent-id"?: string;
  "x-langfuse-session-id"?: string;
}

// Outgoing request: inject headers
function injectTraceHeaders(headers: Headers) {
  const context = getTraceContext();
  if (context) {
    headers.set("x-langfuse-trace-id", context.traceId);
    if (context.parentSpanId) {
      headers.set("x-langfuse-parent-id", context.parentSpanId);
    }
    if (context.sessionId) {
      headers.set("x-langfuse-session-id", context.sessionId);
    }
  }
}

// Incoming request: extract headers and continue trace
function extractTraceContext(request: Request): TraceContext | null {
  const traceId = request.headers.get("x-langfuse-trace-id");
  if (!traceId) return null;

  return {
    traceId,
    parentSpanId: request.headers.get("x-langfuse-parent-id") || undefined,
    sessionId: request.headers.get("x-langfuse-session-id") || undefined,
  };
}

// Create linked trace in downstream service
function createLinkedTrace(parentContext: TraceContext, name: string) {
  return langfuse.trace({
    name,
    sessionId: parentContext.sessionId,
    metadata: {
      parentTraceId: parentContext.traceId,
      parentSpanId: parentContext.parentSpanId,
    },
  });
}
```

## Output
- Singleton SDK pattern with graceful shutdown
- Trace context propagation
- Queue-based async ingestion
- Multi-environment configuration
- Service mesh integration

## Architecture Decision Matrix

| Pattern | Use Case | Complexity | Scale |
|---------|----------|------------|-------|
| Basic Cloud | Small apps | Low | 100K traces/day |
| Self-Hosted | Data privacy | Medium | 1M traces/day |
| Queue-Based | High volume | High | 10M+ traces/day |

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Multiple instances | No singleton | Use singleton pattern |
| Lost traces | No shutdown | Register shutdown handlers |
| Cross-service gaps | No propagation | Implement header injection |
| Scale issues | Direct ingestion | Add message queue buffer |

## Resources
- [Langfuse Self-Hosting](https://langfuse.com/docs/deployment/self-host)
- [Langfuse Architecture](https://langfuse.com/docs)
- [OpenTelemetry Context](https://opentelemetry.io/docs/concepts/context-propagation/)

## Next Steps
For multi-environment setup, see `langfuse-multi-env-setup`.
