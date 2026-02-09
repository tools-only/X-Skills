---
name: openevidence-rate-limits
description: |
  Implement OpenEvidence rate limiting, backoff, and request optimization.
  Use when handling rate limit errors, implementing retry logic,
  or optimizing API request throughput for clinical queries.
  Trigger with phrases like "openevidence rate limit", "openevidence throttling",
  "openevidence 429", "openevidence retry", "openevidence backoff".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Rate Limits

## Overview
Handle OpenEvidence rate limits gracefully with exponential backoff, request queuing, and clinical priority management.

## Prerequisites
- OpenEvidence SDK installed
- Understanding of async/await patterns
- Access to rate limit headers

## Rate Limit Tiers

| Tier | Clinical Queries/min | DeepConsult/hour | Burst | Use Case |
|------|---------------------|------------------|-------|----------|
| Standard | 60 | 5 | 10 | Small practices |
| Professional | 300 | 20 | 50 | Clinics, small hospitals |
| Enterprise | 1,000 | 100 | 200 | Health systems |
| Unlimited | Custom | Custom | Custom | Large integrations |

## Rate Limit Headers

| Header | Description | Example |
|--------|-------------|---------|
| `X-RateLimit-Limit` | Max requests per window | `60` |
| `X-RateLimit-Remaining` | Requests remaining | `45` |
| `X-RateLimit-Reset` | Unix timestamp for reset | `1706295600` |
| `Retry-After` | Seconds to wait (on 429) | `30` |
| `X-DeepConsult-Remaining` | DeepConsult quota | `5` |

## Instructions

### Step 1: Implement Exponential Backoff with Jitter
```typescript
// src/openevidence/rate-limit.ts
interface RetryConfig {
  maxRetries: number;
  baseDelayMs: number;
  maxDelayMs: number;
  jitterMs: number;
}

const DEFAULT_CONFIG: RetryConfig = {
  maxRetries: 5,
  baseDelayMs: 1000,
  maxDelayMs: 60000,
  jitterMs: 500,
};

export async function withExponentialBackoff<T>(
  operation: () => Promise<T>,
  config: Partial<RetryConfig> = {}
): Promise<T> {
  const cfg = { ...DEFAULT_CONFIG, ...config };

  for (let attempt = 0; attempt <= cfg.maxRetries; attempt++) {
    try {
      return await operation();
    } catch (error: any) {
      if (attempt === cfg.maxRetries) throw error;

      const status = error.status || error.response?.status;
      if (status !== 429 && (status < 500 || status >= 600)) throw error;

      // Use Retry-After header if provided
      let delay: number;
      const retryAfter = error.response?.headers?.['retry-after'];
      if (retryAfter) {
        delay = parseInt(retryAfter) * 1000;
      } else {
        // Exponential backoff with jitter
        const exponentialDelay = cfg.baseDelayMs * Math.pow(2, attempt);
        const jitter = Math.random() * cfg.jitterMs;
        delay = Math.min(exponentialDelay + jitter, cfg.maxDelayMs);
      }

      console.log(`Rate limited. Attempt ${attempt + 1}/${cfg.maxRetries}. Retrying in ${delay}ms...`);
      await new Promise(r => setTimeout(r, delay));
    }
  }
  throw new Error('Unreachable');
}
```

### Step 2: Rate Limit Monitor
```typescript
// src/openevidence/rate-monitor.ts
export class RateLimitMonitor {
  private remaining: number = 60;
  private limit: number = 60;
  private resetAt: Date = new Date();
  private deepConsultRemaining: number = 5;

  updateFromHeaders(headers: Headers | Record<string, string>): void {
    const get = (key: string) =>
      headers instanceof Headers ? headers.get(key) : headers[key.toLowerCase()];

    const remaining = get('x-ratelimit-remaining');
    if (remaining) this.remaining = parseInt(remaining);

    const limit = get('x-ratelimit-limit');
    if (limit) this.limit = parseInt(limit);

    const resetTimestamp = get('x-ratelimit-reset');
    if (resetTimestamp) {
      this.resetAt = new Date(parseInt(resetTimestamp) * 1000);
    }

    const deepConsult = get('x-deepconsult-remaining');
    if (deepConsult) this.deepConsultRemaining = parseInt(deepConsult);
  }

  shouldThrottle(): boolean {
    return this.remaining < 5 && new Date() < this.resetAt;
  }

  getWaitTime(): number {
    return Math.max(0, this.resetAt.getTime() - Date.now());
  }

  canDeepConsult(): boolean {
    return this.deepConsultRemaining > 0;
  }

  getStatus(): RateLimitStatus {
    return {
      remaining: this.remaining,
      limit: this.limit,
      resetAt: this.resetAt.toISOString(),
      usagePercent: ((this.limit - this.remaining) / this.limit) * 100,
      deepConsultRemaining: this.deepConsultRemaining,
    };
  }
}

interface RateLimitStatus {
  remaining: number;
  limit: number;
  resetAt: string;
  usagePercent: number;
  deepConsultRemaining: number;
}
```

### Step 3: Priority-Based Request Queue
```typescript
// src/openevidence/request-queue.ts
import PQueue from 'p-queue';

type ClinicalPriority = 'stat' | 'urgent' | 'routine' | 'research';

interface QueuedRequest<T> {
  operation: () => Promise<T>;
  priority: ClinicalPriority;
  resolve: (value: T) => void;
  reject: (error: any) => void;
}

export class ClinicalRequestQueue {
  private queues: Record<ClinicalPriority, PQueue>;
  private monitor: RateLimitMonitor;

  constructor(monitor: RateLimitMonitor) {
    this.monitor = monitor;

    // Priority-based concurrency
    this.queues = {
      stat: new PQueue({ concurrency: 5, interval: 1000, intervalCap: 10 }),
      urgent: new PQueue({ concurrency: 3, interval: 1000, intervalCap: 6 }),
      routine: new PQueue({ concurrency: 2, interval: 1000, intervalCap: 4 }),
      research: new PQueue({ concurrency: 1, interval: 1000, intervalCap: 2 }),
    };
  }

  async enqueue<T>(
    operation: () => Promise<T>,
    priority: ClinicalPriority = 'routine'
  ): Promise<T> {
    const queue = this.queues[priority];

    return queue.add(async () => {
      // Check if we should throttle
      if (this.monitor.shouldThrottle()) {
        const waitTime = this.monitor.getWaitTime();
        console.log(`Throttling ${priority} request for ${waitTime}ms`);
        await new Promise(r => setTimeout(r, waitTime));
      }

      // Execute with retry logic
      return withExponentialBackoff(operation);
    });
  }

  getQueueStatus(): Record<ClinicalPriority, { pending: number; size: number }> {
    return {
      stat: { pending: this.queues.stat.pending, size: this.queues.stat.size },
      urgent: { pending: this.queues.urgent.pending, size: this.queues.urgent.size },
      routine: { pending: this.queues.routine.pending, size: this.queues.routine.size },
      research: { pending: this.queues.research.pending, size: this.queues.research.size },
    };
  }

  async pause(priority?: ClinicalPriority): Promise<void> {
    if (priority) {
      this.queues[priority].pause();
    } else {
      Object.values(this.queues).forEach(q => q.pause());
    }
  }

  async resume(priority?: ClinicalPriority): Promise<void> {
    if (priority) {
      this.queues[priority].start();
    } else {
      Object.values(this.queues).forEach(q => q.start());
    }
  }
}
```

### Step 4: Adaptive Rate Limiting
```typescript
// src/openevidence/adaptive-limiter.ts
export class AdaptiveRateLimiter {
  private windowMs = 60000; // 1 minute window
  private requestTimes: number[] = [];
  private targetUsagePercent = 80; // Stay below 80% of limit
  private currentLimit = 60;

  recordRequest(): void {
    const now = Date.now();
    this.requestTimes.push(now);

    // Clean old entries
    this.requestTimes = this.requestTimes.filter(
      t => t > now - this.windowMs
    );
  }

  updateLimit(limit: number): void {
    this.currentLimit = limit;
  }

  shouldDelay(): { delay: boolean; waitMs: number } {
    const requestsInWindow = this.requestTimes.length;
    const targetMax = Math.floor(this.currentLimit * (this.targetUsagePercent / 100));

    if (requestsInWindow >= targetMax) {
      // Calculate time until oldest request exits window
      const oldestRequest = Math.min(...this.requestTimes);
      const waitMs = Math.max(0, oldestRequest + this.windowMs - Date.now());

      return { delay: true, waitMs };
    }

    return { delay: false, waitMs: 0 };
  }

  getMetrics(): {
    requestsInWindow: number;
    usagePercent: number;
    headroom: number;
  } {
    const requestsInWindow = this.requestTimes.length;
    return {
      requestsInWindow,
      usagePercent: (requestsInWindow / this.currentLimit) * 100,
      headroom: this.currentLimit - requestsInWindow,
    };
  }
}
```

## Output
- Reliable API calls with automatic retry
- Priority-based request queuing
- Rate limit monitoring and alerting
- Adaptive throttling based on usage

## Error Handling
| Header | Description | Action |
|--------|-------------|--------|
| `X-RateLimit-Limit` | Max requests | Adjust queue concurrency |
| `X-RateLimit-Remaining` | Remaining requests | Throttle if low |
| `X-RateLimit-Reset` | Reset timestamp | Wait until reset |
| `Retry-After` | Seconds to wait | Honor exactly |

## Examples

### Complete Rate-Limited Client
```typescript
// src/services/rate-limited-client.ts
import { OpenEvidenceClient } from '@openevidence/sdk';
import { RateLimitMonitor } from '../openevidence/rate-monitor';
import { ClinicalRequestQueue } from '../openevidence/request-queue';

const monitor = new RateLimitMonitor();
const queue = new ClinicalRequestQueue(monitor);

const rawClient = new OpenEvidenceClient({
  apiKey: process.env.OPENEVIDENCE_API_KEY,
  orgId: process.env.OPENEVIDENCE_ORG_ID,
});

export async function clinicalQuery(
  question: string,
  priority: 'stat' | 'urgent' | 'routine' = 'routine'
) {
  return queue.enqueue(
    async () => {
      const response = await rawClient.query({
        question,
        context: { specialty: 'internal-medicine', urgency: priority },
      });

      // Update monitor from response headers
      monitor.updateFromHeaders(response.headers);

      return response.data;
    },
    priority
  );
}

// Usage
await clinicalQuery('What is the dose of amoxicillin for sinusitis?', 'routine');
await clinicalQuery('Contraindications for tPA in stroke patient?', 'stat');
```

### Monitor Dashboard
```typescript
// Expose rate limit metrics for monitoring
app.get('/health/openevidence', (req, res) => {
  const status = monitor.getStatus();
  const queueStatus = queue.getQueueStatus();

  res.json({
    rateLimits: status,
    queues: queueStatus,
    healthy: status.remaining > 5,
  });
});
```

## Resources
- [OpenEvidence API Terms](https://www.openevidence.com/policies/api)
- [p-queue Documentation](https://github.com/sindresorhus/p-queue)

## Next Steps
For security configuration, see `openevidence-security-basics`.
