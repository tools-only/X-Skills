---
name: mistral-rate-limits
description: |
  Implement Mistral AI rate limiting, backoff, and request management.
  Use when handling rate limit errors, implementing retry logic,
  or optimizing API request throughput for Mistral AI.
  Trigger with phrases like "mistral rate limit", "mistral throttling",
  "mistral 429", "mistral retry", "mistral backoff".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Mistral AI Rate Limits

## Overview
Handle Mistral AI rate limits gracefully with exponential backoff and request management.

## Prerequisites
- Mistral AI SDK installed
- Understanding of async/await patterns
- Access to rate limit headers

## Instructions

### Step 1: Understand Rate Limit Tiers

| Tier | Requests/min | Tokens/min | Tokens/month |
|------|-------------|------------|--------------|
| Free | 2 | 500K | 1B |
| Production | 120 | 1M | 10B |
| Enterprise | Custom | Custom | Custom |

**Note:** Limits vary by model and are subject to change. Check console.mistral.ai for current limits.

### Step 2: Implement Exponential Backoff with Jitter

```typescript
interface RetryConfig {
  maxRetries: number;
  baseDelayMs: number;
  maxDelayMs: number;
  jitterMs: number;
}

async function withExponentialBackoff<T>(
  operation: () => Promise<T>,
  config: RetryConfig = {
    maxRetries: 5,
    baseDelayMs: 1000,
    maxDelayMs: 60000,
    jitterMs: 500
  }
): Promise<T> {
  for (let attempt = 0; attempt <= config.maxRetries; attempt++) {
    try {
      return await operation();
    } catch (error: any) {
      if (attempt === config.maxRetries) throw error;

      const status = error.status;

      // Only retry on rate limits (429) or server errors (5xx)
      if (status !== 429 && (status < 500 || status >= 600)) throw error;

      // Check for Retry-After header
      const retryAfter = error.headers?.['retry-after'];
      let delay: number;

      if (retryAfter) {
        delay = parseInt(retryAfter) * 1000;
      } else {
        // Exponential delay with jitter to prevent thundering herd
        const exponentialDelay = config.baseDelayMs * Math.pow(2, attempt);
        const jitter = Math.random() * config.jitterMs;
        delay = Math.min(exponentialDelay + jitter, config.maxDelayMs);
      }

      console.log(`Attempt ${attempt + 1} failed (${status}). Retrying in ${delay.toFixed(0)}ms...`);
      await new Promise(r => setTimeout(r, delay));
    }
  }
  throw new Error('Unreachable');
}

// Usage
const response = await withExponentialBackoff(() =>
  client.chat.complete({
    model: 'mistral-small-latest',
    messages: [{ role: 'user', content: 'Hello!' }],
  })
);
```

### Step 3: Token-Based Rate Limiting

```typescript
class TokenRateLimiter {
  private tokensUsed = 0;
  private windowStart = Date.now();
  private readonly tokensPerMinute: number;
  private readonly windowMs = 60000; // 1 minute

  constructor(tokensPerMinute = 500000) {
    this.tokensPerMinute = tokensPerMinute;
  }

  async waitForCapacity(estimatedTokens: number): Promise<void> {
    const now = Date.now();
    const elapsed = now - this.windowStart;

    // Reset window if needed
    if (elapsed >= this.windowMs) {
      this.tokensUsed = 0;
      this.windowStart = now;
    }

    // Check if we need to wait
    if (this.tokensUsed + estimatedTokens > this.tokensPerMinute) {
      const waitTime = this.windowMs - elapsed;
      console.log(`Token limit approaching. Waiting ${waitTime}ms...`);
      await new Promise(r => setTimeout(r, waitTime));
      this.tokensUsed = 0;
      this.windowStart = Date.now();
    }
  }

  recordUsage(tokensUsed: number): void {
    this.tokensUsed += tokensUsed;
  }
}

// Usage
const rateLimiter = new TokenRateLimiter(500000);

async function rateLimitedChat(messages: Message[]): Promise<string> {
  // Estimate tokens (rough: 4 chars per token)
  const estimatedTokens = JSON.stringify(messages).length / 4;

  await rateLimiter.waitForCapacity(estimatedTokens + 500); // +500 for response

  const response = await client.chat.complete({
    model: 'mistral-small-latest',
    messages,
  });

  if (response.usage) {
    rateLimiter.recordUsage(response.usage.totalTokens || 0);
  }

  return response.choices?.[0]?.message?.content ?? '';
}
```

### Step 4: Request Queue with Concurrency Control

```typescript
import PQueue from 'p-queue';

const requestQueue = new PQueue({
  concurrency: 5,        // Max concurrent requests
  interval: 1000,        // 1 second interval
  intervalCap: 10,       // Max 10 requests per interval
});

async function queuedRequest<T>(operation: () => Promise<T>): Promise<T> {
  return requestQueue.add(async () => {
    return withExponentialBackoff(operation);
  });
}

// Usage
const results = await Promise.all(
  prompts.map(prompt =>
    queuedRequest(() =>
      client.chat.complete({
        model: 'mistral-small-latest',
        messages: [{ role: 'user', content: prompt }],
      })
    )
  )
);
```

### Step 5: Rate Limit Monitor

```typescript
class RateLimitMonitor {
  private requestCount = 0;
  private lastReset = Date.now();
  private readonly alertThreshold: number;

  constructor(alertThreshold = 0.8) {
    this.alertThreshold = alertThreshold;
  }

  recordRequest(): void {
    const now = Date.now();
    if (now - this.lastReset >= 60000) {
      this.requestCount = 0;
      this.lastReset = now;
    }
    this.requestCount++;
  }

  checkThreshold(maxRequests: number): void {
    if (this.requestCount / maxRequests > this.alertThreshold) {
      console.warn(`Rate limit warning: ${this.requestCount}/${maxRequests} requests used`);
    }
  }

  getStats(): { requestCount: number; windowRemaining: number } {
    return {
      requestCount: this.requestCount,
      windowRemaining: 60000 - (Date.now() - this.lastReset),
    };
  }
}
```

## Output
- Reliable API calls with automatic retry
- Token-based rate limiting
- Request queue with concurrency control
- Rate limit monitoring

## Error Handling
| Header | Description | Action |
|--------|-------------|--------|
| Retry-After | Seconds to wait | Honor this value |
| X-RateLimit-Limit | Max requests | Monitor usage |
| X-RateLimit-Remaining | Remaining requests | Throttle if low |
| X-RateLimit-Reset | Reset timestamp | Wait until reset |

## Examples

### Python Rate Limiting
```python
import time
import asyncio
from mistralai import Mistral

async def with_retry(client, max_retries=3):
    for attempt in range(max_retries):
        try:
            return await client.chat.complete(
                model="mistral-small-latest",
                messages=[{"role": "user", "content": "Hello"}]
            )
        except Exception as e:
            if hasattr(e, 'status') and e.status == 429:
                delay = (2 ** attempt) + (random.random() * 0.5)
                print(f"Rate limited. Waiting {delay:.1f}s...")
                await asyncio.sleep(delay)
            else:
                raise
    raise Exception("Max retries exceeded")
```

### Batch Processing with Rate Limiting
```typescript
async function processBatch<T, R>(
  items: T[],
  processor: (item: T) => Promise<R>,
  batchSize = 5,
  delayMs = 1000
): Promise<R[]> {
  const results: R[] = [];

  for (let i = 0; i < items.length; i += batchSize) {
    const batch = items.slice(i, i + batchSize);

    const batchResults = await Promise.all(
      batch.map(item => withExponentialBackoff(() => processor(item)))
    );

    results.push(...batchResults);

    // Delay between batches
    if (i + batchSize < items.length) {
      await new Promise(r => setTimeout(r, delayMs));
    }
  }

  return results;
}
```

## Resources
- [Mistral AI Rate Limits](https://docs.mistral.ai/api/#rate-limits)
- [Mistral AI Console](https://console.mistral.ai/)
- [p-queue Documentation](https://github.com/sindresorhus/p-queue)

## Next Steps
For security configuration, see `mistral-security-basics`.
