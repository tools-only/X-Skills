---
name: mistral-performance-tuning
description: |
  Optimize Mistral AI performance with caching, batching, and latency reduction.
  Use when experiencing slow API responses, implementing caching strategies,
  or optimizing request throughput for Mistral AI integrations.
  Trigger with phrases like "mistral performance", "optimize mistral",
  "mistral latency", "mistral caching", "mistral slow", "mistral batch".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Mistral AI Performance Tuning

## Overview
Optimize Mistral AI API performance with caching, batching, and latency reduction techniques.

## Prerequisites
- Mistral AI SDK installed
- Understanding of async patterns
- Redis or in-memory cache available (optional)
- Performance monitoring in place

## Latency Benchmarks

| Model | P50 | P95 | P99 | Use Case |
|-------|-----|-----|-----|----------|
| mistral-small-latest | 200ms | 500ms | 1s | Fast responses |
| mistral-large-latest | 500ms | 1.5s | 3s | Complex reasoning |
| mistral-embed | 50ms | 150ms | 300ms | Embeddings |

## Instructions

### Step 1: Response Caching

```typescript
import { LRUCache } from 'lru-cache';
import crypto from 'crypto';

const cache = new LRUCache<string, any>({
  max: 1000,
  ttl: 5 * 60 * 1000, // 5 minutes
  updateAgeOnGet: true,
});

function getCacheKey(messages: any[], model: string, options?: any): string {
  const data = JSON.stringify({ messages, model, options });
  return crypto.createHash('sha256').update(data).digest('hex');
}

async function cachedChat(
  client: Mistral,
  messages: any[],
  model: string,
  options?: { temperature?: number; maxTokens?: number }
): Promise<string> {
  // Only cache deterministic requests (temperature = 0)
  const isCacheable = (options?.temperature ?? 0.7) === 0;

  if (isCacheable) {
    const key = getCacheKey(messages, model, options);
    const cached = cache.get(key);
    if (cached) {
      console.log('Cache hit');
      return cached;
    }
  }

  const response = await client.chat.complete({
    model,
    messages,
    ...options,
  });

  const content = response.choices?.[0]?.message?.content ?? '';

  if (isCacheable) {
    const key = getCacheKey(messages, model, options);
    cache.set(key, content);
  }

  return content;
}
```

### Step 2: Redis Distributed Caching

```typescript
import Redis from 'ioredis';
import crypto from 'crypto';

const redis = new Redis(process.env.REDIS_URL);

async function cachedWithRedis<T>(
  key: string,
  fetcher: () => Promise<T>,
  ttlSeconds = 300
): Promise<T> {
  const cached = await redis.get(key);
  if (cached) {
    return JSON.parse(cached);
  }

  const result = await fetcher();
  await redis.setex(key, ttlSeconds, JSON.stringify(result));
  return result;
}

// Semantic cache for similar queries
async function semanticCache(
  client: Mistral,
  query: string,
  threshold = 0.95
): Promise<string | null> {
  // Get embedding for query
  const queryEmbed = await client.embeddings.create({
    model: 'mistral-embed',
    inputs: [query],
  });
  const queryVector = queryEmbed.data[0].embedding;

  // Check cache for similar queries
  const cachedQueries = await redis.keys('semantic:*');

  for (const key of cachedQueries) {
    const cached = JSON.parse(await redis.get(key) || '{}');
    const similarity = cosineSimilarity(queryVector, cached.embedding);

    if (similarity >= threshold) {
      console.log(`Semantic cache hit (similarity: ${similarity.toFixed(3)})`);
      return cached.response;
    }
  }

  return null;
}
```

### Step 3: Request Batching

```typescript
import DataLoader from 'dataloader';

// Batch embedding requests
const embeddingLoader = new DataLoader<string, number[]>(
  async (texts) => {
    const response = await client.embeddings.create({
      model: 'mistral-embed',
      inputs: texts as string[],
    });
    return response.data.map(d => d.embedding);
  },
  {
    maxBatchSize: 100, // Mistral limit
    batchScheduleFn: callback => setTimeout(callback, 10), // 10ms window
  }
);

// Usage - automatically batched
const [embed1, embed2, embed3] = await Promise.all([
  embeddingLoader.load('Text 1'),
  embeddingLoader.load('Text 2'),
  embeddingLoader.load('Text 3'),
]);
```

### Step 4: Connection Optimization

```typescript
import { Agent } from 'https';
import Mistral from '@mistralai/mistralai';

// Keep-alive connection pooling
const agent = new Agent({
  keepAlive: true,
  maxSockets: 10,
  maxFreeSockets: 5,
  timeout: 60000,
});

// Note: Check if Mistral client supports custom agents
// If not, connection pooling happens at the HTTP level
```

### Step 5: Streaming for Perceived Performance

```typescript
// Streaming reduces Time to First Token (TTFT)
async function* streamWithMetrics(
  client: Mistral,
  messages: any[],
  model: string
): AsyncGenerator<{ content: string; metrics: any }> {
  const startTime = Date.now();
  let firstTokenTime: number | null = null;
  let tokenCount = 0;

  const stream = await client.chat.stream({ model, messages });

  for await (const event of stream) {
    const content = event.data?.choices?.[0]?.delta?.content;
    if (content) {
      if (!firstTokenTime) {
        firstTokenTime = Date.now();
      }
      tokenCount++;
      yield {
        content,
        metrics: {
          ttft: firstTokenTime - startTime,
          tokensPerSecond: tokenCount / ((Date.now() - startTime) / 1000),
        },
      };
    }
  }
}

// Usage
let fullResponse = '';
for await (const { content, metrics } of streamWithMetrics(client, messages, 'mistral-small-latest')) {
  fullResponse += content;
  process.stdout.write(content);
}
console.log(`\nTTFT: ${metrics.ttft}ms, Speed: ${metrics.tokensPerSecond.toFixed(1)} tok/s`);
```

### Step 6: Model Selection for Speed

```typescript
type SpeedTier = 'fastest' | 'balanced' | 'quality';

function selectModelForSpeed(tier: SpeedTier, taskComplexity: 'low' | 'medium' | 'high'): string {
  const matrix = {
    fastest: {
      low: 'mistral-small-latest',
      medium: 'mistral-small-latest',
      high: 'mistral-small-latest',
    },
    balanced: {
      low: 'mistral-small-latest',
      medium: 'mistral-small-latest',
      high: 'mistral-large-latest',
    },
    quality: {
      low: 'mistral-small-latest',
      medium: 'mistral-large-latest',
      high: 'mistral-large-latest',
    },
  };

  return matrix[tier][taskComplexity];
}
```

### Step 7: Performance Monitoring

```typescript
interface PerformanceMetrics {
  model: string;
  latencyMs: number;
  ttftMs?: number;
  tokensPerSecond?: number;
  inputTokens: number;
  outputTokens: number;
  cached: boolean;
}

async function measurePerformance(
  operation: () => Promise<any>,
  metadata: Partial<PerformanceMetrics>
): Promise<{ result: any; metrics: PerformanceMetrics }> {
  const start = Date.now();

  const result = await operation();

  const metrics: PerformanceMetrics = {
    model: metadata.model || 'unknown',
    latencyMs: Date.now() - start,
    inputTokens: result.usage?.promptTokens || 0,
    outputTokens: result.usage?.completionTokens || 0,
    cached: metadata.cached || false,
    ...metadata,
  };

  // Log to monitoring system
  console.log('[PERF]', JSON.stringify(metrics));

  return { result, metrics };
}

// Usage
const { result, metrics } = await measurePerformance(
  () => client.chat.complete({ model, messages }),
  { model, cached: false }
);
```

## Output
- Reduced API latency
- Caching layer implemented
- Request batching enabled
- Performance monitoring active

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Cache miss storm | TTL expired | Use stale-while-revalidate |
| Batch timeout | Too many items | Reduce batch size |
| Memory pressure | Cache too large | Set max cache entries |
| Slow TTFT | Large prompts | Reduce prompt size or use smaller model |

## Examples

### Quick Performance Wrapper
```typescript
const withPerformance = async <T>(
  name: string,
  fn: () => Promise<T>
): Promise<T> => {
  const start = Date.now();
  const result = await fn();
  console.log(`[${name}] ${Date.now() - start}ms`);
  return result;
};

// Usage
const response = await withPerformance('chat', () =>
  client.chat.complete({ model, messages })
);
```

### Parallel Requests with Concurrency Limit
```typescript
import pLimit from 'p-limit';

const limit = pLimit(5); // Max 5 concurrent requests

const results = await Promise.all(
  prompts.map(prompt =>
    limit(() => client.chat.complete({
      model: 'mistral-small-latest',
      messages: [{ role: 'user', content: prompt }],
    }))
  )
);
```

## Resources
- [Mistral AI Models](https://docs.mistral.ai/getting-started/models/)
- [LRU Cache Documentation](https://github.com/isaacs/node-lru-cache)
- [DataLoader Documentation](https://github.com/graphql/dataloader)

## Next Steps
For cost optimization, see `mistral-cost-tuning`.
