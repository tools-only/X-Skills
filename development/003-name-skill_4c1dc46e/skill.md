---
name: lokalise-performance-tuning
description: |
  Optimize Lokalise API performance with caching, pagination, and bulk operations.
  Use when experiencing slow API responses, implementing caching strategies,
  or optimizing request throughput for Lokalise integrations.
  Trigger with phrases like "lokalise performance", "optimize lokalise",
  "lokalise latency", "lokalise caching", "lokalise slow", "lokalise batch".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Lokalise Performance Tuning

## Overview
Optimize Lokalise API performance with caching, pagination, and bulk operations.

## Prerequisites
- Lokalise SDK installed
- Understanding of async patterns
- Redis or in-memory cache available (optional)
- Performance monitoring in place

## Latency Benchmarks

| Operation | Typical P50 | Typical P95 | Max Items |
|-----------|-------------|-------------|-----------|
| List projects | 100ms | 300ms | 100 |
| List keys | 150ms | 500ms | 500 |
| Create key | 200ms | 600ms | 1 |
| Bulk create keys | 300ms | 1000ms | 500 |
| Download files | 500ms | 2000ms | All |
| Upload file | 1000ms | 5000ms | 1 |

## Instructions

### Step 1: Enable Compression
```typescript
const client = new LokaliseApi({
  apiKey: process.env.LOKALISE_API_TOKEN!,
  enableCompression: true,  // Enable gzip for large responses
});
```

### Step 2: Implement Response Caching
```typescript
import { LRUCache } from "lru-cache";

const cache = new LRUCache<string, any>({
  max: 1000,
  ttl: 60000,  // 1 minute default TTL
  updateAgeOnGet: true,
});

async function cachedRequest<T>(
  key: string,
  fetcher: () => Promise<T>,
  ttl?: number
): Promise<T> {
  const cached = cache.get(key);
  if (cached !== undefined) {
    return cached as T;
  }

  const result = await fetcher();
  cache.set(key, result, { ttl });
  return result;
}

// Usage
async function getProject(projectId: string) {
  return cachedRequest(
    `project:${projectId}`,
    () => client.projects().get(projectId),
    300000  // 5 minutes
  );
}

async function getKeys(projectId: string) {
  return cachedRequest(
    `keys:${projectId}`,
    () => client.keys().list({ project_id: projectId, limit: 500 }),
    60000  // 1 minute
  );
}
```

### Step 3: Use Cursor Pagination
```typescript
// Cursor pagination is more efficient for large datasets
async function* iterateAllKeys(projectId: string) {
  let cursor: string | undefined;

  do {
    const result = await client.keys().list({
      project_id: projectId,
      limit: 500,  // Max allowed
      pagination: "cursor",
      cursor,
    });

    for (const key of result.items) {
      yield key;
    }

    cursor = result.hasNextCursor() ? result.nextCursor : undefined;
  } while (cursor);
}

// Usage
async function getAllKeys(projectId: string) {
  const keys = [];
  for await (const key of iterateAllKeys(projectId)) {
    keys.push(key);
  }
  return keys;
}
```

### Step 4: Batch Operations
```typescript
// Batch creates are much faster than individual creates
async function createKeysBatched(
  projectId: string,
  keys: any[],
  batchSize = 100
): Promise<any[]> {
  const results: any[] = [];

  for (let i = 0; i < keys.length; i += batchSize) {
    const batch = keys.slice(i, i + batchSize);

    const result = await client.keys().create({
      project_id: projectId,
      keys: batch,
    });

    results.push(...result.items);

    // Small delay to respect rate limits
    await new Promise(r => setTimeout(r, 200));
  }

  return results;
}

// DataLoader for automatic batching
import DataLoader from "dataloader";

const keyLoader = new DataLoader<string, any>(
  async (keyIds) => {
    // Batch fetch keys
    const result = await client.keys().list({
      project_id: projectId,
      filter_key_ids: keyIds.join(","),
    });

    // Return in same order as requested
    return keyIds.map(id =>
      result.items.find(k => k.key_id.toString() === id) || null
    );
  },
  {
    maxBatchSize: 100,
    batchScheduleFn: callback => setTimeout(callback, 10),
  }
);
```

### Step 5: Parallel Downloads with Rate Limiting
```typescript
import PQueue from "p-queue";

const queue = new PQueue({
  concurrency: 5,
  interval: 1000,
  intervalCap: 5,
});

async function downloadMultipleProjects(projectIds: string[]) {
  return Promise.all(
    projectIds.map(id =>
      queue.add(() =>
        client.files().download(id, {
          format: "json",
          original_filenames: false,
        })
      )
    )
  );
}
```

## Output
- Compression enabled for faster transfers
- Response caching implemented
- Cursor pagination for large datasets
- Batch operations for bulk changes

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Cache miss storm | TTL expired | Use stale-while-revalidate |
| Memory pressure | Cache too large | Set max entries, use Redis |
| Slow pagination | Offset pagination | Switch to cursor pagination |
| Rate limit hit | Too many parallel requests | Use request queue |

## Examples

### Redis Caching (Distributed)
```typescript
import Redis from "ioredis";

const redis = new Redis(process.env.REDIS_URL);

async function cachedWithRedis<T>(
  key: string,
  fetcher: () => Promise<T>,
  ttlSeconds = 60
): Promise<T> {
  const cached = await redis.get(key);
  if (cached) {
    return JSON.parse(cached);
  }

  const result = await fetcher();
  await redis.setex(key, ttlSeconds, JSON.stringify(result));
  return result;
}

// Stale-while-revalidate pattern
async function staleWhileRevalidate<T>(
  key: string,
  fetcher: () => Promise<T>,
  staleTtl = 60,
  maxTtl = 3600
): Promise<T> {
  const cached = await redis.get(key);

  if (cached) {
    const { data, timestamp } = JSON.parse(cached);
    const age = (Date.now() - timestamp) / 1000;

    if (age < staleTtl) {
      return data;  // Fresh
    }

    if (age < maxTtl) {
      // Stale but usable - revalidate in background
      revalidate(key, fetcher, staleTtl, maxTtl);
      return data;
    }
  }

  // Cache miss or expired - fetch fresh
  return revalidate(key, fetcher, staleTtl, maxTtl);
}

async function revalidate<T>(
  key: string,
  fetcher: () => Promise<T>,
  staleTtl: number,
  maxTtl: number
): Promise<T> {
  const data = await fetcher();
  await redis.setex(key, maxTtl, JSON.stringify({
    data,
    timestamp: Date.now(),
  }));
  return data;
}
```

### Performance Monitoring
```typescript
async function measuredRequest<T>(
  operation: string,
  fn: () => Promise<T>
): Promise<T> {
  const start = performance.now();

  try {
    const result = await fn();
    const duration = performance.now() - start;

    console.log({
      operation,
      duration: `${duration.toFixed(2)}ms`,
      status: "success",
    });

    // Report to metrics
    await reportMetric("lokalise_request_duration", duration, {
      operation,
      status: "success",
    });

    return result;
  } catch (error: any) {
    const duration = performance.now() - start;

    console.error({
      operation,
      duration: `${duration.toFixed(2)}ms`,
      status: "error",
      error: error.message,
    });

    await reportMetric("lokalise_request_duration", duration, {
      operation,
      status: "error",
    });

    throw error;
  }
}
```

### Preloading Translations
```typescript
// Preload translations during app startup
async function preloadTranslations(
  projectId: string,
  locales: string[]
): Promise<Map<string, any>> {
  const translations = new Map();

  // Download all in parallel
  const downloads = await Promise.all(
    locales.map(async (locale) => {
      const result = await cachedWithRedis(
        `translations:${projectId}:${locale}`,
        () => fetchTranslationsForLocale(projectId, locale),
        3600  // 1 hour cache
      );
      return { locale, translations: result };
    })
  );

  for (const { locale, translations: trans } of downloads) {
    translations.set(locale, trans);
  }

  return translations;
}
```

## Resources
- [Lokalise API Rate Limits](https://developers.lokalise.com/reference/api-rate-limits)
- [Node SDK Documentation](https://lokalise.github.io/node-lokalise-api/)
- [DataLoader Documentation](https://github.com/graphql/dataloader)

## Next Steps
For cost optimization, see `lokalise-cost-tuning`.
