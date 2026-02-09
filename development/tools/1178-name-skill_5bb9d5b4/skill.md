---
name: openevidence-performance-tuning
description: |
  Optimize OpenEvidence clinical query performance and response times.
  Use when improving latency, optimizing query efficiency,
  or tuning caching for clinical AI applications.
  Trigger with phrases like "openevidence performance", "openevidence slow",
  "optimize openevidence", "openevidence latency", "speed up clinical queries".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Performance Tuning

## Overview
Optimize OpenEvidence clinical query performance for point-of-care response times.

## Prerequisites
- OpenEvidence integration running
- Monitoring configured (see `openevidence-observability`)
- Understanding of caching strategies
- Access to performance metrics

## Performance Targets

| Metric | Target | Critical |
|--------|--------|----------|
| Clinical Query P50 | < 3s | > 10s |
| Clinical Query P95 | < 8s | > 15s |
| Clinical Query P99 | < 15s | > 30s |
| DeepConsult Start | < 5s | > 15s |
| Cache Hit Rate | > 70% | < 50% |

## Instructions

### Step 1: Query Optimization
```typescript
// src/services/optimized-query.ts

// 1. Optimize question structure
function optimizeQuestion(question: string): string {
  // Remove unnecessary words
  let optimized = question
    .replace(/\b(please|can you|could you|I want to know)\b/gi, '')
    .replace(/\s+/g, ' ')
    .trim();

  // Ensure question ends with ?
  if (!optimized.endsWith('?')) {
    optimized += '?';
  }

  return optimized;
}

// 2. Optimize context to reduce processing
function optimizeContext(context: ClinicalContext): ClinicalContext {
  return {
    specialty: context.specialty,
    urgency: context.urgency,
    // Only include relevant patient context
    ...(context.patientAge && { patientAge: context.patientAge }),
    ...(context.patientSex && { patientSex: context.patientSex }),
    // Limit conditions to top 5 most relevant
    ...(context.relevantConditions && {
      relevantConditions: context.relevantConditions.slice(0, 5),
    }),
    // Limit medications to top 10
    ...(context.currentMedications && {
      currentMedications: context.currentMedications.slice(0, 10),
    }),
  };
}

// 3. Optimize query options
const OPTIMIZED_OPTIONS = {
  maxCitations: 5, // More citations = slower response
  includeGuidelines: true, // High value, low overhead
  includeDrugInfo: false, // Only when needed
  prioritizeRecent: true, // Faster index access
};
```

### Step 2: Intelligent Caching Layer
```typescript
// src/cache/clinical-cache.ts
import { createHash } from 'crypto';
import Redis from 'ioredis';

interface CacheConfig {
  redis: Redis;
  defaultTTL: number;
  maxEntries: number;
}

interface CachedResponse {
  data: any;
  cachedAt: number;
  ttl: number;
  hitCount: number;
}

export class ClinicalQueryCache {
  private redis: Redis;
  private defaultTTL: number;
  private prefix = 'oe:query:';

  constructor(config: CacheConfig) {
    this.redis = config.redis;
    this.defaultTTL = config.defaultTTL;
  }

  // Generate cache key from query
  private generateKey(question: string, context: ClinicalContext): string {
    const normalized = JSON.stringify({
      q: question.toLowerCase().trim(),
      s: context.specialty,
      u: context.urgency,
    });
    return this.prefix + createHash('sha256').update(normalized).digest('hex').substring(0, 16);
  }

  // Determine TTL based on query type
  private getTTL(question: string, context: ClinicalContext): number {
    const lowerQuestion = question.toLowerCase();

    // Short TTL for time-sensitive queries
    if (context.urgency === 'stat' || context.urgency === 'urgent') {
      return 300; // 5 minutes
    }

    // Longer TTL for stable medical facts
    if (lowerQuestion.includes('mechanism') ||
        lowerQuestion.includes('pharmacokinetics') ||
        lowerQuestion.includes('half-life')) {
      return 86400; // 24 hours
    }

    // Medium TTL for treatment guidelines (may update)
    if (lowerQuestion.includes('guideline') ||
        lowerQuestion.includes('first-line') ||
        lowerQuestion.includes('recommended')) {
      return 3600; // 1 hour
    }

    return this.defaultTTL;
  }

  async get(question: string, context: ClinicalContext): Promise<any | null> {
    const key = this.generateKey(question, context);
    const cached = await this.redis.get(key);

    if (!cached) return null;

    const entry: CachedResponse = JSON.parse(cached);

    // Increment hit count
    entry.hitCount++;
    await this.redis.setex(key, entry.ttl, JSON.stringify(entry));

    return entry.data;
  }

  async set(
    question: string,
    context: ClinicalContext,
    data: any
  ): Promise<void> {
    const key = this.generateKey(question, context);
    const ttl = this.getTTL(question, context);

    const entry: CachedResponse = {
      data,
      cachedAt: Date.now(),
      ttl,
      hitCount: 0,
    };

    await this.redis.setex(key, ttl, JSON.stringify(entry));
  }

  async invalidate(pattern?: string): Promise<number> {
    const keys = await this.redis.keys(this.prefix + (pattern || '*'));
    if (keys.length === 0) return 0;
    return this.redis.del(...keys);
  }

  async getStats(): Promise<{
    totalKeys: number;
    memoryUsage: string;
    hitRate: number;
  }> {
    const keys = await this.redis.keys(this.prefix + '*');
    const info = await this.redis.info('memory');

    let totalHits = 0;
    let totalEntries = 0;

    for (const key of keys.slice(0, 100)) {
      const cached = await this.redis.get(key);
      if (cached) {
        const entry: CachedResponse = JSON.parse(cached);
        totalHits += entry.hitCount;
        totalEntries++;
      }
    }

    return {
      totalKeys: keys.length,
      memoryUsage: info.match(/used_memory_human:(\S+)/)?.[1] || 'unknown',
      hitRate: totalEntries > 0 ? totalHits / totalEntries : 0,
    };
  }
}
```

### Step 3: Connection Pooling & Keep-Alive
```typescript
// src/openevidence/optimized-client.ts
import { OpenEvidenceClient } from '@openevidence/sdk';
import { Agent } from 'https';

// Reuse connections with keep-alive
const httpsAgent = new Agent({
  keepAlive: true,
  keepAliveMsecs: 30000,
  maxSockets: 10,
  maxFreeSockets: 5,
});

export const optimizedClient = new OpenEvidenceClient({
  apiKey: process.env.OPENEVIDENCE_API_KEY!,
  orgId: process.env.OPENEVIDENCE_ORG_ID!,
  httpAgent: httpsAgent,
  timeout: 30000,
});

// Pre-warm connections on startup
export async function warmupConnections(): Promise<void> {
  console.log('[Performance] Warming up OpenEvidence connections...');

  try {
    // Make a lightweight request to establish connection
    await optimizedClient.health.check();
    console.log('[Performance] Connection pool warmed up');
  } catch (error) {
    console.error('[Performance] Warmup failed:', error);
  }
}
```

### Step 4: Request Batching
```typescript
// src/services/batch-query.ts
import DataLoader from 'dataloader';

// Batch multiple clinical queries
const queryBatcher = new DataLoader<
  { question: string; context: ClinicalContext },
  ClinicalQueryResponse
>(
  async (queries) => {
    // OpenEvidence may support batch API in future
    // For now, parallelize within rate limits
    const results = await Promise.all(
      queries.map(q => optimizedClient.query({
        question: q.question,
        context: q.context,
      }))
    );
    return results;
  },
  {
    maxBatchSize: 5, // Respect rate limits
    batchScheduleFn: (callback) => setTimeout(callback, 50), // 50ms window
    cacheKeyFn: (query) => `${query.question}:${query.context.specialty}`,
  }
);

// Usage in service
export async function batchedQuery(
  question: string,
  context: ClinicalContext
): Promise<ClinicalQueryResponse> {
  return queryBatcher.load({ question, context });
}
```

### Step 5: Response Streaming (When Available)
```typescript
// src/services/streaming-query.ts
// Use streaming for faster perceived performance

export async function streamingClinicalQuery(
  question: string,
  context: ClinicalContext,
  onPartialResponse: (partial: string) => void
): Promise<ClinicalQueryResponse> {
  // Check if streaming is supported
  if (optimizedClient.supportsStreaming?.()) {
    const stream = await optimizedClient.query.stream({
      question,
      context,
    });

    let fullAnswer = '';

    for await (const chunk of stream) {
      fullAnswer += chunk.text;
      onPartialResponse(fullAnswer);
    }

    return stream.finalResponse;
  }

  // Fallback to regular query
  const response = await optimizedClient.query({ question, context });
  onPartialResponse(response.answer);
  return response;
}
```

### Step 6: Performance Monitoring
```typescript
// src/monitoring/performance-metrics.ts
import { Histogram, Counter, Gauge } from 'prom-client';

// Latency histogram
const queryLatency = new Histogram({
  name: 'openevidence_query_duration_seconds',
  help: 'Clinical query latency',
  labelNames: ['specialty', 'urgency', 'cached'],
  buckets: [0.5, 1, 2, 5, 10, 15, 30],
});

// Cache metrics
const cacheHits = new Counter({
  name: 'openevidence_cache_hits_total',
  help: 'Cache hit count',
});

const cacheMisses = new Counter({
  name: 'openevidence_cache_misses_total',
  help: 'Cache miss count',
});

// Current queue size
const queueSize = new Gauge({
  name: 'openevidence_queue_size',
  help: 'Current request queue size',
  labelNames: ['priority'],
});

// Instrumented query function
export async function instrumentedQuery(
  question: string,
  context: ClinicalContext,
  cache: ClinicalQueryCache
): Promise<ClinicalQueryResponse> {
  const timer = queryLatency.startTimer({
    specialty: context.specialty,
    urgency: context.urgency,
  });

  // Check cache first
  const cached = await cache.get(question, context);
  if (cached) {
    cacheHits.inc();
    timer({ cached: 'true' });
    return cached;
  }
  cacheMisses.inc();

  // Query API
  const response = await optimizedClient.query({ question, context });

  // Cache response
  await cache.set(question, context, response);

  timer({ cached: 'false' });
  return response;
}
```

## Performance Checklist
- [ ] Query optimization applied
- [ ] Caching layer configured
- [ ] Connection pooling enabled
- [ ] Keep-alive connections
- [ ] Request batching where applicable
- [ ] Monitoring metrics collecting
- [ ] Performance baselines established

## Output
- Optimized query construction
- Intelligent caching with TTL management
- Connection pooling for efficiency
- Performance metrics and monitoring

## Error Handling
| Performance Issue | Detection | Resolution |
|-------------------|-----------|------------|
| High P95 latency | Metrics alert | Check cache hit rate, optimize queries |
| Low cache hit rate | < 50% hits | Review TTL strategy, increase cache size |
| Connection timeouts | Timeout errors | Check keep-alive, increase pool size |
| Memory pressure | Redis alerts | Implement LRU eviction |

## Resources
- [Redis Caching Best Practices](https://redis.io/docs/manual/client-side-caching/)
- [DataLoader Documentation](https://github.com/graphql/dataloader)

## Next Steps
For cost optimization, see `openevidence-cost-tuning`.
