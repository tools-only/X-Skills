---
name: speak-rate-limits
description: |
  Implement Speak rate limiting, backoff, and idempotency patterns for language learning APIs.
  Use when handling rate limit errors, implementing retry logic,
  or optimizing API request throughput for Speak integrations.
  Trigger with phrases like "speak rate limit", "speak throttling",
  "speak 429", "speak retry", "speak backoff".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Speak Rate Limits

## Overview
Handle Speak rate limits gracefully with exponential backoff and session management for language learning applications.

## Prerequisites
- Speak SDK installed
- Understanding of async/await patterns
- Access to rate limit headers

## Instructions

### Step 1: Understand Rate Limit Tiers

| Tier | Sessions/hour | API Calls/min | Audio Processing/min | Burst |
|------|---------------|---------------|---------------------|-------|
| Free | 10 | 60 | 20 | 5 |
| Premium | 100 | 300 | 100 | 25 |
| Education | 500 | 1000 | 300 | 50 |
| Enterprise | Unlimited | Custom | Custom | Custom |

### Step 2: Implement Exponential Backoff with Jitter

```typescript
interface BackoffConfig {
  maxRetries: number;
  baseDelayMs: number;
  maxDelayMs: number;
  jitterMs: number;
}

async function withExponentialBackoff<T>(
  operation: () => Promise<T>,
  config: BackoffConfig = {
    maxRetries: 5,
    baseDelayMs: 1000,
    maxDelayMs: 32000,
    jitterMs: 500,
  }
): Promise<T> {
  for (let attempt = 0; attempt <= config.maxRetries; attempt++) {
    try {
      return await operation();
    } catch (error: any) {
      if (attempt === config.maxRetries) throw error;

      const status = error.status || error.response?.status;
      const isRetryable = status === 429 || (status >= 500 && status < 600);

      if (!isRetryable) throw error;

      // Check for Retry-After header
      const retryAfter = error.headers?.['retry-after'];
      let delay: number;

      if (retryAfter) {
        delay = parseInt(retryAfter) * 1000;
      } else {
        // Exponential backoff with jitter
        const exponentialDelay = config.baseDelayMs * Math.pow(2, attempt);
        const jitter = Math.random() * config.jitterMs;
        delay = Math.min(exponentialDelay + jitter, config.maxDelayMs);
      }

      console.log(`Rate limited. Retrying in ${delay.toFixed(0)}ms...`);
      await new Promise(r => setTimeout(r, delay));
    }
  }
  throw new Error('Unreachable');
}
```

### Step 3: Session-Aware Rate Limiting

```typescript
// Speak uses session-based rate limits for lesson interactions
class SessionRateLimiter {
  private sessionCounts: Map<string, { count: number; windowStart: Date }> = new Map();
  private readonly maxSessionsPerHour: number;
  private readonly maxResponsesPerSession: number;

  constructor(tier: 'free' | 'premium' | 'education' = 'premium') {
    const limits = {
      free: { sessions: 10, responses: 50 },
      premium: { sessions: 100, responses: 200 },
      education: { sessions: 500, responses: 500 },
    };
    this.maxSessionsPerHour = limits[tier].sessions;
    this.maxResponsesPerSession = limits[tier].responses;
  }

  canStartSession(userId: string): boolean {
    const now = new Date();
    const windowMs = 60 * 60 * 1000; // 1 hour

    const userLimits = this.sessionCounts.get(userId);
    if (!userLimits) return true;

    // Check if window has reset
    if (now.getTime() - userLimits.windowStart.getTime() > windowMs) {
      this.sessionCounts.set(userId, { count: 0, windowStart: now });
      return true;
    }

    return userLimits.count < this.maxSessionsPerHour;
  }

  recordSessionStart(userId: string): void {
    const now = new Date();
    const current = this.sessionCounts.get(userId);

    if (!current || now.getTime() - current.windowStart.getTime() > 60 * 60 * 1000) {
      this.sessionCounts.set(userId, { count: 1, windowStart: now });
    } else {
      current.count++;
    }
  }

  getRemainingSession(userId: string): number {
    const current = this.sessionCounts.get(userId);
    if (!current) return this.maxSessionsPerHour;
    return Math.max(0, this.maxSessionsPerHour - current.count);
  }
}
```

### Step 4: Audio Processing Queue

```typescript
import PQueue from 'p-queue';

// Audio processing has stricter limits
const audioQueue = new PQueue({
  concurrency: 2, // Max 2 concurrent audio requests
  interval: 1000, // Per second
  intervalCap: 3, // Max 3 per second
});

async function queuedAudioProcessing(
  client: SpeakClient,
  audioData: ArrayBuffer,
  options: AudioProcessingOptions
): Promise<PronunciationResult> {
  return audioQueue.add(async () => {
    return client.speech.analyze(audioData, options);
  });
}

// Monitor queue status
audioQueue.on('idle', () => {
  console.log('Audio queue idle');
});

audioQueue.on('add', () => {
  console.log(`Audio queue size: ${audioQueue.size}, pending: ${audioQueue.pending}`);
});
```

### Step 5: Rate Limit Monitor

```typescript
class SpeakRateLimitMonitor {
  private apiRemaining: number = 60;
  private audioRemaining: number = 20;
  private apiResetAt: Date = new Date();
  private audioResetAt: Date = new Date();

  updateFromHeaders(headers: Headers): void {
    // API limits
    if (headers.has('X-RateLimit-Remaining')) {
      this.apiRemaining = parseInt(headers.get('X-RateLimit-Remaining')!);
    }
    if (headers.has('X-RateLimit-Reset')) {
      this.apiResetAt = new Date(parseInt(headers.get('X-RateLimit-Reset')!) * 1000);
    }

    // Audio-specific limits
    if (headers.has('X-Audio-RateLimit-Remaining')) {
      this.audioRemaining = parseInt(headers.get('X-Audio-RateLimit-Remaining')!);
    }
    if (headers.has('X-Audio-RateLimit-Reset')) {
      this.audioResetAt = new Date(parseInt(headers.get('X-Audio-RateLimit-Reset')!) * 1000);
    }
  }

  shouldThrottleApi(): boolean {
    return this.apiRemaining < 5 && new Date() < this.apiResetAt;
  }

  shouldThrottleAudio(): boolean {
    return this.audioRemaining < 3 && new Date() < this.audioResetAt;
  }

  getApiWaitTime(): number {
    return Math.max(0, this.apiResetAt.getTime() - Date.now());
  }

  getAudioWaitTime(): number {
    return Math.max(0, this.audioResetAt.getTime() - Date.now());
  }

  getStatus(): RateLimitStatus {
    return {
      api: {
        remaining: this.apiRemaining,
        resetAt: this.apiResetAt,
        throttled: this.shouldThrottleApi(),
      },
      audio: {
        remaining: this.audioRemaining,
        resetAt: this.audioResetAt,
        throttled: this.shouldThrottleAudio(),
      },
    };
  }
}
```

### Step 6: Idempotency for Lesson Responses

```typescript
import crypto from 'crypto';

// Generate deterministic key for lesson responses
function generateIdempotencyKey(
  sessionId: string,
  responseText: string,
  timestamp: number
): string {
  const data = JSON.stringify({ sessionId, responseText, timestamp });
  return crypto.createHash('sha256').update(data).digest('hex').slice(0, 32);
}

async function idempotentSubmitResponse(
  session: LessonSession,
  response: LessonResponse,
  idempotencyKey?: string
): Promise<TutorFeedback> {
  const key = idempotencyKey || generateIdempotencyKey(
    session.id,
    response.text,
    Math.floor(Date.now() / 10000) // 10-second buckets
  );

  return session.submitResponse(response, {
    headers: { 'Idempotency-Key': key },
  });
}
```

## Output
- Reliable API calls with automatic retry
- Session-aware rate limiting
- Audio processing queue
- Rate limit monitoring
- Idempotent lesson responses

## Error Handling
| Header | Description | Action |
|--------|-------------|--------|
| X-RateLimit-Limit | Max requests | Monitor usage |
| X-RateLimit-Remaining | Remaining requests | Throttle if low |
| X-RateLimit-Reset | Reset timestamp | Wait until reset |
| X-Audio-RateLimit-Remaining | Audio requests left | Queue audio processing |
| Retry-After | Seconds to wait | Honor this value |

## Examples

### Pre-flight Rate Check
```typescript
async function safeSpeakCall<T>(
  monitor: SpeakRateLimitMonitor,
  operation: () => Promise<T>,
  isAudioOperation: boolean = false
): Promise<T> {
  // Check if we should wait
  if (isAudioOperation && monitor.shouldThrottleAudio()) {
    const waitTime = monitor.getAudioWaitTime();
    console.log(`Audio rate limit - waiting ${waitTime}ms`);
    await new Promise(r => setTimeout(r, waitTime));
  } else if (monitor.shouldThrottleApi()) {
    const waitTime = monitor.getApiWaitTime();
    console.log(`API rate limit - waiting ${waitTime}ms`);
    await new Promise(r => setTimeout(r, waitTime));
  }

  return withExponentialBackoff(operation);
}
```

### Batch Lesson Processing
```typescript
async function batchProcessLessons(
  client: SpeakClient,
  lessons: LessonConfig[],
  rateLimiter: SessionRateLimiter
): Promise<LessonResult[]> {
  const results: LessonResult[] = [];

  for (const lesson of lessons) {
    // Check rate limit before starting
    if (!rateLimiter.canStartSession(lesson.userId)) {
      console.log(`Rate limited for user ${lesson.userId}, skipping`);
      results.push({ success: false, error: 'RATE_LIMITED' });
      continue;
    }

    rateLimiter.recordSessionStart(lesson.userId);

    try {
      const result = await withExponentialBackoff(() =>
        client.tutor.runLesson(lesson)
      );
      results.push({ success: true, data: result });
    } catch (error) {
      results.push({ success: false, error });
    }

    // Small delay between lessons
    await new Promise(r => setTimeout(r, 100));
  }

  return results;
}
```

## Resources
- [Speak Rate Limits](https://developer.speak.com/docs/rate-limits)
- [p-queue Documentation](https://github.com/sindresorhus/p-queue)
- [Speak API Headers](https://developer.speak.com/docs/headers)

## Next Steps
For security configuration, see `speak-security-basics`.
