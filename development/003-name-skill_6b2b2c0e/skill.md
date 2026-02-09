---
name: speak-sdk-patterns
description: |
  Apply production-ready Speak SDK patterns for TypeScript and Python.
  Use when implementing Speak integrations, refactoring SDK usage,
  or establishing team coding standards for language learning features.
  Trigger with phrases like "speak SDK patterns", "speak best practices",
  "speak code patterns", "idiomatic speak".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Speak SDK Patterns

## Overview
Production-ready patterns for Speak SDK usage in TypeScript and Python for AI-powered language learning.

## Prerequisites
- Completed `speak-install-auth` setup
- Familiarity with async/await patterns
- Understanding of error handling best practices

## Instructions

### Step 1: Implement Singleton Pattern for Client
```typescript
// src/speak/client.ts
import { SpeakClient } from '@speak/language-sdk';

let instance: SpeakClient | null = null;

export interface SpeakConfig {
  apiKey: string;
  appId: string;
  language: string;
  nativeLanguage?: string;
}

export function getSpeakClient(config?: Partial<SpeakConfig>): SpeakClient {
  if (!instance) {
    instance = new SpeakClient({
      apiKey: config?.apiKey || process.env.SPEAK_API_KEY!,
      appId: config?.appId || process.env.SPEAK_APP_ID!,
      language: config?.language || process.env.SPEAK_TARGET_LANGUAGE || 'es',
      nativeLanguage: config?.nativeLanguage || 'en',
    });
  }
  return instance;
}

export function resetSpeakClient(): void {
  instance = null;
}
```

### Step 2: Add Error Handling Wrapper
```typescript
import { SpeakError, RateLimitError, AuthError } from '@speak/language-sdk';

interface SpeakResult<T> {
  data: T | null;
  error: SpeakServiceError | null;
}

class SpeakServiceError extends Error {
  constructor(
    message: string,
    public readonly code: string,
    public readonly retryable: boolean,
    public readonly originalError?: Error
  ) {
    super(message);
    this.name = 'SpeakServiceError';
  }
}

async function safeSpeakCall<T>(
  operation: () => Promise<T>
): Promise<SpeakResult<T>> {
  try {
    const data = await operation();
    return { data, error: null };
  } catch (err) {
    if (err instanceof RateLimitError) {
      return {
        data: null,
        error: new SpeakServiceError(
          'Rate limit exceeded',
          'RATE_LIMITED',
          true,
          err
        ),
      };
    }
    if (err instanceof AuthError) {
      return {
        data: null,
        error: new SpeakServiceError(
          'Authentication failed',
          'AUTH_ERROR',
          false,
          err
        ),
      };
    }
    if (err instanceof SpeakError) {
      console.error({
        code: err.code,
        message: err.message,
        requestId: err.requestId,
      });
    }
    return {
      data: null,
      error: new SpeakServiceError(
        err instanceof Error ? err.message : 'Unknown error',
        'UNKNOWN',
        true,
        err instanceof Error ? err : undefined
      ),
    };
  }
}
```

### Step 3: Implement Lesson Session Manager
```typescript
// src/speak/session-manager.ts
import { LessonSession, SessionConfig } from '@speak/language-sdk';

interface ActiveSession {
  session: LessonSession;
  startedAt: Date;
  topic: string;
  exchanges: number;
}

class LessonSessionManager {
  private activeSessions: Map<string, ActiveSession> = new Map();
  private client: SpeakClient;

  constructor(client: SpeakClient) {
    this.client = client;
  }

  async startSession(
    userId: string,
    config: SessionConfig
  ): Promise<LessonSession> {
    // End any existing session for this user
    await this.endSession(userId);

    const session = await this.client.tutor.startSession(config);

    this.activeSessions.set(userId, {
      session,
      startedAt: new Date(),
      topic: config.topic,
      exchanges: 0,
    });

    return session;
  }

  async submitResponse(
    userId: string,
    response: { text: string; audioData?: ArrayBuffer }
  ) {
    const active = this.activeSessions.get(userId);
    if (!active) {
      throw new Error('No active session for user');
    }

    active.exchanges++;
    return active.session.submitResponse(response);
  }

  async endSession(userId: string): Promise<SessionSummary | null> {
    const active = this.activeSessions.get(userId);
    if (!active) return null;

    const summary = await active.session.end();
    this.activeSessions.delete(userId);

    return {
      ...summary,
      duration: Date.now() - active.startedAt.getTime(),
      totalExchanges: active.exchanges,
    };
  }

  getActiveSession(userId: string): ActiveSession | undefined {
    return this.activeSessions.get(userId);
  }
}
```

### Step 4: Implement Retry Logic with Backoff
```typescript
async function withRetry<T>(
  operation: () => Promise<T>,
  config = { maxRetries: 3, baseDelayMs: 1000, maxDelayMs: 10000 }
): Promise<T> {
  for (let attempt = 1; attempt <= config.maxRetries; attempt++) {
    try {
      return await operation();
    } catch (err: any) {
      if (attempt === config.maxRetries) throw err;

      // Only retry on transient errors
      const isRetryable =
        err.code === 'RATE_LIMITED' ||
        err.code === 'TIMEOUT' ||
        (err.status >= 500 && err.status < 600);

      if (!isRetryable) throw err;

      const delay = Math.min(
        config.baseDelayMs * Math.pow(2, attempt - 1),
        config.maxDelayMs
      );
      console.log(`Speak API retry ${attempt}/${config.maxRetries} in ${delay}ms`);
      await new Promise(r => setTimeout(r, delay));
    }
  }
  throw new Error('Unreachable');
}
```

### Step 5: Speech Recognition Wrapper
```typescript
// src/speak/speech.ts
import { SpeechRecognizer, PronunciationResult } from '@speak/language-sdk';

interface SpeechConfig {
  language: string;
  continuous?: boolean;
  interimResults?: boolean;
}

class SpeakSpeechService {
  private recognizer: SpeechRecognizer | null = null;
  private client: SpeakClient;

  constructor(client: SpeakClient) {
    this.client = client;
  }

  async initialize(config: SpeechConfig): Promise<void> {
    this.recognizer = new SpeechRecognizer(this.client, {
      language: config.language,
      continuous: config.continuous ?? false,
      interimResults: config.interimResults ?? true,
    });
  }

  async recognizeAudio(audioData: ArrayBuffer): Promise<PronunciationResult> {
    if (!this.recognizer) {
      throw new Error('Speech recognizer not initialized');
    }

    return this.recognizer.recognize(audioData);
  }

  async scorePronunciation(
    audioData: ArrayBuffer,
    expectedText: string
  ): Promise<PronunciationResult> {
    return this.client.speech.score({
      audio: audioData,
      expectedText,
      detailed: true,
    });
  }
}
```

## Output
- Type-safe client singleton
- Robust error handling with structured logging
- Session management for concurrent users
- Automatic retry with exponential backoff
- Speech recognition wrapper with scoring

## Error Handling
| Pattern | Use Case | Benefit |
|---------|----------|---------|
| Safe wrapper | All API calls | Prevents uncaught exceptions |
| Retry logic | Transient failures | Improves reliability |
| Session manager | Multi-user apps | Clean session lifecycle |
| Speech service | Audio processing | Centralized recognition |

## Examples

### Factory Pattern (Multi-tenant Language Schools)
```typescript
const clients = new Map<string, SpeakClient>();

export function getClientForSchool(schoolId: string): SpeakClient {
  if (!clients.has(schoolId)) {
    const config = getSchoolConfig(schoolId);
    clients.set(schoolId, new SpeakClient({
      apiKey: config.apiKey,
      appId: config.appId,
      language: config.defaultLanguage,
    }));
  }
  return clients.get(schoolId)!;
}
```

### Python Context Manager
```python
from contextlib import asynccontextmanager
from speak_sdk import SpeakClient, LessonSession

@asynccontextmanager
async def lesson_session(user_id: str, topic: str):
    client = SpeakClient()
    session = await client.tutor.start_session(topic=topic)
    try:
        yield session
    finally:
        summary = await session.end()
        await log_session_summary(user_id, summary)
```

### Zod Validation for Lesson Results
```typescript
import { z } from 'zod';

const pronunciationResultSchema = z.object({
  overall: z.number().min(0).max(100),
  fluency: z.number().min(0).max(100),
  accuracy: z.number().min(0).max(100),
  pronunciation: z.number().min(0).max(100),
  words: z.array(z.object({
    word: z.string(),
    score: z.number(),
    phonemes: z.array(z.object({
      phoneme: z.string(),
      score: z.number(),
    })).optional(),
  })),
});

const lessonFeedbackSchema = z.object({
  text: z.string(),
  audioUrl: z.string().url().optional(),
  pronunciationScore: z.number().min(0).max(100),
  grammarCorrections: z.array(z.object({
    original: z.string(),
    corrected: z.string(),
    explanation: z.string(),
  })),
  suggestions: z.array(z.string()),
});
```

## Resources
- [Speak SDK Reference](https://developer.speak.com/sdk)
- [Speak API Types](https://developer.speak.com/types)
- [Zod Documentation](https://zod.dev/)

## Next Steps
Apply patterns in `speak-core-workflow-a` for conversation practice implementation.
