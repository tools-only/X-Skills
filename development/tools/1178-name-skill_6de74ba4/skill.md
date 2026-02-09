---
name: mistral-sdk-patterns
description: |
  Apply production-ready Mistral AI SDK patterns for TypeScript and Python.
  Use when implementing Mistral integrations, refactoring SDK usage,
  or establishing team coding standards for Mistral AI.
  Trigger with phrases like "mistral SDK patterns", "mistral best practices",
  "mistral code patterns", "idiomatic mistral".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Mistral AI SDK Patterns

## Overview
Production-ready patterns for Mistral AI SDK usage in TypeScript and Python.

## Prerequisites
- Completed `mistral-install-auth` setup
- Familiarity with async/await patterns
- Understanding of error handling best practices

## Instructions

### Step 1: Implement Singleton Pattern (Recommended)

**TypeScript**
```typescript
// src/mistral/client.ts
import Mistral from '@mistralai/mistralai';

let instance: Mistral | null = null;

export interface MistralConfig {
  apiKey: string;
  timeout?: number;
  maxRetries?: number;
}

export function getMistralClient(config?: Partial<MistralConfig>): Mistral {
  if (!instance) {
    const apiKey = config?.apiKey || process.env.MISTRAL_API_KEY;
    if (!apiKey) {
      throw new Error('MISTRAL_API_KEY is required');
    }
    instance = new Mistral({
      apiKey,
      timeout: config?.timeout ?? 30000,
    });
  }
  return instance;
}

// For testing - reset singleton
export function resetMistralClient(): void {
  instance = null;
}
```

**Python**
```python
# src/mistral/client.py
import os
from mistralai import Mistral
from typing import Optional

_instance: Optional[Mistral] = None

def get_mistral_client(api_key: Optional[str] = None) -> Mistral:
    global _instance
    if _instance is None:
        key = api_key or os.environ.get("MISTRAL_API_KEY")
        if not key:
            raise ValueError("MISTRAL_API_KEY is required")
        _instance = Mistral(api_key=key)
    return _instance

def reset_mistral_client() -> None:
    global _instance
    _instance = None
```

### Step 2: Add Error Handling Wrapper

**TypeScript**
```typescript
import Mistral from '@mistralai/mistralai';

interface MistralResult<T> {
  data: T | null;
  error: Error | null;
  usage?: { promptTokens: number; completionTokens: number; totalTokens: number };
}

async function safeMistralCall<T>(
  operation: () => Promise<T>
): Promise<MistralResult<T>> {
  try {
    const data = await operation();
    return { data, error: null };
  } catch (err) {
    const error = err as Error;
    console.error({
      name: error.name,
      message: error.message,
      // Extract status code if available
      status: (err as any).status,
    });
    return { data: null, error };
  }
}

// Usage
const result = await safeMistralCall(() =>
  client.chat.complete({
    model: 'mistral-small-latest',
    messages: [{ role: 'user', content: 'Hello!' }],
  })
);

if (result.error) {
  console.error('Chat failed:', result.error.message);
} else {
  console.log(result.data?.choices?.[0]?.message?.content);
}
```

### Step 3: Implement Retry Logic with Exponential Backoff

```typescript
interface RetryConfig {
  maxRetries: number;
  baseDelayMs: number;
  maxDelayMs: number;
  jitterMs: number;
}

async function withRetry<T>(
  operation: () => Promise<T>,
  config: RetryConfig = { maxRetries: 3, baseDelayMs: 1000, maxDelayMs: 32000, jitterMs: 500 }
): Promise<T> {
  for (let attempt = 0; attempt <= config.maxRetries; attempt++) {
    try {
      return await operation();
    } catch (err) {
      if (attempt === config.maxRetries) throw err;

      const status = (err as any).status;
      // Only retry on rate limits (429) or server errors (5xx)
      if (status !== 429 && (status < 500 || status >= 600)) throw err;

      const exponentialDelay = config.baseDelayMs * Math.pow(2, attempt);
      const jitter = Math.random() * config.jitterMs;
      const delay = Math.min(exponentialDelay + jitter, config.maxDelayMs);

      console.log(`Attempt ${attempt + 1} failed. Retrying in ${delay.toFixed(0)}ms...`);
      await new Promise(r => setTimeout(r, delay));
    }
  }
  throw new Error('Unreachable');
}
```

### Step 4: Type-Safe Chat Wrapper

```typescript
import Mistral from '@mistralai/mistralai';

type MistralModel =
  | 'mistral-large-latest'
  | 'mistral-medium-latest'
  | 'mistral-small-latest'
  | 'open-mistral-7b'
  | 'open-mixtral-8x7b';

interface ChatOptions {
  model?: MistralModel;
  temperature?: number;
  maxTokens?: number;
  topP?: number;
  stream?: boolean;
}

interface Message {
  role: 'system' | 'user' | 'assistant';
  content: string;
}

export async function chat(
  client: Mistral,
  messages: Message[],
  options: ChatOptions = {}
): Promise<string> {
  const response = await client.chat.complete({
    model: options.model ?? 'mistral-small-latest',
    messages,
    temperature: options.temperature,
    maxTokens: options.maxTokens,
    topP: options.topP,
  });

  return response.choices?.[0]?.message?.content ?? '';
}
```

### Step 5: Streaming Helper

```typescript
export async function* streamChat(
  client: Mistral,
  messages: Message[],
  options: ChatOptions = {}
): AsyncGenerator<string, void, unknown> {
  const stream = await client.chat.stream({
    model: options.model ?? 'mistral-small-latest',
    messages,
    temperature: options.temperature,
    maxTokens: options.maxTokens,
  });

  for await (const event of stream) {
    const content = event.data?.choices?.[0]?.delta?.content;
    if (content) {
      yield content;
    }
  }
}

// Usage
for await (const chunk of streamChat(client, messages)) {
  process.stdout.write(chunk);
}
```

## Output
- Type-safe client singleton
- Robust error handling with structured logging
- Automatic retry with exponential backoff
- Streaming support with async generators

## Error Handling
| Pattern | Use Case | Benefit |
|---------|----------|---------|
| Safe wrapper | All API calls | Prevents uncaught exceptions |
| Retry logic | Transient failures | Improves reliability |
| Type guards | Response validation | Catches API changes |
| Logging | All operations | Debugging and monitoring |

## Examples

### Factory Pattern (Multi-tenant)
```typescript
const clients = new Map<string, Mistral>();

export function getClientForTenant(tenantId: string): Mistral {
  if (!clients.has(tenantId)) {
    const apiKey = getTenantApiKey(tenantId);
    clients.set(tenantId, new Mistral({ apiKey }));
  }
  return clients.get(tenantId)!;
}
```

### Python Context Manager
```python
from contextlib import contextmanager
from mistralai import Mistral

@contextmanager
def mistral_client():
    client = Mistral(api_key=os.environ.get("MISTRAL_API_KEY"))
    try:
        yield client
    finally:
        pass  # Cleanup if needed

# Usage
with mistral_client() as client:
    response = client.chat.complete(
        model="mistral-small-latest",
        messages=[{"role": "user", "content": "Hello!"}]
    )
```

### Zod Response Validation
```typescript
import { z } from 'zod';

const mistralResponseSchema = z.object({
  id: z.string(),
  object: z.literal('chat.completion'),
  model: z.string(),
  choices: z.array(z.object({
    index: z.number(),
    message: z.object({
      role: z.enum(['assistant']),
      content: z.string(),
    }),
    finishReason: z.string(),
  })),
  usage: z.object({
    promptTokens: z.number(),
    completionTokens: z.number(),
    totalTokens: z.number(),
  }),
});

function validateResponse(response: unknown) {
  return mistralResponseSchema.parse(response);
}
```

## Resources
- [Mistral AI API Reference](https://docs.mistral.ai/api/)
- [Mistral AI Client Libraries](https://docs.mistral.ai/getting-started/clients/)
- [Zod Documentation](https://zod.dev/)

## Next Steps
Apply patterns in `mistral-core-workflow-a` for chat completions.
