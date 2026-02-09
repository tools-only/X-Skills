---
name: mistral-reference-architecture
description: |
  Implement Mistral AI reference architecture with best-practice project layout.
  Use when designing new Mistral AI integrations, reviewing project structure,
  or establishing architecture standards for Mistral AI applications.
  Trigger with phrases like "mistral architecture", "mistral best practices",
  "mistral project structure", "how to organize mistral", "mistral layout".
allowed-tools: Read, Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Mistral AI Reference Architecture

## Overview
Production-ready architecture patterns for Mistral AI integrations.

## Prerequisites
- Understanding of layered architecture
- Mistral AI SDK knowledge
- TypeScript project setup
- Testing framework configured

## Project Structure

```
my-mistral-project/
├── src/
│   ├── mistral/
│   │   ├── client.ts           # Singleton client wrapper
│   │   ├── config.ts           # Environment configuration
│   │   ├── types.ts            # TypeScript types
│   │   ├── errors.ts           # Custom error classes
│   │   ├── prompts/            # Prompt templates
│   │   │   ├── system.ts       # System prompts
│   │   │   └── templates.ts    # Reusable templates
│   │   └── handlers/
│   │       ├── chat.ts         # Chat completion handlers
│   │       ├── embeddings.ts   # Embedding handlers
│   │       └── tools.ts        # Function calling handlers
│   ├── services/
│   │   └── ai/
│   │       ├── index.ts        # Service facade
│   │       ├── chat.ts         # Chat service
│   │       ├── rag.ts          # RAG implementation
│   │       └── cache.ts        # Caching layer
│   ├── api/
│   │   └── chat/
│   │       ├── route.ts        # API endpoint
│   │       └── stream.ts       # Streaming endpoint
│   └── jobs/
│       └── ai/
│           └── batch.ts        # Background processing
├── tests/
│   ├── unit/
│   │   └── mistral/
│   └── integration/
│       └── mistral/
├── config/
│   ├── mistral.development.json
│   ├── mistral.staging.json
│   └── mistral.production.json
└── docs/
    └── mistral/
        ├── SETUP.md
        └── RUNBOOK.md
```

## Layer Architecture

```
┌─────────────────────────────────────────┐
│             API Layer                    │
│   (Routes, Controllers, Middleware)      │
├─────────────────────────────────────────┤
│           Service Layer                  │
│  (Business Logic, Orchestration)         │
├─────────────────────────────────────────┤
│          Mistral Layer                   │
│   (Client, Prompts, Error Handling)      │
├─────────────────────────────────────────┤
│         Infrastructure Layer             │
│    (Cache, Queue, Monitoring)            │
└─────────────────────────────────────────┘
```

## Key Components

### Step 1: Client Wrapper

```typescript
// src/mistral/client.ts
import Mistral from '@mistralai/mistralai';
import { getMistralConfig } from './config';

let instance: Mistral | null = null;

export function getMistralClient(): Mistral {
  if (!instance) {
    const config = getMistralConfig();
    instance = new Mistral({
      apiKey: config.apiKey,
      timeout: config.timeout,
    });
  }
  return instance;
}

export function resetMistralClient(): void {
  instance = null;
}
```

### Step 2: Configuration Management

```typescript
// src/mistral/config.ts
import { z } from 'zod';

const configSchema = z.object({
  apiKey: z.string().min(1),
  model: z.string().default('mistral-small-latest'),
  timeout: z.number().default(30000),
  maxRetries: z.number().default(3),
  cache: z.object({
    enabled: z.boolean().default(true),
    ttlSeconds: z.number().default(300),
  }).default({}),
});

export type MistralConfig = z.infer<typeof configSchema>;

export function getMistralConfig(): MistralConfig {
  const env = process.env.NODE_ENV || 'development';

  const baseConfig = {
    apiKey: process.env.MISTRAL_API_KEY,
    model: process.env.MISTRAL_MODEL,
    timeout: parseInt(process.env.MISTRAL_TIMEOUT || '30000'),
    maxRetries: parseInt(process.env.MISTRAL_MAX_RETRIES || '3'),
  };

  return configSchema.parse(baseConfig);
}
```

### Step 3: Error Handling

```typescript
// src/mistral/errors.ts
export class MistralServiceError extends Error {
  constructor(
    message: string,
    public readonly code: string,
    public readonly status?: number,
    public readonly retryable: boolean = false,
    public readonly originalError?: Error
  ) {
    super(message);
    this.name = 'MistralServiceError';
  }
}

export function wrapMistralError(error: unknown): MistralServiceError {
  if (error instanceof MistralServiceError) return error;

  const err = error as any;

  // Rate limit
  if (err.status === 429) {
    return new MistralServiceError(
      'Rate limit exceeded',
      'RATE_LIMIT',
      429,
      true,
      err
    );
  }

  // Auth error
  if (err.status === 401) {
    return new MistralServiceError(
      'Authentication failed',
      'AUTH_ERROR',
      401,
      false,
      err
    );
  }

  // Server error
  if (err.status >= 500) {
    return new MistralServiceError(
      'Mistral service error',
      'SERVICE_ERROR',
      err.status,
      true,
      err
    );
  }

  return new MistralServiceError(
    err.message || 'Unknown error',
    'UNKNOWN',
    err.status,
    false,
    err
  );
}
```

### Step 4: Service Layer

```typescript
// src/services/ai/chat.ts
import { getMistralClient } from '../../mistral/client';
import { wrapMistralError } from '../../mistral/errors';
import { withRetry } from '../../utils/retry';
import { cache } from './cache';

export interface ChatOptions {
  model?: string;
  temperature?: number;
  maxTokens?: number;
  useCache?: boolean;
}

export class ChatService {
  async complete(
    messages: Array<{ role: string; content: string }>,
    options: ChatOptions = {}
  ): Promise<string> {
    const client = getMistralClient();
    const model = options.model || 'mistral-small-latest';

    // Check cache for deterministic requests
    if (options.useCache && options.temperature === 0) {
      const cached = await cache.get(messages, model);
      if (cached) return cached;
    }

    try {
      const response = await withRetry(() =>
        client.chat.complete({
          model,
          messages,
          temperature: options.temperature,
          maxTokens: options.maxTokens,
        })
      );

      const content = response.choices?.[0]?.message?.content ?? '';

      // Cache deterministic responses
      if (options.useCache && options.temperature === 0) {
        await cache.set(messages, model, content);
      }

      return content;
    } catch (error) {
      throw wrapMistralError(error);
    }
  }

  async *stream(
    messages: Array<{ role: string; content: string }>,
    options: ChatOptions = {}
  ): AsyncGenerator<string> {
    const client = getMistralClient();
    const model = options.model || 'mistral-small-latest';

    try {
      const stream = await client.chat.stream({
        model,
        messages,
        temperature: options.temperature,
        maxTokens: options.maxTokens,
      });

      for await (const event of stream) {
        const content = event.data?.choices?.[0]?.delta?.content;
        if (content) yield content;
      }
    } catch (error) {
      throw wrapMistralError(error);
    }
  }
}

export const chatService = new ChatService();
```

### Step 5: Health Check

```typescript
// src/mistral/health.ts
import { getMistralClient } from './client';

export interface HealthStatus {
  status: 'healthy' | 'degraded' | 'unhealthy';
  latencyMs: number;
  error?: string;
}

export async function checkMistralHealth(): Promise<HealthStatus> {
  const client = getMistralClient();
  const start = Date.now();

  try {
    await client.models.list();
    return {
      status: 'healthy',
      latencyMs: Date.now() - start,
    };
  } catch (error: any) {
    return {
      status: 'unhealthy',
      latencyMs: Date.now() - start,
      error: error.message,
    };
  }
}
```

### Step 6: Prompt Templates

```typescript
// src/mistral/prompts/templates.ts
export interface PromptTemplate {
  system: string;
  user: (vars: Record<string, string>) => string;
}

export const templates: Record<string, PromptTemplate> = {
  summarize: {
    system: 'You are a helpful assistant that creates concise summaries.',
    user: ({ text, maxWords }) =>
      `Summarize the following text in ${maxWords || '100'} words or less:\n\n${text}`,
  },

  classify: {
    system: 'You are a classifier. Respond with only the category name.',
    user: ({ text, categories }) =>
      `Classify the following text into one of these categories: ${categories}\n\nText: ${text}`,
  },

  codeReview: {
    system: 'You are an expert code reviewer. Be concise and actionable.',
    user: ({ code, language }) =>
      `Review this ${language} code and suggest improvements:\n\n\`\`\`${language}\n${code}\n\`\`\``,
  },
};

export function buildPrompt(
  templateName: keyof typeof templates,
  vars: Record<string, string>
): Array<{ role: string; content: string }> {
  const template = templates[templateName];
  return [
    { role: 'system', content: template.system },
    { role: 'user', content: template.user(vars) },
  ];
}
```

## Data Flow Diagram

```
User Request
     │
     ▼
┌─────────────┐
│   API       │
│   Gateway   │
└──────┬──────┘
       │
       ▼
┌─────────────┐    ┌─────────────┐
│   Service   │───▶│   Cache     │
│   Layer     │    │   (Redis)   │
└──────┬──────┘    └─────────────┘
       │
       ▼
┌─────────────┐
│   Mistral   │
│   Client    │
└──────┬──────┘
       │
       ▼
┌─────────────┐
│  Mistral    │
│   API       │
└─────────────┘
```

## Instructions

### Step 1: Create Directory Structure
```bash
mkdir -p src/mistral/{handlers,prompts}
mkdir -p src/services/ai
mkdir -p src/api/chat
mkdir -p tests/{unit,integration}/mistral
mkdir -p config docs/mistral
```

### Step 2: Implement Core Components
Create client wrapper, config, and error handling.

### Step 3: Build Service Layer
Implement chat service with caching and retry.

### Step 4: Add Health Checks
Configure health endpoint for monitoring.

## Output
- Structured project layout
- Client wrapper with retry
- Error handling implemented
- Health checks configured

## Resources
- [Mistral AI API Reference](https://docs.mistral.ai/api/)
- [Clean Architecture](https://blog.cleancoder.com/uncle-bob/2012/08/13/the-clean-architecture.html)

## Next Steps
For multi-environment setup, see `mistral-multi-env-setup`.
