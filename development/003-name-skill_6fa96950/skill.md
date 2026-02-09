---
name: mistral-migration-deep-dive
description: |
  Execute Mistral AI major migrations and re-architecture strategies.
  Use when migrating to Mistral AI from another provider, performing major refactoring,
  or re-platforming existing AI integrations to Mistral AI.
  Trigger with phrases like "migrate to mistral", "mistral migration",
  "switch to mistral", "mistral replatform", "openai to mistral".
allowed-tools: Read, Write, Edit, Bash(npm:*), Bash(node:*), Bash(kubectl:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Mistral AI Migration Deep Dive

## Overview
Comprehensive guide for migrating to Mistral AI from other providers or major version upgrades.

## Prerequisites
- Current system documentation
- Mistral AI SDK installed
- Feature flag infrastructure
- Rollback strategy tested

## Migration Types

| Type | Complexity | Duration | Risk |
|------|-----------|----------|------|
| Fresh install | Low | Days | Low |
| OpenAI to Mistral | Medium | Weeks | Medium |
| Multi-provider | Medium | Weeks | Medium |
| Full replatform | High | Months | High |

## Instructions

### Step 1: Pre-Migration Assessment

```bash
# Document current implementation
echo "=== Pre-Migration Assessment ==="

# Find all AI-related code
find . -name "*.ts" -o -name "*.py" | xargs grep -l "openai\|anthropic\|ai" > ai-files.txt
echo "Files with AI code: $(wc -l < ai-files.txt)"

# Count integration points
grep -r "chat.completions\|createChatCompletion" src/ --include="*.ts" | wc -l

# Check current SDK versions
npm list openai @anthropic-ai/sdk 2>/dev/null || echo "No existing AI SDKs"
```

```typescript
// scripts/assess-migration.ts
interface MigrationAssessment {
  currentProvider: string;
  integrationPoints: number;
  features: string[];
  estimatedEffort: 'low' | 'medium' | 'high';
  risks: string[];
}

async function assessMigration(): Promise<MigrationAssessment> {
  const files = await glob('src/**/*.{ts,js}');
  const features = new Set<string>();
  let integrationPoints = 0;

  for (const file of files) {
    const content = await fs.readFile(file, 'utf-8');

    // Detect features
    if (/chat\.completions|createChatCompletion/i.test(content)) {
      features.add('chat');
      integrationPoints++;
    }
    if (/embeddings\.create|createEmbedding/i.test(content)) {
      features.add('embeddings');
      integrationPoints++;
    }
    if (/function_call|tools/i.test(content)) {
      features.add('function_calling');
      integrationPoints++;
    }
    if (/stream/i.test(content)) {
      features.add('streaming');
      integrationPoints++;
    }
  }

  return {
    currentProvider: detectProvider(files),
    integrationPoints,
    features: Array.from(features),
    estimatedEffort: integrationPoints > 10 ? 'high' : integrationPoints > 3 ? 'medium' : 'low',
    risks: identifyRisks(features),
  };
}
```

### Step 2: Create Adapter Layer

```typescript
// src/ai/adapter.ts
// Provider-agnostic interface

export interface Message {
  role: 'system' | 'user' | 'assistant';
  content: string;
}

export interface ChatOptions {
  model?: string;
  temperature?: number;
  maxTokens?: number;
  stream?: boolean;
}

export interface ChatResponse {
  content: string;
  usage?: {
    inputTokens: number;
    outputTokens: number;
  };
}

export interface AIAdapter {
  chat(messages: Message[], options?: ChatOptions): Promise<ChatResponse>;
  chatStream(messages: Message[], options?: ChatOptions): AsyncGenerator<string>;
  embed(text: string | string[]): Promise<number[][]>;
}
```

### Step 3: Implement OpenAI Adapter (Current)

```typescript
// src/ai/adapters/openai.ts
import OpenAI from 'openai';
import { AIAdapter, Message, ChatOptions, ChatResponse } from '../adapter';

export class OpenAIAdapter implements AIAdapter {
  private client: OpenAI;

  constructor() {
    this.client = new OpenAI({ apiKey: process.env.OPENAI_API_KEY });
  }

  async chat(messages: Message[], options?: ChatOptions): Promise<ChatResponse> {
    const response = await this.client.chat.completions.create({
      model: options?.model || 'gpt-3.5-turbo',
      messages,
      temperature: options?.temperature,
      max_tokens: options?.maxTokens,
    });

    return {
      content: response.choices[0]?.message?.content || '',
      usage: response.usage ? {
        inputTokens: response.usage.prompt_tokens,
        outputTokens: response.usage.completion_tokens,
      } : undefined,
    };
  }

  async *chatStream(messages: Message[], options?: ChatOptions): AsyncGenerator<string> {
    const stream = await this.client.chat.completions.create({
      model: options?.model || 'gpt-3.5-turbo',
      messages,
      stream: true,
    });

    for await (const chunk of stream) {
      const content = chunk.choices[0]?.delta?.content;
      if (content) yield content;
    }
  }

  async embed(text: string | string[]): Promise<number[][]> {
    const input = Array.isArray(text) ? text : [text];
    const response = await this.client.embeddings.create({
      model: 'text-embedding-ada-002',
      input,
    });
    return response.data.map(d => d.embedding);
  }
}
```

### Step 4: Implement Mistral Adapter (Target)

```typescript
// src/ai/adapters/mistral.ts
import Mistral from '@mistralai/mistralai';
import { AIAdapter, Message, ChatOptions, ChatResponse } from '../adapter';

export class MistralAdapter implements AIAdapter {
  private client: Mistral;

  constructor() {
    this.client = new Mistral({ apiKey: process.env.MISTRAL_API_KEY });
  }

  async chat(messages: Message[], options?: ChatOptions): Promise<ChatResponse> {
    const response = await this.client.chat.complete({
      model: options?.model || 'mistral-small-latest',
      messages,
      temperature: options?.temperature,
      maxTokens: options?.maxTokens,
    });

    return {
      content: response.choices?.[0]?.message?.content || '',
      usage: response.usage ? {
        inputTokens: response.usage.promptTokens || 0,
        outputTokens: response.usage.completionTokens || 0,
      } : undefined,
    };
  }

  async *chatStream(messages: Message[], options?: ChatOptions): AsyncGenerator<string> {
    const stream = await this.client.chat.stream({
      model: options?.model || 'mistral-small-latest',
      messages,
      temperature: options?.temperature,
      maxTokens: options?.maxTokens,
    });

    for await (const event of stream) {
      const content = event.data?.choices?.[0]?.delta?.content;
      if (content) yield content;
    }
  }

  async embed(text: string | string[]): Promise<number[][]> {
    const input = Array.isArray(text) ? text : [text];
    const response = await this.client.embeddings.create({
      model: 'mistral-embed',
      inputs: input,
    });
    return response.data.map(d => d.embedding);
  }
}
```

### Step 5: Feature Flag Controlled Migration

```typescript
// src/ai/factory.ts
import { AIAdapter } from './adapter';
import { OpenAIAdapter } from './adapters/openai';
import { MistralAdapter } from './adapters/mistral';

type Provider = 'openai' | 'mistral';

function getAIProvider(): Provider {
  // Feature flag based migration
  const mistralPercentage = parseInt(process.env.MISTRAL_ROLLOUT_PERCENT || '0');
  const random = Math.random() * 100;

  if (random < mistralPercentage) {
    return 'mistral';
  }

  return 'openai';
}

export function createAIAdapter(): AIAdapter {
  const provider = getAIProvider();

  switch (provider) {
    case 'mistral':
      console.log('[AI] Using Mistral adapter');
      return new MistralAdapter();
    case 'openai':
    default:
      console.log('[AI] Using OpenAI adapter');
      return new OpenAIAdapter();
  }
}
```

### Step 6: Gradual Rollout

```bash
# Phase 1: 0% Mistral (validation)
export MISTRAL_ROLLOUT_PERCENT=0
# Run tests, verify adapter works

# Phase 2: 5% Mistral (canary)
export MISTRAL_ROLLOUT_PERCENT=5
# Monitor for errors, compare latency

# Phase 3: 25% Mistral
export MISTRAL_ROLLOUT_PERCENT=25
# Monitor for 24-48 hours

# Phase 4: 50% Mistral
export MISTRAL_ROLLOUT_PERCENT=50
# Monitor for 24-48 hours

# Phase 5: 100% Mistral
export MISTRAL_ROLLOUT_PERCENT=100
# Full migration complete
```

### Step 7: Model Mapping

```typescript
// src/ai/model-mapping.ts

interface ModelMapping {
  openai: string;
  mistral: string;
  notes: string;
}

const MODEL_MAPPINGS: ModelMapping[] = [
  {
    openai: 'gpt-3.5-turbo',
    mistral: 'mistral-small-latest',
    notes: 'Fast, cost-effective',
  },
  {
    openai: 'gpt-4',
    mistral: 'mistral-large-latest',
    notes: 'Complex reasoning',
  },
  {
    openai: 'gpt-4-turbo',
    mistral: 'mistral-large-latest',
    notes: 'Best available',
  },
  {
    openai: 'text-embedding-ada-002',
    mistral: 'mistral-embed',
    notes: '1024 dimensions',
  },
];

export function mapModel(openaiModel: string): string {
  const mapping = MODEL_MAPPINGS.find(m => m.openai === openaiModel);
  return mapping?.mistral || 'mistral-small-latest';
}
```

### Step 8: Validation & Testing

```typescript
// tests/migration/compare-outputs.test.ts
import { describe, it, expect } from 'vitest';
import { OpenAIAdapter } from '../../src/ai/adapters/openai';
import { MistralAdapter } from '../../src/ai/adapters/mistral';

describe('Migration Validation', () => {
  const openai = new OpenAIAdapter();
  const mistral = new MistralAdapter();

  const testCases = [
    { name: 'Simple greeting', messages: [{ role: 'user', content: 'Hello' }] },
    { name: 'Math question', messages: [{ role: 'user', content: 'What is 2+2?' }] },
    { name: 'Code generation', messages: [{ role: 'user', content: 'Write hello world in Python' }] },
  ];

  for (const testCase of testCases) {
    it(`should produce similar output: ${testCase.name}`, async () => {
      const [openaiResult, mistralResult] = await Promise.all([
        openai.chat(testCase.messages, { temperature: 0 }),
        mistral.chat(testCase.messages, { temperature: 0 }),
      ]);

      // Both should return non-empty content
      expect(openaiResult.content.length).toBeGreaterThan(0);
      expect(mistralResult.content.length).toBeGreaterThan(0);

      // Log for manual review
      console.log(`[${testCase.name}]`);
      console.log('OpenAI:', openaiResult.content);
      console.log('Mistral:', mistralResult.content);
    });
  }
});
```

### Step 9: Rollback Plan

```bash
#!/bin/bash
# rollback-to-openai.sh

echo "=== Rolling back to OpenAI ==="

# 1. Set rollout percentage to 0
kubectl set env deployment/ai-service MISTRAL_ROLLOUT_PERCENT=0

# 2. Verify rollback
kubectl rollout status deployment/ai-service

# 3. Check health
curl -sf https://api.yourapp.com/health | jq '.services.ai'

# 4. Alert team
echo "Rollback complete. Mistral disabled."
```

## Output
- Migration assessment complete
- Adapter layer implemented
- Gradual rollout in progress
- Rollback procedure ready

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Different output format | API differences | Normalize in adapter |
| Missing feature | Not supported | Implement fallback |
| Performance difference | Model characteristics | Adjust timeouts |
| Cost increase | Token differences | Monitor and optimize |

## Examples

### Quick A/B Comparison
```typescript
const [openaiResponse, mistralResponse] = await Promise.all([
  openaiAdapter.chat(messages, { temperature: 0 }),
  mistralAdapter.chat(messages, { temperature: 0 }),
]);

console.log('OpenAI tokens:', openaiResponse.usage);
console.log('Mistral tokens:', mistralResponse.usage);
```

## Resources
- [Mistral AI Documentation](https://docs.mistral.ai/)
- [Strangler Fig Pattern](https://martinfowler.com/bliki/StranglerFigApplication.html)
- [Feature Flags Best Practices](https://www.martinfowler.com/articles/feature-toggles.html)

## Completion
Congratulations! You've completed the Mistral AI skill pack. For ongoing support, visit docs.mistral.ai.
