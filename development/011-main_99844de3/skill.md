# AI-Integration Skill

Best Practices für Google Gemini + Anthropic Claude Integration.

## Hybrid AI-Ansatz (ManufacturingInsideAnalyzer)

```ts
// Gemini 2.0 Flash: Schnelle Analysen
import { GoogleGenerativeAI } from '@google/genai';

const genai = new GoogleGenerativeAI(apiKey);
const model = genai.getGenerativeModel({
  model: 'gemini-2.0-flash-exp',
  systemInstruction: 'Du bist ein Manufacturing-Experte...'
});

// Claude 3.5 Haiku: Chatbot
import Anthropic from '@anthropic-ai/sdk';

const anthropic = new Anthropic({ apiKey });
const message = await anthropic.messages.create({
  model: 'claude-3-5-haiku-20241022',
  max_tokens: 1024,
  messages: [{ role: 'user', content: prompt }]
});
```

## Structured Output (Gemini)

```ts
// ✅ JSON Schema für predictable Responses
const schema = {
  type: 'object',
  properties: {
    kpis: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          name: { type: 'string' },
          value: { type: 'number' }
        }
      }
    }
  },
  required: ['kpis']
};

const result = await model.generateContent({
  contents: [{ role: 'user', parts: [{ text: prompt }] }],
  generationConfig: {
    responseMimeType: 'application/json',
    responseSchema: schema
  }
});
```

## Timeout & Error Handling

```ts
// ✅ 10-Sekunden Vercel-Limit beachten
const controller = new AbortController();
const timeout = setTimeout(() => controller.abort(), 9000);

try {
  const response = await fetch(apiUrl, {
    signal: controller.signal,
    body: JSON.stringify(prompt)
  });
} catch (err) {
  if (err.name === 'AbortError') {
    return { error: 'Timeout', fallback: true };
  }
} finally {
  clearTimeout(timeout);
}
```

## Rate Limit Handling

```ts
// Exponential Backoff
async function callWithRetry(fn, maxRetries = 3) {
  for (let i = 0; i < maxRetries; i++) {
    try {
      return await fn();
    } catch (err) {
      if (err.status === 429 && i < maxRetries - 1) {
        await sleep(Math.pow(2, i) * 1000);
        continue;
      }
      throw err;
    }
  }
}
```

## Kosten-Monitoring

```ts
// Track Tokens für Budget-Kontrolle
const usage = result.usageMetadata;
await logUsage({
  model: 'gemini-2.0-flash',
  inputTokens: usage.promptTokenCount,
  outputTokens: usage.candidatesTokenCount,
  cost: calculateCost(usage)
});
```

## digitalTwin (Gemini 2.5 Flash Exp)

```ts
// Pose-Estimation via Gemini Vision
const imageParts = [{
  inlineData: {
    mimeType: 'image/jpeg',
    data: base64ImageData
  }
}];

const result = await model.generateContent({
  contents: [{
    role: 'user',
    parts: [
      { text: 'Analyze ergonomic posture and return REBA score' },
      ...imageParts
    ]
  }]
});
```

## Sicherheit

- ✅ API-Keys nur in `.env` (NIEMALS in Code)
- ✅ Input Validation (max length, format)
- ✅ Output Sanitization (HTML-escape)
- ✅ Medical Claims Filtering (ManufacturingInsideAnalyzer)
