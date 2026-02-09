---
name: mistral-deploy-integration
description: |
  Deploy Mistral AI integrations to Vercel, Fly.io, and Cloud Run platforms.
  Use when deploying Mistral AI-powered applications to production,
  configuring platform-specific secrets, or setting up deployment pipelines.
  Trigger with phrases like "deploy mistral", "mistral Vercel",
  "mistral production deploy", "mistral Cloud Run", "mistral Fly.io".
allowed-tools: Read, Write, Edit, Bash(vercel:*), Bash(fly:*), Bash(gcloud:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Mistral AI Deploy Integration

## Overview
Deploy Mistral AI-powered applications to popular platforms with proper secrets management.

## Prerequisites
- Mistral AI API keys for production environment
- Platform CLI installed (vercel, fly, or gcloud)
- Application code ready for deployment
- Environment variables documented

## Vercel Deployment

### Environment Setup

```bash
# Add Mistral secrets to Vercel
vercel env add MISTRAL_API_KEY production
# Enter your API key when prompted

# Or use CLI with secret
vercel secrets add mistral-api-key "your-api-key"

# Link to project
vercel link

# Deploy preview
vercel

# Deploy production
vercel --prod
```

### vercel.json Configuration

```json
{
  "env": {
    "MISTRAL_API_KEY": "@mistral-api-key"
  },
  "functions": {
    "api/**/*.ts": {
      "maxDuration": 60
    }
  },
  "regions": ["iad1", "sfo1"]
}
```

### Vercel Edge Function Example

```typescript
// api/chat/route.ts (Next.js App Router)
import { NextRequest } from 'next/server';
import Mistral from '@mistralai/mistralai';

const client = new Mistral({
  apiKey: process.env.MISTRAL_API_KEY!,
});

export async function POST(request: NextRequest) {
  const { messages } = await request.json();

  const stream = await client.chat.stream({
    model: 'mistral-small-latest',
    messages,
  });

  const encoder = new TextEncoder();
  const readable = new ReadableStream({
    async start(controller) {
      for await (const event of stream) {
        const content = event.data?.choices?.[0]?.delta?.content;
        if (content) {
          controller.enqueue(encoder.encode(`data: ${JSON.stringify({ content })}\n\n`));
        }
      }
      controller.enqueue(encoder.encode('data: [DONE]\n\n'));
      controller.close();
    },
  });

  return new Response(readable, {
    headers: {
      'Content-Type': 'text/event-stream',
      'Cache-Control': 'no-cache',
      'Connection': 'keep-alive',
    },
  });
}
```

## Fly.io Deployment

### fly.toml Configuration

```toml
app = "my-mistral-app"
primary_region = "iad"

[env]
  NODE_ENV = "production"
  LOG_LEVEL = "info"

[http_service]
  internal_port = 3000
  force_https = true
  auto_stop_machines = true
  auto_start_machines = true
  min_machines_running = 1

[http_service.concurrency]
  type = "connections"
  hard_limit = 100
  soft_limit = 80

[[services.http_checks]]
  interval = 10000
  timeout = 2000
  path = "/health"
```

### Deploy to Fly.io

```bash
# Set Mistral secrets
fly secrets set MISTRAL_API_KEY="your-api-key"

# Deploy
fly deploy

# View logs
fly logs

# Scale
fly scale count 2 --region iad
```

## Google Cloud Run

### Dockerfile

```dockerfile
FROM node:20-slim AS builder
WORKDIR /app
COPY package*.json ./
RUN npm ci --only=production
COPY . .
RUN npm run build

FROM node:20-slim
WORKDIR /app
COPY --from=builder /app/dist ./dist
COPY --from=builder /app/node_modules ./node_modules
COPY --from=builder /app/package.json ./

ENV NODE_ENV=production
EXPOSE 8080
CMD ["node", "dist/index.js"]
```

### Deploy Script

```bash
#!/bin/bash
# deploy-cloud-run.sh

PROJECT_ID="${GOOGLE_CLOUD_PROJECT}"
SERVICE_NAME="mistral-service"
REGION="us-central1"

# Create secret in Secret Manager
echo -n "your-api-key" | gcloud secrets create mistral-api-key --data-file=-

# Grant Cloud Run access to secret
gcloud secrets add-iam-policy-binding mistral-api-key \
  --member="serviceAccount:${PROJECT_ID}-compute@developer.gserviceaccount.com" \
  --role="roles/secretmanager.secretAccessor"

# Build and push image
gcloud builds submit --tag gcr.io/$PROJECT_ID/$SERVICE_NAME

# Deploy to Cloud Run
gcloud run deploy $SERVICE_NAME \
  --image gcr.io/$PROJECT_ID/$SERVICE_NAME \
  --region $REGION \
  --platform managed \
  --allow-unauthenticated \
  --set-secrets="MISTRAL_API_KEY=mistral-api-key:latest" \
  --memory 512Mi \
  --cpu 1 \
  --timeout 60s \
  --concurrency 80
```

## Environment Configuration Pattern

```typescript
// config/mistral.ts
interface MistralConfig {
  apiKey: string;
  model: string;
  timeout: number;
  maxRetries: number;
}

export function getMistralConfig(): MistralConfig {
  const env = process.env.NODE_ENV || 'development';

  return {
    apiKey: process.env.MISTRAL_API_KEY!,
    model: process.env.MISTRAL_MODEL || 'mistral-small-latest',
    timeout: parseInt(process.env.MISTRAL_TIMEOUT || '30000'),
    maxRetries: parseInt(process.env.MISTRAL_MAX_RETRIES || '3'),
  };
}
```

## Health Check Endpoint

```typescript
// api/health.ts
import Mistral from '@mistralai/mistralai';

export async function GET() {
  const start = Date.now();

  try {
    const client = new Mistral({ apiKey: process.env.MISTRAL_API_KEY! });
    await client.models.list();

    return Response.json({
      status: 'healthy',
      services: {
        mistral: {
          connected: true,
          latencyMs: Date.now() - start,
        },
      },
      timestamp: new Date().toISOString(),
    });
  } catch (error) {
    return Response.json({
      status: 'degraded',
      services: {
        mistral: {
          connected: false,
          latencyMs: Date.now() - start,
          error: 'Connection failed',
        },
      },
      timestamp: new Date().toISOString(),
    }, { status: 503 });
  }
}
```

## Instructions

### Step 1: Choose Deployment Platform
Select the platform that best fits your infrastructure needs.

### Step 2: Configure Secrets
Store Mistral API keys securely using the platform's secrets management.

### Step 3: Deploy Application
Use the platform CLI to deploy your application.

### Step 4: Verify Health
Test the health check endpoint to confirm Mistral connectivity.

## Output
- Application deployed to production
- Mistral secrets securely configured
- Health check endpoint functional
- Environment-specific configuration in place

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Secret not found | Missing configuration | Add secret via platform CLI |
| Deploy timeout | Large build | Increase build timeout |
| Health check fails | Wrong API key | Verify environment variable |
| Cold start slow | No warm-up | Configure minimum instances |

## Examples

### Quick Deploy Script

```bash
#!/bin/bash
# Platform-agnostic deploy helper
case "$1" in
  vercel)
    vercel env add MISTRAL_API_KEY production <<< "$MISTRAL_API_KEY"
    vercel --prod
    ;;
  fly)
    fly secrets set MISTRAL_API_KEY="$MISTRAL_API_KEY"
    fly deploy
    ;;
  cloudrun)
    gcloud run deploy mistral-app \
      --set-env-vars="MISTRAL_API_KEY=$MISTRAL_API_KEY"
    ;;
  *)
    echo "Usage: $0 {vercel|fly|cloudrun}"
    ;;
esac
```

## Resources
- [Vercel Documentation](https://vercel.com/docs)
- [Fly.io Documentation](https://fly.io/docs)
- [Cloud Run Documentation](https://cloud.google.com/run/docs)
- [Mistral AI API Reference](https://docs.mistral.ai/api/)

## Next Steps
For event handling, see `mistral-webhooks-events`.
