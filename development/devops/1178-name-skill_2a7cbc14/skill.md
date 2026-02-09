---
name: langfuse-deploy-integration
description: |
  Deploy Langfuse with your application across different platforms.
  Use when deploying Langfuse to Vercel, AWS, GCP, or Docker,
  or integrating Langfuse into your deployment pipeline.
  Trigger with phrases like "deploy langfuse", "langfuse Vercel",
  "langfuse AWS", "langfuse Docker", "langfuse production deploy".
allowed-tools: Read, Write, Edit, Bash(docker:*), Bash(vercel:*), Bash(gcloud:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Deploy Integration

## Overview
Deploy Langfuse observability with your application across various platforms.

## Prerequisites
- Langfuse production API keys
- Deployment platform access (Vercel, AWS, GCP, etc.)
- Application ready for deployment

## Instructions

### Step 1: Vercel Deployment

```bash
# Set environment variables
vercel env add LANGFUSE_PUBLIC_KEY production
vercel env add LANGFUSE_SECRET_KEY production
vercel env add LANGFUSE_HOST production
```

```typescript
// vercel.json
{
  "env": {
    "LANGFUSE_HOST": "https://cloud.langfuse.com"
  }
}
```

```typescript
// Next.js API route with Langfuse
// app/api/chat/route.ts
import { Langfuse } from "langfuse";
import { NextResponse } from "next/server";

// Initialize outside handler for connection reuse
const langfuse = new Langfuse({
  publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
  secretKey: process.env.LANGFUSE_SECRET_KEY!,
  flushAt: 1, // Flush immediately in serverless
});

export async function POST(request: Request) {
  const { message } = await request.json();

  const trace = langfuse.trace({
    name: "api/chat",
    input: { message },
    metadata: {
      platform: "vercel",
      region: process.env.VERCEL_REGION,
    },
  });

  try {
    const response = await processChat(message, trace);

    trace.update({ output: response });

    // Flush before response in serverless
    await langfuse.flushAsync();

    return NextResponse.json(response);
  } catch (error) {
    trace.update({ level: "ERROR", statusMessage: String(error) });
    await langfuse.flushAsync();
    throw error;
  }
}
```

### Step 2: AWS Lambda Deployment

```typescript
// handler.ts
import { Langfuse } from "langfuse";
import { APIGatewayEvent, Context } from "aws-lambda";

// Initialize outside handler for warm start reuse
const langfuse = new Langfuse({
  publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
  secretKey: process.env.LANGFUSE_SECRET_KEY!,
  flushAt: 1, // Immediate flush for Lambda
});

export async function handler(event: APIGatewayEvent, context: Context) {
  // Don't wait for event loop to empty
  context.callbackWaitsForEmptyEventLoop = false;

  const trace = langfuse.trace({
    name: "lambda/handler",
    input: JSON.parse(event.body || "{}"),
    metadata: {
      requestId: context.awsRequestId,
      functionName: context.functionName,
      coldStart: !global.isWarm,
    },
  });

  global.isWarm = true;

  try {
    const result = await processRequest(event, trace);

    trace.update({ output: result });

    // CRITICAL: Flush before Lambda freezes
    await langfuse.flushAsync();

    return {
      statusCode: 200,
      body: JSON.stringify(result),
    };
  } catch (error) {
    trace.update({ level: "ERROR", statusMessage: String(error) });
    await langfuse.flushAsync();

    return {
      statusCode: 500,
      body: JSON.stringify({ error: "Internal error" }),
    };
  }
}
```

```yaml
# serverless.yml
service: my-llm-app

provider:
  name: aws
  runtime: nodejs20.x
  environment:
    LANGFUSE_PUBLIC_KEY: ${ssm:/langfuse/public-key}
    LANGFUSE_SECRET_KEY: ${ssm:/langfuse/secret-key}
    LANGFUSE_HOST: https://cloud.langfuse.com

functions:
  chat:
    handler: handler.handler
    timeout: 30
    events:
      - http:
          path: /chat
          method: post
```

### Step 3: Google Cloud Run Deployment

```dockerfile
# Dockerfile
FROM node:20-slim

WORKDIR /app

COPY package*.json ./
RUN npm ci --production

COPY . .
RUN npm run build

# Cloud Run sets PORT
ENV PORT=8080

CMD ["node", "dist/server.js"]
```

```typescript
// server.ts
import express from "express";
import { Langfuse } from "langfuse";

const app = express();
const langfuse = new Langfuse({
  publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
  secretKey: process.env.LANGFUSE_SECRET_KEY!,
});

// Graceful shutdown for Cloud Run
process.on("SIGTERM", async () => {
  console.log("SIGTERM received, flushing Langfuse...");
  await langfuse.shutdownAsync();
  process.exit(0);
});

app.post("/chat", async (req, res) => {
  const trace = langfuse.trace({
    name: "cloudrun/chat",
    input: req.body,
    metadata: {
      service: process.env.K_SERVICE,
      revision: process.env.K_REVISION,
    },
  });

  // ... handler logic
});

const port = process.env.PORT || 8080;
app.listen(port, () => {
  console.log(`Server listening on port ${port}`);
});
```

```bash
# Deploy to Cloud Run
gcloud run deploy my-llm-app \
  --source . \
  --region us-central1 \
  --set-env-vars "LANGFUSE_PUBLIC_KEY=$LANGFUSE_PUBLIC_KEY" \
  --set-secrets "LANGFUSE_SECRET_KEY=langfuse-secret-key:latest" \
  --allow-unauthenticated
```

### Step 4: Docker Compose Deployment

```yaml
# docker-compose.yml
version: '3.8'

services:
  app:
    build: .
    ports:
      - "3000:3000"
    environment:
      - NODE_ENV=production
      - LANGFUSE_PUBLIC_KEY=${LANGFUSE_PUBLIC_KEY}
      - LANGFUSE_SECRET_KEY=${LANGFUSE_SECRET_KEY}
      - LANGFUSE_HOST=https://cloud.langfuse.com
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:3000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
    deploy:
      resources:
        limits:
          memory: 512M
    restart: unless-stopped

  # Optional: Local Langfuse for development
  langfuse:
    image: langfuse/langfuse:latest
    ports:
      - "3001:3000"
    environment:
      - DATABASE_URL=postgresql://postgres:postgres@db:5432/langfuse
      - NEXTAUTH_SECRET=${NEXTAUTH_SECRET}
      - NEXTAUTH_URL=http://localhost:3001
      - SALT=${LANGFUSE_SALT}
    depends_on:
      - db
    profiles:
      - local

  db:
    image: postgres:15-alpine
    environment:
      - POSTGRES_USER=postgres
      - POSTGRES_PASSWORD=postgres
      - POSTGRES_DB=langfuse
    volumes:
      - langfuse-db:/var/lib/postgresql/data
    profiles:
      - local

volumes:
  langfuse-db:
```

### Step 5: Kubernetes Deployment

```yaml
# k8s/deployment.yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: llm-app
spec:
  replicas: 3
  selector:
    matchLabels:
      app: llm-app
  template:
    metadata:
      labels:
        app: llm-app
    spec:
      containers:
        - name: app
          image: myregistry/llm-app:latest
          ports:
            - containerPort: 3000
          env:
            - name: LANGFUSE_PUBLIC_KEY
              valueFrom:
                secretKeyRef:
                  name: langfuse-secrets
                  key: public-key
            - name: LANGFUSE_SECRET_KEY
              valueFrom:
                secretKeyRef:
                  name: langfuse-secrets
                  key: secret-key
            - name: LANGFUSE_HOST
              value: "https://cloud.langfuse.com"
          resources:
            requests:
              memory: "256Mi"
              cpu: "250m"
            limits:
              memory: "512Mi"
              cpu: "500m"
          livenessProbe:
            httpGet:
              path: /health
              port: 3000
            initialDelaySeconds: 10
            periodSeconds: 30
          lifecycle:
            preStop:
              exec:
                command: ["/bin/sh", "-c", "sleep 10"]
      terminationGracePeriodSeconds: 30
---
apiVersion: v1
kind: Secret
metadata:
  name: langfuse-secrets
type: Opaque
stringData:
  public-key: pk-lf-...
  secret-key: sk-lf-...
```

## Output
- Platform-specific deployment configurations
- Serverless optimizations (immediate flush)
- Graceful shutdown handling
- Secret management integration
- Health checks configured

## Serverless Considerations

| Platform | Flush Strategy | Shutdown Handling |
|----------|----------------|-------------------|
| Vercel | `flushAt: 1` | Flush before response |
| Lambda | `flushAt: 1` | Flush before return |
| Cloud Run | Default batching | SIGTERM handler |
| Kubernetes | Default batching | preStop hook |

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Lost traces (serverless) | No flush before freeze | Add `flushAsync()` |
| Slow cold starts | Large SDK init | Move init outside handler |
| Secrets exposed | Env vars in logs | Use secret managers |
| Traces delayed | Default batching | Reduce `flushInterval` |

## Resources
- [Vercel Environment Variables](https://vercel.com/docs/environment-variables)
- [AWS Lambda Best Practices](https://docs.aws.amazon.com/lambda/latest/dg/best-practices.html)
- [Cloud Run Secrets](https://cloud.google.com/run/docs/configuring/secrets)
- [Langfuse Self-Hosting](https://langfuse.com/docs/deployment/self-host)

## Next Steps
For webhooks and events, see `langfuse-webhooks-events`.
