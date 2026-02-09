---
name: openevidence-deploy-integration
description: |
  Deploy OpenEvidence integrations to healthcare production environments.
  Use when deploying to production, setting up staging environments,
  or configuring cloud deployments for clinical AI applications.
  Trigger with phrases like "deploy openevidence", "openevidence staging",
  "openevidence production deploy", "release openevidence".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Deploy Integration

## Overview
Deploy OpenEvidence clinical AI integrations to production healthcare environments with proper security, compliance, and rollback capabilities.

## Prerequisites
- CI pipeline configured (see `openevidence-ci-integration`)
- Production credentials from OpenEvidence
- Cloud platform configured (GCP, AWS, Azure)
- Signed BAA for production environment

## Instructions

### Step 1: Multi-Environment Configuration
```typescript
// src/config/environments.ts
interface OpenEvidenceEnvConfig {
  baseUrl: string;
  timeout: number;
  retries: number;
  secretPath: string;
  features: {
    deepConsult: boolean;
    webhooks: boolean;
    caching: boolean;
  };
}

export const environments: Record<string, OpenEvidenceEnvConfig> = {
  development: {
    baseUrl: 'https://api.sandbox.openevidence.com',
    timeout: 60000,
    retries: 1,
    secretPath: 'local', // Uses .env
    features: {
      deepConsult: true,
      webhooks: false,
      caching: false,
    },
  },
  staging: {
    baseUrl: 'https://api.sandbox.openevidence.com',
    timeout: 45000,
    retries: 2,
    secretPath: 'projects/staging/secrets/openevidence',
    features: {
      deepConsult: true,
      webhooks: true,
      caching: true,
    },
  },
  production: {
    baseUrl: 'https://api.openevidence.com',
    timeout: 30000,
    retries: 3,
    secretPath: 'projects/production/secrets/openevidence',
    features: {
      deepConsult: true,
      webhooks: true,
      caching: true,
    },
  },
};

export function getConfig(): OpenEvidenceEnvConfig {
  const env = process.env.NODE_ENV || 'development';
  return environments[env] || environments.development;
}
```

### Step 2: GitHub Actions Deployment Workflow
```yaml
# .github/workflows/deploy.yml
name: Deploy OpenEvidence Integration

on:
  push:
    branches: [main]
    tags: ['v*']
  workflow_dispatch:
    inputs:
      environment:
        description: 'Deployment environment'
        required: true
        default: 'staging'
        type: choice
        options:
          - staging
          - production

env:
  PROJECT_ID: ${{ secrets.GCP_PROJECT_ID }}
  SERVICE_NAME: clinical-evidence-api
  REGION: us-central1

jobs:
  test:
    uses: ./.github/workflows/openevidence-ci.yml
    secrets: inherit

  build:
    needs: test
    runs-on: ubuntu-latest
    outputs:
      image: ${{ steps.build.outputs.image }}
    steps:
      - uses: actions/checkout@v4

      - name: Setup Cloud SDK
        uses: google-github-actions/setup-gcloud@v2
        with:
          project_id: ${{ env.PROJECT_ID }}

      - name: Configure Docker
        run: gcloud auth configure-docker ${{ env.REGION }}-docker.pkg.dev

      - name: Build and push
        id: build
        run: |
          IMAGE="${{ env.REGION }}-docker.pkg.dev/${{ env.PROJECT_ID }}/images/${{ env.SERVICE_NAME }}:${{ github.sha }}"
          docker build -t $IMAGE .
          docker push $IMAGE
          echo "image=$IMAGE" >> $GITHUB_OUTPUT

  deploy-staging:
    needs: build
    runs-on: ubuntu-latest
    environment: staging
    if: github.ref == 'refs/heads/main'
    steps:
      - uses: actions/checkout@v4

      - name: Setup Cloud SDK
        uses: google-github-actions/setup-gcloud@v2
        with:
          project_id: ${{ env.PROJECT_ID }}

      - name: Deploy to Cloud Run (Staging)
        run: |
          gcloud run deploy ${{ env.SERVICE_NAME }}-staging \
            --image ${{ needs.build.outputs.image }} \
            --region ${{ env.REGION }} \
            --set-env-vars NODE_ENV=staging \
            --set-secrets OPENEVIDENCE_API_KEY=openevidence-staging-api-key:latest \
            --set-secrets OPENEVIDENCE_ORG_ID=openevidence-staging-org-id:latest \
            --min-instances 1 \
            --max-instances 10 \
            --memory 512Mi \
            --cpu 1 \
            --timeout 60s \
            --allow-unauthenticated=false

      - name: Run smoke tests
        env:
          STAGING_URL: ${{ secrets.STAGING_URL }}
          STAGING_TOKEN: ${{ secrets.STAGING_TOKEN }}
        run: |
          curl -sf -H "Authorization: Bearer $STAGING_TOKEN" \
            "$STAGING_URL/health/openevidence" | jq .

  deploy-production:
    needs: [build, deploy-staging]
    runs-on: ubuntu-latest
    environment: production
    if: startsWith(github.ref, 'refs/tags/v')
    steps:
      - uses: actions/checkout@v4

      - name: Setup Cloud SDK
        uses: google-github-actions/setup-gcloud@v2
        with:
          project_id: ${{ env.PROJECT_ID }}

      - name: Deploy to Cloud Run (Production) - Canary
        run: |
          # Deploy new revision with no traffic
          gcloud run deploy ${{ env.SERVICE_NAME }} \
            --image ${{ needs.build.outputs.image }} \
            --region ${{ env.REGION }} \
            --set-env-vars NODE_ENV=production \
            --set-secrets OPENEVIDENCE_API_KEY=openevidence-prod-api-key:latest \
            --set-secrets OPENEVIDENCE_ORG_ID=openevidence-prod-org-id:latest \
            --min-instances 2 \
            --max-instances 50 \
            --memory 1Gi \
            --cpu 2 \
            --timeout 30s \
            --no-traffic

      - name: Canary rollout - 10%
        run: |
          gcloud run services update-traffic ${{ env.SERVICE_NAME }} \
            --region ${{ env.REGION }} \
            --to-revisions LATEST=10

      - name: Monitor canary (5 minutes)
        run: |
          echo "Monitoring canary deployment..."
          sleep 300

          # Check error rate
          ERROR_RATE=$(gcloud logging read \
            "resource.type=cloud_run_revision severity>=ERROR" \
            --limit 100 --format json | jq length)

          if [ "$ERROR_RATE" -gt 5 ]; then
            echo "High error rate detected: $ERROR_RATE"
            exit 1
          fi

      - name: Full rollout
        run: |
          gcloud run services update-traffic ${{ env.SERVICE_NAME }} \
            --region ${{ env.REGION }} \
            --to-latest

      - name: Create release
        uses: softprops/action-gh-release@v1
        with:
          generate_release_notes: true
```

### Step 3: Kubernetes Deployment
```yaml
# k8s/openevidence-deployment.yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: clinical-evidence-api
  labels:
    app: clinical-evidence-api
spec:
  replicas: 3
  selector:
    matchLabels:
      app: clinical-evidence-api
  template:
    metadata:
      labels:
        app: clinical-evidence-api
    spec:
      serviceAccountName: clinical-evidence-api
      containers:
        - name: api
          image: gcr.io/PROJECT_ID/clinical-evidence-api:VERSION
          ports:
            - containerPort: 8080
          env:
            - name: NODE_ENV
              value: production
            - name: OPENEVIDENCE_API_KEY
              valueFrom:
                secretKeyRef:
                  name: openevidence-secrets
                  key: api-key
            - name: OPENEVIDENCE_ORG_ID
              valueFrom:
                secretKeyRef:
                  name: openevidence-secrets
                  key: org-id
          resources:
            requests:
              memory: "512Mi"
              cpu: "500m"
            limits:
              memory: "1Gi"
              cpu: "1000m"
          livenessProbe:
            httpGet:
              path: /health
              port: 8080
            initialDelaySeconds: 10
            periodSeconds: 30
          readinessProbe:
            httpGet:
              path: /health/openevidence
              port: 8080
            initialDelaySeconds: 5
            periodSeconds: 10
---
apiVersion: v1
kind: Service
metadata:
  name: clinical-evidence-api
spec:
  selector:
    app: clinical-evidence-api
  ports:
    - port: 80
      targetPort: 8080
  type: ClusterIP
---
apiVersion: networking.k8s.io/v1
kind: Ingress
metadata:
  name: clinical-evidence-api
  annotations:
    kubernetes.io/ingress.class: nginx
    cert-manager.io/cluster-issuer: letsencrypt-prod
spec:
  tls:
    - hosts:
        - api.yourhealthcare.com
      secretName: api-tls
  rules:
    - host: api.yourhealthcare.com
      http:
        paths:
          - path: /clinical
            pathType: Prefix
            backend:
              service:
                name: clinical-evidence-api
                port:
                  number: 80
```

### Step 4: Health Check Endpoint
```typescript
// src/routes/health.ts
import { Router } from 'express';
import { OpenEvidenceClient } from '@openevidence/sdk';
import { getConfig } from '../config/environments';

const router = Router();

router.get('/health', (req, res) => {
  res.json({
    status: 'healthy',
    timestamp: new Date().toISOString(),
    version: process.env.npm_package_version,
  });
});

router.get('/health/openevidence', async (req, res) => {
  const startTime = Date.now();

  try {
    const client = new OpenEvidenceClient({
      apiKey: process.env.OPENEVIDENCE_API_KEY!,
      orgId: process.env.OPENEVIDENCE_ORG_ID!,
      baseUrl: getConfig().baseUrl,
    });

    await client.health.check();

    res.json({
      status: 'healthy',
      latencyMs: Date.now() - startTime,
      environment: process.env.NODE_ENV,
      baseUrl: getConfig().baseUrl,
    });
  } catch (error: any) {
    res.status(503).json({
      status: 'unhealthy',
      latencyMs: Date.now() - startTime,
      error: error.message,
    });
  }
});

export default router;
```

### Step 5: Rollback Procedure
```bash
#!/bin/bash
# scripts/rollback.sh

set -e

SERVICE_NAME="clinical-evidence-api"
REGION="us-central1"

echo "=== OpenEvidence Integration Rollback ==="

# List recent revisions
echo "Recent revisions:"
gcloud run revisions list --service $SERVICE_NAME --region $REGION --limit 5

# Get previous revision
CURRENT=$(gcloud run services describe $SERVICE_NAME --region $REGION \
  --format 'value(status.traffic[0].revisionName)')
PREVIOUS=$(gcloud run revisions list --service $SERVICE_NAME --region $REGION \
  --format 'value(metadata.name)' --limit 2 | tail -1)

echo "Current: $CURRENT"
echo "Rolling back to: $PREVIOUS"

read -p "Proceed with rollback? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
  gcloud run services update-traffic $SERVICE_NAME \
    --region $REGION \
    --to-revisions $PREVIOUS=100

  echo "Rollback complete. Monitoring..."
  sleep 30

  # Verify health
  gcloud run services describe $SERVICE_NAME --region $REGION \
    --format 'value(status.url)' | xargs -I {} curl -sf {}/health/openevidence
fi
```

## Output
- Multi-environment deployment configuration
- CI/CD pipeline with canary releases
- Kubernetes deployment manifests
- Health check endpoints
- Rollback procedures

## Deployment Checklist
- [ ] Staging environment tested
- [ ] Production secrets configured
- [ ] Health checks passing
- [ ] Monitoring alerts configured
- [ ] Rollback procedure tested
- [ ] On-call notified of deployment
- [ ] Clinical team notified of changes

## Error Handling
| Deployment Issue | Detection | Resolution |
|------------------|-----------|------------|
| Health check fails | Readiness probe | Rollback to previous revision |
| High error rate | Monitoring alert | Rollback, investigate logs |
| Secret missing | Container fails to start | Check secret manager |
| Rate limit hit | 429 errors spike | Scale down, check quotas |

## Resources
- [Cloud Run Documentation](https://cloud.google.com/run/docs)
- [Kubernetes Deployments](https://kubernetes.io/docs/concepts/workloads/controllers/deployment/)

## Next Steps
For webhook configuration, see `openevidence-webhooks-events`.
