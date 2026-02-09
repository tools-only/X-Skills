---
name: openevidence-observability
description: |
  Set up comprehensive observability for OpenEvidence integrations with metrics, traces, and alerts.
  Use when implementing monitoring for clinical AI operations, setting up dashboards,
  or configuring alerting for healthcare application health.
  Trigger with phrases like "openevidence monitoring", "openevidence metrics",
  "openevidence observability", "monitor openevidence", "openevidence alerts".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Observability

## Overview
Set up comprehensive observability for OpenEvidence clinical AI integrations in healthcare environments.

## Prerequisites
- Prometheus or compatible metrics backend
- OpenTelemetry SDK installed
- Grafana or similar dashboarding tool
- AlertManager or PagerDuty configured

## Key Metrics

| Metric | Type | Description | Alert Threshold |
|--------|------|-------------|-----------------|
| `openevidence_requests_total` | Counter | Total API requests | N/A |
| `openevidence_request_duration_seconds` | Histogram | Request latency | P95 > 15s |
| `openevidence_errors_total` | Counter | Error count by type | > 5% rate |
| `openevidence_cache_hits_total` | Counter | Cache hit count | < 50% rate |
| `openevidence_rate_limit_remaining` | Gauge | Rate limit headroom | < 10% |
| `openevidence_deepconsult_active` | Gauge | Active DeepConsults | > 50 |

## Instructions

### Step 1: Prometheus Metrics
```typescript
// src/monitoring/metrics.ts
import { Registry, Counter, Histogram, Gauge, Summary } from 'prom-client';

const registry = new Registry();

// Request metrics
export const requestCounter = new Counter({
  name: 'openevidence_requests_total',
  help: 'Total OpenEvidence API requests',
  labelNames: ['method', 'specialty', 'status'],
  registers: [registry],
});

export const requestDuration = new Histogram({
  name: 'openevidence_request_duration_seconds',
  help: 'OpenEvidence request duration',
  labelNames: ['method', 'specialty'],
  buckets: [0.5, 1, 2, 5, 10, 15, 30, 60],
  registers: [registry],
});

// Error metrics
export const errorCounter = new Counter({
  name: 'openevidence_errors_total',
  help: 'OpenEvidence errors by type',
  labelNames: ['error_type', 'specialty'],
  registers: [registry],
});

// Cache metrics
export const cacheHits = new Counter({
  name: 'openevidence_cache_hits_total',
  help: 'Cache hit count',
  registers: [registry],
});

export const cacheMisses = new Counter({
  name: 'openevidence_cache_misses_total',
  help: 'Cache miss count',
  registers: [registry],
});

// Rate limit metrics
export const rateLimitRemaining = new Gauge({
  name: 'openevidence_rate_limit_remaining',
  help: 'Remaining rate limit',
  registers: [registry],
});

export const rateLimitPercent = new Gauge({
  name: 'openevidence_rate_limit_usage_percent',
  help: 'Rate limit usage percentage',
  registers: [registry],
});

// DeepConsult metrics
export const deepConsultActive = new Gauge({
  name: 'openevidence_deepconsult_active',
  help: 'Active DeepConsult requests',
  registers: [registry],
});

export const deepConsultDuration = new Histogram({
  name: 'openevidence_deepconsult_duration_seconds',
  help: 'DeepConsult completion time',
  labelNames: ['specialty', 'status'],
  buckets: [30, 60, 120, 180, 300, 600],
  registers: [registry],
});

// Confidence metrics
export const confidenceScore = new Summary({
  name: 'openevidence_confidence_score',
  help: 'Response confidence score distribution',
  labelNames: ['specialty'],
  percentiles: [0.5, 0.9, 0.95, 0.99],
  registers: [registry],
});

export { registry };
```

### Step 2: Instrumented Client Wrapper
```typescript
// src/monitoring/instrumented-client.ts
import { OpenEvidenceClient } from '@openevidence/sdk';
import {
  requestCounter,
  requestDuration,
  errorCounter,
  cacheHits,
  cacheMisses,
  rateLimitRemaining,
  rateLimitPercent,
  confidenceScore,
} from './metrics';
import { trace, SpanStatusCode, Span } from '@opentelemetry/api';

const tracer = trace.getTracer('openevidence-client');

export class InstrumentedOpenEvidenceClient {
  private client: OpenEvidenceClient;
  private cache: ClinicalQueryCache;

  constructor(config: any, cache: ClinicalQueryCache) {
    this.client = new OpenEvidenceClient(config);
    this.cache = cache;
  }

  async query(request: ClinicalQueryRequest): Promise<ClinicalQueryResponse> {
    const { specialty } = request.context;

    return tracer.startActiveSpan('openevidence.query', async (span: Span) => {
      const timer = requestDuration.startTimer({ method: 'query', specialty });

      try {
        // Add span attributes
        span.setAttribute('openevidence.specialty', specialty);
        span.setAttribute('openevidence.urgency', request.context.urgency);

        // Check cache first
        const cached = await this.cache.get(request.question, request.context);
        if (cached) {
          cacheHits.inc();
          span.setAttribute('cache.hit', true);
          timer({ cached: 'true' });
          requestCounter.inc({ method: 'query', specialty, status: 'cache_hit' });
          return cached;
        }
        cacheMisses.inc();
        span.setAttribute('cache.hit', false);

        // Make API request
        const response = await this.client.query(request);

        // Update rate limit metrics from headers
        this.updateRateLimitMetrics(response.headers);

        // Record confidence score
        confidenceScore.observe({ specialty }, response.confidence);

        // Cache response
        await this.cache.set(request.question, request.context, response);

        // Record success
        timer({ cached: 'false' });
        requestCounter.inc({ method: 'query', specialty, status: 'success' });
        span.setStatus({ code: SpanStatusCode.OK });

        return response;
      } catch (error: any) {
        // Record error
        timer({ cached: 'false' });
        requestCounter.inc({ method: 'query', specialty, status: 'error' });
        errorCounter.inc({ error_type: error.code || 'unknown', specialty });

        span.setStatus({ code: SpanStatusCode.ERROR, message: error.message });
        span.recordException(error);

        throw error;
      } finally {
        span.end();
      }
    });
  }

  private updateRateLimitMetrics(headers: any): void {
    const remaining = parseInt(headers['x-ratelimit-remaining'] || '0');
    const limit = parseInt(headers['x-ratelimit-limit'] || '100');

    rateLimitRemaining.set(remaining);
    rateLimitPercent.set(((limit - remaining) / limit) * 100);
  }
}
```

### Step 3: Distributed Tracing Setup
```typescript
// src/monitoring/tracing.ts
import { NodeSDK } from '@opentelemetry/sdk-node';
import { Resource } from '@opentelemetry/resources';
import { SemanticResourceAttributes } from '@opentelemetry/semantic-conventions';
import { TraceExporter } from '@google-cloud/opentelemetry-cloud-trace-exporter';
import { HttpInstrumentation } from '@opentelemetry/instrumentation-http';
import { ExpressInstrumentation } from '@opentelemetry/instrumentation-express';

export function initTracing(): void {
  const sdk = new NodeSDK({
    resource: new Resource({
      [SemanticResourceAttributes.SERVICE_NAME]: 'clinical-evidence-api',
      [SemanticResourceAttributes.SERVICE_VERSION]: process.env.npm_package_version,
      [SemanticResourceAttributes.DEPLOYMENT_ENVIRONMENT]: process.env.NODE_ENV,
    }),
    traceExporter: new TraceExporter(),
    instrumentations: [
      new HttpInstrumentation({
        requestHook: (span, request) => {
          // Don't trace health checks
          if (request.url?.includes('/health')) {
            span.setAttribute('http.skip', true);
          }
        },
      }),
      new ExpressInstrumentation(),
    ],
  });

  sdk.start();

  process.on('SIGTERM', () => {
    sdk.shutdown().catch(console.error);
  });
}
```

### Step 4: Structured Logging
```typescript
// src/monitoring/logger.ts
import pino from 'pino';

export const logger = pino({
  name: 'clinical-evidence-api',
  level: process.env.LOG_LEVEL || 'info',
  formatters: {
    level: (label) => ({ level: label }),
  },
  // Redact PHI fields
  redact: {
    paths: ['patient.*', 'patientId', 'mrn', '*.ssn'],
    censor: '[REDACTED]',
  },
});

// OpenEvidence-specific logger
export const oeLogger = logger.child({ service: 'openevidence' });

export function logClinicalQuery(
  queryId: string,
  specialty: string,
  duration: number,
  cached: boolean,
  confidence?: number
): void {
  oeLogger.info({
    event: 'clinical_query',
    queryId,
    specialty,
    durationMs: duration,
    cached,
    confidence,
  });
}

export function logDeepConsult(
  consultId: string,
  specialty: string,
  status: string,
  duration?: number
): void {
  oeLogger.info({
    event: 'deep_consult',
    consultId,
    specialty,
    status,
    durationMs: duration,
  });
}

export function logError(
  error: Error,
  context: Record<string, any>
): void {
  oeLogger.error({
    event: 'error',
    errorType: error.name,
    errorMessage: error.message,
    ...context,
  });
}
```

### Step 5: Alert Rules
```yaml
# prometheus/openevidence_alerts.yaml
groups:
  - name: openevidence_alerts
    rules:
      # High error rate
      - alert: OpenEvidenceHighErrorRate
        expr: |
          rate(openevidence_errors_total[5m]) /
          rate(openevidence_requests_total[5m]) > 0.05
        for: 5m
        labels:
          severity: warning
          service: openevidence
        annotations:
          summary: "OpenEvidence error rate > 5%"
          description: "Error rate is {{ $value | humanizePercentage }}"

      # Critical error rate
      - alert: OpenEvidenceCriticalErrorRate
        expr: |
          rate(openevidence_errors_total[5m]) /
          rate(openevidence_requests_total[5m]) > 0.20
        for: 2m
        labels:
          severity: critical
          service: openevidence
        annotations:
          summary: "OpenEvidence error rate > 20%"
          description: "Critical error rate - possible outage"

      # High latency
      - alert: OpenEvidenceHighLatency
        expr: |
          histogram_quantile(0.95,
            rate(openevidence_request_duration_seconds_bucket[5m])
          ) > 15
        for: 5m
        labels:
          severity: warning
          service: openevidence
        annotations:
          summary: "OpenEvidence P95 latency > 15s"
          description: "95th percentile latency is {{ $value | humanizeDuration }}"

      # Low cache hit rate
      - alert: OpenEvidenceLowCacheHitRate
        expr: |
          rate(openevidence_cache_hits_total[1h]) /
          (rate(openevidence_cache_hits_total[1h]) + rate(openevidence_cache_misses_total[1h])) < 0.5
        for: 30m
        labels:
          severity: warning
          service: openevidence
        annotations:
          summary: "OpenEvidence cache hit rate < 50%"
          description: "Consider reviewing cache configuration"

      # Rate limit warning
      - alert: OpenEvidenceRateLimitWarning
        expr: openevidence_rate_limit_remaining < 10
        for: 1m
        labels:
          severity: warning
          service: openevidence
        annotations:
          summary: "OpenEvidence rate limit < 10 remaining"
          description: "Consider throttling requests"

      # Service down
      - alert: OpenEvidenceDown
        expr: up{job="openevidence"} == 0
        for: 1m
        labels:
          severity: critical
          service: openevidence
        annotations:
          summary: "OpenEvidence integration is down"
          description: "No successful health checks for 1 minute"
```

### Step 6: Grafana Dashboard
```json
{
  "dashboard": {
    "title": "OpenEvidence Clinical AI",
    "panels": [
      {
        "title": "Request Rate",
        "type": "stat",
        "targets": [{
          "expr": "sum(rate(openevidence_requests_total[5m]))",
          "legendFormat": "req/s"
        }]
      },
      {
        "title": "Error Rate",
        "type": "gauge",
        "targets": [{
          "expr": "sum(rate(openevidence_errors_total[5m])) / sum(rate(openevidence_requests_total[5m])) * 100"
        }],
        "fieldConfig": {
          "defaults": {
            "max": 100,
            "thresholds": {
              "mode": "absolute",
              "steps": [
                {"color": "green", "value": 0},
                {"color": "yellow", "value": 5},
                {"color": "red", "value": 20}
              ]
            }
          }
        }
      },
      {
        "title": "Latency Distribution",
        "type": "heatmap",
        "targets": [{
          "expr": "sum(rate(openevidence_request_duration_seconds_bucket[5m])) by (le)"
        }]
      },
      {
        "title": "Cache Hit Rate",
        "type": "timeseries",
        "targets": [{
          "expr": "sum(rate(openevidence_cache_hits_total[5m])) / (sum(rate(openevidence_cache_hits_total[5m])) + sum(rate(openevidence_cache_misses_total[5m])))",
          "legendFormat": "Hit Rate"
        }]
      },
      {
        "title": "Rate Limit Usage",
        "type": "gauge",
        "targets": [{
          "expr": "openevidence_rate_limit_usage_percent"
        }],
        "fieldConfig": {
          "defaults": {
            "max": 100,
            "thresholds": {
              "mode": "absolute",
              "steps": [
                {"color": "green", "value": 0},
                {"color": "yellow", "value": 70},
                {"color": "red", "value": 90}
              ]
            }
          }
        }
      },
      {
        "title": "Confidence Score Distribution",
        "type": "histogram",
        "targets": [{
          "expr": "openevidence_confidence_score"
        }]
      }
    ]
  }
}
```

## Metrics Endpoint
```typescript
// src/routes/metrics.ts
import { Router } from 'express';
import { registry } from '../monitoring/metrics';

const router = Router();

router.get('/metrics', async (req, res) => {
  res.set('Content-Type', registry.contentType);
  res.send(await registry.metrics());
});

export default router;
```

## Output
- Prometheus metrics collection
- Distributed tracing with OpenTelemetry
- Structured logging with PHI redaction
- Alert rules for critical conditions
- Grafana dashboard

## Observability Checklist
- [ ] Metrics endpoint exposed
- [ ] Request counters and histograms
- [ ] Error counters with labels
- [ ] Cache hit/miss tracking
- [ ] Rate limit monitoring
- [ ] Distributed tracing enabled
- [ ] Structured logging configured
- [ ] Alert rules deployed
- [ ] Dashboard created

## Resources
- [Prometheus Best Practices](https://prometheus.io/docs/practices/naming/)
- [OpenTelemetry Documentation](https://opentelemetry.io/docs/)
- [Grafana Dashboards](https://grafana.com/docs/grafana/latest/dashboards/)

## Next Steps
For incident response, see `openevidence-incident-runbook`.
