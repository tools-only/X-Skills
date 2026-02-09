---
name: twinmind-observability
description: |
  Set up comprehensive observability for TwinMind integrations with metrics, traces, and alerts.
  Use when implementing monitoring for TwinMind operations, setting up dashboards,
  or configuring alerting for meeting AI integration health.
  Trigger with phrases like "twinmind monitoring", "twinmind metrics",
  "twinmind observability", "monitor twinmind", "twinmind alerts", "twinmind tracing".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Observability

## Overview
Set up comprehensive observability for TwinMind integrations including metrics, distributed tracing, logging, and alerting.

## Prerequisites
- Prometheus or compatible metrics backend
- OpenTelemetry SDK installed
- Grafana or similar dashboarding tool
- AlertManager configured

## Metrics Collection

### Key Metrics

| Metric | Type | Description |
|--------|------|-------------|
| `twinmind_transcriptions_total` | Counter | Total transcription requests |
| `twinmind_transcription_duration_seconds` | Histogram | Transcription processing time |
| `twinmind_transcription_audio_hours` | Counter | Total audio hours processed |
| `twinmind_errors_total` | Counter | Error count by type |
| `twinmind_rate_limit_remaining` | Gauge | Rate limit headroom |
| `twinmind_summary_requests_total` | Counter | Summary generation requests |
| `twinmind_action_items_extracted` | Counter | Action items extracted |
| `twinmind_ai_tokens_used` | Counter | AI tokens consumed |

### Prometheus Metrics Implementation

```typescript
// src/observability/metrics.ts
import { Registry, Counter, Histogram, Gauge, Summary } from 'prom-client';

const registry = new Registry();

// Transcription metrics
export const transcriptionCounter = new Counter({
  name: 'twinmind_transcriptions_total',
  help: 'Total TwinMind transcription requests',
  labelNames: ['status', 'model', 'language'],
  registers: [registry],
});

export const transcriptionDuration = new Histogram({
  name: 'twinmind_transcription_duration_seconds',
  help: 'TwinMind transcription processing duration',
  labelNames: ['model'],
  buckets: [1, 5, 10, 30, 60, 120, 300, 600],  // Up to 10 minutes
  registers: [registry],
});

export const audioHoursProcessed = new Counter({
  name: 'twinmind_transcription_audio_hours',
  help: 'Total audio hours processed',
  labelNames: ['model'],
  registers: [registry],
});

// API metrics
export const apiRequestCounter = new Counter({
  name: 'twinmind_api_requests_total',
  help: 'Total TwinMind API requests',
  labelNames: ['method', 'endpoint', 'status'],
  registers: [registry],
});

export const apiLatency = new Histogram({
  name: 'twinmind_api_latency_seconds',
  help: 'TwinMind API request latency',
  labelNames: ['method', 'endpoint'],
  buckets: [0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10],
  registers: [registry],
});

// Error metrics
export const errorCounter = new Counter({
  name: 'twinmind_errors_total',
  help: 'TwinMind errors by type',
  labelNames: ['error_type', 'operation'],
  registers: [registry],
});

// Rate limit metrics
export const rateLimitRemaining = new Gauge({
  name: 'twinmind_rate_limit_remaining',
  help: 'Remaining rate limit quota',
  labelNames: ['endpoint'],
  registers: [registry],
});

// AI usage metrics
export const aiTokensUsed = new Counter({
  name: 'twinmind_ai_tokens_used',
  help: 'AI tokens consumed',
  labelNames: ['operation'],
  registers: [registry],
});

// Summary generation
export const summaryCounter = new Counter({
  name: 'twinmind_summary_requests_total',
  help: 'Summary generation requests',
  labelNames: ['status', 'format'],
  registers: [registry],
});

export const actionItemsExtracted = new Counter({
  name: 'twinmind_action_items_extracted',
  help: 'Action items extracted from meetings',
  registers: [registry],
});

export { registry };
```

### Instrumented Client

```typescript
// src/twinmind/instrumented-client.ts
import {
  transcriptionCounter,
  transcriptionDuration,
  audioHoursProcessed,
  apiRequestCounter,
  apiLatency,
  errorCounter,
  rateLimitRemaining,
} from '../observability/metrics';

export class InstrumentedTwinMindClient {
  private client: TwinMindClient;

  constructor(config: TwinMindConfig) {
    this.client = new TwinMindClient(config);
  }

  async transcribe(audioUrl: string, options?: TranscriptionOptions): Promise<Transcript> {
    const timer = transcriptionDuration.startTimer({ model: options?.model || 'ear-3' });

    try {
      const result = await this.client.transcribe(audioUrl, options);

      // Record success metrics
      transcriptionCounter.inc({
        status: 'success',
        model: options?.model || 'ear-3',
        language: result.language,
      });

      // Record audio hours
      audioHoursProcessed.inc(
        { model: options?.model || 'ear-3' },
        result.duration_seconds / 3600
      );

      return result;
    } catch (error: any) {
      transcriptionCounter.inc({
        status: 'error',
        model: options?.model || 'ear-3',
        language: 'unknown',
      });

      errorCounter.inc({
        error_type: error.code || 'unknown',
        operation: 'transcribe',
      });

      throw error;
    } finally {
      timer();
    }
  }

  async request<T>(method: string, endpoint: string, data?: any): Promise<T> {
    const timer = apiLatency.startTimer({ method, endpoint });

    try {
      const response = await this.client.request(method, endpoint, data);

      apiRequestCounter.inc({ method, endpoint, status: 'success' });

      // Update rate limit gauge from response headers
      if (response.headers?.['x-ratelimit-remaining']) {
        rateLimitRemaining.set(
          { endpoint },
          parseInt(response.headers['x-ratelimit-remaining'])
        );
      }

      return response.data;
    } catch (error: any) {
      apiRequestCounter.inc({
        method,
        endpoint,
        status: error.response?.status || 'error',
      });

      throw error;
    } finally {
      timer();
    }
  }
}
```

## Distributed Tracing

### OpenTelemetry Setup

```typescript
// src/observability/tracing.ts
import { NodeSDK } from '@opentelemetry/sdk-node';
import { getNodeAutoInstrumentations } from '@opentelemetry/auto-instrumentations-node';
import { OTLPTraceExporter } from '@opentelemetry/exporter-trace-otlp-http';
import { Resource } from '@opentelemetry/resources';
import { SemanticResourceAttributes } from '@opentelemetry/semantic-conventions';
import { trace, SpanStatusCode, SpanKind } from '@opentelemetry/api';

// Initialize SDK
const sdk = new NodeSDK({
  resource: new Resource({
    [SemanticResourceAttributes.SERVICE_NAME]: 'twinmind-integration',
    [SemanticResourceAttributes.SERVICE_VERSION]: process.env.npm_package_version,
    [SemanticResourceAttributes.DEPLOYMENT_ENVIRONMENT]: process.env.NODE_ENV,
  }),
  traceExporter: new OTLPTraceExporter({
    url: process.env.OTEL_EXPORTER_OTLP_ENDPOINT || 'http://localhost:4318/v1/traces',
  }),
  instrumentations: [getNodeAutoInstrumentations()],
});

sdk.start();

// Custom tracer for TwinMind operations
export const tracer = trace.getTracer('twinmind-client');

// Traced operation wrapper
export async function tracedOperation<T>(
  operationName: string,
  operation: () => Promise<T>,
  attributes?: Record<string, string | number>
): Promise<T> {
  return tracer.startActiveSpan(
    `twinmind.${operationName}`,
    { kind: SpanKind.CLIENT, attributes },
    async (span) => {
      try {
        const result = await operation();
        span.setStatus({ code: SpanStatusCode.OK });
        return result;
      } catch (error: any) {
        span.setStatus({
          code: SpanStatusCode.ERROR,
          message: error.message,
        });
        span.recordException(error);
        throw error;
      } finally {
        span.end();
      }
    }
  );
}

// Traced TwinMind client
export class TracedTwinMindClient {
  private client: TwinMindClient;

  async transcribe(audioUrl: string, options?: TranscriptionOptions): Promise<Transcript> {
    return tracedOperation(
      'transcribe',
      () => this.client.transcribe(audioUrl, options),
      {
        'twinmind.audio_url': audioUrl,
        'twinmind.model': options?.model || 'ear-3',
        'twinmind.diarization': options?.diarization ? 'true' : 'false',
      }
    );
  }

  async summarize(transcriptId: string): Promise<Summary> {
    return tracedOperation(
      'summarize',
      () => this.client.summarize(transcriptId),
      { 'twinmind.transcript_id': transcriptId }
    );
  }
}
```

## Structured Logging

```typescript
// src/observability/logging.ts
import pino from 'pino';

export const logger = pino({
  name: 'twinmind',
  level: process.env.LOG_LEVEL || 'info',
  formatters: {
    level: (label) => ({ level: label }),
  },
  redact: ['apiKey', 'authorization', 'password'],
});

// Operation-specific logging
export function logTwinMindOperation(
  operation: string,
  data: Record<string, any>,
  duration?: number
): void {
  logger.info({
    service: 'twinmind',
    operation,
    duration_ms: duration,
    ...data,
  });
}

// Error logging
export function logTwinMindError(
  operation: string,
  error: Error,
  context?: Record<string, any>
): void {
  logger.error({
    service: 'twinmind',
    operation,
    error: {
      name: error.name,
      message: error.message,
      stack: error.stack,
    },
    ...context,
  });
}
```

## Alert Configuration

### Prometheus AlertManager Rules

```yaml
# alerts/twinmind_alerts.yaml
groups:
  - name: twinmind_alerts
    rules:
      - alert: TwinMindHighErrorRate
        expr: |
          rate(twinmind_errors_total[5m]) /
          rate(twinmind_api_requests_total[5m]) > 0.05
        for: 5m
        labels:
          severity: warning
          service: twinmind
        annotations:
          summary: "TwinMind error rate > 5%"
          description: "Error rate is {{ $value | humanizePercentage }}"

      - alert: TwinMindHighLatency
        expr: |
          histogram_quantile(0.95,
            rate(twinmind_api_latency_seconds_bucket[5m])
          ) > 5
        for: 5m
        labels:
          severity: warning
          service: twinmind
        annotations:
          summary: "TwinMind P95 latency > 5s"
          description: "P95 latency is {{ $value | humanizeDuration }}"

      - alert: TwinMindTranscriptionFailures
        expr: |
          increase(twinmind_transcriptions_total{status="error"}[15m]) > 5
        for: 5m
        labels:
          severity: warning
          service: twinmind
        annotations:
          summary: "TwinMind transcription failures increasing"
          description: "{{ $value }} failures in last 15 minutes"

      - alert: TwinMindRateLimitApproaching
        expr: |
          twinmind_rate_limit_remaining < 10
        for: 1m
        labels:
          severity: warning
          service: twinmind
        annotations:
          summary: "TwinMind rate limit approaching"
          description: "Only {{ $value }} requests remaining"

      - alert: TwinMindRateLimitExceeded
        expr: |
          increase(twinmind_errors_total{error_type="RATE_LIMITED"}[5m]) > 0
        for: 1m
        labels:
          severity: critical
          service: twinmind
        annotations:
          summary: "TwinMind rate limit exceeded"
          description: "Rate limit errors detected"

      - alert: TwinMindAPIDown
        expr: |
          up{job="twinmind"} == 0
        for: 1m
        labels:
          severity: critical
          service: twinmind
        annotations:
          summary: "TwinMind integration is down"
          description: "Health check failing"

      - alert: TwinMindHighTokenUsage
        expr: |
          increase(twinmind_ai_tokens_used[24h]) > 1500000
        for: 5m
        labels:
          severity: warning
          service: twinmind
        annotations:
          summary: "High TwinMind token usage"
          description: "{{ $value }} tokens used in 24h (limit: 2M)"
```

## Grafana Dashboard

```json
{
  "dashboard": {
    "title": "TwinMind Integration",
    "panels": [
      {
        "title": "Transcription Rate",
        "type": "graph",
        "targets": [{
          "expr": "rate(twinmind_transcriptions_total[5m])",
          "legendFormat": "{{status}}"
        }]
      },
      {
        "title": "API Latency (P50/P95/P99)",
        "type": "graph",
        "targets": [
          {
            "expr": "histogram_quantile(0.5, rate(twinmind_api_latency_seconds_bucket[5m]))",
            "legendFormat": "P50"
          },
          {
            "expr": "histogram_quantile(0.95, rate(twinmind_api_latency_seconds_bucket[5m]))",
            "legendFormat": "P95"
          },
          {
            "expr": "histogram_quantile(0.99, rate(twinmind_api_latency_seconds_bucket[5m]))",
            "legendFormat": "P99"
          }
        ]
      },
      {
        "title": "Error Rate",
        "type": "graph",
        "targets": [{
          "expr": "rate(twinmind_errors_total[5m])",
          "legendFormat": "{{error_type}}"
        }]
      },
      {
        "title": "Rate Limit Remaining",
        "type": "gauge",
        "targets": [{
          "expr": "twinmind_rate_limit_remaining"
        }]
      },
      {
        "title": "Audio Hours Processed (24h)",
        "type": "stat",
        "targets": [{
          "expr": "increase(twinmind_transcription_audio_hours[24h])"
        }]
      },
      {
        "title": "AI Tokens Used (24h)",
        "type": "stat",
        "targets": [{
          "expr": "increase(twinmind_ai_tokens_used[24h])"
        }]
      }
    ]
  }
}
```

## Metrics Endpoint

```typescript
// src/api/metrics.ts
import express from 'express';
import { registry } from '../observability/metrics';

const router = express.Router();

router.get('/metrics', async (req, res) => {
  res.set('Content-Type', registry.contentType);
  res.send(await registry.metrics());
});

export default router;
```

## Output
- Prometheus metrics implementation
- Distributed tracing with OpenTelemetry
- Structured logging with Pino
- AlertManager rules
- Grafana dashboard configuration
- Metrics endpoint

## Error Handling

| Issue | Cause | Solution |
|-------|-------|----------|
| Missing metrics | No instrumentation | Wrap client calls |
| Trace gaps | Missing propagation | Check context headers |
| Alert storms | Wrong thresholds | Tune alert rules |
| High cardinality | Too many labels | Reduce label values |

## Resources
- [Prometheus Best Practices](https://prometheus.io/docs/practices/naming/)
- [OpenTelemetry Documentation](https://opentelemetry.io/docs/)
- [Grafana Dashboard Examples](https://grafana.com/grafana/dashboards/)

## Next Steps
For incident response, see `twinmind-incident-runbook`.
