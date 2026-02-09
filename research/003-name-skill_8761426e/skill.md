---
name: openevidence-reference-architecture
description: |
  Implement OpenEvidence reference architecture with best-practice project layout.
  Use when designing new clinical AI integrations, reviewing project structure,
  or establishing architecture standards for healthcare applications.
  Trigger with phrases like "openevidence architecture", "openevidence best practices",
  "openevidence project structure", "clinical ai architecture".
allowed-tools: Read, Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Reference Architecture

## Overview
Production-ready architecture patterns for OpenEvidence clinical AI integrations in healthcare environments.

## Prerequisites
- Understanding of layered architecture
- OpenEvidence SDK knowledge
- Healthcare compliance requirements
- TypeScript/Node.js project setup

## Project Structure

```
clinical-evidence-api/
├── src/
│   ├── openevidence/
│   │   ├── client.ts           # Singleton client wrapper
│   │   ├── config.ts           # Environment configuration
│   │   ├── types.ts            # TypeScript types
│   │   ├── errors.ts           # Custom error classes
│   │   ├── cache.ts            # Caching layer
│   │   └── handlers/
│   │       ├── webhooks.ts     # Webhook handlers
│   │       └── events.ts       # Event processing
│   ├── services/
│   │   └── clinical/
│   │       ├── index.ts        # Service facade
│   │       ├── query.ts        # Clinical query service
│   │       ├── deepconsult.ts  # DeepConsult service
│   │       └── drug-info.ts    # Drug information service
│   ├── api/
│   │   ├── routes/
│   │   │   ├── clinical.ts     # Clinical query endpoints
│   │   │   └── webhooks.ts     # Webhook endpoints
│   │   └── middleware/
│   │       ├── auth.ts         # Authentication
│   │       ├── audit.ts        # HIPAA audit logging
│   │       └── rate-limit.ts   # Rate limiting
│   ├── integrations/
│   │   ├── ehr/
│   │   │   ├── fhir.ts         # FHIR integration
│   │   │   ├── epic.ts         # Epic EHR hooks
│   │   │   └── cerner.ts       # Cerner integration
│   │   └── notifications/
│   │       ├── email.ts        # Email notifications
│   │       └── push.ts         # Push notifications
│   ├── compliance/
│   │   ├── hipaa/
│   │   │   ├── audit-log.ts    # HIPAA audit logging
│   │   │   ├── phi-handler.ts  # PHI sanitization
│   │   │   └── encryption.ts   # Data encryption
│   │   └── retention.ts        # Data retention policies
│   └── jobs/
│       ├── cleanup.ts          # Data cleanup job
│       └── reporting.ts        # Usage reporting job
├── tests/
│   ├── unit/
│   │   └── services/
│   ├── integration/
│   │   └── openevidence/
│   └── clinical-validation/
│       └── known-answers/
├── config/
│   ├── default.json
│   ├── development.json
│   ├── staging.json
│   └── production.json
├── docs/
│   ├── architecture.md
│   ├── runbook.md
│   └── api-reference.md
└── scripts/
    ├── migrate.sh
    └── deploy.sh
```

## Layer Architecture

```
┌─────────────────────────────────────────────────────┐
│                  API Layer                           │
│   (Controllers, Routes, Webhooks, Middleware)        │
├─────────────────────────────────────────────────────┤
│                Service Layer                         │
│    (Business Logic, Orchestration, Validation)       │
├─────────────────────────────────────────────────────┤
│              OpenEvidence Layer                      │
│     (Client, Types, Caching, Error Handling)         │
├─────────────────────────────────────────────────────┤
│             Compliance Layer                         │
│       (HIPAA Audit, PHI Handling, Encryption)        │
├─────────────────────────────────────────────────────┤
│            Integration Layer                         │
│    (EHR/FHIR, Notifications, External Services)      │
├─────────────────────────────────────────────────────┤
│            Infrastructure Layer                      │
│      (Database, Cache, Queue, Monitoring)            │
└─────────────────────────────────────────────────────┘
```

## Key Components

### Step 1: Client Wrapper with Caching & Monitoring
```typescript
// src/openevidence/client.ts
import { OpenEvidenceClient } from '@openevidence/sdk';
import { ClinicalQueryCache } from './cache';
import { MetricsCollector } from '../monitoring/metrics';
import { HIPAAAuditLogger } from '../compliance/hipaa/audit-log';

export interface OpenEvidenceServiceConfig {
  apiKey: string;
  orgId: string;
  baseUrl: string;
  timeout: number;
  cache: ClinicalQueryCache;
  metrics: MetricsCollector;
  auditLogger: HIPAAAuditLogger;
}

export class OpenEvidenceService {
  private client: OpenEvidenceClient;
  private cache: ClinicalQueryCache;
  private metrics: MetricsCollector;
  private auditLogger: HIPAAAuditLogger;

  constructor(config: OpenEvidenceServiceConfig) {
    this.client = new OpenEvidenceClient({
      apiKey: config.apiKey,
      orgId: config.orgId,
      baseUrl: config.baseUrl,
      timeout: config.timeout,
    });
    this.cache = config.cache;
    this.metrics = config.metrics;
    this.auditLogger = config.auditLogger;
  }

  async query(
    request: ClinicalQueryRequest,
    context: RequestContext
  ): Promise<ClinicalQueryResponse> {
    const timer = this.metrics.startTimer('clinical_query');

    try {
      // Check cache
      const cached = await this.cache.get(request.question, request.context);
      if (cached) {
        this.metrics.incrementCounter('cache_hit');
        timer.stop({ cached: 'true' });
        return cached;
      }
      this.metrics.incrementCounter('cache_miss');

      // Query API
      const response = await this.client.query(request);

      // Cache response
      await this.cache.set(request.question, request.context, response);

      // Audit log (without PHI)
      await this.auditLogger.logQuery(context.userId, context.userRole, response.id, true);

      timer.stop({ cached: 'false' });
      return response;
    } catch (error: any) {
      this.metrics.incrementCounter('query_error', { error: error.code });
      await this.auditLogger.logQuery(context.userId, context.userRole, null, false);
      timer.stop({ error: 'true' });
      throw error;
    }
  }

  async deepConsult(
    request: DeepConsultRequest,
    context: RequestContext
  ): Promise<string> {
    // Return consultId for async tracking
    const consultId = await this.client.deepConsult.create(request);
    await this.auditLogger.logDeepConsult(context.userId, context.userRole, consultId, true);
    return consultId;
  }

  async health(): Promise<HealthStatus> {
    try {
      const start = Date.now();
      await this.client.health.check();
      return {
        status: 'healthy',
        latencyMs: Date.now() - start,
      };
    } catch (error: any) {
      return { status: 'unhealthy', error: error.message };
    }
  }
}
```

### Step 2: Service Facade
```typescript
// src/services/clinical/index.ts
import { OpenEvidenceService } from '../../openevidence/client';
import { PHIHandler } from '../../compliance/hipaa/phi-handler';

export class ClinicalEvidenceService {
  constructor(
    private openEvidence: OpenEvidenceService,
    private phiHandler: PHIHandler
  ) {}

  async queryClinicalEvidence(
    question: string,
    patientContext: PatientContext | undefined,
    context: RequestContext
  ): Promise<FormattedClinicalAnswer> {
    // 1. Sanitize input (remove PHI)
    const sanitizedRequest = this.phiHandler.sanitizeQuery(question, patientContext);

    // 2. Query OpenEvidence
    const response = await this.openEvidence.query(sanitizedRequest, context);

    // 3. Format response for clinical use
    return this.formatClinicalAnswer(response);
  }

  async requestDeepConsult(
    question: string,
    options: DeepConsultOptions,
    context: RequestContext
  ): Promise<{ consultId: string; estimatedTime: number }> {
    const sanitizedQuestion = this.phiHandler.sanitizeText(question);

    const consultId = await this.openEvidence.deepConsult({
      question: sanitizedQuestion,
      context: {
        specialty: options.specialty,
        researchFocus: options.researchFocus,
        timeframe: options.timeframe,
      },
    }, context);

    return {
      consultId,
      estimatedTime: 180, // 3 minutes typical
    };
  }

  async checkDrugInteraction(
    drugs: string[],
    context: RequestContext
  ): Promise<DrugInteractionResult> {
    const question = `What are the drug interactions between ${drugs.join(' and ')}?`;

    const response = await this.openEvidence.query({
      question,
      context: {
        specialty: 'pharmacology',
        urgency: 'urgent',
      },
    }, context);

    return this.parseDrugInteractionResponse(response);
  }

  private formatClinicalAnswer(response: ClinicalQueryResponse): FormattedClinicalAnswer {
    return {
      summary: response.answer.split('.')[0] + '.',
      detailedAnswer: response.answer,
      citations: response.citations.map(c => ({
        source: c.source,
        title: c.title,
        year: c.year,
      })),
      confidence: {
        score: response.confidence,
        level: response.confidence > 0.9 ? 'high' : response.confidence > 0.7 ? 'moderate' : 'low',
      },
      disclaimer: 'For clinical decision support only. Verify with current guidelines.',
    };
  }

  private parseDrugInteractionResponse(response: ClinicalQueryResponse): DrugInteractionResult {
    // Parse response to structured format
    const answerLower = response.answer.toLowerCase();

    return {
      hasInteraction: !answerLower.includes('no significant interaction'),
      severity: this.determineSeverity(answerLower),
      details: response.answer,
      citations: response.citations,
    };
  }

  private determineSeverity(text: string): 'major' | 'moderate' | 'minor' | 'none' {
    if (text.includes('contraindicated') || text.includes('major')) return 'major';
    if (text.includes('moderate') || text.includes('caution')) return 'moderate';
    if (text.includes('minor')) return 'minor';
    return 'none';
  }
}
```

### Step 3: EHR Integration Layer
```typescript
// src/integrations/ehr/fhir.ts
import { ClinicalEvidenceService } from '../../services/clinical';

// HL7 FHIR CDS Hooks implementation
interface CDSHooksRequest {
  hook: string;
  hookInstance: string;
  context: {
    patientId: string;
    encounterId?: string;
    medications?: FHIRMedication[];
    conditions?: FHIRCondition[];
  };
  prefetch?: {
    patient?: FHIRPatient;
    medications?: FHIRBundle<FHIRMedication>;
  };
}

interface CDSHooksResponse {
  cards: CDSCard[];
}

export class FHIRIntegration {
  constructor(private clinicalService: ClinicalEvidenceService) {}

  async handleCDSHook(request: CDSHooksRequest): Promise<CDSHooksResponse> {
    switch (request.hook) {
      case 'medication-prescribe':
        return this.handleMedicationPrescribe(request);
      case 'order-sign':
        return this.handleOrderSign(request);
      default:
        return { cards: [] };
    }
  }

  private async handleMedicationPrescribe(
    request: CDSHooksRequest
  ): Promise<CDSHooksResponse> {
    const medications = request.context.medications?.map(
      m => m.medicationCodeableConcept?.text
    ).filter(Boolean) as string[];

    if (medications.length < 2) {
      return { cards: [] };
    }

    // Check for drug interactions
    const interaction = await this.clinicalService.checkDrugInteraction(
      medications,
      { userId: 'system', userRole: 'cds-hook' }
    );

    if (!interaction.hasInteraction) {
      return { cards: [] };
    }

    return {
      cards: [{
        uuid: crypto.randomUUID(),
        summary: `Drug Interaction Alert: ${medications.join(', ')}`,
        detail: interaction.details,
        indicator: interaction.severity === 'major' ? 'critical' : 'warning',
        source: {
          label: 'OpenEvidence Clinical Decision Support',
          url: 'https://openevidence.com',
        },
        suggestions: interaction.severity === 'major' ? [{
          label: 'Review alternatives',
          uuid: crypto.randomUUID(),
        }] : undefined,
      }],
    };
  }

  private async handleOrderSign(request: CDSHooksRequest): Promise<CDSHooksResponse> {
    // Similar implementation for order signing
    return { cards: [] };
  }
}
```

## Data Flow Diagram

```
User/EHR Request
       │
       ▼
┌─────────────────┐
│   API Gateway   │
│  (Auth, Rate)   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐    ┌─────────────────┐
│  PHI Sanitizer  │───▶│   HIPAA Audit   │
└────────┬────────┘    └─────────────────┘
         │
         ▼
┌─────────────────┐    ┌─────────────────┐
│  Clinical       │───▶│   Cache Layer   │
│  Query Service  │    │   (Redis)       │
└────────┬────────┘    └────────┬────────┘
         │                      │
         ▼ (cache miss)         │
┌─────────────────┐             │
│  OpenEvidence   │◄────────────┘
│  API Client     │   (cache hit)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  OpenEvidence   │
│  Cloud API      │
└─────────────────┘
```

## Configuration Management

```typescript
// src/config/index.ts
import convict from 'convict';

const config = convict({
  env: {
    doc: 'Application environment',
    format: ['development', 'staging', 'production'],
    default: 'development',
    env: 'NODE_ENV',
  },
  openevidence: {
    apiKey: {
      doc: 'OpenEvidence API key',
      format: String,
      default: '',
      env: 'OPENEVIDENCE_API_KEY',
      sensitive: true,
    },
    orgId: {
      doc: 'OpenEvidence organization ID',
      format: String,
      default: '',
      env: 'OPENEVIDENCE_ORG_ID',
    },
    baseUrl: {
      doc: 'OpenEvidence API base URL',
      format: 'url',
      default: 'https://api.openevidence.com',
      env: 'OPENEVIDENCE_BASE_URL',
    },
    timeout: {
      doc: 'Request timeout in milliseconds',
      format: 'int',
      default: 30000,
      env: 'OPENEVIDENCE_TIMEOUT',
    },
  },
  cache: {
    enabled: {
      doc: 'Enable caching',
      format: Boolean,
      default: true,
    },
    ttlSeconds: {
      doc: 'Default cache TTL',
      format: 'int',
      default: 3600,
    },
  },
  hipaa: {
    auditLogRetentionDays: {
      doc: 'HIPAA audit log retention',
      format: 'int',
      default: 2190, // 6 years
    },
  },
});

config.loadFile(`./config/${config.get('env')}.json`);
config.validate({ allowed: 'strict' });

export default config;
```

## Output
- Layered architecture with separation of concerns
- HIPAA-compliant data handling
- EHR integration ready
- Comprehensive caching
- Health checks and monitoring

## Architecture Checklist
- [ ] Layered architecture implemented
- [ ] PHI sanitization in place
- [ ] HIPAA audit logging
- [ ] Caching strategy configured
- [ ] EHR integration patterns
- [ ] Health checks implemented
- [ ] Configuration management
- [ ] Error boundaries at each layer

## Resources
- [HL7 FHIR](https://www.hl7.org/fhir/)
- [CDS Hooks](https://cds-hooks.hl7.org/)
- [SMART on FHIR](https://smarthealthit.org/)

## Flagship Skills
For multi-environment setup, see `openevidence-multi-env-setup`.
