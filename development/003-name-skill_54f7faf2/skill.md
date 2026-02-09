---
name: openevidence-security-basics
description: |
  Apply OpenEvidence security best practices for HIPAA compliance and PHI protection.
  Use when securing API keys, implementing PHI handling,
  or auditing OpenEvidence security configuration.
  Trigger with phrases like "openevidence security", "openevidence hipaa",
  "openevidence phi", "secure openevidence", "openevidence compliance".
allowed-tools: Read, Write, Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Security Basics

## Overview
Security best practices for OpenEvidence integrations handling Protected Health Information (PHI) in compliance with HIPAA regulations.

## Prerequisites
- OpenEvidence SDK installed
- Understanding of HIPAA requirements
- Signed Business Associate Agreement (BAA)
- Access to organization security policies

## OpenEvidence Security Certifications
- SOC 2 Type II certified
- HIPAA compliant
- AES-256 encryption at rest
- TLS 1.2+ encryption in transit
- Google Cloud Platform hosted

## Instructions

### Step 1: Secure Credential Management
```bash
# .env (NEVER commit to git)
OPENEVIDENCE_API_KEY=oe_live_***
OPENEVIDENCE_ORG_ID=org_***
OPENEVIDENCE_WEBHOOK_SECRET=whsec_***

# .gitignore (REQUIRED)
.env
.env.local
.env.*.local
*.pem
*.key
credentials*.json
```

```typescript
// src/config/secrets.ts
// Use secret manager in production (AWS Secrets Manager, GCP Secret Manager, HashiCorp Vault)
import { SecretManagerServiceClient } from '@google-cloud/secret-manager';

const client = new SecretManagerServiceClient();

export async function getOpenEvidenceCredentials(): Promise<{
  apiKey: string;
  orgId: string;
}> {
  const [apiKeyVersion] = await client.accessSecretVersion({
    name: 'projects/my-project/secrets/openevidence-api-key/versions/latest',
  });
  const [orgIdVersion] = await client.accessSecretVersion({
    name: 'projects/my-project/secrets/openevidence-org-id/versions/latest',
  });

  return {
    apiKey: apiKeyVersion.payload!.data!.toString(),
    orgId: orgIdVersion.payload!.data!.toString(),
  };
}
```

### Step 2: PHI Handling - Input Sanitization
```typescript
// src/openevidence/phi-handler.ts
// IMPORTANT: OpenEvidence should NOT receive direct PHI in queries

interface SanitizedQuery {
  question: string;
  context: {
    ageRange?: string; // "65-74" not exact age
    sex?: string;
    conditionCategories?: string[]; // "cardiovascular" not specific diagnosis
  };
}

export function sanitizeQueryForOpenEvidence(
  question: string,
  patientContext?: PatientContext
): SanitizedQuery {
  // Remove any PHI that might have leaked into the question
  let sanitized = question;

  // Remove potential names
  sanitized = sanitized.replace(/\b(Mr\.|Mrs\.|Ms\.|Dr\.)\s+[A-Z][a-z]+\b/g, '[PATIENT]');

  // Remove dates of birth
  sanitized = sanitized.replace(/\b\d{1,2}[\/\-]\d{1,2}[\/\-]\d{2,4}\b/g, '[DATE]');

  // Remove MRNs
  sanitized = sanitized.replace(/\bMRN[:\s]*\d+\b/gi, '[MRN]');

  // Remove SSNs
  sanitized = sanitized.replace(/\b\d{3}-\d{2}-\d{4}\b/g, '[SSN]');

  return {
    question: sanitized,
    context: patientContext ? {
      ageRange: getAgeRange(patientContext.age),
      sex: patientContext.sex,
      conditionCategories: categorizeConditions(patientContext.conditions),
    } : undefined,
  };
}

function getAgeRange(age?: number): string | undefined {
  if (!age) return undefined;
  if (age < 18) return 'pediatric';
  if (age < 40) return '18-39';
  if (age < 65) return '40-64';
  if (age < 75) return '65-74';
  return '75+';
}

function categorizeConditions(conditions?: string[]): string[] | undefined {
  // Map specific conditions to categories
  const categoryMap: Record<string, string> = {
    'hypertension': 'cardiovascular',
    'diabetes': 'endocrine',
    'asthma': 'pulmonary',
    // ... more mappings
  };

  return conditions?.map(c => categoryMap[c.toLowerCase()] || 'other');
}
```

### Step 3: Audit Logging (HIPAA Required)
```typescript
// src/openevidence/audit-log.ts
interface AuditEntry {
  timestamp: Date;
  eventType: 'query' | 'deepconsult' | 'access' | 'export';
  userId: string;
  userRole: string;
  action: string;
  resourceType: 'clinical_answer' | 'research_report';
  resourceId?: string;
  ipAddress: string;
  userAgent: string;
  success: boolean;
  errorCode?: string;
  // Never log the actual query content or response - may contain PHI
}

export class HIPAAAuditLogger {
  private logStore: AuditLogStore;

  constructor(store: AuditLogStore) {
    this.logStore = store;
  }

  async logClinicalQuery(
    userId: string,
    userRole: string,
    queryId: string,
    success: boolean,
    request: Request
  ): Promise<void> {
    await this.logStore.insert({
      timestamp: new Date(),
      eventType: 'query',
      userId,
      userRole,
      action: 'clinical_query',
      resourceType: 'clinical_answer',
      resourceId: queryId,
      ipAddress: this.getClientIP(request),
      userAgent: request.headers.get('user-agent') || 'unknown',
      success,
    });
  }

  async logDeepConsult(
    userId: string,
    userRole: string,
    consultId: string,
    success: boolean,
    request: Request
  ): Promise<void> {
    await this.logStore.insert({
      timestamp: new Date(),
      eventType: 'deepconsult',
      userId,
      userRole,
      action: 'research_synthesis',
      resourceType: 'research_report',
      resourceId: consultId,
      ipAddress: this.getClientIP(request),
      userAgent: request.headers.get('user-agent') || 'unknown',
      success,
    });
  }

  private getClientIP(request: Request): string {
    // Handle proxies
    const forwarded = request.headers.get('x-forwarded-for');
    if (forwarded) return forwarded.split(',')[0].trim();
    return 'unknown';
  }
}
```

### Step 4: Webhook Signature Verification
```typescript
// src/openevidence/webhook-security.ts
import crypto from 'crypto';

export function verifyWebhookSignature(
  payload: string | Buffer,
  signature: string,
  secret: string,
  tolerance: number = 300 // 5 minute tolerance
): boolean {
  // Parse signature header: t=timestamp,v1=signature
  const parts = signature.split(',').reduce((acc, part) => {
    const [key, value] = part.split('=');
    acc[key] = value;
    return acc;
  }, {} as Record<string, string>);

  const timestamp = parseInt(parts['t']);
  const providedSignature = parts['v1'];

  // Check timestamp to prevent replay attacks
  const now = Math.floor(Date.now() / 1000);
  if (Math.abs(now - timestamp) > tolerance) {
    console.error('Webhook timestamp too old or in future');
    return false;
  }

  // Compute expected signature
  const signedPayload = `${timestamp}.${payload}`;
  const expectedSignature = crypto
    .createHmac('sha256', secret)
    .update(signedPayload)
    .digest('hex');

  // Timing-safe comparison
  try {
    return crypto.timingSafeEqual(
      Buffer.from(providedSignature),
      Buffer.from(expectedSignature)
    );
  } catch {
    return false;
  }
}
```

### Step 5: Data Retention Compliance
```typescript
// src/openevidence/data-retention.ts
// HIPAA requires minimum 6-year retention of access logs

interface RetentionPolicy {
  auditLogs: number;      // days
  queryResults: number;   // days
  deepConsultReports: number; // days
}

const HIPAA_RETENTION: RetentionPolicy = {
  auditLogs: 2190,        // 6 years
  queryResults: 0,        // Do not store (re-query as needed)
  deepConsultReports: 365, // 1 year (user can re-run)
};

export class DataRetentionManager {
  constructor(private policy: RetentionPolicy = HIPAA_RETENTION) {}

  async cleanupExpiredData(): Promise<CleanupReport> {
    const cutoffs = {
      auditLogs: new Date(Date.now() - this.policy.auditLogs * 24 * 60 * 60 * 1000),
      deepConsultReports: new Date(Date.now() - this.policy.deepConsultReports * 24 * 60 * 60 * 1000),
    };

    // Note: Audit logs should be archived, not deleted
    const archivedLogs = await db.auditLogs.archive({
      where: { timestamp: { lt: cutoffs.auditLogs } },
    });

    const deletedReports = await db.deepConsultReports.deleteMany({
      where: { createdAt: { lt: cutoffs.deepConsultReports } },
    });

    return {
      auditLogsArchived: archivedLogs.count,
      reportsDeleted: deletedReports.count,
      executedAt: new Date(),
    };
  }
}
```

## Security Checklist
- [ ] API keys stored in secret manager (not env vars in production)
- [ ] BAA signed with OpenEvidence
- [ ] PHI sanitization before sending queries
- [ ] Audit logging enabled for all access
- [ ] Webhook signatures verified
- [ ] TLS 1.2+ enforced
- [ ] IP allowlist configured (if available)
- [ ] Different keys for dev/staging/prod
- [ ] Key rotation schedule established
- [ ] Data retention policy implemented

## Output
- Secure credential management
- PHI sanitization layer
- HIPAA-compliant audit logging
- Webhook signature verification
- Data retention compliance

## Error Handling
| Security Issue | Detection | Mitigation |
|----------------|-----------|------------|
| Exposed API key | Git scanning alert | Rotate immediately, audit access |
| PHI in query | Log pattern matching | Block request, alert compliance |
| Failed signature | Webhook verification | Reject webhook, alert security |
| Unauthorized access | Audit log review | Revoke access, investigate |

## Examples

### Secure Query Service
```typescript
const auditLogger = new HIPAAAuditLogger(auditStore);

export async function secureClinicaQuery(
  question: string,
  patientContext: PatientContext | undefined,
  user: AuthenticatedUser,
  request: Request
): Promise<ClinicalResponse> {
  // 1. Sanitize input
  const sanitized = sanitizeQueryForOpenEvidence(question, patientContext);

  // 2. Make query
  const response = await client.query({
    question: sanitized.question,
    context: {
      specialty: 'internal-medicine',
      urgency: 'routine',
      ...sanitized.context,
    },
  });

  // 3. Audit log
  await auditLogger.logClinicalQuery(
    user.id,
    user.role,
    response.id,
    true,
    request
  );

  return response;
}
```

## Resources
- [OpenEvidence Security](https://www.openevidence.com/security)
- [HIPAA Security Rule](https://www.hhs.gov/hipaa/for-professionals/security/index.html)
- [NIST Cybersecurity Framework](https://www.nist.gov/cyberframework)

## Next Steps
For production deployment, see `openevidence-prod-checklist`.
