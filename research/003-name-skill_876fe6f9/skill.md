---
name: openevidence-migration-deep-dive
description: |
  Execute complex OpenEvidence migrations including EHR integration, data migration, and system transitions.
  Use when migrating from legacy clinical decision support systems, integrating with new EHRs,
  or performing major platform transitions.
  Trigger with phrases like "openevidence migration", "ehr integration",
  "migrate to openevidence", "clinical ai migration", "legacy cds migration".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Migration Deep Dive

## Overview
Execute complex migrations including EHR integration, legacy system transitions, and multi-site deployments for OpenEvidence clinical AI.

## Prerequisites
- Understanding of source and target systems
- Access to both environments during migration
- Downtime window approved (if required)
- Rollback plan documented
- Clinical staff communication plan

## Migration Scenarios

| Scenario | Complexity | Downtime | Duration |
|----------|------------|----------|----------|
| Legacy CDS to OpenEvidence | High | 4-8 hours | 2-4 weeks |
| Add EHR Integration | Medium | None | 1-2 weeks |
| Multi-site Expansion | Medium | None | 1-2 weeks |
| API Version Upgrade | Low | None | 1-3 days |

## Instructions

### Scenario 1: Legacy CDS to OpenEvidence Migration

#### Step 1: Assessment and Planning
```typescript
// src/migration/assessment.ts

interface LegacySystemAssessment {
  systemName: string;
  queryVolume: {
    dailyQueries: number;
    peakQueriesPerMinute: number;
    uniqueUsers: number;
  };
  features: {
    clinicalQuery: boolean;
    drugInfo: boolean;
    guidelines: boolean;
    deepResearch: boolean;
  };
  integrations: {
    ehr: string[];
    sso: string[];
    audit: string[];
  };
  dataToMigrate: {
    savedSearches: boolean;
    userPreferences: boolean;
    auditLogs: boolean;
  };
}

async function assessLegacySystem(): Promise<LegacySystemAssessment> {
  // Analyze legacy system usage
  const queryLogs = await legacyDb.queryLogs.aggregate({
    _count: { daily: true },
    _max: { queriesPerMinute: true },
  });

  const users = await legacyDb.users.count();

  return {
    systemName: 'LegacyCDS',
    queryVolume: {
      dailyQueries: queryLogs._count.daily,
      peakQueriesPerMinute: queryLogs._max.queriesPerMinute,
      uniqueUsers: users,
    },
    features: {
      clinicalQuery: true,
      drugInfo: true,
      guidelines: true,
      deepResearch: false, // OpenEvidence advantage
    },
    integrations: {
      ehr: ['Epic', 'Cerner'],
      sso: ['ADFS'],
      audit: ['Splunk'],
    },
    dataToMigrate: {
      savedSearches: true,
      userPreferences: true,
      auditLogs: true, // Required for HIPAA
    },
  };
}
```

#### Step 2: Parallel Running Period
```typescript
// src/migration/parallel-runner.ts

interface ParallelConfig {
  legacyClient: LegacyClient;
  openEvidenceClient: OpenEvidenceClient;
  comparisonRatio: number; // 0-1, how many queries to compare
  autoFallback: boolean;
}

export class ParallelRunner {
  constructor(private config: ParallelConfig) {}

  async query(question: string, context: any): Promise<{
    result: any;
    source: 'legacy' | 'openevidence';
    comparison?: ComparisonResult;
  }> {
    const shouldCompare = Math.random() < this.config.comparisonRatio;

    // Always query OpenEvidence (new system)
    const oePromise = this.config.openEvidenceClient.query({
      question,
      context,
    });

    // Optionally query legacy for comparison
    let legacyPromise: Promise<any> | null = null;
    if (shouldCompare) {
      legacyPromise = this.config.legacyClient.query(question);
    }

    try {
      const oeResult = await oePromise;

      // Log comparison if available
      if (legacyPromise) {
        const legacyResult = await legacyPromise.catch(() => null);
        if (legacyResult) {
          await this.logComparison(question, oeResult, legacyResult);
        }
      }

      return {
        result: oeResult,
        source: 'openevidence',
      };
    } catch (error) {
      if (this.config.autoFallback && legacyPromise) {
        // Fallback to legacy if OpenEvidence fails
        console.warn('[Migration] OpenEvidence failed, falling back to legacy');
        const legacyResult = await legacyPromise;
        return {
          result: legacyResult,
          source: 'legacy',
        };
      }
      throw error;
    }
  }

  private async logComparison(
    question: string,
    oeResult: any,
    legacyResult: any
  ): Promise<void> {
    const comparison: ComparisonResult = {
      timestamp: new Date(),
      questionHash: hashQuestion(question), // Don't log PHI
      openEvidence: {
        responseTime: oeResult.responseTime,
        confidence: oeResult.confidence,
        citationCount: oeResult.citations.length,
      },
      legacy: {
        responseTime: legacyResult.responseTime,
        sourceCount: legacyResult.sources?.length || 0,
      },
    };

    await db.migrationComparisons.create({ data: comparison });
  }
}
```

#### Step 3: Data Migration
```typescript
// src/migration/data-migration.ts

interface MigrationJob {
  id: string;
  type: 'users' | 'preferences' | 'audit_logs';
  status: 'pending' | 'running' | 'completed' | 'failed';
  progress: number;
  errors: string[];
}

export class DataMigrator {
  async migrateUsers(): Promise<MigrationJob> {
    const job = await this.createJob('users');

    try {
      const legacyUsers = await legacyDb.users.findMany();
      let processed = 0;

      for (const user of legacyUsers) {
        try {
          await db.users.upsert({
            where: { email: user.email },
            create: {
              email: user.email,
              firstName: user.firstName,
              lastName: user.lastName,
              role: this.mapRole(user.role),
              migratedFrom: 'legacy',
              migratedAt: new Date(),
            },
            update: {
              migratedFrom: 'legacy',
              migratedAt: new Date(),
            },
          });
        } catch (error: any) {
          job.errors.push(`User ${user.email}: ${error.message}`);
        }

        processed++;
        job.progress = (processed / legacyUsers.length) * 100;
        await this.updateJob(job);
      }

      job.status = 'completed';
      await this.updateJob(job);
      return job;
    } catch (error: any) {
      job.status = 'failed';
      job.errors.push(error.message);
      await this.updateJob(job);
      throw error;
    }
  }

  async migrateAuditLogs(): Promise<MigrationJob> {
    // HIPAA requires preserving audit logs
    const job = await this.createJob('audit_logs');

    const batchSize = 10000;
    let offset = 0;
    let hasMore = true;

    while (hasMore) {
      const logs = await legacyDb.auditLogs.findMany({
        skip: offset,
        take: batchSize,
      });

      if (logs.length === 0) {
        hasMore = false;
        continue;
      }

      // Transform and insert
      const transformed = logs.map(log => ({
        timestamp: log.timestamp,
        userId: log.userId,
        action: this.mapAction(log.action),
        resourceType: 'legacy_query',
        resourceId: log.queryId,
        ipAddress: log.ipAddress,
        migratedFrom: 'legacy',
      }));

      await db.auditLogs.createMany({ data: transformed });

      offset += batchSize;
      job.progress = Math.min(99, job.progress + 10);
      await this.updateJob(job);
    }

    job.status = 'completed';
    job.progress = 100;
    await this.updateJob(job);
    return job;
  }

  private mapRole(legacyRole: string): ClinicalRole {
    const mapping: Record<string, ClinicalRole> = {
      'doctor': ClinicalRole.Physician,
      'nurse': ClinicalRole.Nurse,
      'pharmacist': ClinicalRole.Pharmacist,
      'admin': ClinicalRole.Admin,
    };
    return mapping[legacyRole] || ClinicalRole.Nurse;
  }

  private mapAction(legacyAction: string): string {
    const mapping: Record<string, string> = {
      'search': 'clinical_query',
      'view': 'view_result',
      'export': 'export',
    };
    return mapping[legacyAction] || legacyAction;
  }
}
```

### Scenario 2: EHR Integration (Epic)

#### Step 1: SMART on FHIR Configuration
```typescript
// src/integrations/epic/smart-config.ts

interface SMARTConfig {
  clientId: string;
  clientSecret: string;
  redirectUri: string;
  scope: string[];
  epicEndpoints: {
    authorize: string;
    token: string;
    fhir: string;
  };
}

export const epicConfig: SMARTConfig = {
  clientId: process.env.EPIC_CLIENT_ID!,
  clientSecret: process.env.EPIC_CLIENT_SECRET!,
  redirectUri: `${process.env.APP_URL}/auth/epic/callback`,
  scope: [
    'openid',
    'fhirUser',
    'launch/patient',
    'patient/Patient.read',
    'patient/MedicationRequest.read',
    'patient/Condition.read',
  ],
  epicEndpoints: {
    authorize: 'https://fhir.epic.com/interconnect-fhir-oauth/oauth2/authorize',
    token: 'https://fhir.epic.com/interconnect-fhir-oauth/oauth2/token',
    fhir: 'https://fhir.epic.com/interconnect-fhir-oauth/api/FHIR/R4',
  },
};
```

#### Step 2: CDS Hooks Integration
```typescript
// src/integrations/epic/cds-hooks.ts
import { Router } from 'express';
import { ClinicalEvidenceService } from '../../services/clinical';

const router = Router();

// CDS Hooks Discovery Endpoint
router.get('/cds-services', (req, res) => {
  res.json({
    services: [
      {
        id: 'openevidence-clinical-decision-support',
        hook: 'patient-view',
        title: 'OpenEvidence Clinical Decision Support',
        description: 'AI-powered clinical evidence at the point of care',
        prefetch: {
          patient: 'Patient/{{context.patientId}}',
          conditions: 'Condition?patient={{context.patientId}}',
          medications: 'MedicationRequest?patient={{context.patientId}}&status=active',
        },
      },
      {
        id: 'openevidence-drug-interaction',
        hook: 'medication-prescribe',
        title: 'OpenEvidence Drug Interaction Check',
        description: 'Check for drug interactions using AI clinical evidence',
        prefetch: {
          patient: 'Patient/{{context.patientId}}',
          medications: 'MedicationRequest?patient={{context.patientId}}&status=active',
        },
      },
    ],
  });
});

// CDS Hooks Service Endpoint
router.post('/cds-services/:serviceId', async (req, res) => {
  const { serviceId } = req.params;
  const { hook, context, prefetch } = req.body;

  try {
    let cards: any[] = [];

    switch (serviceId) {
      case 'openevidence-clinical-decision-support':
        cards = await handlePatientView(context, prefetch);
        break;
      case 'openevidence-drug-interaction':
        cards = await handleMedicationPrescribe(context, prefetch);
        break;
    }

    res.json({ cards });
  } catch (error) {
    console.error('[CDS Hooks] Error:', error);
    res.json({ cards: [] });
  }
});

async function handlePatientView(context: any, prefetch: any): Promise<any[]> {
  // Extract relevant conditions
  const conditions = prefetch.conditions?.entry?.map(
    (e: any) => e.resource.code?.text
  ).filter(Boolean) || [];

  if (conditions.length === 0) {
    return [];
  }

  // Query OpenEvidence for relevant clinical information
  const evidence = await clinicalService.queryClinicalEvidence(
    `What are the latest treatment guidelines for ${conditions.join(', ')}?`,
    {
      specialty: 'internal-medicine',
      urgency: 'routine',
    },
    { userId: 'epic-cds', userRole: 'integration' }
  );

  return [{
    uuid: crypto.randomUUID(),
    summary: 'Clinical Evidence Available',
    detail: evidence.summary,
    indicator: 'info',
    source: {
      label: 'OpenEvidence',
      url: 'https://openevidence.com',
      icon: 'https://openevidence.com/icon.png',
    },
    links: [{
      label: 'View Full Evidence',
      url: `${process.env.APP_URL}/evidence?q=${encodeURIComponent(conditions.join(','))}`,
      type: 'absolute',
    }],
  }];
}

async function handleMedicationPrescribe(context: any, prefetch: any): Promise<any[]> {
  const newMedication = context.draftOrders?.[0]?.resource?.medicationCodeableConcept?.text;
  const currentMedications = prefetch.medications?.entry?.map(
    (e: any) => e.resource.medicationCodeableConcept?.text
  ).filter(Boolean) || [];

  if (!newMedication || currentMedications.length === 0) {
    return [];
  }

  // Check for drug interactions
  const interaction = await clinicalService.checkDrugInteraction(
    [newMedication, ...currentMedications],
    { userId: 'epic-cds', userRole: 'integration' }
  );

  if (!interaction.hasInteraction) {
    return [];
  }

  return [{
    uuid: crypto.randomUUID(),
    summary: `Drug Interaction Alert: ${newMedication}`,
    detail: interaction.details,
    indicator: interaction.severity === 'major' ? 'critical' : 'warning',
    source: {
      label: 'OpenEvidence',
      url: 'https://openevidence.com',
    },
    suggestions: interaction.severity === 'major' ? [{
      label: 'Review Alternatives',
      uuid: crypto.randomUUID(),
      isRecommended: true,
    }] : [],
  }];
}

export default router;
```

### Scenario 3: Multi-Site Expansion

#### Step 1: Site Configuration
```typescript
// src/config/multi-site.ts

interface SiteConfig {
  siteId: string;
  siteName: string;
  region: string;
  openEvidenceOrgId: string;
  ehrSystem: 'epic' | 'cerner' | 'meditech' | 'none';
  features: {
    deepConsult: boolean;
    webhooks: boolean;
  };
  quotas: {
    dailyQueries: number;
    monthlyDeepConsults: number;
  };
}

const SITE_CONFIGS: Record<string, SiteConfig> = {
  'main-hospital': {
    siteId: 'main-hospital',
    siteName: 'Main Hospital',
    region: 'us-central',
    openEvidenceOrgId: 'org_main_***',
    ehrSystem: 'epic',
    features: { deepConsult: true, webhooks: true },
    quotas: { dailyQueries: 10000, monthlyDeepConsults: 500 },
  },
  'clinic-network': {
    siteId: 'clinic-network',
    siteName: 'Clinic Network',
    region: 'us-east',
    openEvidenceOrgId: 'org_clinic_***',
    ehrSystem: 'cerner',
    features: { deepConsult: false, webhooks: true },
    quotas: { dailyQueries: 5000, monthlyDeepConsults: 0 },
  },
  'urgent-care': {
    siteId: 'urgent-care',
    siteName: 'Urgent Care Centers',
    region: 'us-west',
    openEvidenceOrgId: 'org_uc_***',
    ehrSystem: 'none',
    features: { deepConsult: false, webhooks: false },
    quotas: { dailyQueries: 2000, monthlyDeepConsults: 0 },
  },
};

export function getSiteConfig(siteId: string): SiteConfig {
  const config = SITE_CONFIGS[siteId];
  if (!config) {
    throw new Error(`Unknown site: ${siteId}`);
  }
  return config;
}
```

## Migration Checklist

### Pre-Migration
- [ ] Stakeholder approval obtained
- [ ] Clinical staff notified
- [ ] Downtime window scheduled (if needed)
- [ ] Rollback plan documented
- [ ] Data backup completed
- [ ] Test environment validated

### During Migration
- [ ] Parallel running period completed
- [ ] Comparison metrics reviewed
- [ ] User data migrated
- [ ] Audit logs migrated (HIPAA)
- [ ] EHR integration tested
- [ ] SSO integration tested

### Post-Migration
- [ ] Legacy system decommissioned
- [ ] All users transitioned
- [ ] Performance baseline established
- [ ] Clinical staff trained
- [ ] Support documentation updated
- [ ] Lessons learned documented

## Output
- Legacy system assessment
- Parallel running configuration
- Data migration scripts
- EHR integration (SMART on FHIR, CDS Hooks)
- Multi-site configuration

## Error Handling
| Migration Issue | Detection | Resolution |
|-----------------|-----------|------------|
| Data corruption | Validation checks | Restore from backup |
| SSO failure | Auth errors | Verify IdP config |
| EHR integration fails | CDS Hooks errors | Check FHIR endpoints |
| Performance regression | Metrics comparison | Optimize or rollback |

## Resources
- [SMART on FHIR](https://smarthealthit.org/)
- [CDS Hooks](https://cds-hooks.hl7.org/)
- [Epic App Orchard](https://appmarket.epic.com/)
- [Cerner Code](https://code.cerner.com/)

## Next Steps
Congratulations! You have completed the OpenEvidence skill pack. For ongoing support, contact support@openevidence.com or visit https://openevidence.com.
