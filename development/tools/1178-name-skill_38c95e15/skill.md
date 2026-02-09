---
name: openevidence-data-handling
description: |
  Implement HIPAA-compliant PHI data handling for OpenEvidence integrations.
  Use when implementing data protection, configuring retention policies,
  or ensuring compliance for clinical AI data flows.
  Trigger with phrases like "openevidence phi", "openevidence data",
  "openevidence hipaa data", "clinical data handling", "patient data protection".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Data Handling

## Overview
Implement HIPAA-compliant Protected Health Information (PHI) handling for OpenEvidence clinical AI integrations.

## Prerequisites
- Signed BAA with OpenEvidence
- Understanding of HIPAA regulations
- Data classification policy
- Encryption infrastructure

## HIPAA Data Categories

| Category | Examples | Handling |
|----------|----------|----------|
| PHI Identifiers | Name, DOB, SSN, MRN | Never send to OpenEvidence |
| Clinical Data | Conditions, medications | May send de-identified |
| Query Results | Answers, citations | Cache with encryption, audit access |
| Audit Logs | User actions, timestamps | Retain 6 years, encrypt |

## 18 HIPAA Identifiers (Never Send to OpenEvidence)

1. Names
2. Geographic data (smaller than state)
3. Dates (except year) related to individual
4. Phone numbers
5. Fax numbers
6. Email addresses
7. Social Security numbers
8. Medical record numbers
9. Health plan beneficiary numbers
10. Account numbers
11. Certificate/license numbers
12. Vehicle identifiers and serial numbers
13. Device identifiers and serial numbers
14. Web URLs
15. IP addresses
16. Biometric identifiers
17. Full-face photographs
18. Any other unique identifying characteristic

## Instructions

### Step 1: PHI Detection and Removal
```typescript
// src/compliance/phi-detector.ts

// Patterns for common PHI
const PHI_PATTERNS = {
  ssn: /\b\d{3}[-\s]?\d{2}[-\s]?\d{4}\b/g,
  mrn: /\b(MRN|Medical Record)[\s:#]*\d{6,12}\b/gi,
  phone: /\b(\+?1[-.\s]?)?\(?\d{3}\)?[-.\s]?\d{3}[-.\s]?\d{4}\b/g,
  email: /\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b/g,
  dob: /\b(DOB|Date of Birth|Born)[\s:]*\d{1,2}[\/\-]\d{1,2}[\/\-]\d{2,4}\b/gi,
  date: /\b\d{1,2}[\/\-]\d{1,2}[\/\-]\d{2,4}\b/g,
  name: /\b(Mr\.|Mrs\.|Ms\.|Dr\.|Patient)\s+[A-Z][a-z]+\s+[A-Z][a-z]+\b/g,
  address: /\b\d{1,5}\s+[A-Z][a-z]+\s+(Street|St|Avenue|Ave|Road|Rd|Boulevard|Blvd)\b/gi,
};

export interface PHIDetectionResult {
  containsPHI: boolean;
  detectedPatterns: string[];
  sanitizedText: string;
}

export function detectPHI(text: string): PHIDetectionResult {
  const detectedPatterns: string[] = [];
  let sanitized = text;

  for (const [pattern, regex] of Object.entries(PHI_PATTERNS)) {
    if (regex.test(text)) {
      detectedPatterns.push(pattern);
      // Reset regex lastIndex for global patterns
      regex.lastIndex = 0;
      sanitized = sanitized.replace(regex, `[${pattern.toUpperCase()}_REDACTED]`);
    }
  }

  return {
    containsPHI: detectedPatterns.length > 0,
    detectedPatterns,
    sanitizedText: sanitized,
  };
}

export function sanitizeForOpenEvidence(text: string): string {
  const result = detectPHI(text);

  if (result.containsPHI) {
    console.warn('[PHI] Detected PHI in query, sanitizing:', result.detectedPatterns);
  }

  return result.sanitizedText;
}
```

### Step 2: Patient Context De-identification
```typescript
// src/compliance/deidentify.ts

interface IdentifiedPatientContext {
  patientId?: string;
  name?: string;
  dateOfBirth?: Date;
  age?: number;
  sex?: 'male' | 'female' | 'other';
  conditions?: string[];
  medications?: string[];
  allergies?: string[];
}

interface DeidentifiedPatientContext {
  ageRange?: string;
  sex?: string;
  conditionCategories?: string[];
  medicationClasses?: string[];
  hasAllergies?: boolean;
}

// Age ranges for de-identification
function getAgeRange(age?: number): string | undefined {
  if (!age) return undefined;
  if (age < 1) return 'infant';
  if (age < 5) return 'toddler';
  if (age < 13) return 'child';
  if (age < 18) return 'adolescent';
  if (age < 30) return 'young-adult';
  if (age < 50) return 'adult';
  if (age < 65) return 'middle-aged';
  if (age < 80) return 'elderly';
  return 'elderly-80+';
}

// Condition category mapping
const CONDITION_CATEGORIES: Record<string, string> = {
  // Cardiovascular
  'hypertension': 'cardiovascular',
  'coronary artery disease': 'cardiovascular',
  'heart failure': 'cardiovascular',
  'atrial fibrillation': 'cardiovascular',

  // Metabolic
  'diabetes': 'metabolic',
  'type 2 diabetes': 'metabolic',
  'hyperlipidemia': 'metabolic',
  'obesity': 'metabolic',

  // Respiratory
  'asthma': 'respiratory',
  'copd': 'respiratory',
  'pneumonia': 'respiratory',

  // Mental Health
  'depression': 'mental-health',
  'anxiety': 'mental-health',

  // Add more as needed
};

function categorizeConditions(conditions?: string[]): string[] | undefined {
  if (!conditions) return undefined;

  const categories = new Set<string>();
  for (const condition of conditions) {
    const category = CONDITION_CATEGORIES[condition.toLowerCase()];
    if (category) {
      categories.add(category);
    } else {
      categories.add('other');
    }
  }

  return Array.from(categories);
}

// Drug class mapping
const DRUG_CLASSES: Record<string, string> = {
  'metformin': 'antidiabetic',
  'lisinopril': 'ace-inhibitor',
  'atorvastatin': 'statin',
  'amlodipine': 'calcium-channel-blocker',
  'metoprolol': 'beta-blocker',
  // Add more as needed
};

function classifyMedications(medications?: string[]): string[] | undefined {
  if (!medications) return undefined;

  const classes = new Set<string>();
  for (const med of medications) {
    const drugClass = DRUG_CLASSES[med.toLowerCase()];
    if (drugClass) {
      classes.add(drugClass);
    } else {
      classes.add('other-medication');
    }
  }

  return Array.from(classes);
}

export function deidentifyPatientContext(
  context: IdentifiedPatientContext
): DeidentifiedPatientContext {
  return {
    ageRange: context.age ? getAgeRange(context.age) : undefined,
    sex: context.sex,
    conditionCategories: categorizeConditions(context.conditions),
    medicationClasses: classifyMedications(context.medications),
    hasAllergies: context.allergies ? context.allergies.length > 0 : undefined,
  };
}
```

### Step 3: Encrypted Storage
```typescript
// src/compliance/encryption.ts
import crypto from 'crypto';

const ALGORITHM = 'aes-256-gcm';
const IV_LENGTH = 16;
const AUTH_TAG_LENGTH = 16;

export class EncryptionService {
  private key: Buffer;

  constructor(encryptionKey: string) {
    // Key should be 32 bytes for AES-256
    this.key = crypto.scryptSync(encryptionKey, 'salt', 32);
  }

  encrypt(plaintext: string): string {
    const iv = crypto.randomBytes(IV_LENGTH);
    const cipher = crypto.createCipheriv(ALGORITHM, this.key, iv);

    let encrypted = cipher.update(plaintext, 'utf8', 'hex');
    encrypted += cipher.final('hex');

    const authTag = cipher.getAuthTag();

    // Return IV + AuthTag + Encrypted data
    return iv.toString('hex') + authTag.toString('hex') + encrypted;
  }

  decrypt(ciphertext: string): string {
    const iv = Buffer.from(ciphertext.slice(0, IV_LENGTH * 2), 'hex');
    const authTag = Buffer.from(
      ciphertext.slice(IV_LENGTH * 2, IV_LENGTH * 2 + AUTH_TAG_LENGTH * 2),
      'hex'
    );
    const encrypted = ciphertext.slice(IV_LENGTH * 2 + AUTH_TAG_LENGTH * 2);

    const decipher = crypto.createDecipheriv(ALGORITHM, this.key, iv);
    decipher.setAuthTag(authTag);

    let decrypted = decipher.update(encrypted, 'hex', 'utf8');
    decrypted += decipher.final('utf8');

    return decrypted;
  }
}

// Usage for caching OpenEvidence responses
export class EncryptedCacheStore {
  private encryption: EncryptionService;
  private redis: Redis;

  constructor(redis: Redis, encryptionKey: string) {
    this.redis = redis;
    this.encryption = new EncryptionService(encryptionKey);
  }

  async set(key: string, value: any, ttlSeconds: number): Promise<void> {
    const encrypted = this.encryption.encrypt(JSON.stringify(value));
    await this.redis.setex(key, ttlSeconds, encrypted);
  }

  async get<T>(key: string): Promise<T | null> {
    const encrypted = await this.redis.get(key);
    if (!encrypted) return null;

    try {
      const decrypted = this.encryption.decrypt(encrypted);
      return JSON.parse(decrypted) as T;
    } catch (error) {
      console.error('Decryption failed, cache entry may be corrupted');
      await this.redis.del(key);
      return null;
    }
  }
}
```

### Step 4: Data Retention Policies
```typescript
// src/compliance/retention.ts

interface RetentionPolicy {
  dataType: string;
  retentionDays: number;
  archiveAfterDays?: number;
  encryptionRequired: boolean;
}

const HIPAA_RETENTION_POLICIES: RetentionPolicy[] = [
  {
    dataType: 'audit_logs',
    retentionDays: 2190, // 6 years (HIPAA minimum)
    archiveAfterDays: 365, // Archive after 1 year
    encryptionRequired: true,
  },
  {
    dataType: 'query_cache',
    retentionDays: 1, // 24 hours max
    encryptionRequired: true,
  },
  {
    dataType: 'deepconsult_reports',
    retentionDays: 365, // 1 year
    archiveAfterDays: 90,
    encryptionRequired: true,
  },
  {
    dataType: 'error_logs',
    retentionDays: 90, // 90 days
    encryptionRequired: false, // No PHI in error logs
  },
];

export class DataRetentionManager {
  constructor(private db: Database, private archive: ArchiveStorage) {}

  async enforceRetentionPolicies(): Promise<RetentionReport> {
    const report: RetentionReport = {
      timestamp: new Date(),
      actions: [],
    };

    for (const policy of HIPAA_RETENTION_POLICIES) {
      const result = await this.enforcePolicy(policy);
      report.actions.push(result);
    }

    return report;
  }

  private async enforcePolicy(policy: RetentionPolicy): Promise<RetentionAction> {
    const cutoffDate = new Date();
    cutoffDate.setDate(cutoffDate.getDate() - policy.retentionDays);

    // Delete expired data
    const deletedCount = await this.db.deleteOlderThan(
      policy.dataType,
      cutoffDate
    );

    // Archive if applicable
    let archivedCount = 0;
    if (policy.archiveAfterDays) {
      const archiveCutoff = new Date();
      archiveCutoff.setDate(archiveCutoff.getDate() - policy.archiveAfterDays);

      const toArchive = await this.db.findOlderThan(policy.dataType, archiveCutoff);
      if (toArchive.length > 0) {
        await this.archive.store(policy.dataType, toArchive);
        archivedCount = toArchive.length;
      }
    }

    return {
      dataType: policy.dataType,
      deletedCount,
      archivedCount,
      cutoffDate,
    };
  }
}

// Run daily via cron job
export async function runDailyRetentionJob(): Promise<void> {
  const manager = new DataRetentionManager(db, archiveStorage);
  const report = await manager.enforceRetentionPolicies();

  console.log('[Retention] Daily job completed:', JSON.stringify(report));

  // Alert if significant data deleted
  const totalDeleted = report.actions.reduce((sum, a) => sum + a.deletedCount, 0);
  if (totalDeleted > 10000) {
    await alertOps(`Data retention: ${totalDeleted} records deleted`);
  }
}
```

### Step 5: Audit Trail for PHI Access
```typescript
// src/compliance/audit-trail.ts

interface PHIAccessLog {
  timestamp: Date;
  userId: string;
  userRole: string;
  action: 'query' | 'view_result' | 'export' | 'deepconsult';
  resourceType: string;
  resourceId: string;
  accessReason?: string;
  ipAddress: string;
  userAgent: string;
  // Never log: actual query content, patient identifiers, clinical data
}

export class PHIAuditLogger {
  private db: Database;
  private encryption: EncryptionService;

  constructor(db: Database, encryptionKey: string) {
    this.db = db;
    this.encryption = new EncryptionService(encryptionKey);
  }

  async logAccess(entry: Omit<PHIAccessLog, 'timestamp'>): Promise<void> {
    const log: PHIAccessLog = {
      ...entry,
      timestamp: new Date(),
    };

    // Encrypt sensitive fields
    const encryptedLog = {
      ...log,
      ipAddress: this.encryption.encrypt(log.ipAddress),
      userAgent: this.encryption.encrypt(log.userAgent),
    };

    await this.db.auditLogs.create({ data: encryptedLog });
  }

  async queryLogs(
    filters: {
      userId?: string;
      startDate?: Date;
      endDate?: Date;
      action?: string;
    },
    pagination: { page: number; pageSize: number }
  ): Promise<{ logs: PHIAccessLog[]; total: number }> {
    // Require authorization to query audit logs
    const result = await this.db.auditLogs.findMany({
      where: {
        ...(filters.userId && { userId: filters.userId }),
        ...(filters.startDate && { timestamp: { gte: filters.startDate } }),
        ...(filters.endDate && { timestamp: { lte: filters.endDate } }),
        ...(filters.action && { action: filters.action }),
      },
      skip: (pagination.page - 1) * pagination.pageSize,
      take: pagination.pageSize,
      orderBy: { timestamp: 'desc' },
    });

    // Decrypt sensitive fields for authorized viewer
    const decryptedLogs = result.map(log => ({
      ...log,
      ipAddress: this.encryption.decrypt(log.ipAddress),
      userAgent: this.encryption.decrypt(log.userAgent),
    }));

    const total = await this.db.auditLogs.count({ where: filters });

    return { logs: decryptedLogs, total };
  }
}
```

## Data Flow Diagram (HIPAA Compliant)

```
Clinical User
      │
      ▼
┌─────────────────┐
│ PHI Detection   │──── Block if PHI detected
│ & Sanitization  │
└────────┬────────┘
         │ (De-identified query)
         ▼
┌─────────────────┐    ┌─────────────────┐
│ Audit Logger    │───▶│ Encrypted       │
│ (Access Log)    │    │ Audit Storage   │
└────────┬────────┘    └─────────────────┘
         │
         ▼
┌─────────────────┐    ┌─────────────────┐
│ OpenEvidence    │───▶│ Encrypted       │
│ API Call        │    │ Cache           │
└────────┬────────┘    └─────────────────┘
         │
         ▼
┌─────────────────┐
│ Response to     │
│ Clinical User   │
└─────────────────┘
```

## Output
- PHI detection and sanitization
- Patient context de-identification
- Encrypted storage for cached data
- HIPAA-compliant retention policies
- Comprehensive audit trail

## Data Handling Checklist
- [ ] PHI patterns defined and detected
- [ ] Patient context de-identified before queries
- [ ] All cached data encrypted (AES-256)
- [ ] Audit logging for all PHI access
- [ ] 6-year retention for audit logs
- [ ] Daily retention job running
- [ ] Encryption keys managed securely
- [ ] BAA signed with OpenEvidence

## Resources
- [HIPAA Privacy Rule](https://www.hhs.gov/hipaa/for-professionals/privacy/index.html)
- [HIPAA Security Rule](https://www.hhs.gov/hipaa/for-professionals/security/index.html)
- [OpenEvidence Security](https://www.openevidence.com/security)

## Next Steps
For enterprise access control, see `openevidence-enterprise-rbac`.
