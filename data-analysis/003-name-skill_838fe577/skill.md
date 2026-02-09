---
name: langfuse-data-handling
description: |
  Manage Langfuse data export, retention, and compliance requirements.
  Use when exporting trace data, configuring retention policies,
  or implementing data compliance for LLM observability.
  Trigger with phrases like "langfuse data export", "langfuse retention",
  "langfuse GDPR", "langfuse compliance", "export langfuse traces".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Data Handling

## Overview
Manage Langfuse data lifecycle including export, retention, and compliance.

## Prerequisites
- Langfuse account with data access
- Understanding of data retention requirements
- Compliance framework knowledge (GDPR, SOC2, etc.)

## Data Categories

| Category | Examples | Retention | Sensitivity |
|----------|----------|-----------|-------------|
| Traces | Operation metadata | 90 days | Medium |
| Generations | LLM inputs/outputs | 30-90 days | High |
| Scores | Evaluation results | 1 year | Low |
| Sessions | User interactions | 90 days | High |
| Metrics | Aggregated stats | 2 years | Low |

## Instructions

### Step 1: Export Trace Data via API

```typescript
// lib/langfuse/export.ts
import { Langfuse } from "langfuse";
import fs from "fs/promises";

interface ExportOptions {
  fromDate: Date;
  toDate: Date;
  format: "json" | "csv";
  includeInputs: boolean;
  includeOutputs: boolean;
}

async function exportTraces(options: ExportOptions) {
  const langfuse = new Langfuse();

  let page = 1;
  const pageSize = 100;
  const allTraces: any[] = [];

  while (true) {
    const response = await langfuse.fetchTraces({
      fromTimestamp: options.fromDate,
      toTimestamp: options.toDate,
      limit: pageSize,
      page,
    });

    if (response.data.length === 0) break;

    // Process each trace
    for (const trace of response.data) {
      const exportedTrace = {
        id: trace.id,
        name: trace.name,
        timestamp: trace.timestamp,
        userId: trace.userId,
        sessionId: trace.sessionId,
        metadata: trace.metadata,
        level: trace.level,
        ...(options.includeInputs && { input: trace.input }),
        ...(options.includeOutputs && { output: trace.output }),
      };

      allTraces.push(exportedTrace);
    }

    page++;

    // Rate limiting
    await new Promise((r) => setTimeout(r, 100));
  }

  // Write to file
  const filename = `langfuse-export-${options.fromDate.toISOString().split("T")[0]}-${options.toDate.toISOString().split("T")[0]}`;

  if (options.format === "json") {
    await fs.writeFile(`${filename}.json`, JSON.stringify(allTraces, null, 2));
  } else {
    const csv = convertToCSV(allTraces);
    await fs.writeFile(`${filename}.csv`, csv);
  }

  return {
    filename,
    count: allTraces.length,
  };
}

function convertToCSV(data: any[]): string {
  if (data.length === 0) return "";

  const headers = Object.keys(data[0]);
  const rows = data.map((item) =>
    headers.map((h) => JSON.stringify(item[h] ?? "")).join(",")
  );

  return [headers.join(","), ...rows].join("\n");
}
```

### Step 2: Implement Data Retention Policy

```typescript
// lib/langfuse/retention.ts
import { Langfuse } from "langfuse";

interface RetentionPolicy {
  defaultRetentionDays: number;
  traceRetentionDays: number;
  generationRetentionDays: number;
  scoreRetentionDays: number;
  piiRetentionDays: number;
}

const DEFAULT_POLICY: RetentionPolicy = {
  defaultRetentionDays: 90,
  traceRetentionDays: 90,
  generationRetentionDays: 30,
  scoreRetentionDays: 365,
  piiRetentionDays: 30,
};

// Note: Data deletion requires self-hosted Langfuse with DB access
// or using Langfuse's data retention settings in dashboard

async function applyRetentionPolicy(policy: RetentionPolicy = DEFAULT_POLICY) {
  const langfuse = new Langfuse();

  // Get traces older than retention period
  const cutoffDate = new Date();
  cutoffDate.setDate(cutoffDate.getDate() - policy.traceRetentionDays);

  const oldTraces = await langfuse.fetchTraces({
    toTimestamp: cutoffDate,
    limit: 1000,
  });

  console.log(`Found ${oldTraces.data.length} traces older than ${policy.traceRetentionDays} days`);

  // For self-hosted: Delete via direct DB access
  // For cloud: Configure retention in dashboard settings

  return {
    tracesIdentified: oldTraces.data.length,
    cutoffDate,
  };
}
```

### Step 3: Implement GDPR Data Subject Requests

```typescript
// lib/langfuse/gdpr.ts

interface DataSubjectRequest {
  type: "access" | "deletion" | "rectification";
  userId: string;
  requestId: string;
  requestedAt: Date;
}

async function handleDataAccessRequest(userId: string): Promise<any> {
  const langfuse = new Langfuse();

  // Fetch all user data
  const traces = await langfuse.fetchTraces({
    userId,
    limit: 10000,
  });

  const userData = {
    requestType: "GDPR Data Access",
    userId,
    generatedAt: new Date().toISOString(),
    data: {
      traces: traces.data.map((trace) => ({
        id: trace.id,
        timestamp: trace.timestamp,
        name: trace.name,
        // Include relevant data fields
        sessionId: trace.sessionId,
        // Redact sensitive internal data
        input: "[Available on request]",
        output: "[Available on request]",
      })),
      totalTraces: traces.data.length,
    },
  };

  return userData;
}

async function handleDataDeletionRequest(userId: string): Promise<any> {
  const langfuse = new Langfuse();

  // Fetch all user traces
  const traces = await langfuse.fetchTraces({
    userId,
    limit: 10000,
  });

  // For self-hosted: Delete from database
  // For cloud: Submit deletion request to Langfuse

  const deletionRecord = {
    requestType: "GDPR Deletion",
    userId,
    processedAt: new Date().toISOString(),
    tracesMarkedForDeletion: traces.data.length,
    status: "pending", // or "completed" after actual deletion
  };

  // Log for audit trail
  console.log("Data deletion request processed:", deletionRecord);

  return deletionRecord;
}

async function processDataSubjectRequest(
  request: DataSubjectRequest
): Promise<any> {
  const auditLog = {
    requestId: request.requestId,
    type: request.type,
    userId: request.userId,
    receivedAt: request.requestedAt.toISOString(),
    processedAt: new Date().toISOString(),
  };

  switch (request.type) {
    case "access":
      return {
        ...auditLog,
        result: await handleDataAccessRequest(request.userId),
      };

    case "deletion":
      return {
        ...auditLog,
        result: await handleDataDeletionRequest(request.userId),
      };

    case "rectification":
      return {
        ...auditLog,
        status: "manual_review_required",
        message: "Rectification requests require manual review",
      };
  }
}
```

### Step 4: Implement Data Anonymization

```typescript
// lib/langfuse/anonymize.ts

interface AnonymizationConfig {
  hashUserId: boolean;
  removeInputs: boolean;
  removeOutputs: boolean;
  removeMetadataFields: string[];
}

function anonymizeTrace(
  trace: any,
  config: AnonymizationConfig = {
    hashUserId: true,
    removeInputs: false,
    removeOutputs: false,
    removeMetadataFields: ["email", "name", "phone"],
  }
): any {
  const anonymized = { ...trace };

  // Hash user ID
  if (config.hashUserId && anonymized.userId) {
    anonymized.userId = hashString(anonymized.userId);
    anonymized.sessionId = anonymized.sessionId
      ? hashString(anonymized.sessionId)
      : null;
  }

  // Remove sensitive content
  if (config.removeInputs) {
    anonymized.input = "[REDACTED]";
  }

  if (config.removeOutputs) {
    anonymized.output = "[REDACTED]";
  }

  // Clean metadata
  if (anonymized.metadata) {
    for (const field of config.removeMetadataFields) {
      delete anonymized.metadata[field];
    }
  }

  return anonymized;
}

function hashString(str: string): string {
  const crypto = require("crypto");
  return crypto.createHash("sha256").update(str).digest("hex").slice(0, 16);
}

// Export anonymized data for analytics
async function exportAnonymizedData(
  fromDate: Date,
  toDate: Date
): Promise<any[]> {
  const langfuse = new Langfuse();

  const traces = await langfuse.fetchTraces({
    fromTimestamp: fromDate,
    toTimestamp: toDate,
    limit: 10000,
  });

  return traces.data.map((trace) =>
    anonymizeTrace(trace, {
      hashUserId: true,
      removeInputs: true,
      removeOutputs: true,
      removeMetadataFields: ["email", "name", "phone", "ip"],
    })
  );
}
```

### Step 5: Create Audit Trail

```typescript
// lib/langfuse/audit.ts

interface AuditEvent {
  timestamp: Date;
  action: string;
  actor: string;
  resource: string;
  details: Record<string, any>;
  result: "success" | "failure";
}

class AuditLogger {
  private events: AuditEvent[] = [];

  log(event: Omit<AuditEvent, "timestamp">) {
    const auditEvent: AuditEvent = {
      ...event,
      timestamp: new Date(),
    };

    this.events.push(auditEvent);

    // Also write to permanent storage
    this.persist(auditEvent);
  }

  private async persist(event: AuditEvent) {
    // Write to audit log file
    const logLine = JSON.stringify(event) + "\n";
    await fs.appendFile("audit.log", logLine);

    // Or send to audit service
    // await auditService.log(event);
  }

  // Query audit trail
  async query(options: {
    from?: Date;
    to?: Date;
    actor?: string;
    action?: string;
  }): Promise<AuditEvent[]> {
    // Read from audit log
    const content = await fs.readFile("audit.log", "utf-8");
    let events = content
      .split("\n")
      .filter(Boolean)
      .map((line) => JSON.parse(line));

    if (options.from) {
      events = events.filter((e) => new Date(e.timestamp) >= options.from!);
    }
    if (options.to) {
      events = events.filter((e) => new Date(e.timestamp) <= options.to!);
    }
    if (options.actor) {
      events = events.filter((e) => e.actor === options.actor);
    }
    if (options.action) {
      events = events.filter((e) => e.action === options.action);
    }

    return events;
  }
}

export const auditLogger = new AuditLogger();

// Usage
auditLogger.log({
  action: "data_export",
  actor: "admin@company.com",
  resource: "traces",
  details: { count: 1500, dateRange: "2024-01-01 to 2024-01-31" },
  result: "success",
});
```

## Output
- Data export functionality (JSON/CSV)
- Retention policy implementation
- GDPR compliance handlers
- Data anonymization
- Audit trail logging

## Compliance Checklist

| Requirement | Implementation | Status |
|-------------|----------------|--------|
| Data inventory | Track what data is collected | [ ] |
| Retention limits | Auto-delete old data | [ ] |
| Access requests | Export user data | [ ] |
| Deletion requests | Remove user data | [ ] |
| Anonymization | Hash/redact PII | [ ] |
| Audit trail | Log all data access | [ ] |
| Encryption | TLS in transit, at rest | [ ] |

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Export timeout | Too much data | Use pagination |
| Missing user data | Wrong userId format | Verify identifier |
| Deletion failed | No DB access | Contact Langfuse support |
| Audit gaps | Async logging | Use sync logging for critical |

## Resources
- [Langfuse Privacy Policy](https://langfuse.com/privacy)
- [GDPR Guidelines](https://gdpr.eu/)
- [Langfuse Self-Hosting](https://langfuse.com/docs/deployment/self-host)
- [Data Protection](https://langfuse.com/docs/security)

## Next Steps
For enterprise RBAC, see `langfuse-enterprise-rbac`.
