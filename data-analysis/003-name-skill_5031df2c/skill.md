---
name: langfuse-migration-deep-dive
description: |
  Execute complex Langfuse migrations including data migration and platform changes.
  Use when migrating from other observability platforms, moving between Langfuse instances,
  or performing major infrastructure migrations.
  Trigger with phrases like "langfuse migration", "migrate to langfuse",
  "langfuse data migration", "langfuse platform migration", "switch to langfuse".
allowed-tools: Read, Write, Edit, Bash(npm:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Migration Deep Dive

## Overview
Comprehensive guide for complex migrations to or between Langfuse instances.

## Prerequisites
- Understanding of source and target systems
- Database access (for data migrations)
- Downtime window planned (if needed)
- Rollback plan prepared

## Migration Scenarios

| Scenario | Complexity | Downtime | Data Loss Risk |
|----------|------------|----------|----------------|
| Cloud to Cloud (same) | Low | None | None |
| Self-hosted to Cloud | Medium | Minutes | Low |
| Cloud to Self-hosted | Medium | Minutes | Low |
| Other platform to Langfuse | High | Hours | Medium |
| SDK version upgrade | Low | None | None |

## Instructions

### Scenario 1: Migrate from Cloud to Self-Hosted

#### Step 1: Export Data from Cloud

```typescript
// scripts/export-cloud-data.ts
import { Langfuse } from "langfuse";
import fs from "fs/promises";

async function exportCloudData() {
  const langfuse = new Langfuse({
    publicKey: process.env.LANGFUSE_CLOUD_PUBLIC_KEY!,
    secretKey: process.env.LANGFUSE_CLOUD_SECRET_KEY!,
    baseUrl: "https://cloud.langfuse.com",
  });

  const exportDir = `export-${Date.now()}`;
  await fs.mkdir(exportDir, { recursive: true });

  // Export traces
  console.log("Exporting traces...");
  let page = 1;
  let totalTraces = 0;

  while (true) {
    const traces = await langfuse.fetchTraces({
      limit: 100,
      page,
    });

    if (traces.data.length === 0) break;

    await fs.writeFile(
      `${exportDir}/traces-${page}.json`,
      JSON.stringify(traces.data, null, 2)
    );

    totalTraces += traces.data.length;
    page++;

    // Rate limiting
    await new Promise((r) => setTimeout(r, 200));
  }

  console.log(`Exported ${totalTraces} traces`);

  // Export scores
  console.log("Exporting scores...");
  // Similar pagination for scores

  // Export datasets
  console.log("Exporting datasets...");
  const datasets = await langfuse.fetchDatasets({});
  await fs.writeFile(
    `${exportDir}/datasets.json`,
    JSON.stringify(datasets.data, null, 2)
  );

  return exportDir;
}
```

#### Step 2: Set Up Self-Hosted Instance

```yaml
# docker-compose.self-hosted.yml
version: "3.8"

services:
  langfuse:
    image: langfuse/langfuse:latest
    ports:
      - "3000:3000"
    environment:
      - DATABASE_URL=postgresql://postgres:${DB_PASSWORD}@db:5432/langfuse
      - NEXTAUTH_SECRET=${NEXTAUTH_SECRET}
      - NEXTAUTH_URL=${LANGFUSE_URL}
      - SALT=${LANGFUSE_SALT}
      - ENCRYPTION_KEY=${ENCRYPTION_KEY}
    depends_on:
      db:
        condition: service_healthy

  db:
    image: postgres:15-alpine
    environment:
      - POSTGRES_USER=postgres
      - POSTGRES_PASSWORD=${DB_PASSWORD}
      - POSTGRES_DB=langfuse
    volumes:
      - langfuse-db:/var/lib/postgresql/data
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U postgres"]
      interval: 5s
      timeout: 5s
      retries: 5

volumes:
  langfuse-db:
```

#### Step 3: Import Data to Self-Hosted

```typescript
// scripts/import-to-self-hosted.ts
import { Langfuse } from "langfuse";
import fs from "fs/promises";
import path from "path";

async function importToSelfHosted(exportDir: string) {
  const langfuse = new Langfuse({
    publicKey: process.env.LANGFUSE_SELF_HOSTED_PUBLIC_KEY!,
    secretKey: process.env.LANGFUSE_SELF_HOSTED_SECRET_KEY!,
    baseUrl: process.env.LANGFUSE_SELF_HOSTED_URL!,
  });

  // Get all trace files
  const files = await fs.readdir(exportDir);
  const traceFiles = files.filter((f) => f.startsWith("traces-"));

  for (const file of traceFiles) {
    const content = await fs.readFile(path.join(exportDir, file), "utf-8");
    const traces = JSON.parse(content);

    for (const traceData of traces) {
      // Recreate trace
      const trace = langfuse.trace({
        name: traceData.name,
        userId: traceData.userId,
        sessionId: traceData.sessionId,
        input: traceData.input,
        output: traceData.output,
        metadata: {
          ...traceData.metadata,
          migratedFrom: "cloud",
          originalId: traceData.id,
        },
      });

      // Note: Observations/generations need to be recreated separately
      // if you have access to the detailed data
    }

    await langfuse.flushAsync();
    console.log(`Imported ${file}`);
  }
}
```

### Scenario 2: Migrate from LangSmith to Langfuse

#### Step 1: Create Adapter for LangSmith Data

```typescript
// lib/migration/langsmith-adapter.ts

interface LangSmithRun {
  id: string;
  name: string;
  run_type: "chain" | "llm" | "tool";
  inputs: any;
  outputs: any;
  start_time: string;
  end_time: string;
  extra?: {
    metadata?: Record<string, any>;
  };
  child_runs?: LangSmithRun[];
}

interface LangfuseTraceData {
  name: string;
  input: any;
  output: any;
  metadata: Record<string, any>;
  spans: LangfuseSpanData[];
  generations: LangfuseGenerationData[];
}

interface LangfuseSpanData {
  name: string;
  input: any;
  output: any;
  startTime: Date;
  endTime: Date;
}

interface LangfuseGenerationData {
  name: string;
  model: string;
  input: any;
  output: any;
  startTime: Date;
  endTime: Date;
  usage?: { promptTokens: number; completionTokens: number };
}

function convertLangSmithRun(run: LangSmithRun): LangfuseTraceData {
  const spans: LangfuseSpanData[] = [];
  const generations: LangfuseGenerationData[] = [];

  // Process child runs
  function processRun(r: LangSmithRun, depth: number = 0) {
    if (r.run_type === "llm") {
      generations.push({
        name: r.name,
        model: r.extra?.metadata?.model || "unknown",
        input: r.inputs,
        output: r.outputs,
        startTime: new Date(r.start_time),
        endTime: new Date(r.end_time),
        usage: r.extra?.metadata?.usage,
      });
    } else {
      spans.push({
        name: r.name,
        input: r.inputs,
        output: r.outputs,
        startTime: new Date(r.start_time),
        endTime: new Date(r.end_time),
      });
    }

    for (const child of r.child_runs || []) {
      processRun(child, depth + 1);
    }
  }

  processRun(run);

  return {
    name: run.name,
    input: run.inputs,
    output: run.outputs,
    metadata: {
      ...run.extra?.metadata,
      migratedFrom: "langsmith",
      originalId: run.id,
    },
    spans,
    generations,
  };
}
```

#### Step 2: Bulk Import Converted Data

```typescript
// scripts/migrate-from-langsmith.ts
import { Langfuse } from "langfuse";

async function migrateFromLangSmith(runs: LangSmithRun[]) {
  const langfuse = new Langfuse();
  let migrated = 0;
  let failed = 0;

  for (const run of runs) {
    try {
      const converted = convertLangSmithRun(run);

      const trace = langfuse.trace({
        name: converted.name,
        input: converted.input,
        output: converted.output,
        metadata: converted.metadata,
      });

      // Create spans
      for (const span of converted.spans) {
        const s = trace.span({
          name: span.name,
          input: span.input,
          startTime: span.startTime,
        });
        s.end({
          output: span.output,
          endTime: span.endTime,
        });
      }

      // Create generations
      for (const gen of converted.generations) {
        const g = trace.generation({
          name: gen.name,
          model: gen.model,
          input: gen.input,
          startTime: gen.startTime,
        });
        g.end({
          output: gen.output,
          endTime: gen.endTime,
          usage: gen.usage,
        });
      }

      migrated++;
    } catch (error) {
      console.error(`Failed to migrate run ${run.id}:`, error);
      failed++;
    }

    // Periodic flush
    if (migrated % 100 === 0) {
      await langfuse.flushAsync();
      console.log(`Progress: ${migrated} migrated, ${failed} failed`);
    }
  }

  await langfuse.flushAsync();
  return { migrated, failed };
}
```

### Scenario 3: Zero-Downtime SDK Migration

```typescript
// lib/migration/dual-write.ts

// During migration, write to both old and new systems
class DualWriteLangfuse {
  private oldClient: any; // Old observability system
  private newClient: Langfuse;
  private writeToOld: boolean = true;
  private writeToNew: boolean = true;

  constructor(config: {
    oldConfig: any;
    newConfig: ConstructorParameters<typeof Langfuse>[0];
  }) {
    this.oldClient = createOldClient(config.oldConfig);
    this.newClient = new Langfuse(config.newConfig);
  }

  trace(params: any) {
    const traces: any[] = [];

    if (this.writeToOld) {
      try {
        traces.push({ type: "old", trace: this.oldClient.trace(params) });
      } catch (error) {
        console.error("Old client trace failed:", error);
      }
    }

    if (this.writeToNew) {
      try {
        traces.push({ type: "new", trace: this.newClient.trace(params) });
      } catch (error) {
        console.error("New client trace failed:", error);
      }
    }

    // Return proxy that writes to both
    return new Proxy(traces[0]?.trace || {}, {
      get: (target, prop) => {
        if (prop === "span" || prop === "generation") {
          return (...args: any[]) => {
            for (const { trace } of traces) {
              trace[prop]?.(...args);
            }
            return target[prop]?.(...args);
          };
        }
        return target[prop];
      },
    });
  }

  // Gradual cutover
  setWriteMode(options: { old: boolean; new: boolean }) {
    this.writeToOld = options.old;
    this.writeToNew = options.new;
  }

  async flush() {
    await Promise.all([
      this.writeToOld && this.oldClient.flush?.(),
      this.writeToNew && this.newClient.flushAsync(),
    ]);
  }
}

// Migration timeline
// Week 1: writeToOld: true, writeToNew: true (dual write)
// Week 2: writeToOld: true, writeToNew: true (verify new system)
// Week 3: writeToOld: false, writeToNew: true (cutover)
// Week 4: Remove old client code
```

### Step 4: Validation and Verification

```typescript
// scripts/validate-migration.ts

interface MigrationValidation {
  source: {
    traceCount: number;
    generationCount: number;
    dateRange: { start: Date; end: Date };
  };
  target: {
    traceCount: number;
    generationCount: number;
    dateRange: { start: Date; end: Date };
  };
  discrepancies: string[];
}

async function validateMigration(
  sourceClient: any,
  targetClient: Langfuse
): Promise<MigrationValidation> {
  const discrepancies: string[] = [];

  // Count source data
  const sourceTraces = await sourceClient.fetchTraces({ limit: 1 });
  const sourceCount = sourceTraces.totalCount || sourceTraces.data.length;

  // Count target data
  const targetTraces = await targetClient.fetchTraces({ limit: 1 });
  const targetCount = targetTraces.totalCount || targetTraces.data.length;

  // Compare counts
  const countDiff = Math.abs(sourceCount - targetCount);
  if (countDiff > sourceCount * 0.01) {
    discrepancies.push(
      `Trace count mismatch: source=${sourceCount}, target=${targetCount}`
    );
  }

  // Spot check random traces
  const sampleTraces = await sourceClient.fetchTraces({ limit: 10 });
  for (const trace of sampleTraces.data) {
    const targetTrace = await targetClient.fetchTraces({
      filter: { metadata: { originalId: trace.id } },
    });

    if (targetTrace.data.length === 0) {
      discrepancies.push(`Missing trace: ${trace.id}`);
    }
  }

  return {
    source: {
      traceCount: sourceCount,
      generationCount: 0, // Calculate as needed
      dateRange: { start: new Date(), end: new Date() },
    },
    target: {
      traceCount: targetCount,
      generationCount: 0,
      dateRange: { start: new Date(), end: new Date() },
    },
    discrepancies,
  };
}
```

## Migration Checklist

| Phase | Task | Status |
|-------|------|--------|
| **Planning** | | |
| | Inventory source data | [ ] |
| | Plan downtime window | [ ] |
| | Create rollback plan | [ ] |
| | Notify stakeholders | [ ] |
| **Execution** | | |
| | Export source data | [ ] |
| | Set up target system | [ ] |
| | Import data | [ ] |
| | Update application config | [ ] |
| **Validation** | | |
| | Verify trace counts | [ ] |
| | Spot check data quality | [ ] |
| | Test new integration | [ ] |
| | Monitor for errors | [ ] |
| **Cleanup** | | |
| | Remove old config | [ ] |
| | Archive source data | [ ] |
| | Update documentation | [ ] |

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Data loss | Incomplete export | Re-run with pagination |
| Duplicate traces | Re-import | Dedupe by originalId |
| Missing metadata | Format mismatch | Update adapter |
| Performance issues | Large import | Use batch processing |

## Rollback Plan

```typescript
// If migration fails, rollback to old system
async function rollback() {
  // 1. Update environment variables
  process.env.LANGFUSE_PUBLIC_KEY = process.env.OLD_LANGFUSE_PUBLIC_KEY;
  process.env.LANGFUSE_SECRET_KEY = process.env.OLD_LANGFUSE_SECRET_KEY;
  process.env.LANGFUSE_HOST = process.env.OLD_LANGFUSE_HOST;

  // 2. Restart application
  console.log("Rollback complete. Restart application to apply.");

  // 3. Notify team
  await sendSlackNotification("Langfuse migration rolled back");
}
```

## Resources
- [Langfuse Self-Hosting](https://langfuse.com/docs/deployment/self-host)
- [Langfuse API Reference](https://langfuse.com/docs/api-reference)
- [Data Migration Best Practices](https://langfuse.com/docs/deployment/migration)

## Next Steps
Migration complete. Return to `langfuse-install-auth` for new project setup.
