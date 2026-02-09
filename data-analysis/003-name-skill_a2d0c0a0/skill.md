---
name: twinmind-migration-deep-dive
description: |
  Comprehensive migration guide from other meeting AI tools to TwinMind.
  Use when migrating from Otter.ai, Fireflies, Rev, or other transcription services.
  Trigger with phrases like "migrate to twinmind", "switch from otter",
  "twinmind migration", "move to twinmind", "replace fireflies".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Migration Deep Dive

## Overview
Comprehensive guide for migrating from other meeting AI services to TwinMind.

## Prerequisites
- TwinMind Pro/Enterprise account
- Export access on source platform
- Understanding of data formats
- Sufficient storage for migration

## Supported Migration Sources

| Source | Export Format | API Available | Difficulty |
|--------|--------------|---------------|------------|
| Otter.ai | TXT, JSON | Limited | Medium |
| Fireflies.ai | JSON, CSV | Yes | Easy |
| Rev.ai | JSON | Yes | Easy |
| Trint | DOCX, JSON | Yes | Medium |
| Descript | JSON | Yes | Medium |
| Zoom transcripts | VTT, TXT | Yes | Easy |
| Google Meet | SRT | Limited | Medium |
| Microsoft Teams | DOCX, VTT | Yes | Medium |
| Grain | JSON | Yes | Easy |
| Chorus.ai | JSON | Limited | Hard |

## Instructions

### Step 1: Assess Current Data

```typescript
// scripts/migration-assessment.ts
import * as fs from 'fs';
import * as path from 'path';

interface MigrationAssessment {
  source: string;
  totalTranscripts: number;
  totalDurationHours: number;
  dateRange: { earliest: Date; latest: Date };
  dataFormats: string[];
  estimatedMigrationTime: string;
  potentialIssues: string[];
}

async function assessMigration(
  source: 'otter' | 'fireflies' | 'rev' | 'zoom' | 'teams',
  dataPath: string
): Promise<MigrationAssessment> {
  const files = fs.readdirSync(dataPath);
  const issues: string[] = [];

  let totalTranscripts = 0;
  let totalDuration = 0;
  let earliest = new Date();
  let latest = new Date(0);
  const formats = new Set<string>();

  for (const file of files) {
    const ext = path.extname(file).toLowerCase();
    formats.add(ext);

    const filePath = path.join(dataPath, file);
    const stat = fs.statSync(filePath);

    if (ext === '.json') {
      try {
        const content = JSON.parse(fs.readFileSync(filePath, 'utf-8'));
        totalTranscripts++;

        // Extract duration if available
        if (content.duration) {
          totalDuration += content.duration;
        }

        // Extract date
        const date = new Date(content.created_at || content.date || stat.mtime);
        if (date < earliest) earliest = date;
        if (date > latest) latest = date;
      } catch (e) {
        issues.push(`Failed to parse: ${file}`);
      }
    }
  }

  // Estimate migration time (rough: 1 transcript per second + processing)
  const estimatedMinutes = Math.ceil(totalTranscripts / 60) + Math.ceil(totalDuration / 3600 * 10);

  return {
    source,
    totalTranscripts,
    totalDurationHours: totalDuration / 3600,
    dateRange: { earliest, latest },
    dataFormats: Array.from(formats),
    estimatedMigrationTime: `${estimatedMinutes} minutes`,
    potentialIssues: issues,
  };
}

// Run assessment
const assessment = await assessMigration('otter', './exports/otter');
console.log('Migration Assessment:');
console.log(JSON.stringify(assessment, null, 2));
```

### Step 2: Export Data from Source

#### Otter.ai Export

```typescript
// scripts/exporters/otter.ts
import puppeteer from 'puppeteer';

async function exportFromOtter(credentials: { email: string; password: string }) {
  const browser = await puppeteer.launch({ headless: false });
  const page = await browser.newPage();

  // Login
  await page.goto('https://otter.ai/login');
  await page.type('input[name="email"]', credentials.email);
  await page.type('input[name="password"]', credentials.password);
  await page.click('button[type="submit"]');
  await page.waitForNavigation();

  // Navigate to conversations
  await page.goto('https://otter.ai/my-notes');
  await page.waitForSelector('.conversation-list');

  // Get all conversation IDs
  const conversationIds = await page.$$eval(
    '.conversation-item',
    items => items.map(item => item.getAttribute('data-id'))
  );

  console.log(`Found ${conversationIds.length} conversations`);

  // Export each conversation
  const exports = [];
  for (const id of conversationIds) {
    const data = await page.evaluate(async (convId) => {
      const response = await fetch(`/api/conversation/${convId}/export`);
      return response.json();
    }, id);

    exports.push(data);
  }

  await browser.close();
  return exports;
}
```

#### Fireflies.ai Export

```typescript
// scripts/exporters/fireflies.ts
const FIREFLIES_API_URL = 'https://api.fireflies.ai/graphql';

async function exportFromFireflies(apiKey: string): Promise<any[]> {
  const query = `
    query GetAllTranscripts {
      transcripts {
        id
        title
        date
        duration
        transcript_text
        summary
        action_items {
          text
          assignee
        }
        speakers {
          name
          id
        }
        sentences {
          text
          speaker_id
          start_time
          end_time
        }
      }
    }
  `;

  const response = await fetch(FIREFLIES_API_URL, {
    method: 'POST',
    headers: {
      'Authorization': `Bearer ${apiKey}`,
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ query }),
  });

  const data = await response.json();
  return data.data.transcripts;
}
```

#### Rev.ai Export

```typescript
// scripts/exporters/rev.ts
async function exportFromRev(apiKey: string): Promise<any[]> {
  // Get all jobs
  const jobsResponse = await fetch('https://api.rev.ai/speechtotext/v1/jobs', {
    headers: { 'Authorization': `Bearer ${apiKey}` },
  });

  const jobs = await jobsResponse.json();
  const transcripts = [];

  for (const job of jobs) {
    if (job.status === 'transcribed') {
      // Get transcript
      const transcriptResponse = await fetch(
        `https://api.rev.ai/speechtotext/v1/jobs/${job.id}/transcript`,
        {
          headers: {
            'Authorization': `Bearer ${apiKey}`,
            'Accept': 'application/vnd.rev.transcript.v1.0+json',
          },
        }
      );

      const transcript = await transcriptResponse.json();
      transcripts.push({
        id: job.id,
        name: job.name,
        created_on: job.created_on,
        duration_seconds: job.duration_seconds,
        transcript,
      });
    }
  }

  return transcripts;
}
```

### Step 3: Transform Data to TwinMind Format

```typescript
// scripts/migration/transform.ts
import { Transcript, Segment, Speaker } from '../../src/twinmind/types';

interface TransformOptions {
  preserveTimestamps: boolean;
  preserveSpeakers: boolean;
  generateNewIds: boolean;
}

// Otter.ai to TwinMind
export function transformOtterToTwinMind(otterData: any): Transcript {
  return {
    id: `tm_imported_${otterData.id}`,
    text: otterData.text || otterData.transcript_text,
    duration_seconds: otterData.duration || 0,
    language: otterData.language || 'en',
    created_at: otterData.created_at || new Date().toISOString(),
    metadata: {
      imported_from: 'otter.ai',
      original_id: otterData.id,
    },
    speakers: otterData.speakers?.map((s: any, i: number) => ({
      id: `spk_${i}`,
      name: s.name || s.speaker_name || `Speaker ${i + 1}`,
    })),
    segments: otterData.segments?.map((seg: any) => ({
      start: seg.start_time || seg.start,
      end: seg.end_time || seg.end,
      text: seg.text || seg.transcript,
      speaker_id: seg.speaker_id,
      confidence: seg.confidence || 0.95,
    })),
  };
}

// Fireflies.ai to TwinMind
export function transformFirefliesToTwinMind(ffData: any): Transcript {
  return {
    id: `tm_imported_${ffData.id}`,
    title: ffData.title,
    text: ffData.transcript_text,
    duration_seconds: ffData.duration,
    language: 'en',  // Fireflies doesn't export language
    created_at: ffData.date,
    metadata: {
      imported_from: 'fireflies.ai',
      original_id: ffData.id,
      summary: ffData.summary,
      action_items: ffData.action_items,
    },
    speakers: ffData.speakers?.map((s: any) => ({
      id: s.id,
      name: s.name,
    })),
    segments: ffData.sentences?.map((sent: any) => ({
      start: sent.start_time,
      end: sent.end_time,
      text: sent.text,
      speaker_id: sent.speaker_id,
      confidence: 0.95,
    })),
  };
}

// Rev.ai to TwinMind
export function transformRevToTwinMind(revData: any): Transcript {
  const monologues = revData.transcript.monologues || [];

  const segments: Segment[] = [];
  const speakers = new Map<number, string>();

  for (const monologue of monologues) {
    const speakerId = monologue.speaker;
    if (!speakers.has(speakerId)) {
      speakers.set(speakerId, `Speaker ${speakerId + 1}`);
    }

    for (const element of monologue.elements) {
      if (element.type === 'text') {
        segments.push({
          start: element.ts,
          end: element.end_ts,
          text: element.value,
          speaker_id: `spk_${speakerId}`,
          confidence: element.confidence || 0.95,
        });
      }
    }
  }

  return {
    id: `tm_imported_${revData.id}`,
    title: revData.name,
    text: segments.map(s => s.text).join(' '),
    duration_seconds: revData.duration_seconds,
    language: 'en',
    created_at: revData.created_on,
    metadata: {
      imported_from: 'rev.ai',
      original_id: revData.id,
    },
    speakers: Array.from(speakers.entries()).map(([id, name]) => ({
      id: `spk_${id}`,
      name,
    })),
    segments,
  };
}

// Zoom VTT to TwinMind
export function transformZoomVTTToTwinMind(
  vttContent: string,
  metadata: { title: string; date: string; duration: number }
): Transcript {
  const segments: Segment[] = [];
  const lines = vttContent.split('\n');

  let currentSegment: Partial<Segment> = {};

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();

    // Parse timestamp line: 00:00:00.000 --> 00:00:02.000
    if (line.includes('-->')) {
      const [start, end] = line.split('-->').map(t => parseVTTTimestamp(t.trim()));
      currentSegment.start = start;
      currentSegment.end = end;
    }
    // Parse speaker and text
    else if (line && !line.match(/^\d+$/) && currentSegment.start !== undefined) {
      // Format: "Speaker Name: Text" or just "Text"
      const speakerMatch = line.match(/^(.+?):\s*(.*)$/);

      if (speakerMatch) {
        currentSegment.speaker_id = speakerMatch[1];
        currentSegment.text = speakerMatch[2];
      } else {
        currentSegment.text = line;
      }

      if (currentSegment.text) {
        segments.push({
          start: currentSegment.start,
          end: currentSegment.end || currentSegment.start + 1,
          text: currentSegment.text,
          speaker_id: currentSegment.speaker_id,
          confidence: 0.9,
        } as Segment);
        currentSegment = {};
      }
    }
  }

  // Extract unique speakers
  const speakerIds = [...new Set(segments.map(s => s.speaker_id).filter(Boolean))];

  return {
    id: `tm_imported_zoom_${Date.now()}`,
    title: metadata.title,
    text: segments.map(s => s.text).join(' '),
    duration_seconds: metadata.duration,
    language: 'en',
    created_at: metadata.date,
    metadata: {
      imported_from: 'zoom',
    },
    speakers: speakerIds.map((id, i) => ({ id: `spk_${i}`, name: id as string })),
    segments,
  };
}

function parseVTTTimestamp(timestamp: string): number {
  const parts = timestamp.split(':');
  const seconds = parts.pop()!;
  const minutes = parts.pop() || '0';
  const hours = parts.pop() || '0';

  return (
    parseFloat(hours) * 3600 +
    parseFloat(minutes) * 60 +
    parseFloat(seconds)
  );
}
```

### Step 4: Import to TwinMind

```typescript
// scripts/migration/import.ts
import { getTwinMindClient } from '../../src/twinmind/client';
import { Transcript } from '../../src/twinmind/types';

interface ImportResult {
  successful: number;
  failed: number;
  errors: Array<{ id: string; error: string }>;
}

export async function importToTwinMind(
  transcripts: Transcript[],
  options: {
    batchSize?: number;
    delayMs?: number;
    dryRun?: boolean;
  } = {}
): Promise<ImportResult> {
  const client = getTwinMindClient();
  const batchSize = options.batchSize || 10;
  const delayMs = options.delayMs || 1000;

  const result: ImportResult = {
    successful: 0,
    failed: 0,
    errors: [],
  };

  console.log(`Starting import of ${transcripts.length} transcripts...`);

  for (let i = 0; i < transcripts.length; i += batchSize) {
    const batch = transcripts.slice(i, i + batchSize);

    for (const transcript of batch) {
      try {
        if (options.dryRun) {
          console.log(`[DRY RUN] Would import: ${transcript.title || transcript.id}`);
          result.successful++;
          continue;
        }

        // Check if already imported
        const existing = await client.get('/transcripts', {
          params: {
            'metadata.original_id': transcript.metadata?.original_id,
          },
        });

        if (existing.data.length > 0) {
          console.log(`Skipping duplicate: ${transcript.title || transcript.id}`);
          continue;
        }

        // Import transcript
        await client.post('/transcripts/import', {
          text: transcript.text,
          title: transcript.title,
          duration_seconds: transcript.duration_seconds,
          language: transcript.language,
          created_at: transcript.created_at,
          metadata: transcript.metadata,
          speakers: transcript.speakers,
          segments: transcript.segments,
        });

        result.successful++;
        console.log(`Imported: ${transcript.title || transcript.id}`);
      } catch (error: any) {
        result.failed++;
        result.errors.push({
          id: transcript.id,
          error: error.message,
        });
        console.error(`Failed: ${transcript.id} - ${error.message}`);
      }
    }

    // Rate limit delay between batches
    if (i + batchSize < transcripts.length) {
      await new Promise(r => setTimeout(r, delayMs));
    }
  }

  return result;
}
```

### Step 5: Migration Verification

```typescript
// scripts/migration/verify.ts
interface VerificationResult {
  sourceCount: number;
  importedCount: number;
  matchRate: number;
  missingTranscripts: string[];
  contentMismatches: Array<{
    id: string;
    issue: string;
  }>;
}

export async function verifyMigration(
  sourceTranscripts: Transcript[],
  importedIds: string[]
): Promise<VerificationResult> {
  const client = getTwinMindClient();
  const missingTranscripts: string[] = [];
  const contentMismatches: Array<{ id: string; issue: string }> = [];

  for (const source of sourceTranscripts) {
    // Find in TwinMind
    const imported = await client.get('/transcripts', {
      params: {
        'metadata.original_id': source.metadata?.original_id,
      },
    });

    if (imported.data.length === 0) {
      missingTranscripts.push(source.id);
      continue;
    }

    const importedTranscript = imported.data[0];

    // Verify content
    if (importedTranscript.text !== source.text) {
      const sourceWords = source.text.split(/\s+/).length;
      const importedWords = importedTranscript.text.split(/\s+/).length;

      if (Math.abs(sourceWords - importedWords) > sourceWords * 0.1) {
        contentMismatches.push({
          id: source.id,
          issue: `Word count mismatch: ${sourceWords} vs ${importedWords}`,
        });
      }
    }

    // Verify duration
    if (Math.abs(importedTranscript.duration_seconds - source.duration_seconds) > 1) {
      contentMismatches.push({
        id: source.id,
        issue: `Duration mismatch: ${source.duration_seconds}s vs ${importedTranscript.duration_seconds}s`,
      });
    }
  }

  return {
    sourceCount: sourceTranscripts.length,
    importedCount: sourceTranscripts.length - missingTranscripts.length,
    matchRate: (sourceTranscripts.length - missingTranscripts.length) / sourceTranscripts.length,
    missingTranscripts,
    contentMismatches,
  };
}
```

### Step 6: Post-Migration Checklist

```markdown
## Post-Migration Checklist

### Data Verification
- [ ] All transcripts imported (compare counts)
- [ ] Spot-check 10 random transcripts for content accuracy
- [ ] Verify speaker labels preserved
- [ ] Verify timestamps are accurate
- [ ] Check metadata (dates, durations) are correct

### Feature Verification
- [ ] Search finds imported transcripts
- [ ] Summaries can be generated for imported data
- [ ] Action items can be extracted
- [ ] Calendar integration works with imported meetings

### User Acceptance
- [ ] Key users can find their transcripts
- [ ] Historical data accessible
- [ ] Reports include migrated data

### Cleanup
- [ ] Remove temporary import files
- [ ] Update user documentation
- [ ] Archive source platform export
- [ ] Cancel source platform subscription

### Performance
- [ ] Search performance acceptable with migrated data
- [ ] No increase in error rates
- [ ] Memory usage within limits
```

## Output
- Migration assessment tool
- Platform-specific exporters
- Data transformation functions
- TwinMind import utility
- Verification scripts
- Post-migration checklist

## Platform Comparison

| Feature | TwinMind | Otter.ai | Fireflies |
|---------|----------|----------|-----------|
| WER | 5.26% | ~8% | ~7% |
| Languages | 140+ | 30+ | 60+ |
| On-device | Yes | No | No |
| No audio storage | Yes | No | No |
| API access | Pro+ | Business | Pro+ |
| Custom models | Enterprise | No | No |

## Error Handling

| Issue | Cause | Solution |
|-------|-------|----------|
| Export rate limited | Too many requests | Add delays |
| Data format changed | API version | Update transformer |
| Import duplicate | Already migrated | Skip or overwrite |
| Missing speakers | Source didn't export | Set as unknown |

## Resources
- [TwinMind Import API](https://twinmind.com/docs/import)
- [Data Export Best Practices](https://twinmind.com/docs/migration)
- [Platform Comparison](https://twinmind.com/compare)

## Congratulations!
You have completed the TwinMind skill pack. For questions, visit [TwinMind Support](https://twinmind.com/support).
