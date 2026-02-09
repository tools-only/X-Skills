---
name: obsidian-migration-deep-dive
description: |
  Execute major Obsidian plugin rewrites and migration strategies.
  Use when migrating to or from Obsidian, performing major plugin rewrites,
  or re-platforming existing note systems to Obsidian.
  Trigger with phrases like "migrate to obsidian", "obsidian migration",
  "convert notes to obsidian", "obsidian replatform".
allowed-tools: Read, Write, Edit, Bash(npm:*), Bash(node:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Obsidian Migration Deep Dive

## Overview
Comprehensive guide for migrating to Obsidian from other note-taking apps, or performing major plugin architecture rewrites.

## Prerequisites
- Source data access
- Understanding of Obsidian vault structure
- Node.js for scripted migrations
- Backup of source data

## Migration Types

| Type | Complexity | Duration | Risk |
|------|-----------|----------|------|
| Single app import | Low | Hours | Low |
| Multi-source merge | Medium | Days | Medium |
| Plugin major rewrite | Medium | Weeks | Medium |
| Enterprise migration | High | Months | High |

## Instructions

### Step 1: Pre-Migration Assessment
```typescript
// scripts/migration-assessment.ts
interface MigrationAssessment {
  sourceSystem: string;
  noteCount: number;
  attachmentCount: number;
  totalSize: number;
  linkCount: number;
  tagCount: number;
  uniqueTags: string[];
  folderStructure: string[];
  issues: MigrationIssue[];
}

interface MigrationIssue {
  type: 'encoding' | 'format' | 'link' | 'attachment' | 'metadata';
  severity: 'warning' | 'error';
  description: string;
  affectedFiles: string[];
}

async function assessMigration(sourcePath: string): Promise<MigrationAssessment> {
  const assessment: MigrationAssessment = {
    sourceSystem: 'unknown',
    noteCount: 0,
    attachmentCount: 0,
    totalSize: 0,
    linkCount: 0,
    tagCount: 0,
    uniqueTags: [],
    folderStructure: [],
    issues: [],
  };

  // Scan source directory
  // Count files, measure sizes
  // Identify formats and potential issues

  return assessment;
}

// Generate report
function generateAssessmentReport(assessment: MigrationAssessment): string {
  return `
# Migration Assessment Report

## Source System: ${assessment.sourceSystem}

### Content Summary
- Notes: ${assessment.noteCount}
- Attachments: ${assessment.attachmentCount}
- Total Size: ${(assessment.totalSize / 1024 / 1024).toFixed(2)} MB
- Links: ${assessment.linkCount}
- Tags: ${assessment.tagCount} (${assessment.uniqueTags.length} unique)

### Folder Structure
${assessment.folderStructure.map(f => `- ${f}`).join('\n')}

### Issues Found
${assessment.issues.map(i => `- [${i.severity.toUpperCase()}] ${i.type}: ${i.description}`).join('\n')}

### Recommendations
${assessment.issues.length === 0 ? '- No issues found, proceed with migration' : '- Address issues before migration'}
  `;
}
```

### Step 2: Format Converters
```typescript
// scripts/converters/evernote.ts
import * as fs from 'fs';
import * as path from 'path';
import { parseStringPromise } from 'xml2js';

interface EvernoteNote {
  title: string;
  content: string;
  created: string;
  updated: string;
  tags: string[];
  attachments: EvernoteAttachment[];
}

interface EvernoteAttachment {
  filename: string;
  mime: string;
  data: string; // base64
}

export async function convertEvernoteExport(
  enexPath: string,
  outputPath: string
): Promise<{ notes: number; attachments: number }> {
  const content = fs.readFileSync(enexPath, 'utf-8');
  const parsed = await parseStringPromise(content);
  const notes = parsed['en-export']?.note || [];

  let noteCount = 0;
  let attachmentCount = 0;

  for (const note of notes) {
    const converted = convertEvernoteNote(note);
    const fileName = sanitizeFileName(converted.title) + '.md';
    const filePath = path.join(outputPath, fileName);

    // Convert HTML content to Markdown
    const markdown = convertHtmlToMarkdown(converted.content);

    // Add frontmatter
    const frontmatter = `---
title: ${converted.title}
created: ${converted.created}
updated: ${converted.updated}
tags: [${converted.tags.join(', ')}]
source: evernote
---

`;

    fs.writeFileSync(filePath, frontmatter + markdown);
    noteCount++;

    // Handle attachments
    for (const attachment of converted.attachments) {
      const attachmentPath = path.join(outputPath, 'attachments', attachment.filename);
      const data = Buffer.from(attachment.data, 'base64');
      fs.writeFileSync(attachmentPath, data);
      attachmentCount++;
    }
  }

  return { notes: noteCount, attachments: attachmentCount };
}

function convertEvernoteNote(note: any): EvernoteNote {
  return {
    title: note.title?.[0] || 'Untitled',
    content: note.content?.[0] || '',
    created: formatDate(note.created?.[0]),
    updated: formatDate(note.updated?.[0]),
    tags: note.tag || [],
    attachments: extractAttachments(note.resource || []),
  };
}

// scripts/converters/notion.ts
export async function convertNotionExport(
  notionPath: string,
  outputPath: string
): Promise<{ notes: number; databases: number }> {
  // Notion exports as nested folders with markdown/CSV
  // Walk directory and convert

  let noteCount = 0;
  let databaseCount = 0;

  // Implementation...

  return { notes: noteCount, databases: databaseCount };
}

// scripts/converters/roam.ts
export async function convertRoamExport(
  roamJsonPath: string,
  outputPath: string
): Promise<{ pages: number; blocks: number }> {
  const content = fs.readFileSync(roamJsonPath, 'utf-8');
  const roamData = JSON.parse(content);

  let pageCount = 0;
  let blockCount = 0;

  for (const page of roamData) {
    const markdown = convertRoamPage(page);
    const fileName = sanitizeFileName(page.title) + '.md';

    fs.writeFileSync(path.join(outputPath, fileName), markdown);
    pageCount++;
    blockCount += countBlocks(page);
  }

  return { pages: pageCount, blocks: blockCount };
}

function convertRoamPage(page: any): string {
  const lines: string[] = [`# ${page.title}`, ''];

  if (page.children) {
    for (const block of page.children) {
      lines.push(...convertRoamBlock(block, 0));
    }
  }

  return lines.join('\n');
}

function convertRoamBlock(block: any, depth: number): string[] {
  const lines: string[] = [];
  const indent = '  '.repeat(depth);
  const content = convertRoamSyntax(block.string || '');

  lines.push(`${indent}- ${content}`);

  if (block.children) {
    for (const child of block.children) {
      lines.push(...convertRoamBlock(child, depth + 1));
    }
  }

  return lines;
}

function convertRoamSyntax(text: string): string {
  // Convert Roam-specific syntax to Obsidian
  return text
    .replace(/\[\[([^\]]+)\]\]/g, '[[$1]]') // Links same
    .replace(/\(\(([^)]+)\)\)/g, '^$1') // Block refs to block IDs
    .replace(/#\[\[([^\]]+)\]\]/g, '#$1') // Tag pages to tags
    .replace(/{{embed: \[\[([^\]]+)\]\]}}/g, '![[$ 1]]'); // Embeds
}
```

### Step 3: Link Migration
```typescript
// scripts/migrate-links.ts
import * as fs from 'fs';
import * as path from 'path';
import * as glob from 'glob';

interface LinkMapping {
  original: string;
  converted: string;
  type: 'internal' | 'external' | 'attachment';
}

export class LinkMigrator {
  private linkMappings: Map<string, LinkMapping> = new Map();
  private orphanedLinks: string[] = [];

  async buildLinkIndex(vaultPath: string): Promise<void> {
    const files = glob.sync('**/*.md', { cwd: vaultPath });

    for (const file of files) {
      const baseName = path.basename(file, '.md');
      this.linkMappings.set(baseName.toLowerCase(), {
        original: baseName,
        converted: baseName,
        type: 'internal',
      });
    }
  }

  async migrateLinks(vaultPath: string): Promise<{
    updated: number;
    orphaned: string[];
  }> {
    const files = glob.sync('**/*.md', { cwd: vaultPath });
    let updatedCount = 0;

    for (const file of files) {
      const filePath = path.join(vaultPath, file);
      let content = fs.readFileSync(filePath, 'utf-8');
      let modified = false;

      // Find all wiki-style links
      const linkRegex = /\[\[([^\]|]+)(\|[^\]]+)?\]\]/g;
      let match;

      while ((match = linkRegex.exec(content)) !== null) {
        const originalLink = match[1];
        const alias = match[2] || '';
        const resolvedLink = this.resolveLink(originalLink);

        if (resolvedLink !== originalLink) {
          const newLink = `[[${resolvedLink}${alias}]]`;
          content = content.replace(match[0], newLink);
          modified = true;
        }
      }

      if (modified) {
        fs.writeFileSync(filePath, content);
        updatedCount++;
      }
    }

    return {
      updated: updatedCount,
      orphaned: this.orphanedLinks,
    };
  }

  private resolveLink(link: string): string {
    // Try exact match
    const mapping = this.linkMappings.get(link.toLowerCase());
    if (mapping) {
      return mapping.converted;
    }

    // Try without path
    const baseName = path.basename(link);
    const baseMapping = this.linkMappings.get(baseName.toLowerCase());
    if (baseMapping) {
      return baseMapping.converted;
    }

    // Mark as orphaned
    if (!this.orphanedLinks.includes(link)) {
      this.orphanedLinks.push(link);
    }

    return link;
  }

  async createOrphanedLinksReport(vaultPath: string): Promise<void> {
    const report = `# Orphaned Links Report

These links could not be resolved during migration:

${this.orphanedLinks.map(link => `- [[${link}]]`).join('\n')}

## Actions Needed
- Create missing notes
- Update or remove broken links
- Check for renamed files
`;

    fs.writeFileSync(
      path.join(vaultPath, '_migration', 'orphaned-links.md'),
      report
    );
  }
}
```

### Step 4: Batch Migration Script
```typescript
// scripts/migrate.ts
import * as fs from 'fs';
import * as path from 'path';
import { convertEvernoteExport } from './converters/evernote';
import { convertNotionExport } from './converters/notion';
import { convertRoamExport } from './converters/roam';
import { LinkMigrator } from './migrate-links';

interface MigrationConfig {
  source: {
    type: 'evernote' | 'notion' | 'roam' | 'markdown';
    path: string;
  };
  target: {
    vaultPath: string;
    createBackup: boolean;
  };
  options: {
    preserveFolderStructure: boolean;
    convertTags: boolean;
    migrateAttachments: boolean;
    fixLinks: boolean;
    dryRun: boolean;
  };
}

async function runMigration(config: MigrationConfig): Promise<void> {
  console.log('Starting migration...');
  console.log(`Source: ${config.source.type} from ${config.source.path}`);
  console.log(`Target: ${config.target.vaultPath}`);

  // Create backup if requested
  if (config.target.createBackup && !config.options.dryRun) {
    const backupPath = `${config.target.vaultPath}-backup-${Date.now()}`;
    fs.cpSync(config.target.vaultPath, backupPath, { recursive: true });
    console.log(`Backup created at: ${backupPath}`);
  }

  // Create migration folder for reports
  const migrationFolder = path.join(config.target.vaultPath, '_migration');
  if (!config.options.dryRun) {
    fs.mkdirSync(migrationFolder, { recursive: true });
  }

  // Run appropriate converter
  let result: { notes: number; [key: string]: number };

  switch (config.source.type) {
    case 'evernote':
      result = await convertEvernoteExport(
        config.source.path,
        config.target.vaultPath
      );
      break;
    case 'notion':
      result = await convertNotionExport(
        config.source.path,
        config.target.vaultPath
      );
      break;
    case 'roam':
      result = await convertRoamExport(
        config.source.path,
        config.target.vaultPath
      );
      break;
    default:
      throw new Error(`Unsupported source type: ${config.source.type}`);
  }

  console.log(`Converted ${result.notes} notes`);

  // Fix links if requested
  if (config.options.fixLinks) {
    const linkMigrator = new LinkMigrator();
    await linkMigrator.buildLinkIndex(config.target.vaultPath);
    const linkResult = await linkMigrator.migrateLinks(config.target.vaultPath);

    console.log(`Updated links in ${linkResult.updated} files`);
    console.log(`Found ${linkResult.orphaned.length} orphaned links`);

    if (linkResult.orphaned.length > 0) {
      await linkMigrator.createOrphanedLinksReport(config.target.vaultPath);
    }
  }

  // Generate migration report
  const report = generateMigrationReport(config, result);
  if (!config.options.dryRun) {
    fs.writeFileSync(
      path.join(migrationFolder, 'migration-report.md'),
      report
    );
  }

  console.log('Migration complete!');
}

function generateMigrationReport(
  config: MigrationConfig,
  result: { notes: number; [key: string]: number }
): string {
  return `# Migration Report

## Summary
- **Date:** ${new Date().toISOString()}
- **Source:** ${config.source.type}
- **Notes migrated:** ${result.notes}

## Configuration
\`\`\`json
${JSON.stringify(config, null, 2)}
\`\`\`

## Results
${Object.entries(result).map(([key, value]) => `- ${key}: ${value}`).join('\n')}

## Next Steps
1. Review migrated content
2. Check orphaned links report
3. Test in Obsidian
4. Remove _migration folder when satisfied
`;
}

// Run migration
const config: MigrationConfig = {
  source: {
    type: 'evernote',
    path: '/path/to/export.enex',
  },
  target: {
    vaultPath: '/path/to/obsidian/vault',
    createBackup: true,
  },
  options: {
    preserveFolderStructure: true,
    convertTags: true,
    migrateAttachments: true,
    fixLinks: true,
    dryRun: false,
  },
};

runMigration(config).catch(console.error);
```

### Step 5: Plugin Architecture Migration
```typescript
// For major plugin rewrites

interface PluginMigrationPlan {
  currentVersion: string;
  targetVersion: string;
  phases: MigrationPhase[];
  rollbackPlan: string;
}

interface MigrationPhase {
  name: string;
  description: string;
  changes: string[];
  breakingChanges: string[];
  migrationSteps: string[];
}

const migrationPlan: PluginMigrationPlan = {
  currentVersion: '1.x',
  targetVersion: '2.0',
  phases: [
    {
      name: 'Phase 1: Settings Migration',
      description: 'Migrate settings to new format',
      changes: [
        'New settings schema',
        'Split monolithic settings into categories',
      ],
      breakingChanges: [
        'Old settings format deprecated',
      ],
      migrationSteps: [
        'Load old settings on upgrade',
        'Transform to new format',
        'Save new settings',
        'Backup old settings',
      ],
    },
    {
      name: 'Phase 2: API Changes',
      description: 'Update internal APIs',
      changes: [
        'New service-based architecture',
        'Event-driven communication',
      ],
      breakingChanges: [
        'Direct vault access deprecated',
        'Command IDs changed',
      ],
      migrationSteps: [
        'Update command registrations',
        'Migrate to new services',
        'Update event handlers',
      ],
    },
  ],
  rollbackPlan: 'Install previous version from GitHub releases',
};
```

## Output
- Pre-migration assessment
- Format converters for common apps
- Link migration and fixing
- Batch migration scripts
- Migration reports

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Encoding errors | Non-UTF8 content | Detect and convert encoding |
| Broken links | Renamed/deleted files | Generate orphaned links report |
| Missing attachments | Export incomplete | Re-export with attachments |
| Duplicate files | Same name different folders | Add path prefix |

## Examples

### Command Line Usage
```bash
# Install dependencies
npm install xml2js glob

# Run migration
npx ts-node scripts/migrate.ts

# Dry run first
MIGRATION_DRY_RUN=true npx ts-node scripts/migrate.ts
```

### Post-Migration Checklist
```markdown
## Post-Migration Checklist

- [ ] Open vault in Obsidian
- [ ] Check random sample of notes (10-20)
- [ ] Verify links resolve correctly
- [ ] Check attachments display
- [ ] Verify tags imported
- [ ] Test search functionality
- [ ] Check folder structure
- [ ] Review orphaned links report
- [ ] Delete _migration folder
- [ ] Update any external integrations
```

## Resources
- [Obsidian Import/Export](https://help.obsidian.md/import)
- [Evernote Export Format](https://evernote.com/blog/how-evernotes-xml-export-format-works)
- [Notion Export](https://www.notion.so/help/export-your-content)

## Flagship+ Skills
Migration complete! You now have comprehensive Obsidian plugin development skills.
