---
name: twinmind-core-workflow-b
description: |
  Execute TwinMind secondary workflow: Action item extraction and follow-up automation.
  Use when automating meeting follow-ups, extracting tasks,
  or integrating with project management tools.
  Trigger with phrases like "twinmind action items",
  "meeting follow-up automation", "extract tasks from meeting".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Core Workflow B: Action Items & Follow-ups

## Overview
Secondary workflow for extracting action items, automating follow-ups, and integrating with task management systems.

## Prerequisites
- Completed `twinmind-core-workflow-a` (transcription)
- Valid transcript or summary available
- Integration tokens for external services (optional)

## Instructions

### Step 1: Extract Action Items

```typescript
// src/workflows/action-items.ts
import { getTwinMindClient } from '../twinmind/client';

export interface ActionItem {
  id: string;
  text: string;
  assignee?: string;
  dueDate?: string;
  priority?: 'high' | 'medium' | 'low';
  context?: string;
  status: 'pending' | 'completed';
}

export interface ActionItemExtractionOptions {
  includeContext?: boolean;
  assignFromSpeakers?: boolean;
  inferDueDates?: boolean;
  maxItems?: number;
}

export class ActionItemExtractor {
  private client = getTwinMindClient();

  async extractFromTranscript(
    transcriptId: string,
    options: ActionItemExtractionOptions = {}
  ): Promise<ActionItem[]> {
    const response = await this.client.post('/extract/action-items', {
      transcript_id: transcriptId,
      include_context: options.includeContext ?? true,
      assign_from_speakers: options.assignFromSpeakers ?? true,
      infer_due_dates: options.inferDueDates ?? true,
      max_items: options.maxItems || 20,
    });

    return response.data.action_items.map((item: any) => ({
      id: item.id,
      text: item.text,
      assignee: item.assignee,
      dueDate: item.due_date,
      priority: this.inferPriority(item),
      context: item.context,
      status: 'pending',
    }));
  }

  private inferPriority(item: any): 'high' | 'medium' | 'low' {
    const text = item.text.toLowerCase();
    if (text.includes('urgent') || text.includes('asap') || text.includes('critical')) {
      return 'high';
    }
    if (text.includes('important') || text.includes('soon')) {
      return 'medium';
    }
    return 'low';
  }

  async categorizeItems(items: ActionItem[]): Promise<Map<string, ActionItem[]>> {
    const categories = new Map<string, ActionItem[]>();

    for (const item of items) {
      const category = this.determineCategory(item);
      if (!categories.has(category)) {
        categories.set(category, []);
      }
      categories.get(category)!.push(item);
    }

    return categories;
  }

  private determineCategory(item: ActionItem): string {
    const text = item.text.toLowerCase();
    if (text.includes('review') || text.includes('check')) return 'Review';
    if (text.includes('create') || text.includes('build') || text.includes('develop')) return 'Development';
    if (text.includes('send') || text.includes('email') || text.includes('contact')) return 'Communication';
    if (text.includes('schedule') || text.includes('meeting')) return 'Meetings';
    if (text.includes('document') || text.includes('write')) return 'Documentation';
    return 'General';
  }
}
```

### Step 2: Automate Follow-up Emails

```typescript
// src/workflows/follow-up.ts
export interface FollowUpEmailOptions {
  recipients: string[];
  includeTranscript?: boolean;
  includeSummary?: boolean;
  includeActionItems?: boolean;
  customNote?: string;
  sendImmediately?: boolean;
}

export interface GeneratedEmail {
  to: string[];
  cc?: string[];
  subject: string;
  body: string;
  attachments?: Array<{
    filename: string;
    content: string;
    contentType: string;
  }>;
}

export class FollowUpAutomation {
  private client = getTwinMindClient();

  async generateFollowUp(
    transcriptId: string,
    options: FollowUpEmailOptions
  ): Promise<GeneratedEmail> {
    const response = await this.client.post('/generate/follow-up', {
      transcript_id: transcriptId,
      recipients: options.recipients,
      include_transcript: options.includeTranscript ?? false,
      include_summary: options.includeSummary ?? true,
      include_action_items: options.includeActionItems ?? true,
      custom_note: options.customNote,
    });

    return {
      to: options.recipients,
      subject: response.data.subject,
      body: response.data.body,
      attachments: response.data.attachments,
    };
  }

  async sendFollowUp(email: GeneratedEmail): Promise<{ messageId: string }> {
    // TwinMind can send via connected email accounts
    const response = await this.client.post('/email/send', {
      to: email.to,
      cc: email.cc,
      subject: email.subject,
      body: email.body,
      attachments: email.attachments,
    });

    return { messageId: response.data.message_id };
  }

  async scheduleFollowUp(
    email: GeneratedEmail,
    sendAt: Date
  ): Promise<{ scheduleId: string }> {
    const response = await this.client.post('/email/schedule', {
      ...email,
      send_at: sendAt.toISOString(),
    });

    return { scheduleId: response.data.schedule_id };
  }
}
```

### Step 3: Integrate with Task Management

```typescript
// src/integrations/task-sync.ts
export interface TaskIntegration {
  name: string;
  createTask(item: ActionItem): Promise<string>;
  updateTask(taskId: string, updates: Partial<ActionItem>): Promise<void>;
}

// Asana integration
export class AsanaIntegration implements TaskIntegration {
  name = 'Asana';

  constructor(private accessToken: string, private projectId: string) {}

  async createTask(item: ActionItem): Promise<string> {
    const response = await fetch('https://app.asana.com/api/1.0/tasks', {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${this.accessToken}`,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        data: {
          name: item.text,
          notes: item.context,
          due_on: item.dueDate,
          projects: [this.projectId],
          assignee: item.assignee, // Email or user ID
        },
      }),
    });

    const data = await response.json();
    return data.data.gid;
  }

  async updateTask(taskId: string, updates: Partial<ActionItem>): Promise<void> {
    await fetch(`https://app.asana.com/api/1.0/tasks/${taskId}`, {
      method: 'PUT',
      headers: {
        'Authorization': `Bearer ${this.accessToken}`,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        data: {
          name: updates.text,
          completed: updates.status === 'completed',
        },
      }),
    });
  }
}

// Linear integration
export class LinearIntegration implements TaskIntegration {
  name = 'Linear';

  constructor(private apiKey: string, private teamId: string) {}

  async createTask(item: ActionItem): Promise<string> {
    const response = await fetch('https://api.linear.app/graphql', {
      method: 'POST',
      headers: {
        'Authorization': this.apiKey,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        query: `
          mutation CreateIssue($input: IssueCreateInput!) {
            issueCreate(input: $input) {
              success
              issue { id identifier }
            }
          }
        `,
        variables: {
          input: {
            teamId: this.teamId,
            title: item.text,
            description: item.context,
            dueDate: item.dueDate,
            priority: this.mapPriority(item.priority),
          },
        },
      }),
    });

    const data = await response.json();
    return data.data.issueCreate.issue.id;
  }

  private mapPriority(priority?: string): number {
    switch (priority) {
      case 'high': return 1;
      case 'medium': return 2;
      case 'low': return 3;
      default: return 4;
    }
  }

  async updateTask(taskId: string, updates: Partial<ActionItem>): Promise<void> {
    // Implementation for Linear update
  }
}

// Factory for integrations
export function getTaskIntegration(
  service: 'asana' | 'linear' | 'notion' | 'jira'
): TaskIntegration {
  switch (service) {
    case 'asana':
      return new AsanaIntegration(
        process.env.ASANA_ACCESS_TOKEN!,
        process.env.ASANA_PROJECT_ID!
      );
    case 'linear':
      return new LinearIntegration(
        process.env.LINEAR_API_KEY!,
        process.env.LINEAR_TEAM_ID!
      );
    default:
      throw new Error(`Integration ${service} not supported`);
  }
}
```

### Step 4: Complete Follow-up Workflow

```typescript
// src/workflows/complete-followup.ts
import { ActionItemExtractor, ActionItem } from './action-items';
import { FollowUpAutomation } from './follow-up';
import { getTaskIntegration, TaskIntegration } from '../integrations/task-sync';

export interface FollowUpWorkflowOptions {
  transcriptId: string;
  attendees: string[];
  taskIntegration?: 'asana' | 'linear' | 'notion' | 'jira';
  sendEmail?: boolean;
  emailDelay?: number; // minutes
}

export interface FollowUpResult {
  actionItems: ActionItem[];
  createdTasks: Array<{ item: ActionItem; externalId: string }>;
  emailSent?: { messageId: string };
  emailScheduled?: { scheduleId: string; sendAt: Date };
}

export async function runFollowUpWorkflow(
  options: FollowUpWorkflowOptions
): Promise<FollowUpResult> {
  const extractor = new ActionItemExtractor();
  const followUp = new FollowUpAutomation();

  // Step 1: Extract action items
  console.log('Extracting action items...');
  const actionItems = await extractor.extractFromTranscript(options.transcriptId, {
    includeContext: true,
    assignFromSpeakers: true,
    inferDueDates: true,
  });
  console.log(`Found ${actionItems.length} action items`);

  const result: FollowUpResult = {
    actionItems,
    createdTasks: [],
  };

  // Step 2: Create tasks in external system
  if (options.taskIntegration) {
    console.log(`Creating tasks in ${options.taskIntegration}...`);
    const integration = getTaskIntegration(options.taskIntegration);

    for (const item of actionItems) {
      try {
        const externalId = await integration.createTask(item);
        result.createdTasks.push({ item, externalId });
        console.log(`Created task: ${item.text.substring(0, 50)}...`);
      } catch (error) {
        console.error(`Failed to create task: ${item.text}`, error);
      }
    }
  }

  // Step 3: Send or schedule follow-up email
  if (options.sendEmail && options.attendees.length > 0) {
    const email = await followUp.generateFollowUp(options.transcriptId, {
      recipients: options.attendees,
      includeSummary: true,
      includeActionItems: true,
    });

    if (options.emailDelay && options.emailDelay > 0) {
      const sendAt = new Date(Date.now() + options.emailDelay * 60 * 1000);
      result.emailScheduled = await followUp.scheduleFollowUp(email, sendAt);
      result.emailScheduled.sendAt = sendAt;
      console.log(`Follow-up email scheduled for ${sendAt.toISOString()}`);
    } else {
      result.emailSent = await followUp.sendFollowUp(email);
      console.log(`Follow-up email sent: ${result.emailSent.messageId}`);
    }
  }

  return result;
}

// Example usage
async function main() {
  const result = await runFollowUpWorkflow({
    transcriptId: 'tr_abc123',
    attendees: ['alice@example.com', 'bob@example.com'],
    taskIntegration: 'linear',
    sendEmail: true,
    emailDelay: 30, // Send 30 minutes after meeting
  });

  console.log('\n=== Follow-up Workflow Complete ===');
  console.log(`Action Items: ${result.actionItems.length}`);
  console.log(`Tasks Created: ${result.createdTasks.length}`);
  if (result.emailScheduled) {
    console.log(`Email scheduled for: ${result.emailScheduled.sendAt}`);
  }
}
```

## Output
- Extracted action items with assignees and due dates
- Tasks created in project management tool
- Follow-up email sent or scheduled
- Complete audit trail

Example console output:
```
Extracting action items...
Found 5 action items
Creating tasks in linear...
Created task: Review Q1 budget spreadsheet...
Created task: Schedule design review meeting...
Created task: Update project timeline document...
Follow-up email scheduled for 2025-01-15T15:30:00Z

=== Follow-up Workflow Complete ===
Action Items: 5
Tasks Created: 5
Email scheduled for: 2025-01-15T15:30:00Z
```

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| No action items found | Transcript too vague | Verify meeting had clear action items |
| Task creation failed | Invalid project/team ID | Check integration credentials |
| Email send failed | Invalid recipients | Verify email addresses |
| Assignee not found | Name mismatch | Map speakers to user accounts |
| Due date parsing failed | Ambiguous date | Use explicit date format |

## Supported Integrations

| Service | Task Creation | Due Dates | Assignees | Priority |
|---------|---------------|-----------|-----------|----------|
| Asana | Yes | Yes | Yes | No |
| Linear | Yes | Yes | Yes | Yes |
| Notion | Yes | Yes | Yes | No |
| Jira | Yes | Yes | Yes | Yes |
| Trello | Yes | Yes | Yes | No |
| Monday.com | Yes | Yes | Yes | Yes |

## Resources
- [TwinMind Action Items API](https://twinmind.com/docs/action-items)
- [Follow-up Automation](https://twinmind.com/docs/follow-up)
- [Asana API](https://developers.asana.com)
- [Linear API](https://developers.linear.app)

## Next Steps
For troubleshooting issues, see `twinmind-common-errors`.
