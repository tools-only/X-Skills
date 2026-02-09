---
name: twinmind-core-workflow-a
description: |
  Execute TwinMind primary workflow: Meeting transcription and summary generation.
  Use when implementing meeting capture, building transcription features,
  or automating meeting documentation.
  Trigger with phrases like "twinmind transcription workflow",
  "meeting transcription", "capture meeting with twinmind".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Core Workflow A: Meeting Transcription & Summary

## Overview
Primary workflow for capturing meetings, generating transcripts, and creating AI summaries.

## Prerequisites
- Completed `twinmind-install-auth` setup
- TwinMind Pro/Enterprise for API access
- Valid API credentials configured
- Audio source available (live or file)

## Instructions

### Step 1: Initialize Meeting Capture

```typescript
// src/workflows/meeting-capture.ts
import { getTwinMindClient } from '../twinmind/client';
import { Transcript, Summary } from '../twinmind/types';

interface MeetingOptions {
  title?: string;
  calendarEventId?: string;
  language?: string;
  enableDiarization?: boolean;
}

export class MeetingCapture {
  private client = getTwinMindClient();

  async startLiveCapture(options: MeetingOptions = {}): Promise<string> {
    const response = await this.client.post('/meetings/live/start', {
      title: options.title || `Meeting ${new Date().toISOString()}`,
      calendar_event_id: options.calendarEventId,
      language: options.language || 'auto',
      diarization: options.enableDiarization ?? true,
      model: 'ear-3',
    });

    return response.data.session_id;
  }

  async stopCapture(sessionId: string): Promise<Transcript> {
    const response = await this.client.post(`/meetings/live/${sessionId}/stop`);
    return response.data.transcript;
  }

  async transcribeRecording(audioUrl: string, options: MeetingOptions = {}): Promise<Transcript> {
    const response = await this.client.post('/transcribe', {
      audio_url: audioUrl,
      title: options.title,
      language: options.language || 'auto',
      diarization: options.enableDiarization ?? true,
      model: 'ear-3',
    });

    return this.waitForTranscript(response.data.transcript_id);
  }

  private async waitForTranscript(transcriptId: string, maxWaitMs = 300000): Promise<Transcript> {
    const startTime = Date.now();
    const pollIntervalMs = 2000;

    while (Date.now() - startTime < maxWaitMs) {
      const response = await this.client.get(`/transcripts/${transcriptId}`);

      if (response.data.status === 'completed') {
        return response.data;
      }

      if (response.data.status === 'failed') {
        throw new Error(`Transcription failed: ${response.data.error}`);
      }

      await new Promise(r => setTimeout(r, pollIntervalMs));
    }

    throw new Error('Transcription timeout');
  }
}
```

### Step 2: Generate AI Summary

```typescript
// src/workflows/summary-generation.ts
export interface SummaryOptions {
  format?: 'brief' | 'detailed' | 'bullet-points';
  includeActionItems?: boolean;
  includeKeyPoints?: boolean;
  maxLength?: number;
}

export class SummaryGenerator {
  private client = getTwinMindClient();

  async generateSummary(
    transcriptId: string,
    options: SummaryOptions = {}
  ): Promise<Summary> {
    const response = await this.client.post('/summarize', {
      transcript_id: transcriptId,
      format: options.format || 'detailed',
      include_action_items: options.includeActionItems ?? true,
      include_key_points: options.includeKeyPoints ?? true,
      max_length: options.maxLength || 500,
    });

    return response.data;
  }

  async generateFollowUpEmail(transcriptId: string): Promise<string> {
    const response = await this.client.post('/generate/follow-up-email', {
      transcript_id: transcriptId,
    });

    return response.data.email_content;
  }

  async generateMeetingNotes(transcriptId: string): Promise<string> {
    const response = await this.client.post('/generate/meeting-notes', {
      transcript_id: transcriptId,
      format: 'markdown',
    });

    return response.data.notes;
  }
}
```

### Step 3: Handle Speaker Identification

```typescript
// src/workflows/speaker-handling.ts
export interface Speaker {
  id: string;
  name?: string;
  email?: string;
  speakingTime: number;
  segments: number;
}

export class SpeakerManager {
  async identifySpeakers(transcript: Transcript, attendees?: string[]): Promise<Speaker[]> {
    const speakers = new Map<string, Speaker>();

    for (const segment of transcript.segments) {
      const speakerId = segment.speaker_id || 'unknown';
      const existing = speakers.get(speakerId);

      if (existing) {
        existing.speakingTime += segment.end - segment.start;
        existing.segments += 1;
      } else {
        speakers.set(speakerId, {
          id: speakerId,
          speakingTime: segment.end - segment.start,
          segments: 1,
        });
      }
    }

    // Match speakers to attendees if provided
    if (attendees && attendees.length > 0) {
      // TwinMind can match speakers to calendar attendees
      const matched = await this.matchToAttendees(
        Array.from(speakers.values()),
        attendees
      );
      return matched;
    }

    return Array.from(speakers.values());
  }

  private async matchToAttendees(speakers: Speaker[], attendees: string[]): Promise<Speaker[]> {
    const client = getTwinMindClient();
    const response = await client.post('/speakers/match', {
      speakers: speakers.map(s => s.id),
      attendees,
    });

    return speakers.map((speaker, idx) => ({
      ...speaker,
      name: response.data.matches[idx]?.name,
      email: response.data.matches[idx]?.email,
    }));
  }
}
```

### Step 4: Complete Workflow Orchestration

```typescript
// src/workflows/full-meeting-workflow.ts
import { MeetingCapture } from './meeting-capture';
import { SummaryGenerator } from './summary-generation';
import { SpeakerManager } from './speaker-handling';

export interface MeetingResult {
  transcriptId: string;
  transcript: Transcript;
  summary: Summary;
  speakers: Speaker[];
  followUpEmail?: string;
  meetingNotes?: string;
}

export async function processMeeting(
  audioUrl: string,
  options: {
    title?: string;
    attendees?: string[];
    generateEmail?: boolean;
    generateNotes?: boolean;
  } = {}
): Promise<MeetingResult> {
  const capture = new MeetingCapture();
  const summaryGen = new SummaryGenerator();
  const speakerMgr = new SpeakerManager();

  // Step 1: Transcribe
  console.log('Starting transcription...');
  const transcript = await capture.transcribeRecording(audioUrl, {
    title: options.title,
    enableDiarization: true,
  });
  console.log(`Transcription complete: ${transcript.id}`);

  // Step 2: Generate summary (parallel with speaker identification)
  console.log('Generating summary and identifying speakers...');
  const [summary, speakers] = await Promise.all([
    summaryGen.generateSummary(transcript.id, {
      format: 'detailed',
      includeActionItems: true,
      includeKeyPoints: true,
    }),
    speakerMgr.identifySpeakers(transcript, options.attendees),
  ]);

  const result: MeetingResult = {
    transcriptId: transcript.id,
    transcript,
    summary,
    speakers,
  };

  // Step 3: Optional outputs
  if (options.generateEmail) {
    result.followUpEmail = await summaryGen.generateFollowUpEmail(transcript.id);
  }

  if (options.generateNotes) {
    result.meetingNotes = await summaryGen.generateMeetingNotes(transcript.id);
  }

  return result;
}

// Example usage
async function main() {
  const result = await processMeeting(
    'https://storage.example.com/meetings/standup-2025-01-15.mp3',
    {
      title: 'Daily Standup - Jan 15',
      attendees: ['alice@example.com', 'bob@example.com'],
      generateEmail: true,
      generateNotes: true,
    }
  );

  console.log('Meeting processed successfully!');
  console.log(`Summary: ${result.summary.summary}`);
  console.log(`Action Items: ${result.summary.action_items.length}`);
  console.log(`Speakers: ${result.speakers.map(s => s.name || s.id).join(', ')}`);
}
```

## Output
- Complete meeting transcript with timestamps
- Speaker-labeled segments
- AI-generated summary
- Extracted action items with assignees
- Optional follow-up email draft
- Optional formatted meeting notes

Example console output:
```
Starting transcription...
Transcription complete: tr_abc123
Generating summary and identifying speakers...
Meeting processed successfully!
Summary: Daily standup covering sprint progress and blockers...
Action Items: 3
Speakers: Alice, Bob, Charlie
```

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| Transcription timeout | Large audio file | Increase maxWaitMs or use async callback |
| Speaker match failed | No calendar data | Provide attendees list manually |
| Summary generation failed | Transcript too short | Ensure minimum 30s of audio |
| Audio format unsupported | Wrong codec | Convert to MP3/WAV/M4A |
| Rate limit exceeded | Too many requests | Implement queue-based processing |

## Audio Format Support

| Format | Extension | Supported | Notes |
|--------|-----------|-----------|-------|
| MP3 | .mp3 | Yes | Recommended |
| WAV | .wav | Yes | Best quality |
| M4A | .m4a | Yes | iOS recordings |
| WebM | .webm | Yes | Browser recordings |
| OGG | .ogg | Yes | Open format |
| FLAC | .flac | Yes | Lossless |

## Resources
- [TwinMind Transcription API](https://twinmind.com/docs/transcription)
- [Ear-3 Model Details](https://twinmind.com/ear-3)
- [Audio Format Guide](https://twinmind.com/docs/audio-formats)

## Next Steps
For action item extraction and follow-up automation, see `twinmind-core-workflow-b`.
