---
name: twinmind-performance-tuning
description: |
  Optimize TwinMind transcription accuracy and processing speed.
  Use when improving transcription quality, reducing latency,
  or tuning model parameters for specific use cases.
  Trigger with phrases like "twinmind performance", "improve transcription accuracy",
  "faster twinmind", "optimize twinmind", "transcription quality".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Performance Tuning

## Overview
Optimize TwinMind for better transcription accuracy, faster processing, and improved user experience.

## Prerequisites
- TwinMind Pro/Enterprise account
- Understanding of audio processing concepts
- Access to quality metrics and logs

## Instructions

### Step 1: Understand Performance Metrics

```typescript
// src/twinmind/metrics/performance.ts
export interface TranscriptionMetrics {
  // Accuracy metrics
  wordErrorRate: number;         // WER - lower is better (Ear-3: ~5.26%)
  diarizationErrorRate: number;  // DER - speaker labeling accuracy (~3.8%)
  confidenceScore: number;       // Average confidence (0-1)

  // Timing metrics
  processingTime: number;        // Time to complete (ms)
  realtimeFactor: number;        // Processing time / audio duration
  firstWordLatency: number;      // Time to first result (streaming)

  // Quality metrics
  speakerCount: number;          // Detected speakers
  languageDetected: string;      // Detected language
  noiseLevel: string;            // low/medium/high
}

export async function getTranscriptionMetrics(transcriptId: string): Promise<TranscriptionMetrics> {
  const client = getTwinMindClient();
  const response = await client.get(`/transcripts/${transcriptId}/metrics`);
  return response.data;
}

// Log and analyze metrics
export function analyzePerformance(metrics: TranscriptionMetrics): string[] {
  const recommendations: string[] = [];

  if (metrics.wordErrorRate > 0.10) {
    recommendations.push('High WER - consider audio quality improvements');
  }

  if (metrics.diarizationErrorRate > 0.05) {
    recommendations.push('Speaker labeling issues - ensure clear audio separation');
  }

  if (metrics.realtimeFactor > 0.5) {
    recommendations.push('Slow processing - consider model optimization');
  }

  if (metrics.noiseLevel === 'high') {
    recommendations.push('High background noise - recommend noise reduction');
  }

  return recommendations;
}
```

### Step 2: Audio Quality Optimization

```typescript
// src/twinmind/preprocessing/audio.ts
import ffmpeg from 'fluent-ffmpeg';

interface AudioPreprocessOptions {
  targetSampleRate?: number;   // Default: 16000 Hz
  channels?: number;           // Default: 1 (mono)
  noiseReduction?: boolean;    // Enable noise reduction
  normalization?: boolean;     // Normalize volume
  format?: string;             // Output format (mp3, wav)
}

const defaultOptions: AudioPreprocessOptions = {
  targetSampleRate: 16000,
  channels: 1,
  noiseReduction: true,
  normalization: true,
  format: 'mp3',
};

export async function preprocessAudio(
  inputPath: string,
  outputPath: string,
  options: AudioPreprocessOptions = {}
): Promise<void> {
  const opts = { ...defaultOptions, ...options };

  return new Promise((resolve, reject) => {
    let command = ffmpeg(inputPath)
      .audioFrequency(opts.targetSampleRate!)
      .audioChannels(opts.channels!);

    // Noise reduction using highpass and lowpass filters
    if (opts.noiseReduction) {
      command = command.audioFilters([
        'highpass=f=200',     // Remove low frequency noise
        'lowpass=f=3000',     // Remove high frequency noise
        'afftdn=nf=-25',      // FFT denoiser
      ]);
    }

    // Volume normalization
    if (opts.normalization) {
      command = command.audioFilters([
        'loudnorm=I=-16:TP=-1.5:LRA=11', // EBU R128 loudness normalization
      ]);
    }

    command
      .toFormat(opts.format!)
      .on('end', () => resolve())
      .on('error', reject)
      .save(outputPath);
  });
}

// Audio quality assessment
export async function assessAudioQuality(filePath: string): Promise<{
  quality: 'excellent' | 'good' | 'fair' | 'poor';
  issues: string[];
  recommendations: string[];
}> {
  const issues: string[] = [];
  const recommendations: string[] = [];

  // Get audio metadata
  const metadata = await getAudioMetadata(filePath);

  // Check sample rate
  if (metadata.sampleRate < 16000) {
    issues.push(`Low sample rate: ${metadata.sampleRate} Hz`);
    recommendations.push('Use 16 kHz or higher sample rate');
  }

  // Check bit depth
  if (metadata.bitDepth && metadata.bitDepth < 16) {
    issues.push(`Low bit depth: ${metadata.bitDepth} bits`);
    recommendations.push('Use 16-bit or higher audio');
  }

  // Check for clipping
  if (metadata.peakLevel && metadata.peakLevel > -1) {
    issues.push('Audio clipping detected');
    recommendations.push('Reduce recording volume to prevent distortion');
  }

  // Check noise floor
  if (metadata.noiseFloor && metadata.noiseFloor > -40) {
    issues.push(`High noise floor: ${metadata.noiseFloor} dB`);
    recommendations.push('Use noise reduction or improve recording environment');
  }

  const quality = issues.length === 0 ? 'excellent' :
                  issues.length <= 1 ? 'good' :
                  issues.length <= 2 ? 'fair' : 'poor';

  return { quality, issues, recommendations };
}
```

### Step 3: Model Selection and Configuration

```typescript
// src/twinmind/models/config.ts
export interface ModelConfig {
  model: 'ear-3' | 'ear-2' | 'ear-3-custom';
  language?: string;           // 'auto' or ISO code (en, es, fr, etc.)
  diarization: boolean;        // Speaker labeling
  punctuation: boolean;        // Auto-punctuation
  profanityFilter: boolean;    // Filter explicit content
  vocabulary?: string[];       // Custom vocabulary boost
  speakerHints?: string[];     // Expected speaker names
}

// Optimized configs for different scenarios
export const modelConfigs: Record<string, ModelConfig> = {
  // Standard meeting transcription
  meeting: {
    model: 'ear-3',
    language: 'auto',
    diarization: true,
    punctuation: true,
    profanityFilter: false,
  },

  // Technical presentation with jargon
  technical: {
    model: 'ear-3',
    language: 'en',
    diarization: true,
    punctuation: true,
    profanityFilter: false,
    vocabulary: [
      'API', 'SDK', 'microservice', 'Kubernetes', 'Docker',
      'CI/CD', 'serverless', 'GraphQL', 'PostgreSQL'
    ],
  },

  // Call center / customer support
  callCenter: {
    model: 'ear-3',
    language: 'auto',
    diarization: true,  // Important for customer vs agent
    punctuation: true,
    profanityFilter: true,
  },

  // Medical / healthcare
  medical: {
    model: 'ear-3-custom',  // Enterprise custom model
    language: 'en',
    diarization: true,
    punctuation: true,
    profanityFilter: false,
    vocabulary: [
      // Medical terminology
      'diagnosis', 'prognosis', 'contraindication',
      'hematology', 'cardiology', 'oncology',
    ],
  },

  // Lecture / educational
  lecture: {
    model: 'ear-3',
    language: 'auto',
    diarization: false,  // Usually single speaker
    punctuation: true,
    profanityFilter: false,
  },

  // Podcast / interview
  podcast: {
    model: 'ear-3',
    language: 'auto',
    diarization: true,
    punctuation: true,
    profanityFilter: false,
  },
};

export function getOptimalConfig(useCase: string): ModelConfig {
  return modelConfigs[useCase] || modelConfigs.meeting;
}
```

### Step 4: Streaming Optimization

```typescript
// src/twinmind/streaming/optimized.ts
export interface StreamingConfig {
  chunkDurationMs: number;     // Audio chunk size
  overlapMs: number;           // Overlap for continuity
  maxBufferMs: number;         // Maximum buffer before processing
  interimResults: boolean;     // Return partial results
  endpointDetection: boolean;  // Auto-detect speech endpoints
}

const defaultStreamingConfig: StreamingConfig = {
  chunkDurationMs: 100,        // 100ms chunks for low latency
  overlapMs: 50,               // 50ms overlap
  maxBufferMs: 5000,           // 5 second max buffer
  interimResults: true,
  endpointDetection: true,
};

export class OptimizedStreamingClient {
  private config: StreamingConfig;
  private buffer: Float32Array[] = [];
  private lastInterimResult = '';

  constructor(config: Partial<StreamingConfig> = {}) {
    this.config = { ...defaultStreamingConfig, ...config };
  }

  async processChunk(audioChunk: Float32Array): Promise<{
    interim?: string;
    final?: string;
    confidence: number;
  }> {
    this.buffer.push(audioChunk);

    // Check if we have enough data
    const totalMs = this.buffer.length * (this.config.chunkDurationMs);

    if (totalMs >= this.config.chunkDurationMs * 3) {
      // Process accumulated buffer
      const result = await this.sendToApi(this.concatenateBuffer());

      if (result.isFinal) {
        this.buffer = [];  // Clear buffer on final result
        return {
          final: result.text,
          confidence: result.confidence,
        };
      } else {
        // Interim result - keep some buffer for context
        this.trimBuffer();
        return {
          interim: result.text,
          confidence: result.confidence,
        };
      }
    }

    return { confidence: 0 };
  }

  private concatenateBuffer(): Float32Array {
    const totalLength = this.buffer.reduce((sum, arr) => sum + arr.length, 0);
    const result = new Float32Array(totalLength);
    let offset = 0;
    for (const arr of this.buffer) {
      result.set(arr, offset);
      offset += arr.length;
    }
    return result;
  }

  private trimBuffer(): void {
    // Keep last 2 chunks for context overlap
    if (this.buffer.length > 2) {
      this.buffer = this.buffer.slice(-2);
    }
  }

  private async sendToApi(audio: Float32Array): Promise<{
    text: string;
    confidence: number;
    isFinal: boolean;
  }> {
    // Implementation to send to TwinMind streaming API
    const client = getTwinMindClient();
    const response = await client.post('/stream/process', {
      audio: Buffer.from(audio.buffer).toString('base64'),
      interim_results: this.config.interimResults,
    });
    return response.data;
  }
}
```

### Step 5: Caching and Deduplication

```typescript
// src/twinmind/optimization/cache.ts
import crypto from 'crypto';

interface CachedTranscript {
  hash: string;
  transcriptId: string;
  createdAt: Date;
  expiresAt: Date;
}

class TranscriptCache {
  private cache = new Map<string, CachedTranscript>();
  private ttlMs = 24 * 60 * 60 * 1000; // 24 hours

  // Generate hash of audio content for deduplication
  async hashAudio(audioUrl: string): Promise<string> {
    const response = await fetch(audioUrl);
    const buffer = await response.arrayBuffer();

    return crypto
      .createHash('sha256')
      .update(Buffer.from(buffer))
      .digest('hex');
  }

  get(hash: string): string | null {
    const cached = this.cache.get(hash);

    if (!cached) return null;

    if (new Date() > cached.expiresAt) {
      this.cache.delete(hash);
      return null;
    }

    return cached.transcriptId;
  }

  set(hash: string, transcriptId: string): void {
    this.cache.set(hash, {
      hash,
      transcriptId,
      createdAt: new Date(),
      expiresAt: new Date(Date.now() + this.ttlMs),
    });
  }

  // Transcribe with deduplication
  async transcribeWithCache(audioUrl: string): Promise<string> {
    const hash = await this.hashAudio(audioUrl);

    // Check cache
    const cachedId = this.get(hash);
    if (cachedId) {
      console.log(`Cache hit for audio: ${hash.substring(0, 8)}...`);
      return cachedId;
    }

    // Process new audio
    console.log(`Cache miss for audio: ${hash.substring(0, 8)}...`);
    const client = getTwinMindClient();
    const result = await client.transcribe(audioUrl);

    // Store in cache
    this.set(hash, result.id);

    return result.id;
  }
}

export const transcriptCache = new TranscriptCache();
```

## Output
- Performance metrics tracking
- Audio preprocessing pipeline
- Model configuration for use cases
- Streaming optimization
- Caching and deduplication

## Performance Benchmarks

| Metric | Target | Ear-3 Actual |
|--------|--------|--------------|
| Word Error Rate | < 10% | ~5.26% |
| Diarization Error Rate | < 5% | ~3.8% |
| Real-time Factor | < 0.5x | ~0.3x |
| First Word Latency | < 500ms | ~300ms |
| Languages Supported | 100+ | 140+ |

## Error Handling

| Issue | Cause | Solution |
|-------|-------|----------|
| High WER | Poor audio quality | Apply preprocessing |
| Slow processing | Large file | Use streaming API |
| Wrong language | Auto-detect failed | Specify language explicitly |
| Missing speakers | Low audio separation | Improve microphone setup |

## Resources
- [TwinMind Ear-3 Model](https://twinmind.com/ear-3)
- [Audio Best Practices](https://twinmind.com/docs/audio-quality)
- [Streaming API](https://twinmind.com/docs/streaming)

## Next Steps
For cost optimization, see `twinmind-cost-tuning`.
