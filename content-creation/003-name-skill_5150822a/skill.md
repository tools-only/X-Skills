---
name: share-social
description: >
  This skill should be used when the user asks to "optimize for Instagram",
  "YouTube Shorts format", "make it 9:16", "square video", "TikTok format",
  "Reels format", "prepare for social media", "encode for Twitter",
  "optimize for Facebook", "LinkedIn video", "crop for portrait",
  or mentions any platform-specific video format or upload requirements.
allowed-tools:
  - Bash(ffprobe:*)
  - Bash(ffmpeg:*)
  - AskUserQuestion
---

# Video Social

Prepare video for social media platforms: correct aspect ratios, resolutions, bitrates, and container settings.

## Platform Presets

| Platform | Aspect | Max Resolution | Max Duration | Video Bitrate | Audio |
|----------|--------|----------------|--------------|---------------|-------|
| Instagram Feed | 4:5 portrait or 1:1 | 1080×1350 / 1080×1080 | 60s | 3.5 Mbps | AAC 128k |
| Instagram Reels | 9:16 | 1080×1920 | 90s | 8 Mbps | AAC 192k |
| TikTok | 9:16 | 1080×1920 | 10min | 8 Mbps | AAC 192k |
| YouTube Shorts | 9:16 | 1080×1920 | 60s | 8 Mbps | AAC 192k |
| YouTube Standard | 16:9 | 1920×1080 | unlimited | 8 Mbps (1080p) | AAC 192k |
| Twitter / X | 16:9 or 1:1 | 1920×1200 | 140s | 25 Mbps cap | AAC 128k |
| Facebook | 16:9 or 9:16 | 1920×1080 | 240min | 4 Mbps | AAC 128k |
| LinkedIn | 16:9 | 1920×1080 | 10min | 5 Mbps | AAC 128k |

## Process

### 1. Identify platform

If the request is ambiguous (e.g., "make it vertical"), ask which platform — bitrate and audio requirements differ.

### 2. Probe source

```bash
ffprobe -v quiet -print_format json -show_streams -show_format "$INPUT"
```

Extract: width, height, duration, existing bitrate, audio codec.

### 3. Determine required transforms

Compare probe output against the platform row in the presets table. Apply only the transforms that are actually needed:

- **Aspect ratio mismatch** → crop or pad (see Key Decisions for the choice rule)
- **Resolution too large** → scale down (never upscale; social platforms reject oversized files at upload)
- **Duration exceeds platform limit** → trim; confirm cut point with user first
- **Bitrate over limit** → re-encode with target bitrate (platforms reject or silently degrade over-bitrate uploads)

Exit condition: when all four properties (aspect ratio, resolution, duration, bitrate) are within platform bounds and `-movflags +faststart` will be set, proceed to Step 4. If source already matches all properties, skip to Step 5 with a simple re-encode plan.

### 4. Construct command

**Crop to fit** (preferred when subject is centered):
```bash
# 16:9 source → 9:16 (1080×1920), center crop:
ffmpeg -i "$INPUT" \
  -vf "crop=ih*9/16:ih,scale=1080:1920" \
  -c:v libx264 -b:v 8000k \
  -c:a aac -b:a 192k \
  -movflags +faststart "$OUTPUT"
```

**Pad to fit** (preserves full frame, adds letterbox/pillarbox):
```bash
# 16:9 source → 9:16 with black bars:
ffmpeg -i "$INPUT" \
  -vf "scale=1080:-2,pad=1080:1920:(ow-iw)/2:(oh-ih)/2:black" \
  -c:v libx264 -b:v 8000k \
  -c:a aac -b:a 192k \
  -movflags +faststart "$OUTPUT"
```

**Square (1:1) from 16:9 — center crop:**
```bash
ffmpeg -i "$INPUT" \
  -vf "crop=ih:ih,scale=1080:1080" \
  -c:v libx264 -b:v 3500k \
  -c:a aac -b:a 128k \
  -movflags +faststart "$OUTPUT"
```

### 5. Confirm

State: crop vs. pad choice, any trim, output resolution and bitrate, output path. Wait for approval.

### 6. Run and verify

After encoding, verify bitrate: `ffprobe -v quiet -show_format "$OUTPUT" | grep bit_rate`

## Key Decisions

- Always include `-movflags +faststart` — relocates the moov atom to the file start, enabling progressive playback before full download. Required for all social platforms.
- **Crop vs. pad**: default to crop for talking-head or centered subjects; default to pad for wide scenic shots or text-heavy content. When the user says "don't cut anything off", use pad. When uncertain with off-center subjects, ask before running.
- Never upscale: if source is smaller than target resolution, scale to fit within bounds or keep original dimensions.
