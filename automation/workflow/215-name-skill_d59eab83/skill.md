---
name: youtube-research
description: >
  This skill should be used when the user asks to "search YouTube", "find videos about",
  "get a transcript", "download subtitles", "extract audio from YouTube", "scan a channel",
  "research a topic on YouTube", "get video metadata", "what videos exist about",
  "download YouTube audio", "YouTube research", "summarize this video", "what is this
  video about", "pull captions from", "grab the audio from", or provides a YouTube/Vimeo/video
  URL and wants to extract information from it. Also triggers on "batch download transcripts",
  "analyze a channel", or any multi-video research workflow.
argument-hint: "<url-or-query>"
allowed-tools:
  - Read
  - Grep
  - Glob
  - Bash(python3:*)
  - Bash(yt-dlp:*)
  - Bash(jq:*)
  - AskUserQuestion
  - WebSearch
  - Task
---

# YouTube Research

Multi-platform video research toolkit. Operates in two modes: toolkit (individual operations)
and research (autonomous search-to-synthesis pipeline). All operations use a single CLI at
`${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py`.

Run `python3 ${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py <subcommand> --help`
for full flag reference on any subcommand.

## Toolkit Mode

Invoke individual subcommands for targeted operations. Default mode when the user requests
a specific action (transcript, search, metadata, audio, channel scan).

### Search

Find videos matching a query. Returns structured results with metadata.

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py search "<query>" --count 10
```

Add filters to narrow results: `--min-duration 600` (seconds), `--after 20250101` (YYYYMMDD),
`--min-views 50000`. Filters are applied client-side after fetching, so the tool over-fetches
automatically to compensate. Output is JSON by default; add `-f text` for human-readable.

### Transcript

Download and clean subtitles to LLM-ready text.

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py transcript "<url>"
```

Outputs clean prose to stdout by default. Add `--timestamps` for SRT with timing cues.
Add `--save -t <topic>` to persist to `~/youtube-research/<topic>/`. Use `--lang all` to
list available subtitle languages before downloading. Fallback chain: manual subs then
auto-generated. Exit 4 if no subtitles exist in the requested language.

### Metadata

Extract full video information without downloading.

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py metadata "<url>"
```

Returns: title, description, channel, duration, chapters, view/like counts, tags, available
subtitle languages, thumbnail URL. Add `--playlist` for playlist entry listings.

### Audio

Download audio in the requested format.

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py audio "<url>" --audio-format mp3 -t <topic>
```

Saves to `~/youtube-research/<topic>/audio/`. Supported formats: mp3, m4a, opus, wav.
Always saves to disk (audio cannot go to stdout). Prints the file path on success.

### Channel

Scan a channel's content.

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py channel "<url-or-@handle>" --limit 20
```

Supports tabs: `--tab videos` (default), `shorts`, `streams`, `playlists`. Filter with
`--after`/`--before` (YYYYMMDD). Sort with `--sort views` for most-viewed-first.

### Batch Processing

Any subcommand except `search` accepts `--batch <file>` (or `--batch -` for stdin) to
process multiple URLs. One URL per line; lines starting with `#` or `;` are skipped.

```bash
python3 ${CLAUDE_PLUGIN_ROOT}/skills/youtube-research/scripts/yt_research.py transcript --batch urls.txt --save -t <topic>
```

## Research Mode

Activate when the user asks to "research", "investigate", or "find out about" a topic
using YouTube as a source. This is an autonomous multi-step pipeline:

1. Search for the topic with broad terms, request 15-20 results
2. Evaluate results by reading titles, channels, view counts, and durations. Select 3-5
   videos using these criteria in priority order: depth over breadth (a 30-minute deep-dive
   from a domain expert outweighs five 3-minute overviews), channel authority (subject-matter
   experts and recognized educators over generic aggregators), recency (prefer recent content
   unless the topic is stable), and view count as a tiebreaker (higher views signal broader
   validation, not necessarily quality)
3. Extract metadata for selected videos to confirm relevance (check descriptions, chapters)
4. Download transcripts for the final selection
5. Read each transcript and synthesize findings into a report

Use `WebSearch` to cross-reference claims or fill gaps when YouTube sources disagree or
leave questions unanswered. Use `Task` to spawn parallel transcript downloads when
processing more than 3 videos.

### Research Report Format

Present the synthesis as a structured markdown report:

- Title and one-paragraph summary of the research question
- Key findings (3-7 bullet points of the most important takeaways across all sources)
- Points of agreement between sources (what multiple videos confirm)
- Points of disagreement (where sources contradict, with attribution)
- Unique insights (notable points from individual videos not repeated elsewhere)
- Gaps in coverage (what the sources collectively missed)
- Sources table: video title, channel, duration, and URL for each video used

Attribute specific claims to their source video. Include timestamps when the transcript
preserves them. Save all transcripts under the topic directory for future reference.
See `examples/research-report.md` for a sample report structure.

### Research Composability

See `references/cli-reference.md` for pipeline patterns that chain subcommands with
standard Unix tools (search → jq → batch transcript).

## Error Recovery

| Exit Code | Meaning | Recovery Action |
|-----------|---------|-----------------|
| 0 | Success | — |
| 1 | Usage error | Check `--help` for correct syntax |
| 2 | yt-dlp not found | Tell user to install: `pip install yt-dlp` |
| 3 | Network/download error | Check URL validity; try `--cookies <browser>` for private/restricted content |
| 4 | No results | For transcripts: try `--lang all` to list available languages. For search: broaden query or remove filters |

## Platform Notes

Load `references/platforms.md` when processing a non-YouTube URL or when a yt-dlp command
fails with exit code 3 on an unfamiliar platform. YouTube is the primary platform, but any
yt-dlp-supported URL works (Vimeo, Twitter, Twitch, etc.).

After extracting a transcript, read the output and summarize key points for the user unless
they asked for raw output only.
