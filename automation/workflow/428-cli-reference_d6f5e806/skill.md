# CLI Reference

Load this when you need full flag details or pipeline composition patterns beyond what
`--help` provides.

## Global Flags

These flags work on any subcommand:

| Flag | Purpose |
|------|---------|
| `-f json/text` | Output format (default: json for search/metadata/channel, text for transcript) |
| `-t <topic>` | Topic subdirectory for saved files (default: general) |
| `-d <path>` | Override base dir (default: ~/youtube-research, env: `YT_RESEARCH_DIR`) |
| `--cookies <browser>` | Browser cookie auth for age-restricted/private content (chrome, firefox, safari, brave, edge) |
| `-q` | Suppress progress messages on stderr |
| `--dry-run` | Show yt-dlp command without executing |
| `-b <file>` / `--batch <file>` | Read URLs from file, one per line. Use `-` for stdin. Not available on `search`. |
| `--no-color` | Disable colored output (also auto-disabled when `NO_COLOR` is set or stdout is not a TTY) |

## Environment Variables

| Variable | Equivalent flag |
|----------|-----------------|
| `YT_RESEARCH_DIR` | `--dir` |
| `YT_RESEARCH_TOPIC` | `--topic` |
| `YT_RESEARCH_COOKIES` | `--cookies` |
| `NO_COLOR` | `--no-color` |

Precedence: flags > env vars > defaults.

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Usage error (invalid args, unknown subcommand) |
| 2 | yt-dlp not found or incompatible version |
| 3 | Network/download error (yt-dlp failed, URL unreachable, private video) |
| 4 | No results (search empty, no subtitles available, channel has no videos) |

On batch operations: exit 0 if at least one URL succeeded. Exit 3 only if every URL failed.

## Pipeline Patterns

Chain subcommands with standard Unix tools for multi-step workflows:

```bash
# Search, filter with jq, extract transcripts
python3 ${CLAUDE_PLUGIN_ROOT}/skills/search-youtube/scripts/yt_research.py search "topic" --count 20 \
  | jq -r '.[].url' \
  | python3 ${CLAUDE_PLUGIN_ROOT}/skills/search-youtube/scripts/yt_research.py transcript --batch - --save -t topic

# Filter search results by duration before processing
python3 ${CLAUDE_PLUGIN_ROOT}/skills/search-youtube/scripts/yt_research.py search "topic" --count 20 \
  | jq '[.[] | select(.duration > 300)]' \
  | jq -r '.[].url' \
  | python3 ${CLAUDE_PLUGIN_ROOT}/skills/search-youtube/scripts/yt_research.py transcript --batch - --save -t topic

# Channel audit: get all videos, then deep metadata for each
python3 ${CLAUDE_PLUGIN_ROOT}/skills/search-youtube/scripts/yt_research.py channel "@someone" --limit 100 \
  | jq -r '.[].url' \
  | python3 ${CLAUDE_PLUGIN_ROOT}/skills/search-youtube/scripts/yt_research.py metadata --batch - > channel_metadata.json
```
