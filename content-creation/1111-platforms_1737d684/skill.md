# Platform-Specific Notes

Load this reference when processing non-YouTube URLs or when a YouTube operation
behaves unexpectedly.

## YouTube

Primary platform. All subcommands fully supported.

Quirks:
- `ytsearchdate:` was removed in Feb 2026 (PR #15959). Date-sorted search is no longer
  available; use `--after`/`--before` client-side filtering on `ytsearch:` results instead.
- `--flat-playlist` on channels returns limited per-video metadata (no description, no
  like_count). Follow up with individual `metadata` calls when full detail is needed.
- Auto-generated captions are available for most English-language videos. Other languages
  vary. Always check `--lang all` before assuming availability.
- Age-restricted content requires `--cookies <browser>` with an authenticated browser session.
- Private and unlisted videos require cookies from an account with access.
- YouTube Shorts are accessible via `--tab shorts` on channel scans. They have standard
  video IDs and subtitles work normally despite the short format.
- Live streams in progress: yt-dlp downloads the currently available segment range,
  not a continuous recording. For ongoing live capture, external tooling is needed.

## Vimeo

Supported for metadata, transcript, and audio extraction.

Quirks:
- Subtitle availability is lower than YouTube. Many Vimeo videos have no captions.
- Some videos are password-protected; yt-dlp supports `--video-password <pw>`.
- Private videos on Vimeo use token-based URLs; the full URL with token must be provided.
- Channel/user page scanning works but is less reliable than YouTube.

## Twitter / X

Supported for video extraction from tweets.

Quirks:
- Requires cookies for most content (`--cookies chrome` or similar). Twitter aggressively
  rate-limits unauthenticated requests.
- No subtitle/transcript support — Twitter videos rarely have captions.
- Audio extraction works normally.
- Search is not supported through yt-dlp; use the Twitter API or `twitter` skill for
  discovery, then pass individual tweet URLs to yt_research.py for extraction.

## Podcast Feeds

yt-dlp supports some podcast platforms and direct audio URLs.

Quirks:
- RSS/podcast feed URLs may not work directly. Extract individual episode URLs first.
- Spotify is not supported by yt-dlp (DRM-protected).
- Apple Podcasts pages can sometimes be extracted; results vary.
- For reliable podcast transcript extraction, download audio via `audio` subcommand
  and use Whisper or another transcription tool separately.

## Other Platforms

yt-dlp supports 1000+ extractors. Most video platforms work for metadata and download.
Common ones: Dailymotion, Twitch (VODs and clips), Reddit (video posts), Instagram
(reels/posts with video), TikTok, Bilibili, PeerTube.

General guidance: try the URL directly. If it fails, check `yt-dlp --list-extractors`
for the platform name and `yt-dlp --extractor-descriptions` for details.
