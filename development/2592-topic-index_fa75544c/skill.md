# Topic Index

Maps Gemini API topics to documentation files. All paths relative to `~/Documents/google-gemini-context/`.

## Topic-to-File Mapping

| Topic | JS Doc | Python Doc | Common Doc |
|-------|--------|------------|------------|
| Quickstart / setup | javascript/quickstart.md | python/quickstart.md | -- |
| Client initialisation | javascript/client-setup.md | python/client-setup.md | -- |
| Text generation | javascript/text-generation.md | python/text-generation.md | -- |
| Chat / multi-turn | javascript/chat.md | python/chat.md | -- |
| Streaming | javascript/streaming.md | python/streaming.md | -- |
| Function calling / tools | javascript/function-calling.md | python/function-calling.md | -- |
| Structured output / JSON | javascript/structured-output.md | python/structured-output.md | -- |
| Vision / image understanding | javascript/vision.md | python/vision.md | -- |
| Image generation | javascript/image-generation.md | python/image-generation.md | -- |
| Video understanding | javascript/video.md | python/video.md | -- |
| Video generation | javascript/video-generation.md | python/video-generation.md | -- |
| Audio | javascript/audio.md | python/audio.md | -- |
| Speech generation | javascript/speech-generation.md | python/speech-generation.md | -- |
| Music generation | javascript/music-generation.md | python/music-generation.md | -- |
| Document understanding | javascript/document-understanding.md | python/document-understanding.md | -- |
| Embeddings | javascript/embeddings.md | python/embeddings.md | -- |
| File API / uploads | javascript/files-api.md | python/files-api.md | -- |
| Context caching | javascript/caching.md | python/caching.md | -- |
| Code execution | javascript/code-execution.md | python/code-execution.md | -- |
| Thinking / reasoning | javascript/thinking.md | python/thinking.md | -- |
| Token counting | javascript/token-counting.md | python/token-counting.md | -- |
| Grounding / search | javascript/grounding.md | python/grounding.md | -- |
| URL context | javascript/url-context.md | python/url-context.md | -- |
| Batch processing | javascript/batch.md | python/batch.md | -- |
| Safety settings | -- | -- | common/safety.md |
| Authentication | -- | -- | common/authentication.md |
| Pricing | -- | -- | common/pricing.md |
| Rate limits | -- | -- | common/rate-limits.md |
| Error handling | -- | -- | common/errors.md |
| Regions | -- | -- | common/regions.md |
| OpenAI compatibility | -- | -- | common/openai-compatibility.md |
| Models reference | -- | -- | MODELS.md |
| SDK overview (JS+Python) | -- | -- | googlegenai-gemini-api.md |

## Keyword Mapping

When the user's query doesn't map directly to a topic name:

| User Says | Maps To |
|-----------|---------|
| "JSON output", "schema", "typed response" | Structured output |
| "tool use", "tool calling", "tools" | Function calling |
| "upload file", "send image", "attach" | File API |
| "real-time", "live", "websocket" | Streaming |
| "reasoning", "chain of thought", "deep thinking" | Thinking |
| "search grounding", "Google Search" | Grounding |
| "OCR", "PDF", "read document" | Document understanding |
| "RAG", "semantic search", "vector" | Embeddings |
| "save tokens", "reduce cost" | Context caching |
| "run code", "execute Python" | Code execution |
| "how much does it cost" | Pricing |
| "which model", "model comparison" | Models reference |
| "API key", "credentials" | Authentication |
| "429", "quota", "too many requests" | Rate limits |
| "OpenAI compatible", "drop-in replacement" | OpenAI compatibility |
| "TTS", "text to speech", "voice" | Speech generation |
| "generate music", "audio creation" | Music generation |

## Live Documentation Sources

When local docs seem stale or a topic isn't covered locally:

| Source | URL | Notes |
|--------|-----|-------|
| GitHub JS SDK (codegen_instructions) | `https://raw.githubusercontent.com/googleapis/js-genai/refs/heads/main/codegen_instructions.md` | Always current SDK patterns |
| GitHub JS SDK README | `https://raw.githubusercontent.com/googleapis/js-genai/refs/heads/main/README.md` | Package overview |
| Google AI docs (try markdown) | `https://ai.google.dev/gemini-api/docs/{topic}` | HTML only -- use WebFetch |
| NPM package info | `https://www.npmjs.com/package/@google/genai` | Version tracking |
