---
name: 'llamafile: Fix HuggingFace model URLs (wrong org name + fabricated repos)'
description: 'SourceForge mirror is confirmed legitimate (auto-mirror with disclaimer).
  HuggingFace URLs need two fixes: (1) org name `Mozilla` must be changed to `mozilla-ai`,
  (2) model repo names need verification against actual `mozilla-ai` org repos — current
  names appear fabricated. The `llava-v1.5-7b-llamafile` URL works only via redirect.

  **Citations**:

  - SourceForge mirror: <https://sourceforge.net/projects/llamafile.mirror/files/0.9.3/>
  — disclaimer states "exact mirror of the llamafile project" (accessed 2026-02-21)

  - HuggingFace 404s: `Mozilla/Qwen3-0.6B-gguf`, `Mozilla/Mistral-7B-gguf`, `Mozilla/gemma-3-3b-it-gguf`,
  `Mozilla/Llama-3.1-8B-gguf` (accessed 2026-02-21)

  - HuggingFace redirect: `Mozilla/llava-v1.5-7b-llamafile` redirects to `mozilla-ai/llava-v1.5-7b-llamafile`
  (accessed 2026-02-21)

  - GitHub: `mozilla-ai/llamafile` resolves correctly; `Mozilla-Ocho/llamafile` redirects
  (accessed 2026-02-21)'
metadata:
  topic: llamafile-fix-huggingface-model-urls-wrong-org-name-fabricat
  source: Plugin code review session 2026-02-21
  added: '2026-02-21'
  priority: P2
  type: Feature
  status: open
  issue: '#98'
---