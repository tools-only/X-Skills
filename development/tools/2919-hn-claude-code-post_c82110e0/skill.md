# Show HN: NadirClaw - Open-source router that cuts Claude Code costs 40-70%

I got tired of watching my Claude Code bill climb when half the prompts were things like "read this file" or "what does this variable do." Those don't need a frontier model.

So I built NadirClaw, an open-source proxy that sits between Claude Code and the API. It classifies each prompt using sentence embeddings (takes about 10ms) and routes simple ones to a cheaper model like Gemini Flash. Complex prompts, refactoring, architecture questions, those still go to Claude.

It's an OpenAI-compatible proxy, so setup is just two environment variables:

```
export ANTHROPIC_BASE_URL=http://localhost:8856/v1
export ANTHROPIC_API_KEY=local
nadirclaw serve
claude
```

The classifier uses pre-computed centroid vectors from ~170 seed prompts and all-MiniLM-L6-v2 for embeddings. No training needed, no GPU needed. It also detects agentic patterns (tool use loops, multi-step sessions) and forces those to the complex model automatically.

In my usage, about 60% of Claude Code prompts end up on the cheap model. Streaming works fine. Session persistence keeps you on the same model within a conversation so you don't get jarring switches.

Works with Codex, OpenClaw, and anything that speaks the OpenAI API too. Not just Claude Code.

MIT licensed, Python, runs locally.

GitHub: https://github.com/doramirdor/NadirClaw
