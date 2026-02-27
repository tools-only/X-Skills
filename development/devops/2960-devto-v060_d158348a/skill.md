---
title: "How I Cut My LLM Costs 60% With a Local Router (Open Source)"
published: false
tags: ai, llm, opensource, devops
---

I run AI coding agents all day. Claude Code, Codex, various tools through OpenClaw. My API bill was climbing fast.

The problem: every prompt, no matter how simple, was hitting expensive models. "List the files in this directory" costs the same as "refactor this authentication system" when you're sending everything to Claude Opus.

So I built NadirClaw. It's an open source LLM router that classifies prompts in ~10ms and sends them to the right model. Simple stuff goes to cheap or free models. Complex stuff goes to premium.

## How it works

NadirClaw runs as a local proxy on port 8856. Any tool that speaks the OpenAI API can use it. You just change your base URL.

```bash
pip install nadirclaw
nadirclaw serve
```

Under the hood, a sentence-embedding classifier scores prompt complexity in about 10ms. No LLM call for the classification itself.

## What we just shipped in v0.6.0

**Fallback chains** - This is the one I'm most proud of. When a model fails (rate limit, server error, timeout), NadirClaw cascades through a list of fallback models automatically. Your agent never sees the failure.

```
NADIRCLAW_FALLBACK_CHAIN=gpt-4.1,claude-sonnet-4-5,gemini-2.5-flash
```

**Budget alerts** - Set daily and monthly spend limits. Get warnings at 80%. Every single request logs its cost.

```bash
NADIRCLAW_DAILY_BUDGET=5.00
NADIRCLAW_MONTHLY_BUDGET=50.00
nadirclaw budget  # See spend by model
```

**Prompt caching** - LRU cache with configurable TTL. Identical prompts skip the API entirely. Huge for CI/testing where you repeat prompts.

**Web dashboard** - Visit localhost:8856/dashboard for real-time routing stats. Dark theme, auto-refreshes, zero dependencies.

**Docker** - One command to get NadirClaw + Ollama running locally:

```bash
docker compose up
```

Fully local, zero API costs.

## Real numbers

On a typical coding session with mixed prompts:

- ~40% classified as simple (routed to Gemini Flash or Ollama)
- ~50% classified as complex (routed to Claude/GPT)
- ~10% agentic detection override (tool-calling requests forced to complex)

The 40% simple routing is where the savings come from. If you're paying $3/M tokens for everything and 40% of your traffic could use a $0.15/M model, that's a big difference over a month.

## Who it's for

Anyone running AI tools locally. Works with:

- OpenClaw
- Codex CLI
- Claude Code
- Continue
- Cursor
- Any OpenAI-compatible client

If you're paying for OpenRouter or a similar routing service, NadirClaw does the same thing but runs on your machine with your own keys. No platform dependency.

## Try it

```bash
pip install nadirclaw
nadirclaw setup  # Interactive wizard
nadirclaw serve
```

Or Docker:

```bash
git clone https://github.com/doramirdor/NadirClaw
cd NadirClaw
docker compose up
```

GitHub: [github.com/doramirdor/NadirClaw](https://github.com/doramirdor/NadirClaw)

Stars help visibility. If this is useful, consider dropping one.
