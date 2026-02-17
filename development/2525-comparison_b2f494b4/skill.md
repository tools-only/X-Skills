# NadirClaw vs OpenRouter

Both NadirClaw and [OpenRouter](https://openrouter.ai) give you access to multiple LLM providers through a single API. But they solve different problems and work in fundamentally different ways.

## TL;DR

| | NadirClaw | OpenRouter |
|---|---|---|
| **What it does** | Routes prompts to the right model automatically | Gives you access to many models via one API |
| **Runs where** | Locally on your machine | Cloud service |
| **Routing** | Automatic — classifies every prompt in ~10ms | Manual — you pick the model per request |
| **Cost** | Free (open source) + you pay providers directly | Markup on provider pricing |
| **Privacy** | Prompts never leave your machine (except to your chosen provider) | Prompts pass through OpenRouter's servers |
| **Latency** | Direct to provider (no middleman) | Extra hop through OpenRouter |
| **Subscription support** | Use your existing ChatGPT/Claude/Gemini subscription via OAuth | API keys only |
| **Local models** | Native Ollama support | No local model support |

## The Core Difference

**OpenRouter** is a unified API gateway. You send a request, you choose the model, OpenRouter forwards it. It's a proxy that simplifies multi-provider access.

**NadirClaw** is an intelligent router. You send a request, NadirClaw decides which model should handle it based on the prompt complexity. Simple questions go to cheap/free models, complex tasks go to premium models — automatically.

```
OpenRouter:   You pick the model  ──>  OpenRouter  ──>  Provider
NadirClaw:    You send the prompt ──>  NadirClaw classifies it  ──>  Right model
```

## Why This Matters

Without routing, you're either overpaying or underserving:

<p align="center">
  <img src="images/usage-distribution.png" alt="Typical LLM usage distribution" width="500" />
</p>

Most prompts (55%+) are simple — "what does this error mean?", "rename this variable", "summarize this". Sending these to GPT-4.1 or Claude Opus costs 10-50x more than Gemini Flash or a local Ollama model, with no quality difference.

Without routing, your premium quota burns out fast:

<p align="center">
  <img src="images/quota-comparison.png" alt="Quota usage without routing" width="500" />
</p>

NadirClaw makes sure premium models are reserved for prompts that actually need them.

## Feature Comparison

### Routing Intelligence

| Feature | NadirClaw | OpenRouter |
|---|---|---|
| Auto model selection | Sentence-embedding classifier (~10ms) | No (you pick the model) |
| Agentic task detection | Detects tool use, agent loops, system prompts | No |
| Reasoning detection | Routes chain-of-thought to reasoning models | No |
| Session persistence | Pins model for multi-turn conversations | No |
| Context window filtering | Auto-swaps when conversation exceeds context | No |
| Rate limit fallback | Auto-retries then falls back to other tier | Returns error |

### Cost

| | NadirClaw | OpenRouter |
|---|---|---|
| Software cost | Free (MIT license) | Free tier + paid plans |
| Model pricing | Direct provider pricing (no markup) | Provider price + markup |
| Subscription support | Use ChatGPT/Claude/Gemini subscriptions via OAuth | API keys only |
| Local models | Ollama (free, unlimited) | Not supported |
| Savings mechanism | Routes 55%+ of prompts to free/cheap models | None (you pay for what you pick) |

### Privacy and Control

| | NadirClaw | OpenRouter |
|---|---|---|
| Where it runs | Your machine | OpenRouter's cloud |
| Data routing | Direct to provider | Through OpenRouter servers |
| Open source | Yes (MIT) | No |
| Customizable | Full control over routing logic | Limited to API options |
| Self-hosted | Always | No |

### Provider Support

| | NadirClaw | OpenRouter |
|---|---|---|
| OpenAI | Yes | Yes |
| Anthropic | Yes | Yes |
| Google Gemini | Yes (native SDK, not proxied) | Yes |
| DeepSeek | Yes | Yes |
| Ollama (local) | Yes | No |
| Other providers | Via LiteLLM (100+ providers) | 200+ models |
| OAuth login | OpenAI, Anthropic, Google | No |

## When to Use What

### Use NadirClaw if you:
- Want to **save money automatically** without thinking about which model to use
- Use **coding agents** (Codex, Claude Code, Cursor) that make many LLM calls per session
- Have a **ChatGPT/Claude/Gemini subscription** and want to use it programmatically
- Want to run models **locally with Ollama** alongside cloud models
- Care about **privacy** — your prompts go directly to the provider, not through a third party
- Want **full control** — it's open source, runs locally, and you can customize everything

### Use OpenRouter if you:
- Need access to **many models** and want to pick them manually
- Don't want to run anything locally
- Need a **hosted solution** with no setup
- Want to compare outputs across many models side by side

## Quick Comparison: A Coding Session

A typical 2-hour coding session with an AI agent generates ~150 LLM calls:

| | Without Routing | With NadirClaw |
|---|---|---|
| Simple prompts (55%) | 83 x Claude Sonnet = ~$3.70 | 83 x Gemini Flash = ~$0.04 |
| Complex prompts (30%) | 45 x Claude Sonnet = ~$2.00 | 45 x Claude Sonnet = ~$2.00 |
| Reasoning (15%) | 22 x Claude Sonnet = ~$1.00 | 22 x o3 = ~$0.90 |
| **Total** | **~$6.70** | **~$2.94** |
| **Savings** | — | **~56%** |

With a subscription (OAuth), the savings are even larger since simple prompts go to Gemini's free tier.

## Getting Started

```bash
pip install nadirclaw
nadirclaw setup
nadirclaw serve
```

That's it. Point your AI tool at `http://localhost:8856/v1` and NadirClaw handles the rest.

[Back to README](../README.md)
