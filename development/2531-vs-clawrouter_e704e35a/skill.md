# NadirClaw vs ClawRouter

Both NadirClaw and [ClawRouter](https://github.com/BlockRunAI/ClawRouter) route LLM prompts to different models. But they have very different philosophies on how you should pay for and access those models.

## TL;DR

| | NadirClaw | ClawRouter |
|---|---|---|
| **What it does** | Routes prompts to the right model automatically | Routes prompts across 41+ models using 15-dimension analysis |
| **Runs where** | Locally on your machine | Hosted by BlockRunAI |
| **Payment** | Direct to providers (your own API keys or subscriptions) | USDC on Base blockchain |
| **Setup** | `pip install nadirclaw` | Fund a crypto wallet, deposit USDC |
| **Open source** | Yes (MIT) | Yes (source available) |
| **Stars** | Growing | ~3.4K |
| **Privacy** | Prompts stay between you and your provider | Prompts pass through BlockRunAI infrastructure |
| **Local models** | Ollama support built in | No |
| **Subscription support** | OAuth into ChatGPT/Claude/Gemini | No |

## The Core Difference

**ClawRouter** is a hosted routing service. You fund a crypto wallet with USDC, send prompts to BlockRunAI's servers, and they route to the best model. Payment happens on-chain.

**NadirClaw** runs on your machine. You bring your own API keys (or log in with OAuth to use existing subscriptions). Routing happens locally in ~10ms using a sentence-embedding classifier. No middleman, no crypto, no wallet.

```
ClawRouter:  You fund USDC wallet  ──>  BlockRunAI servers  ──>  Provider
NadirClaw:   You send the prompt   ──>  Local classifier     ──>  Direct to provider
```

## Do You Need a Crypto Wallet?

This is the biggest practical difference.

ClawRouter requires you to hold USDC on the Base network. That means setting up a wallet, bridging funds, and managing on-chain transactions. If you're already in crypto, that's fine. If you're a developer who just wants cheaper LLM calls, it's a lot of friction.

NadirClaw uses whatever you already have. Got OpenAI and Anthropic API keys? Done. Got a ChatGPT Plus subscription? Log in with OAuth and use it programmatically. Want to run models locally for free? Point it at Ollama.

## Feature Comparison

### Routing

| Feature | NadirClaw | ClawRouter |
|---|---|---|
| Routing method | Sentence-embedding classifier (~10ms) | 15-dimension scoring |
| Model count | Any OpenAI-compatible + Ollama | 41+ models |
| Agentic task detection | Yes | Unclear |
| Reasoning detection | Routes chain-of-thought to reasoning models | Part of multi-dimension scoring |
| Session persistence | Pins model across multi-turn conversations | Not documented |
| Context window filtering | Auto-swaps when conversation exceeds limits | Not documented |
| Rate limit fallback | Auto-retries, then falls back | Handled by BlockRunAI |

ClawRouter's 15-dimension routing is more granular on paper. NadirClaw's classifier is simpler but fast, and you can customize the routing logic yourself since it runs locally.

### Cost

| | NadirClaw | ClawRouter |
|---|---|---|
| Software cost | Free (MIT) | Free to use, pay per token |
| How you pay | Direct to providers (API keys) | USDC on Base |
| Markup | None | BlockRunAI takes a cut |
| Subscription support | OAuth into ChatGPT/Claude/Gemini | No |
| Local models | Ollama (free, unlimited) | No |
| Savings mechanism | Routes simple prompts to cheap/free models | Routes to cost-optimal model per task |

### Privacy and Control

| | NadirClaw | ClawRouter |
|---|---|---|
| Runs on | Your machine | BlockRunAI cloud |
| Data path | You to provider (direct) | You to BlockRunAI to provider |
| Source code | MIT license, fully modifiable | Open source, tied to BlockRunAI infra |
| Custom routing rules | Full control | Limited to their dimensions |
| Vendor lock-in | None, works with any OpenAI-compatible tool | Tied to USDC payments and BlockRunAI |

### Ecosystem

| | NadirClaw | ClawRouter |
|---|---|---|
| Works with | Anything OpenAI-compatible (Cursor, Cline, aider, Claude Code, etc.) | OpenClaw ecosystem, OpenAI-compatible clients |
| Provider access | OpenAI, Anthropic, Google, DeepSeek, Ollama, 100+ via LiteLLM | 41+ models via BlockRunAI |
| OAuth login | OpenAI, Anthropic, Google | No |
| Local model support | Ollama | No |

## When to Use What

### Use NadirClaw if you:
- Want routing without crypto wallets or on-chain payments
- Already have API keys or LLM subscriptions you want to use
- Care about privacy and want prompts going directly to providers
- Run local models with Ollama
- Want full control over routing logic (it's your code, on your machine)
- Need something that works with any OpenAI-compatible tool

### Use ClawRouter if you:
- Prefer crypto-native payment (already hold USDC, comfortable with wallets)
- Want a hosted solution with no local setup
- Like the 15-dimension routing approach
- Are already in the OpenClaw/BlockRunAI ecosystem
- Don't mind prompts routing through a third party

## Quick Setup Comparison

**NadirClaw:**
```bash
pip install nadirclaw
nadirclaw setup        # add your API keys or log in with OAuth
nadirclaw serve        # starts on localhost:8856
```

**ClawRouter:**
1. Set up a Base-compatible crypto wallet
2. Bridge USDC to Base network
3. Deposit USDC into ClawRouter
4. Configure your client to point at ClawRouter's endpoint

## The Bottom Line

ClawRouter is a solid project with a large model catalog and sophisticated routing dimensions. If you're comfortable with crypto payments and want a hosted solution, it works.

NadirClaw takes a different approach: run locally, pay directly, no middleman. You keep your API keys, your subscriptions, and your privacy. There's nothing to fund, nothing to bridge, and nothing between you and your providers.

For most developers who just want to cut their LLM costs, `pip install nadirclaw` is the shorter path.

[Back to README](../README.md)
