# How I Cut My Claude Code Bill by 60% with a 200-Line Classifier

I love Claude Code. It's the best AI coding assistant I've used. But after a few weeks of heavy usage, I checked my API costs and nearly fell out of my chair.

$247 in a month.

Most of that was from simple prompts that didn't need Claude Sonnet. Things like "what does this function do?" or "read the file at src/main.py" or "add a test for this function." Basic stuff that any decent LLM can handle.

So I built a router. It sits between Claude Code and the API, classifies each prompt in about 10ms, and sends simple stuff to cheaper models while keeping the complex work on Claude.

After a month of using it, my bill dropped to $98. Same usage pattern. Same quality. 60% less money.

Here's how it works and how you can do the same.

## The Problem: Every Prompt Costs the Same

When you use Claude Code (or most AI coding tools), every single request hits the same expensive model. It doesn't matter if you're asking it to refactor a complex async system or just read a file. You pay full price.

In my case, about 65% of my prompts were simple enough that they didn't need Claude Sonnet. But there's no built-in way to route them differently.

I needed a classifier that could decide in real-time: does this prompt need the expensive model, or can it go to something cheaper?

## The Solution: A Tiny Embedding Classifier

I didn't want to train a big ML model or add a bunch of latency. The classifier needed to be fast (under 20ms), lightweight (no GPU), and accurate enough to not mess up complex prompts.

Here's what I built:

1. **Pre-compute two centroid vectors** (one for simple prompts, one for complex) using a sentence embedding model
2. **For each incoming prompt**, compute its embedding and measure cosine similarity to both centroids
3. **Route based on which centroid is closer**

The entire classifier is about 200 lines of Python. It uses [sentence-transformers](https://github.com/UKPLab/sentence-transformers) with the all-MiniLM-L6-v2 model (80 MB, runs on CPU).

### Code: The Classifier

```python
from sentence_transformers import SentenceTransformer
import numpy as np
from pathlib import Path

class PromptClassifier:
    def __init__(self, threshold=0.06):
        self.encoder = SentenceTransformer('all-MiniLM-L6-v2')
        
        # Load pre-computed centroids (shipped with the package)
        pkg_dir = Path(__file__).parent
        self.simple_centroid = np.load(pkg_dir / "simple_centroid.npy")
        self.complex_centroid = np.load(pkg_dir / "complex_centroid.npy")
        
        self.threshold = threshold
    
    def classify(self, prompt: str) -> dict:
        # Encode the prompt
        embedding = self.encoder.encode([prompt], normalize_embeddings=True)[0]
        
        # Measure cosine similarity to both centroids
        simple_sim = np.dot(embedding, self.simple_centroid)
        complex_sim = np.dot(embedding, self.complex_centroid)
        
        # Normalize to 0-1 range
        score = (complex_sim - simple_sim + 2) / 4
        confidence = abs(complex_sim - simple_sim)
        
        # If confidence is low, default to complex (safer)
        if confidence < self.threshold:
            tier = "complex"
        else:
            tier = "complex" if score > 0.5 else "simple"
        
        return {
            "tier": tier,
            "score": score,
            "confidence": confidence
        }
```

That's the core logic. The magic is in the centroids.

### How I Built the Centroids

I collected about 170 real prompts from my own Claude Code sessions and manually labeled them as simple or complex. Then I computed embeddings for all of them and took the mean of each group:

```python
from sentence_transformers import SentenceTransformer
import numpy as np

# Load prompts (SIMPLE_PROMPTS and COMPLEX_PROMPTS are lists of strings)
encoder = SentenceTransformer('all-MiniLM-L6-v2')

simple_embeddings = encoder.encode(SIMPLE_PROMPTS, normalize_embeddings=True)
complex_embeddings = encoder.encode(COMPLEX_PROMPTS, normalize_embeddings=True)

simple_centroid = np.mean(simple_embeddings, axis=0)
complex_centroid = np.mean(complex_embeddings, axis=0)

# Save them
np.save("simple_centroid.npy", simple_centroid)
np.save("complex_centroid.npy", complex_centroid)
```

Those two .npy files are about 1.5 KB each. I ship them with the package. No training step needed when you install it.

## Wrapping It in a Proxy Server

The classifier alone doesn't help unless you can actually route requests to different models. I built a FastAPI server that exposes an OpenAI-compatible API and routes requests based on the classification:

```python
from fastapi import FastAPI, Request
from fastapi.responses import StreamingResponse
import httpx
from nadirclaw.classifier import PromptClassifier

app = FastAPI()
classifier = PromptClassifier()

# Model routing config
SIMPLE_MODEL = "gemini-2.5-flash"
COMPLEX_MODEL = "claude-sonnet-4-5"

@app.post("/v1/chat/completions")
async def chat_completions(request: Request):
    body = await request.json()
    messages = body.get("messages", [])
    
    # Extract the last user message
    last_user = next((m["content"] for m in reversed(messages) if m["role"] == "user"), "")
    
    # Classify it
    result = classifier.classify(last_user)
    tier = result["tier"]
    
    # Pick the model
    target_model = COMPLEX_MODEL if tier == "complex" else SIMPLE_MODEL
    
    # Route to the appropriate provider
    # (Gemini via Google GenAI SDK, others via LiteLLM)
    response = await dispatch_to_provider(target_model, body)
    
    return response
```

Point Claude Code at `http://localhost:8856/v1` instead of the Anthropic API, and every request flows through the router.

## What Gets Routed Where?

After running this for a month on my real Claude Code usage, here's what the distribution looks like:

**Simple tier (65% of requests):**
- "What does this function do?"
- "Read the file at src/main.py"
- "Add a docstring to this class"
- "Show me the git log for this file"
- "What's the error on line 42?"

**Complex tier (35% of requests):**
- "Refactor this module to use dependency injection"
- "Design a caching layer for this API"
- "Explain why this async operation deadlocks"
- Multi-file changes
- Architecture discussions

**Accuracy:** I spot-checked about 200 routed requests. 94% were routed correctly. The 6% that were wrong were borderline cases that worked fine on the cheaper model anyway.

## Beyond Basic Classification: Smart Overrides

A pure embedding classifier isn't enough. I added a few rules on top:

### 1. Agentic Task Detection

If the request includes tool definitions (like when Claude Code is using shell commands or file operations), it always goes to the complex model. Agents need the premium model to handle multi-step reasoning.

```python
def detect_agentic(request: dict) -> bool:
    # Tool definitions in the request
    if request.get("tools"):
        return True
    
    # Tool-role messages (active execution loop)
    messages = request.get("messages", [])
    if any(m.get("role") == "tool" for m in messages):
        return True
    
    # Agent-like system prompts
    system = next((m["content"] for m in messages if m["role"] == "system"), "")
    if any(marker in system.lower() for marker in ["you are a coding agent", "you can execute"]):
        return True
    
    return False
```

### 2. Reasoning Detection

If the prompt has multiple reasoning markers ("step by step", "prove that", "analyze the tradeoffs"), it goes to the complex model.

```python
REASONING_MARKERS = [
    "step by step", "think through", "chain of thought",
    "prove that", "derive the", "mathematically show",
    "analyze the tradeoffs", "compare and contrast"
]

def detect_reasoning(prompt: str) -> bool:
    lower = prompt.lower()
    return sum(1 for marker in REASONING_MARKERS if marker in lower) >= 2
```

### 3. Session Persistence

Once a conversation is routed to a model, follow-up messages in the same session stick to that model. This prevents jarring mid-conversation switches.

```python
from hashlib import sha256

def get_session_key(messages: list) -> str:
    # Hash the system prompt + first user message
    system = next((m["content"] for m in messages if m["role"] == "system"), "")
    first_user = next((m["content"] for m in messages if m["role"] == "user"), "")
    return sha256((system + first_user).encode()).hexdigest()[:16]

# Cache model choice per session (TTL: 30 minutes)
session_cache = {}

def get_model_for_session(session_key: str, default_tier: str):
    if session_key in session_cache:
        return session_cache[session_key]
    model = pick_model_for_tier(default_tier)
    session_cache[session_key] = model
    return model
```

## Results: 60% Cost Reduction

Here's my actual usage for February 2026:

**Before NadirClaw:**
- Total requests: 1,847
- All to Claude Sonnet 4.5
- Total cost: $247.13

**After NadirClaw:**
- Simple tier (65%): 1,201 requests to Gemini 2.5 Flash
  - Cost: $14.82
- Complex tier (35%): 646 requests to Claude Sonnet 4.5
  - Cost: $83.19
- **Total cost: $98.01**

**Savings: $149.12 (60% reduction)**

No quality loss. Same conversations. Just smarter routing.

## How to Use It Yourself

I open-sourced the whole thing: [NadirClaw on GitHub](https://github.com/doramirdor/NadirClaw)

Install it:

```bash
pip install nadirclaw
```

Run the setup wizard:

```bash
nadirclaw setup
```

Start the router:

```bash
nadirclaw serve
```

Point Claude Code at it:

```bash
export ANTHROPIC_BASE_URL=http://localhost:8856/v1
claude
```

Or use it with any tool that speaks the OpenAI API format (Cursor, Continue, OpenClaw, etc.).

## What I Learned

1. **Most LLM usage doesn't need the premium model.** In my case, 65% of prompts were simple enough for a much cheaper model.

2. **A tiny classifier is enough.** You don't need a big ML model or a GPU. Sentence embeddings + cosine similarity gets you 94% accuracy in under 20ms.

3. **Smart overrides matter.** Pure classification isn't enough. You need rules for agentic tasks, reasoning prompts, and session persistence.

4. **Local control beats platform lock-in.** Running the router locally means your API keys stay on your machine, you control the routing logic, and no one can pull the rug out from under you.

The whole classifier is about 200 lines. The cost savings are real. And it works with any LLM tool that speaks the OpenAI API.

If you're spending serious money on Claude (or any other premium LLM), try routing. It's the easiest 60% cost cut I've ever made.

---

**Follow-up questions? Issues? Want to contribute?**  
GitHub: [doramirdor/NadirClaw](https://github.com/doramirdor/NadirClaw)
