# Just-in-Time Integrations

> **"We don't have 500+ stale integrations. We have infinite fresh ones."**

Traditional frameworks ship with massive libraries of "wrappers" (e.g., `LangChain.tools.HubSpotTool`). These are often:
1.  **Outdated**: APIs change faster than framework maintainers can merge PRs.
2.  **Bloated**: You import a massive library just to update one record.
3.  **Black Box**: You don't know how it handles errors or retries until it crashes.

Lár takes a different approach: **The Integration Builder**.

## The Workflow

Instead of waiting for us to build a tool for `Stripe`, `Linear`, or `Notion`, you generate it yourself in 30 seconds.

### 1. The Integration Prompt
We provide a specialized prompt file: `lar/IDE_MASTER_PROMPT.md`.
> "Read `lar/IDE_MASTER_PROMPT.md` and `lar/examples/` first."
This file teaches your IDE (Cursor, Windsurf, Copilot) the "Gold Standard" for writing robust Lár `ToolNode` wrappers.

### 2. How to Use It

1.  **Open your IDE** (Cursor/Windsurf recommended).
2.  **Drag & Drop** `lar/IDE_MASTER_PROMPT.md` into the chat context.
3.  **Ask**: *"Make me a Lár tool that searches Linear tickets."*

### 3. The Result
Your IDE will read the prompt and generate a **production-ready** file that:
1.  Imports the official SDK (e.g., `linear-sdk`).
2.  Authenticates using `os.getenv`.
3.  Wraps the logic in a deterministic function.
4.  Returns a flat dictionary for `GraphState` merging.
5.  Wraps it all in a `ToolNode`.

**Example Output:**

```python
# Generated in < 30s
from lar import ToolNode
import linear_sdk

def search_linear(state):
    # Authenticate
    # Search
    # Return Dict
    pass

linear_tool = ToolNode(
    tool_function=search_linear,
    input_keys=["query"],
    output_key="tickets",
    next_node=None
)
```

## Why this wins
*   **You Own It**: The code is in your repo, visible and editable.
*   **Always Fresh**: You stick to the latest official SDK, not a wrapper.
*   **Zero Dependencies**: You don't need to install `lar-hubspot` or `lar-stripe`. Just `lar` and `hubspot-client`.

## Proof of Concept
Check out [`examples/patterns/7_integration_test.py`](https://github.com/snath-ai/lar/blob/main/examples/patterns/7_integration_test.py) to see a generated integration for the CoinCap API.
