# Release Notes: Lár v1.4.0

**"The Polish Update"**

This release transforms Lár from a powerful engine into a refined product. It introduces a CLI for instant scaffolding, "Low Code" decorators for cleaner code, and crucial developer experience (DX) improvements.

## 🚀 New Features

### 1. The Lár CLI (`lar`)
Stop manually creating files. Scaffold production-ready agents in seconds.
```bash
pip install lar-engine
lar new agent my-bot
```
This generates a complete project structure with `pyproject.toml`, `.env`, and `agent.py`.

### 2. "Low Code" Nodes (`@node`)
Define nodes as simple Python functions.
```python
from lar import node

@node(output_key="summary")
def summarize(state):
    return llm.generate(state["text"])
```

### 3. Dictionary Access for State
Interact with `GraphState` naturally.
```python
# Old
val = state.get("key")
state.set("key", val)

# New
val = state["key"]
state["key"] = val
```

## 🛠 Improvements & Fixes

- **Robust Concurrency**: Fixed a crash in `BatchNode` when using `functools.partial` or dynamic Callables. Logging is now safe for all object types.
- **Identity Resolution**: The package remains `lar-engine` on PyPI, but the CLI ensures your project is set up correctly to `import lar` without confusion.

## 🔄 Backward Compatibility
**Safe to Upgrade.** Lár v1.4.0 is fully compatible with v1.3 code.
- Existing `LLMNode` / class-based setups work unchanged.
- You can mix "Old" and "New" styles in the same graph.
- See `examples/v1_4_showcase.py` for a demonstration.

## 📦 Installation
```bash
pip install lar-engine==1.4.0
```
