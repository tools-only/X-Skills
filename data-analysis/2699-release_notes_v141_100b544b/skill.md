# Lár v1.4.1 Release Notes
*(Patch Release - February 15, 2026)*

This patch release addresses a critical bug in audit logging and formalizes support for Reasoning Models (System 2 thinking).

## Critical Fixes

### 1. Robust Audit Logging (Fixes "Missing Logs")
**Problem:** In v1.4.0, if a script exited early (e.g., user `break` in a loop, `sys.exit`, or crash), the `GraphExecutor` would bypass the log-saving step.
**Fix:** Wrapped the execution loop in a `try...finally` block.
**Impact:** Logs are now **always saved** to `lar_logs/run_{uuid}.json` (or your custom `log_dir`), guaranteeing a complete audit trail.

## New Features

### 1. Reasoning Model Support (System 2 Thinking)
Lár now automatically detects and captures "Reasoning Traces" (Chain of Thought), keeping your final answer clean while preserving the thought process for auditing.

**Supported Models:**
- **DeepSeek R1** (API & Distilled/Ollama)
- **OpenAI o1** (Preview/Mini)
- **Liquid Thinking** (`ollama/liquid-thinking`)

**How it Works:**
1.  **Metadata Capture:** If the API returns `reasoning_content` (Standard), it saves to `state["__last_run_metadata"]["reasoning_content"]`.
2.  **Robust Regex Fallback:** If the model outputs raw `<think>...</think>` tags (common in local models), Lár extracts them.
    - **Robustness:** Handles missing closing tags (cut-off thoughts) and missing opening tags (hallucinated starts).
    - **Clean State:** The main `output_key` (e.g., `state['answer']`) contains *only* the final response. The reasoning is safely stored in metadata.

### 2. New Examples
Added a dedicated directory `examples/reasoning_models/` with patterns for:
- `1_deepseek_r1.py`: DeepSeek R1 / Generic Ollama.
- `2_openai_o1.py`: OpenAI o1 Series.
- `3_liquid_thinking.py`: Liquid Thinking (Lateral Logic).

## Internal Changes
- `executor.py`: Added `finally` block for log persistence.
- `node.py`: Enhanced `LLMNode` with robust regex for tag extraction.
- `tracker.py`: valid logic confirmed for per-model token tracking.

---
*Lár v1.4.1 is a recommended update. It ensures your logs are safe and your reasoning models work out-of-the-box.*
