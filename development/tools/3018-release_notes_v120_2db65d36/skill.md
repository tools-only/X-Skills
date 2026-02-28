# Lár v1.2.0 Release Notes

**Observability & Debugging Update**

This release focuses on "Glass Box" visibility, incorporating feedback from power users who need granular cost tracking and cleaner debug logs.

## New Features

### 1. Cost Attribution per Model
The execution summary now includes a `tokens_by_model` breakdown. This allows you to audit costs across different providers (e.g., separating Ollama usage from OpenAI usage).

**Example Output:**
```json
"summary": {
  "total_tokens": 1500,
  "tokens_by_model": {
    "ollama/phi4": 500,
    "gemini-pro": 1000
  }
}
```

### 2. Console Noise Reduction
Large data structures (like 50kb state dumps or long prompts) are now automatically truncated in the **console output** to keep your terminal readable.
*   **Console**: Shows `... [truncated, total len: 50000 chars]`
*   **JSON Log**: Preserves the full, untruncated data for auditing.

### 3. Granular Node Logging
Every `LLMNode` now explicitly logs `prompt_tokens`, `output_tokens`, and `model` in the step metadata, enabling precise per-step cost analysis.

## Changes
*   **utils.py**: Added `truncate_for_log` utility.
*   **executor.py**: Updated `GraphExecutor` to aggregate token usage by model.
*   **node.py**: Updated `ToolNode`, `LLMNode`, and `AddValueNode` to use truncation.

## Upgrading
```bash
pip install lar-engine==1.2.0
```
