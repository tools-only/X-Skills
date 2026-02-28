# Reasoning Models (System 2)

**New in v1.4.1**

Lár treats "Reasoning Traces" (also known as Chain of Thought or System 2 Thinking) as a first-class citizen. 

Instead of cluttering your final answer with `<think>` tags or internal monologues, Lár automatically extracts this thought process, saves it to a hidden metadata field for auditing, and delivers only the clean answer to your downstream nodes.

## Supported Models

Lár supports reasoning capture for all major "Thinking" models:

| Provider | Model | Method |
| :--- | :--- | :--- |
| **DeepSeek** | `deepseek-reasoner` (R1) | API Metadata (`reasoning_content`) |
| **OpenAI** | `o1-preview`, `o1-mini` | API Metadata (`reasoning_content`) |
| **Ollama** | `deepseek-r1`, `liquid-thinking` | Regex Fallback (`<think>` tags) |
| **Liquid** | `liquid-thinking` | Regex Fallback (`<think>` tags) |

## How it Works

When `LLMNode` executes, it performs a two-step extraction:

1.  **API Check:** It first checks if the API response contains a specific `reasoning_content` field (Standard for OpenAI/DeepSeek APIs).
2.  **Regex Fallback:** If not found, it scans the text for XML tags like `<think>...</think>`.
    *   It extracts the content inside the tags.
    *   It **removes** the tags and content from the `output_key` (the state).
    *   It saves the extraction to `run_metadata`.

### Robustness (v1.4.1)
The regex engine is robust to malformed outputs common in smaller local models:
*   **Missing Closing Tag:** Captures everything after `<think>` (e.g., if the model hit a token limit).
*   **Missing Opening Tag:** Captures everything before `</think>` (e.g., if the model started thinking "silently").

## Usage

You do not need special configuration. Just use `LLMNode` as usual.

### 1. Using DeepSeek R1 (Ollama)

```python
from lar import LLMNode

node = LLMNode(
    model_name="ollama/deepseek-r1:7b",
    prompt_template="Solve this riddle: {riddle}",
    output_key="answer"
)
```

**The Result:**

*   **State (`state['answer']`)**:
    ```text
    The answer is a sponge.
    ```
*   **Audit Log (`run_metadata`)**:
    ```json
    {
      "model": "ollama/deepseek-r1:7b",
      "total_tokens": 150,
      "reasoning_content": "First I need to analyze the properties... It holds water... It has holes..."
    }
    ```

### 2. Using OpenAI o1

```python
node = LLMNode(
    model_name="openai/o1-mini",
    prompt_template="Write a secure Python function for: {task}",
    output_key="code"
)
```

## Best Practices

1.  **Audit Logs:** Always inspect `lar_logs/` (or your custom log dir) to see the reasoning. It is the only place it is stored.
2.  **Prompting:** For local models (like Liquid Thinking), you may need to explicitly ask for tags in the system prompt if the model is stubborn:
    ```python
    system_instruction="You are a reasoning engine. You MUST enclose your internal monologue in <think> tags."
    ```
