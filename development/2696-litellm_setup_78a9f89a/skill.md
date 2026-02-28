# Universal Model Support (LiteLLM)

Lár is built on top of **LiteLLM**, which means it supports **100+ LLM providers** out of the box.

You do NOT need to install separate SDKs (`openai`, `anthropic`, `google-generativeai`). Lár manages the unified interface for you.

---

## 1. Quick Start: The `.env` File

To switch providers, you simply need to:
1.  Set the correct **Environment Variable** for that provider.
2.  Use the correct **Model String** in your `LLMNode`.

Create a `.env` file in your project root:

```bash
# --- OpenAI ---
OPENAI_API_KEY="sk-..."

# --- Anthropic ---
ANTHROPIC_API_KEY="sk-ant-..."

# --- Google Gemini ---
GEMINI_API_KEY="AIza..."

# --- AWS Bedrock ---
AWS_ACCESS_KEY_ID="AKIA..."
AWS_SECRET_ACCESS_KEY="abc..."
AWS_REGION_NAME="us-east-1"

# --- Azure OpenAI ---
AZURE_API_KEY="..."
AZURE_API_BASE="https://my-resource.openai.azure.com/"
AZURE_API_VERSION="2023-05-15"
```

---

## 2. Supported Model Strings (Cheat Sheet)

Pass these strings to the `model_name` argument of `LLMNode`.

| Provider | Model String | Notes |
| :--- | :--- | :--- |
| **OpenAI** | `gpt-4o` | Standard. |
| **OpenAI** | `gpt-3.5-turbo` | |
| **Anthropic** | `claude-3-opus-20240229` | |
| **Anthropic** | `claude-3-5-sonnet-20240620` | |
| **Google** | `gemini/gemini-1.5-pro` | Prefix with `gemini/`. |
| **Google** | `gemini/gemini-1.5-flash` | |
| **Ollama** | `ollama/phi4` | Prefix with `ollama/`. |
| **Ollama** | `ollama/llama3` | |
| **Bedrock** | `bedrock/anthropic.claude-3-sonnet-20240229-v1:0` | Prefix with `bedrock/`. |
| **Azure** | `azure/gpt-4-turbo` | Needs `AZURE_API_BASE`. |
| **Groq** | `groq/llama3-70b-8192` | Extremely fast. |

---

## 3. How to Use Local Models

Lár treats local models as first-class citizens.

### Option A: Ollama (Easiest)
1.  **Install Ollama**: [ollama.com](https://ollama.com)
2.  **Pull a Model**: `ollama pull phi4`
3.  **Run**: Make sure the app is running (it listens on port 11434).
4.  **Code**:

```python
node = LLMNode(
    model_name="ollama/phi4", 
    prompt_template="Analyze this: {task}",
    output_key="result"
)
```

### Option B: Llama.cpp / LM Studio / LocalAI
Most local servers offer an "OpenAI Compatible Endpoint". You can point Lár to this endpoint using `api_base`.

1.  **Start Server**: e.g., `python -m llama_cpp.server --model model.gguf --port 8080`
2.  **Code**: use `openai/custom` prefix to force the OpenAI protocol.

```python
node = LLMNode(
    model_name="openai/my-local-model",
    prompt_template="...",
    output_key="res",
    generation_config={
        "api_base": "http://localhost:8080/v1",  # URL of your local server
        "api_key": "sk-local-key"                # Required by some servers (can be fake)
    }
)
```

---

## 4. Advanced Configuration

You can tune model parameters using the `generation_config` dictionary. This passes arguments directly to the provider's API.

```python
node = LLMNode(
    model_name="gpt-4o",
    prompt_template="...",
    output_key="res",
    generation_config={
        "temperature": 0.2,       # Lower = Deterministic, Higher = Creative
        "max_tokens": 1000,       # Limit the output length
        "top_p": 0.9,             # Nucleus sampling
        "frequency_penalty": 0.0, # Reduce repetition
        "stop": ["User:", "\n"]   # Stop sequences
    }
)
```

---

## 5. Troubleshooting / FAQ

#### "ProviderNotFoundError: Model not found"
*   **Check Prefix**: Did you forget `gemini/` or `ollama/`?
*   **Check LiteLLM**: Run `pip install -U litellm` to get the latest model lists.

#### "AuthenticationError"
*   **Check .env**: Ensure `load_dotenv()` is called at the TOP of your script.
*   **Check Key Name**: LiteLLM expects specific names (e.g., `ANTHROPIC_API_KEY`, not `CLAUDE_KEY`).

#### "I want to see the raw API request"
Enable verbose logging to debug exactly what is being sent to the LLM.

```python
import os
os.environ["LITELLM_LOG"] = "DEBUG"
```
