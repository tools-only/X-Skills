# OpenAI Codex Response API Guide

## Overview

CCProxy provides access to OpenAI's [Response API](https://platform.openai.com/docs/api-reference/responses) through your ChatGPT Plus subscription. This experimental feature allows programmatic access to ChatGPT models without requiring separate API keys or usage-based billing.

## Prerequisites

- **ChatGPT Plus Subscription**: An active ChatGPT Plus subscription is required
- **Codex CLI** (Optional): If you have the official Codex CLI installed, CCProxy can reuse its credentials
- **OAuth2 Authentication**: Uses the same authentication flow as the official Codex CLI

## Architecture

The Codex integration in CCProxy acts as a reverse proxy to the ChatGPT backend:

```
Client -> CCProxy -> chatgpt.com/backend-api/codex -> ChatGPT Response
```

Key components:

- **OAuth2 PKCE Flow**: Secure authentication without client secrets
- **Token Management**: Automatic token refresh and credential reuse
- **Session Management**: Maintains conversation context across requests
- **Instruction Injection**: Automatically adds required Codex instruction prompt

## Authentication

### Credential Storage

Credentials are stored in `$HOME/.codex/auth.json` with the following structure:

```json
{
  "access_token": "...",
  "refresh_token": "...",
  "id_token": "...",
  "expires_at": 1234567890,
  "account_id": "user-..."
}
```

### Authentication Flow

CCProxy follows this authentication priority:

1. **Check Existing Credentials**: Looks for valid credentials in `$HOME/.codex/auth.json`
2. **Reuse Codex CLI Credentials**: If Codex CLI credentials exist and are valid, uses them
3. **Auto-Refresh**: If access token is expired but refresh token is valid, automatically renews
4. **Manual Login Required**: If no valid credentials exist, user must authenticate

### Login Methods

#### Using CCProxy CLI

```bash
# Enable Codex provider first
ccproxy config codex --enable

# Authenticate (opens browser for OAuth2 flow)
ccproxy auth login-openai

# Verify authentication
ccproxy auth status
```

#### Using Official Codex CLI

```bash
# Install Codex CLI if not already installed
npm install -g @openai/codex-cli

# Authenticate
codex auth login

# CCProxy will automatically detect and use these credentials
```

### OAuth2 Technical Details

The authentication uses OAuth2 PKCE (Proof Key for Code Exchange) flow:

- **Authorization Endpoint**: `https://auth.openai.com/authorize`
- **Token Endpoint**: `https://auth.openai.com/token`
- **Client ID**: Uses the same client ID as Codex CLI
- **Scopes**: Standard OpenAI scopes for ChatGPT access
- **PKCE Challenge**: SHA256 code challenge for secure authorization

## API Usage

### OpenAI-Compatible Chat Completions Endpoint

CCProxy provides OpenAI-compatible endpoints for easier integration with existing tools:

```bash
# Standard OpenAI format
curl -X POST http://localhost:8000/codex/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "model": "gpt-5",
    "messages": [
      {"role": "user", "content": "Hello, how are you?"}
    ]
  }'

# Alternative endpoint format
curl -X POST http://localhost:8000/codex/v1/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "model": "gpt-5",
    "messages": [
      {"role": "user", "content": "Hello, how are you?"}
    ]
  }'
```

**Important Limitations for Chat Completions:**

- **Limited Model Support**: Only certain models work (e.g., `gpt-5` confirmed, others may fail)
- **No Parameter Support**: OpenAI parameters like `temperature`, `top_p`, `frequency_penalty`, etc. are not supported
- **No Tool Calling**: Function calling and tool use are NOT supported (use `/codex/responses` for tool calls)
- **No System Prompts**: System messages and instructions are overridden by the required Codex instruction prompt
- **Reasoning in XML**: Models with reasoning show reasoning content in `<reasoning>...</reasoning>` tags

### Response API (Direct Backend Access)

The Response API provides access to additional ChatGPT backend features:

```bash
curl -X POST http://localhost:8000/codex/responses \
  -H "Content-Type: application/json" \
  -d '{
    "model": "gpt-5",
    "messages": [
      {"role": "user", "content": "Hello, how are you?"}
    ],
    "temperature": 0.7,
    "max_tokens": 150
  }'
```

**Response API Features:**

- **Tool Calling**: Function calling and tool use are supported (main feature advantage)
- **Parameter Support**: More OpenAI parameters may work (temperature, max_tokens, etc.)
- **Session Management**: Full session management capabilities
- **No Custom Instructions**: System prompts/instructions are overridden by the required Codex instruction prompt
- **Backend Dependent**: Features depend on what ChatGPT's backend API actually supports - users need to test individual parameters and capabilities

### Response Format

```json
{
  "id": "chatcmpl-...",
  "object": "chat.completion",
  "created": 1234567890,
  "model": "gpt-5",
  "choices": [
    {
      "index": 0,
      "message": {
        "role": "assistant",
        "content": "I'm doing well, thank you! How can I help you today?"
      },
      "finish_reason": "stop"
    }
  ],
  "usage": {
    "prompt_tokens": 10,
    "completion_tokens": 15,
    "total_tokens": 25
  }
}
```

### Streaming Responses

Enable streaming for real-time responses:

```bash
curl -X POST http://localhost:8000/codex/responses \
  -H "Content-Type: application/json" \
  -d '{
    "model": "gpt-5",
    "messages": [
      {"role": "user", "content": "Write a short story"}
    ],
    "stream": true
  }'
```

Streaming returns Server-Sent Events (SSE):

```
data: {"id":"chatcmpl-...","object":"chat.completion.chunk","created":1234567890,"model":"gpt-5","choices":[{"index":0,"delta":{"content":"Once"},"finish_reason":null}]}

data: {"id":"chatcmpl-...","object":"chat.completion.chunk","created":1234567890,"model":"gpt-5","choices":[{"index":0,"delta":{"content":" upon"},"finish_reason":null}]}

data: [DONE]
```

## Session Management

### Auto-Generated Sessions

Each request to `/codex/responses` creates a new session:

```python
import requests

# Each request gets a new session
response1 = requests.post("http://localhost:8000/codex/responses", json={
    "model": "gpt-5",
    "messages": [{"role": "user", "content": "Hello"}]
})

# This is a completely new conversation
response2 = requests.post("http://localhost:8000/codex/responses", json={
    "model": "gpt-5",
    "messages": [{"role": "user", "content": "Do you remember me?"}]
})
```

### Persistent Sessions

Maintain conversation context using session IDs:

```python
# Start a conversation with a specific session
session_id = "my-conversation-123"

# First message
response1 = requests.post(f"http://localhost:8000/codex/{session_id}/responses", json={
    "model": "gpt-5",
    "messages": [{"role": "user", "content": "My name is Alice"}]
})

# Continue the same conversation
response2 = requests.post(f"http://localhost:8000/codex/{session_id}/responses", json={
    "model": "gpt-5",
    "messages": [{"role": "user", "content": "What's my name?"}]
})
# The model will remember "Alice" from the previous message
```

### Session ID via Headers

You can also provide session IDs via headers:

```bash
curl -X POST http://localhost:8000/codex/responses \
  -H "Content-Type: application/json" \
  -H "session_id: my-session-456" \
  -d '{"model": "gpt-5", "messages": [{"role": "user", "content": "Hello"}]}'
```

## Instruction Prompt Injection

### What is Instruction Injection?

CCProxy automatically injects the Codex instruction prompt into every conversation. This is a **required** component for the ChatGPT backend to function properly.

### How It Works

1. **User sends**: Your original messages
2. **CCProxy injects**: Prepends the Codex instruction prompt
3. **Backend receives**: Combined prompt + your messages
4. **Response generated**: Based on the full context

### Impact on Token Usage

The instruction prompt consumes tokens in every request:

- **Additional tokens**: ~100-200 tokens per request (varies)
- **Cannot be disabled**: Required by the ChatGPT backend
- **Counts against limits**: Reduces available tokens for your content
- **Billing impact**: Uses your ChatGPT Plus quota

### Example

Your request:

```json
{
  "messages": [{ "role": "user", "content": "Hello" }]
}
```

What the backend actually receives:

```json
{
  "messages": [
    { "role": "system", "content": "[Codex instruction prompt...]" },
    { "role": "user", "content": "Hello" }
  ]
}
```

## Model Differences

### Available Models

The Response API uses ChatGPT Plus models, which differ from standard OpenAI API models:

| Response API Model | Equivalent To      | Notes                  |
| ------------------ | ------------------ | ---------------------- |
| `gpt-5`            | ChatGPT Plus GPT-4 | Latest GPT-4 version   |
| `gpt-5-turbo`      | ChatGPT Plus Turbo | Faster, more efficient |
| `gpt-3.5-turbo`    | ChatGPT Free tier  | Basic model            |

### Behavioral Differences

- **Response Style**: Matches ChatGPT web interface behavior
- **Context Window**: Limited by ChatGPT Plus subscription
- **Rate Limits**: Based on ChatGPT Plus terms, not API limits
- **Features**: May include ChatGPT-specific capabilities

## Client Integration Examples

### Using with aichat

Configure aichat to use the Codex endpoint:

```yaml
# ~/.config/aichat/config.yaml
clients:
  - type: claude
    api_base: http://127.0.0.1:8000/codex
```

Usage:

```bash
# Note: Only certain models work (gpt-5 confirmed working)
aichat --model openai:gpt-5 "Hello world"

# Reasoning models show XML tags
aichat --model openai:gpt-5 "Solve this step by step: 2+2*3"
# Output will include: <reasoning>...</reasoning> followed by the answer
```

### OpenAI SDK Example (Response API)

Using the official OpenAI Python SDK with the Response API:

```python
import os
from openai import OpenAI

# Configure to use CCProxy's Codex endpoint
client = OpenAI(
    api_key="dummy-key",  # Required by SDK but not used
    base_url="http://localhost:8000/codex"
)

# Use the Response API with gpt-5
response = client.responses.create(
    model="gpt-5",
    input="How do I check if a Python object is an instance of a class?",
)

print(response.output_text)
```

**Note**: This uses the `/codex/responses` endpoint which supports tool calling and more parameters than the chat completions endpoint.

### Python Client Example

```python
import requests
import json

class CodexClient:
    def __init__(self, base_url="http://localhost:8000"):
        self.base_url = base_url
        self.session_id = None

    def create_chat_completion(self, messages, model="gpt-5", stream=False):
        """Create a chat completion using OpenAI-compatible endpoint."""
        endpoint = f"{self.base_url}/codex/chat/completions"

        payload = {
            "model": model,
            "messages": messages,
            "stream": stream
        }

        response = requests.post(endpoint, json=payload)

        if stream:
            return self._handle_stream(response)
        else:
            return response.json()

    def create_completion(self, messages, model="gpt-5", session_id=None, stream=False):
        """Create a completion with optional session management (Response API)."""

        # Determine endpoint based on session preference
        if session_id:
            endpoint = f"{self.base_url}/codex/{session_id}/responses"
        else:
            endpoint = f"{self.base_url}/codex/responses"

        payload = {
            "model": model,
            "messages": messages,
            "stream": stream
        }

        response = requests.post(endpoint, json=payload)

        if stream:
            return self._handle_stream(response)
        else:
            return response.json()

    def _handle_stream(self, response):
        """Process streaming responses."""
        for line in response.iter_lines():
            if line:
                line = line.decode('utf-8')
                if line.startswith('data: '):
                    data = line[6:]  # Remove 'data: ' prefix
                    if data == '[DONE]':
                        break
                    yield json.loads(data)

# Usage examples
client = CodexClient()

# OpenAI-compatible endpoint (limited functionality but easier integration)
result = client.create_chat_completion([
    {"role": "user", "content": "What is Python?"}
], model="gpt-5")
print(result['choices'][0]['message']['content'])

# Response API (more features but auto-generated sessions)
result = client.create_completion([
    {"role": "user", "content": "What is Python?"}
])
print(result['choices'][0]['message']['content'])

# Streaming with reasoning model
for chunk in client.create_chat_completion(
    [{"role": "user", "content": "Explain quantum computing step by step"}],
    model="gpt-5",
    stream=True
):
    if chunk['choices'][0].get('delta', {}).get('content'):
        print(chunk['choices'][0]['delta']['content'], end='')
```

## Troubleshooting

### Authentication Issues

#### "No valid OpenAI credentials found"

```bash
# Check current status
ccproxy auth status

# Check detailed OpenAI credentials
ccproxy auth openai-info

# Re-authenticate if needed
ccproxy auth login-openai
# or
codex auth login
```

The `openai-info` command shows detailed credential status including:

- ChatGPT Plus subscription status (must show "PLUS")
- Token expiration and time remaining
- Storage location (`$HOME/.codex/auth.json`)
- Refresh token availability

#### "Token refresh failed"

- Your refresh token may have expired
- Re-authenticate using one of the login methods above

#### "ChatGPT Plus subscription required"

- Ensure your OpenAI account has an active ChatGPT Plus subscription
- The Response API is not available for free accounts

### Request Errors

#### "Session not found"

- Session IDs expire after inactivity
- Create a new session or use auto-generated sessions

#### "Model not available"

- Use ChatGPT Plus compatible models (gpt-5, ...)
- Check model availability in your region

#### "Rate limit exceeded"

- ChatGPT Plus has usage limits
- Wait before making additional requests
- Consider implementing exponential backoff

### Connection Issues

#### "Failed to connect to ChatGPT backend"

- Check your internet connection
- Verify ChatGPT service status
- Try again after a few moments

## Best Practices

1. **Session Management**
   - Use persistent sessions for multi-turn conversations
   - Generate new sessions for unrelated queries
   - Store session IDs for conversation continuity

2. **Error Handling**
   - Implement retry logic with exponential backoff
   - Handle both streaming and non-streaming errors
   - Log errors for debugging

3. **Token Optimization**
   - Account for instruction prompt overhead
   - Monitor token usage in responses
   - Implement token counting before requests

4. **Security**
   - Never expose your `$HOME/.codex/auth.json` file

## Limitations

- **ChatGPT Plus Required**: Not available for free OpenAI accounts
- **Instruction Prompt Overhead**: Mandatory prompt injection consumes tokens
- **Rate Limits**: Subject to ChatGPT Plus usage limits
- **Model Availability**: Limited to ChatGPT Plus models
- **Geographic Restrictions**: May not be available in all regions

## References

- [OpenAI Response API Documentation](https://platform.openai.com/docs/api-reference/responses)
- [OAuth2 PKCE Specification](https://datatracker.ietf.org/doc/html/rfc7636)
- [ChatGPT Plus Subscription](https://openai.com/chatgpt/pricing)
