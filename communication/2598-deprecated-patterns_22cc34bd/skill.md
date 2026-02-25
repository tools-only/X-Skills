# Deprecated Patterns

Common mistakes Claude makes when generating Gemini API code. Check this BEFORE writing any Gemini code.

Source: curated from `~/Documents/google-gemini-context/CLAUDE.md`

---

## Package Names

| Language | Wrong (old) | Correct (current) |
|----------|-------------|-------------------|
| JavaScript | `@google/generative-ai` | `@google/genai` |
| Python | `google-generativeai` or `google.generativeai` | `google-genai` |

## Imports

### JavaScript
```javascript
// WRONG
import { GoogleGenerativeAI } from "@google/generative-ai";

// CORRECT
import { GoogleGenAI } from "@google/genai";
```

### Python
```python
# WRONG
import google.generativeai as genai
genai.configure(api_key=os.environ["GEMINI_API_KEY"])

# CORRECT
from google import genai
client = genai.Client()  # auto-reads GEMINI_API_KEY
```

## Initialisation

### JavaScript
```javascript
// WRONG - old pattern
const genAI = new GoogleGenerativeAI(process.env.GEMINI_API_KEY);
const model = genAI.getGenerativeModel({ model: "gemini-pro" });
const result = await model.generateContent("prompt");

// CORRECT - new pattern
const ai = new GoogleGenAI({});  // auto-reads GEMINI_API_KEY
const response = await ai.models.generateContent({
  model: "gemini-2.5-flash",
  contents: "prompt"
});
```

### Python
```python
# WRONG - old pattern
genai.configure(api_key=os.environ["GEMINI_API_KEY"])
model = genai.GenerativeModel("gemini-pro")
response = model.generate_content("prompt")

# CORRECT - new pattern
client = genai.Client()  # auto-reads GEMINI_API_KEY
response = client.models.generate_content(
    model="gemini-2.5-flash",
    contents="prompt"
)
```

## Model IDs

| Wrong (deprecated) | Correct (current) |
|--------------------|-------------------|
| `gemini-pro` | `gemini-2.5-flash` |
| `gemini-pro-vision` | `gemini-2.5-flash` (unified multimodal) |
| `gemini-2.0-flash-exp` | `gemini-2.0-flash` |
| `gemini-1.5-pro-002` | `gemini-2.5-pro` |

Default recommendation: `gemini-2.5-flash` for most tasks, `gemini-2.5-pro` for complex reasoning.

## Chat API

### JavaScript
```javascript
// WRONG - old pattern
const chat = model.startChat();
const result = await chat.sendMessage("Hello");

// CORRECT - new pattern
const chat = await ai.chats.create({ model: "gemini-2.5-flash" });
const result = await chat.send("Hello");
```

### Python
```python
# WRONG - old streaming parameter
response = model.generate_content("prompt", stream=True)

# CORRECT - config-based streaming
response = client.models.generate_content(
    model="gemini-2.5-flash",
    contents="prompt",
    config={"stream": True}
)
```

## Safety Categories

Claude consistently forgets the 5th category and uses the wrong name for the 4th:

```python
# WRONG - 4 categories, wrong name
categories = [
    "HARM_CATEGORY_HARASSMENT",
    "HARM_CATEGORY_HATE_SPEECH",
    "HARM_CATEGORY_SEXUALLY_EXPLICIT",
    "HARM_CATEGORY_DANGEROUS_CONTENT"       # wrong: no _CONTENT
]

# CORRECT - 5 categories
categories = [
    "HARM_CATEGORY_HARASSMENT",
    "HARM_CATEGORY_HATE_SPEECH",
    "HARM_CATEGORY_SEXUALLY_EXPLICIT",
    "HARM_CATEGORY_DANGEROUS",               # no _CONTENT suffix
    "HARM_CATEGORY_CIVIC_INTEGRITY"          # don't forget this one
]
```

## REST API Headers

```bash
# WRONG - Claude often capitalises this
-H "X-Goog-Api-Key: $GEMINI_API_KEY"

# CORRECT - lowercase
-H "x-goog-api-key: $GEMINI_API_KEY"
```

## Parameter Names

| Wrong (old) | Correct (current) |
|-------------|-------------------|
| `generationConfig` | `config` |
| `safetySettings` (as constructor param) | `config.safety_settings` |
| `stream=True` (as method param) | `config={"stream": True}` |

## Rate Limits

- There are NO daily limits -- only per-minute (RPM, TPM)
- Claude sometimes invents "1,500 requests per day" -- this is wrong
- Free tier: per-model RPM/TPM limits only
