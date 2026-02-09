---
description: Control thinking token limits via environment variable
allowedTools: ["Bash"]
---

# Thinking Token Control

```bash
ARG="${ARGUMENTS:-on}"
SETTINGS_FILE=".claude/settings.local.json"

case "$ARG" in
  off)
    TOKEN_VALUE="0"
    MESSAGE="Thinking disabled (0 tokens)"
    ;;
  little)
    TOKEN_VALUE="1025"
    MESSAGE="Minimal thinking enabled (1025 tokens)"
    ;;
  on)
    TOKEN_VALUE="8192"
    MESSAGE="Standard thinking enabled (8192 tokens)"
    ;;
  lot)
    TOKEN_VALUE="31999"
    MESSAGE="Extended thinking enabled (31999 tokens)"
    ;;
  *)
    echo "Invalid argument: '$ARG'"
    echo "Usage: /think [off|little|on|lot]"
    echo "  off    - Disable thinking (0 tokens)"
    echo "  little - Minimal thinking (1025 tokens)"
    echo "  on     - Standard thinking (8192 tokens) [default]"
    echo "  lot    - Extended thinking (31999 tokens)"
    exit 1
    ;;
esac

# Update MAX_THINKING_TOKENS in settings file
if [ -f "$SETTINGS_FILE" ]; then
  sed -i '' "s/\"MAX_THINKING_TOKENS\": \"[^\"]*\"/\"MAX_THINKING_TOKENS\": \"$TOKEN_VALUE\"/" "$SETTINGS_FILE"
  echo "$MESSAGE"
  echo "Updated $SETTINGS_FILE with MAX_THINKING_TOKENS: $TOKEN_VALUE"
else
  echo "Error: $SETTINGS_FILE not found"
  exit 1
fi
```