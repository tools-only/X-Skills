# Configuration Guide

This document covers all configuration options for OpenAkita.

## Environment Variables

OpenAkita is configured primarily through environment variables, typically stored in a `.env` file.

### Required Settings

| Variable | Description | Example |
|----------|-------------|---------|
| `ANTHROPIC_API_KEY` | Your Anthropic API key | `sk-ant-...` |

### API Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `ANTHROPIC_BASE_URL` | `https://api.anthropic.com` | API endpoint URL |
| `DEFAULT_MODEL` | `claude-sonnet-4-20250514` | Model to use |
| `MAX_TOKENS` | `8192` | Maximum response tokens |

### Agent Behavior

| Variable | Default | Description |
|----------|---------|-------------|
| `AGENT_NAME` | `OpenAkita` | Display name |
| `MAX_ITERATIONS` | `100` | Max Ralph loop iterations |
| `AUTO_CONFIRM` | `false` | Auto-confirm dangerous operations |

### Storage

| Variable | Default | Description |
|----------|---------|-------------|
| `DATABASE_PATH` | `data/agent.db` | SQLite database location |
| `LOG_LEVEL` | `INFO` | Logging verbosity |

### GitHub Integration

| Variable | Default | Description |
|----------|---------|-------------|
| `GITHUB_TOKEN` | - | GitHub PAT for skill search |

## IM Channel Configuration

### Telegram

```bash
TELEGRAM_ENABLED=true
TELEGRAM_BOT_TOKEN=your-bot-token
```

To get a bot token:
1. Message [@BotFather](https://t.me/botfather) on Telegram
2. Use `/newbot` command
3. Follow the prompts
4. Copy the token

### DingTalk

```bash
DINGTALK_ENABLED=true
DINGTALK_CLIENT_ID=your-client-id
DINGTALK_CLIENT_SECRET=your-client-secret
```

### Feishu (Lark)

```bash
FEISHU_ENABLED=true
FEISHU_APP_ID=your-app-id
FEISHU_APP_SECRET=your-app-secret
```

### WeCom (WeChat Work)

```bash
WEWORK_ENABLED=true
WEWORK_CORP_ID=your-corp-id
WEWORK_AGENT_ID=your-agent-id
WEWORK_SECRET=your-secret
```

### QQ (OneBot)

```bash
QQ_ENABLED=true
QQ_ONEBOT_URL=ws://127.0.0.1:8080
```

## Configuration File

You can also use a YAML configuration file at `config/agent.yaml`:

```yaml
# Agent settings
agent:
  name: OpenAkita
  max_iterations: 100
  auto_confirm: false

# Model settings
model:
  provider: anthropic
  name: claude-sonnet-4-20250514
  max_tokens: 8192

# Tools configuration
tools:
  shell:
    enabled: true
    timeout: 30
    blocked_commands:
      - rm -rf /
      - format
  file:
    enabled: true
    allowed_paths:
      - ./
      - /tmp
  web:
    enabled: true
    timeout: 30

# Logging
logging:
  level: INFO
  format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
  file: logs/agent.log
```

## Document Configuration

### SOUL.md

The soul document defines core values. Generally should not be modified:

```markdown
# Soul Overview

OpenAkita is a self-evolving AI assistant...

## Core Values
1. Safety and human oversight
2. Ethical behavior
3. Following guidelines
4. Being genuinely helpful
```

### AGENT.md

Behavioral specifications and workflows:

```markdown
# Agent Behavior Specification

## Working Mode
- Ralph Wiggum Mode (never give up)
- Task execution flow
- Validation requirements
```

### USER.md

User-specific preferences (auto-updated):

```markdown
# User Profile

## Preferences
- Language: English
- Technical level: Advanced
- Preferred tools: Python, Git
```

### MEMORY.md

Working memory (auto-managed):

```markdown
# Working Memory

## Current Task
- Description: ...
- Progress: 50%
- Next steps: ...

## Lessons Learned
- Issue X was solved by Y
```

## Command Line Options

```bash
# Override config file
openakita --config /path/to/config.yaml

# Override log level
openakita --log-level DEBUG

# Override model
openakita --model claude-opus-4-0-20250514

# Disable confirmation prompts
openakita --auto-confirm

# Run in specific mode
openakita --mode chat|task|test
```

## Advanced Configuration

### Proxy Settings

For users behind firewalls:

```bash
# HTTP proxy
HTTP_PROXY=http://proxy:8080
HTTPS_PROXY=http://proxy:8080

# Or use custom API endpoint
ANTHROPIC_BASE_URL=https://your-proxy-service.com
```

### Rate Limiting

```bash
# Requests per minute
RATE_LIMIT_RPM=60

# Tokens per minute
RATE_LIMIT_TPM=100000
```

### Resource Limits

```bash
# Memory limit (MB)
MEMORY_LIMIT=2048

# CPU cores
CPU_LIMIT=4

# Disk space (MB)
DISK_LIMIT=10240
```

## Validation

To validate your configuration:

```bash
openakita config validate
```

To show current configuration:

```bash
openakita config show
```

## Best Practices

1. **Never commit `.env`** - Add to `.gitignore`
2. **Use separate configs** for dev/staging/production
3. **Rotate API keys** regularly
4. **Enable logging** in production
5. **Set resource limits** to prevent runaway processes
