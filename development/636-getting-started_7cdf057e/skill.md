# Getting Started

This guide will help you get OpenAkita up and running quickly.

## Prerequisites

Before you begin, ensure you have:

- **Python 3.11+** installed
- An **Anthropic API key** ([get one here](https://console.anthropic.com/))
- **Git** for cloning the repository

## Installation

### Option 1: Install from Source (Recommended)

```bash
# Clone the repository
git clone https://github.com/openakita/openakita.git
cd openakita

# Create and activate virtual environment
python -m venv venv

# On Linux/macOS:
source venv/bin/activate

# On Windows:
venv\Scripts\activate

# Install the package
pip install -e .
```

### Option 2: Install from PyPI (Coming Soon)

```bash
pip install openakita
```

## Configuration

### 1. Create Environment File

```bash
cp .env.example .env
```

### 2. Add Your API Key

Edit `.env` and set your Anthropic API key:

```bash
ANTHROPIC_API_KEY=sk-your-api-key-here
```

### 3. Optional Settings

```bash
# Custom API endpoint (useful for proxies)
ANTHROPIC_BASE_URL=https://api.anthropic.com

# Model selection
DEFAULT_MODEL=claude-sonnet-4-20250514

# Agent behavior
MAX_ITERATIONS=100
AUTO_CONFIRM=false
```

## Your First Run

### Start the CLI

```bash
openakita
```

You should see:

```
╭─────────────────────────────────────────╮
│           OpenAkita v0.5.9                │
│   A Self-Evolving AI Agent              │
╰─────────────────────────────────────────╯

You> 
```

### Try a Simple Task

```
You> Hello, what can you do?
```

OpenAkita will introduce itself and explain its capabilities.

### Try a Complex Task

```
You> Create a Python script that calculates prime numbers up to 100
```

Watch as OpenAkita:
1. Analyzes the task
2. Writes the code
3. Tests it
4. Reports the results

## Common Commands

| Command | Description |
|---------|-------------|
| `openakita` | Start interactive mode |
| `openakita run "task"` | Execute a single task |
| `openakita status` | Show agent status |
| `openakita selfcheck` | Run self-diagnostics |
| `openakita --help` | Show all commands |

## Next Steps

- [Architecture Overview](architecture.md) - Understand how OpenAkita works
- [Configuration Guide](configuration.md) - All configuration options
- [Skills System](skills.md) - Create custom skills
- [IM Channels](im-channels.md) - Set up Telegram, etc.

## Troubleshooting

### "API key not found"

Ensure your `.env` file exists and contains `ANTHROPIC_API_KEY`.

### "Connection timeout"

Check your network connection. If in China, consider using a proxy:

```bash
ANTHROPIC_BASE_URL=https://your-proxy-url
```

### "Python version error"

OpenAkita requires Python 3.11+. Check your version:

```bash
python --version
```

### Need More Help?

- Check [GitHub Issues](https://github.com/openakita/openakita/issues)
- Join [GitHub Discussions](https://github.com/openakita/openakita/discussions)
