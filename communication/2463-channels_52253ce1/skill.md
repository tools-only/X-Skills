# Channel Adapters

Connect PocketPaw to Discord, Slack, WhatsApp, Telegram, and the web dashboard.

## Overview

PocketPaw uses a message bus architecture where channel adapters translate between external platforms and the internal `InboundMessage`/`OutboundMessage` format. Multiple channels can run simultaneously.

## Available Channels

| Channel | Library | Dependency | Streaming |
|---------|---------|-----------|-----------|
| **Web Dashboard** | FastAPI WebSocket | Core | Real-time chunks |
| **Telegram** | python-telegram-bot | Core | Edit-in-place (1.5s rate limit) |
| **Discord** | discord.py | `pocketpaw[discord]` | Edit-in-place (1.5s rate limit) |
| **Slack** | slack-bolt | `pocketpaw[slack]` | Edit-in-place |
| **WhatsApp Business** | httpx | Core | No streaming (sends on completion) |
| **WhatsApp Personal** | neonize | `pocketpaw[whatsapp-personal]` | No streaming (sends on completion) |

## Quick Start

### Web Dashboard (Default)

No flags needed — the dashboard is the default mode:

```bash
uv run pocketpaw
# Dashboard at http://localhost:8888
```

All configured adapters auto-start on launch. Manage channels from the sidebar "Channels" button.

### Headless Mode

Run specific channels without the dashboard:

```bash
# Single channel
uv run pocketpaw --discord
uv run pocketpaw --slack
uv run pocketpaw --whatsapp

# Multiple channels
uv run pocketpaw --discord --slack

# Legacy Telegram mode
uv run pocketpaw --telegram
```

## Discord

### Setup

1. Create a Discord bot at [discord.com/developers](https://discord.com/developers/applications)
2. Enable **Message Content Intent** under Bot settings
3. Generate a bot token
4. Invite the bot to your server with appropriate permissions

### Configuration

```json
{
  "discord_bot_token": "your-bot-token",
  "discord_allowed_guild_ids": [123456789],
  "discord_allowed_user_ids": [987654321]
}
```

Or via environment variables:

```bash
export POCKETCLAW_DISCORD_BOT_TOKEN=your-bot-token
export POCKETCLAW_DISCORD_ALLOWED_GUILD_IDS='[123456789]'
```

### Features

- **Slash command**: `/paw <message>` for agent interaction
- **DM support**: Send direct messages to the bot
- **Mentions**: `@BotName <message>` in any allowed channel
- **Streaming**: Responses are edited in-place as chunks arrive (1.5s rate limit)
- **Access control**: Restrict by guild ID and/or user ID

### Auto-Install

If `discord.py` is not installed, PocketPaw will automatically install it when the adapter starts.

## Slack

### Setup

1. Create a Slack app at [api.slack.com/apps](https://api.slack.com/apps)
2. Enable **Socket Mode** (no public URL needed)
3. Add bot scopes: `app_mentions:read`, `chat:write`, `im:history`, `im:read`, `im:write`
4. Generate a Bot Token (`xoxb-...`) and App-Level Token (`xapp-...`)

### Configuration

```json
{
  "slack_bot_token": "xoxb-your-bot-token",
  "slack_app_token": "xapp-your-app-token",
  "slack_allowed_channel_ids": ["C01234567"]
}
```

### Features

- **Socket Mode**: No public URL or webhooks needed
- **App mentions**: `@PocketPaw <message>` in channels
- **DM support**: Direct messages to the bot
- **Thread support**: Conversations maintain thread context via `thread_ts`
- **Access control**: Restrict by channel ID

### Auto-Install

If `slack-bolt` is not installed, PocketPaw will automatically install it when the adapter starts.

## WhatsApp

PocketPaw supports two WhatsApp modes:

### Personal Mode (QR Scan)

Uses the neonize library (WhatsApp Web multi-device protocol). No Meta Developer account needed.

**How it works:**

1. Start PocketPaw with the web dashboard
2. Go to Channels > WhatsApp > select "Personal (QR Scan)"
3. Click "Start WhatsApp"
4. Scan the QR code with your phone's WhatsApp (Linked Devices)
5. You're connected — send messages from any WhatsApp chat

**Configuration:**

```json
{
  "whatsapp_mode": "personal",
  "whatsapp_neonize_db": ""
}
```

The `whatsapp_neonize_db` field is optional — defaults to `~/.pocketclaw/neonize.sqlite3`. This stores your WhatsApp session credentials.

**Features:**

- QR-code-scan pairing (just like WhatsApp Web)
- No Meta Developer account, webhooks, or tunnels required
- Credentials persist across restarts (re-scan not needed)
- Real-time QR display in dashboard (polls every 2s)
- Connected status indicator

**Auto-Install:** If `neonize` is not installed, PocketPaw will automatically install it when the adapter starts.

### Business Mode (Cloud API)

Uses the WhatsApp Business Cloud API (Meta). Requires a Meta Developer account.

**Configuration:**

```json
{
  "whatsapp_mode": "business",
  "whatsapp_access_token": "your-access-token",
  "whatsapp_phone_number_id": "your-phone-number-id",
  "whatsapp_verify_token": "your-verify-token",
  "whatsapp_allowed_phone_numbers": ["+1234567890"]
}
```

**Requirements:**

- Meta Developer account with WhatsApp Business API access
- Registered business phone number
- Public webhook URL (use PocketPaw's Cloudflare tunnel or your own)

**Webhook routes:**

- `GET /webhook/whatsapp` — Meta verification challenge
- `POST /webhook/whatsapp` — Incoming message delivery

## Dashboard Channel Management

The web dashboard provides a GUI for managing all channels:

### REST API

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/channels/status` | Get status of all channels (configured, running, mode) |
| `POST` | `/api/channels/save` | Save channel configuration |
| `POST` | `/api/channels/toggle` | Start or stop a channel adapter |
| `GET` | `/api/whatsapp/qr` | Get current WhatsApp QR code (personal mode) |

### Channel Status Response

```json
{
  "discord": {"configured": true, "running": true},
  "slack": {"configured": false, "running": false},
  "whatsapp": {"configured": true, "running": true, "mode": "personal"},
  "telegram": {"configured": true, "running": false}
}
```

## Auto-Install System

Optional dependencies (discord.py, slack-bolt, neonize) are automatically installed when an adapter starts and the package is missing. The system:

1. Detects `ImportError` when the adapter tries to import its library
2. Runs `uv pip install pocketpaw[extra]` (falls back to `pip install` if uv is not available)
3. Verifies the import succeeds after installation
4. Raises `RuntimeError` with a clear message if installation fails

This means users don't need to manually install optional dependencies — just configure a channel and start it.

## Implementation

| File | Description |
|------|-------------|
| `bus/adapters/__init__.py` | `BaseChannelAdapter` base class, `auto_install()` utility |
| `bus/adapters/discord_adapter.py` | Discord adapter |
| `bus/adapters/slack_adapter.py` | Slack adapter |
| `bus/adapters/whatsapp_adapter.py` | WhatsApp Business adapter |
| `bus/adapters/neonize_adapter.py` | WhatsApp Personal adapter |
| `bus/adapters/telegram_adapter.py` | Telegram adapter |
| `bus/adapters/websocket_adapter.py` | WebSocket adapter |
| `bus/events.py` | `Channel` enum, message types |
| `bus/queue.py` | `MessageBus` pub/sub |
| `dashboard.py` | Channel management REST API, QR endpoint |
| `frontend/js/features/channels.js` | Channel management UI logic |
| `frontend/templates/components/modals/channels.html` | Channel config modal |

## Tests

- `tests/test_discord_adapter.py` — 11 tests
- `tests/test_slack_adapter.py` — 11 tests
- `tests/test_whatsapp_adapter.py` — 13 tests
- `tests/test_neonize_adapter.py` — 15 tests
- `tests/test_bus.py` — 9 tests (including cross-channel pub/sub)
