# Web Dashboard

The default interface for PocketPaw — a browser-based control panel for chatting, managing channels, and monitoring agent activity.

## Overview

The web dashboard is PocketPaw's primary interface. It starts by default (no flags needed) and provides:

- Real-time chat with the AI agent
- Activity panel showing tool use, thinking, and errors
- Channel management (Discord, Slack, WhatsApp, Telegram)
- Settings configuration
- Memory browser
- Audit log viewer

## Quick Start

```bash
uv run pocketpaw
# Opens http://localhost:8888
```

## Architecture

The dashboard is built with:

- **Backend**: FastAPI with Jinja2 templates
- **Frontend**: Vanilla JS/CSS/HTML (no build step)
- **Real-time**: WebSocket for streaming agent responses and system events
- **UI Framework**: Alpine.js for reactive components

### Key Files

| File | Description |
|------|-------------|
| `dashboard.py` | FastAPI app, REST endpoints, WebSocket handler |
| `web_server.py` | Uvicorn server startup |
| `frontend/js/app.js` | Main application state and WebSocket client |
| `frontend/js/features/channels.js` | Channel management UI |
| `frontend/templates/base.html` | Root HTML template |
| `frontend/templates/components/sidebar.html` | Sidebar navigation |
| `frontend/templates/components/modals/` | Modal dialogs |

## Default Mode Behavior

When started without flags, the dashboard:

1. Starts the FastAPI web server on port 8888
2. Auto-starts all configured channel adapters (Discord, Slack, WhatsApp, Telegram)
3. Serves the web frontend
4. Opens WebSocket connections for real-time streaming

## Channel Management

The sidebar "Channels" button opens a modal with tabs for each channel. A badge shows the count of currently running adapters.

### REST API

| Method | Endpoint | Description |
|--------|----------|-------------|
| `GET` | `/api/channels/status` | Status of all channels |
| `POST` | `/api/channels/save` | Save channel config to `~/.pocketclaw/config.json` |
| `POST` | `/api/channels/toggle` | Start or stop a channel |
| `GET` | `/api/whatsapp/qr` | WhatsApp personal mode QR code |

### Dynamic Start/Stop

Channels can be started and stopped at runtime without restarting PocketPaw. The toggle endpoint:

1. Validates the channel is properly configured
2. Creates the appropriate adapter instance
3. Starts/stops the adapter on the message bus
4. Returns the new status

## WebSocket Protocol

The dashboard communicates with the backend over a single WebSocket connection at `/ws?token=<access_token>`.

### Message Types

**Client → Server:**

```json
{"type": "message", "content": "Hello, agent!"}
{"type": "stop"}
```

**Server → Client:**

```json
{"type": "message", "content": "Hello! How can I help?", "is_stream_chunk": true}
{"type": "message", "content": "", "is_stream_end": true}
{"type": "system_event", "event_type": "tool_start", "data": {"tool": "shell", "input": "ls"}}
{"type": "system_event", "event_type": "tool_result", "data": {"tool": "shell", "output": "..."}}
{"type": "system_event", "event_type": "thinking", "data": {"content": "Let me check..."}}
{"type": "system_event", "event_type": "error", "data": {"message": "API rate limited"}}
```

## Authentication

- **Localhost**: Auth bypassed for `127.0.0.1` connections
- **Remote**: UUID access token required (configured in settings)
- **Exempt paths**: `/webhook/whatsapp`, `/api/whatsapp/qr`, static assets
