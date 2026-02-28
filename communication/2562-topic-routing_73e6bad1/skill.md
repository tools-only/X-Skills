# Topic Routing

Map user problem domains to the most relevant doc slug prefixes. Load when the question domain needs identification or when broad questions require multi-area research.

## Domain Map

| Domain | Slug Prefix | Entry Points |
|--------|-------------|--------------|
| Channels (Telegram, Discord, Slack, WhatsApp, etc.) | `channels/` | `channels/index`, `channels/<name>`, `channels/troubleshooting`, `channels/pairing` |
| Gateway & networking | `gateway/` | `gateway/configuration`, `gateway/security`, `gateway/sandboxing`, `gateway/troubleshooting`, `gateway/remote` |
| Installation & setup | `install/`, `start/` | `start/getting-started`, `install/index`, `start/wizard`, `start/onboarding` |
| Automation (cron, hooks, webhooks) | `automation/` | `automation/cron-jobs`, `automation/webhook`, `automation/hooks`, `automation/troubleshooting` |
| CLI commands | `cli/` | `cli/index`, `cli/<command>` |
| Concepts & architecture | `concepts/` | `concepts/architecture`, `concepts/agent`, `concepts/models`, `concepts/sessions`, `concepts/memory` |
| Tools & skills | `tools/` | `tools/skills`, `tools/subagents`, `tools/index`, `tools/exec`, `tools/browser` |
| Platforms | `platforms/` | `platforms/macos`, `platforms/linux`, `platforms/windows`, `platforms/ios`, `platforms/android` |
| Model providers | `providers/` | `providers/index`, `providers/anthropic`, `providers/openai`, `providers/openrouter` |
| Security | `security/`, `gateway/security` | `gateway/sandboxing`, `gateway/security`, `security/formal-verification` |
| Troubleshooting | `help/` | `help/troubleshooting`, `help/faq`, `help/debugging`, `help/environment` |
| Nodes (audio, camera, voice) | `nodes/` | `nodes/index`, `nodes/audio`, `nodes/talk`, `nodes/voicewake` |
| Web UI | `web/` | `web/dashboard`, `web/webchat`, `web/tui` |
| Agent system | `concepts/agent*`, `cli/agent*` | `concepts/agent`, `concepts/agent-loop`, `concepts/agent-workspace`, `cli/agent`, `cli/agents` |
| Templates & reference | `reference/` | `reference/templates/*`, `reference/token-use`, `reference/rpc` |

## Cross-Cutting Patterns

Some questions span multiple domains. Common combinations:

- **Remote access**: `gateway/remote` + `gateway/tailscale` + `gateway/security`
- **Channel + automation**: `channels/<name>` + `automation/webhook` or `automation/cron-jobs`
- **Multi-agent**: `concepts/multi-agent` + `tools/subagents` + `gateway/sandboxing`
- **Model selection**: `concepts/models` + `providers/<name>` + `concepts/model-failover`
- **Getting started**: `start/getting-started` + `start/wizard` + `install/index`
- **macOS app**: `platforms/macos` + `platforms/mac/*` (10+ sub-pages)

## Discovery

When the topic doesn't map cleanly to a domain:

```bash
clawdocs search "<user's keywords>" --slugs-only
```

Then fetch the top 1-3 results.
