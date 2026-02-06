# iTerm2 + Claude Code Session Awareness Optimization

## Problem

When running multiple Claude Code sessions simultaneously in iTerm2, you lose track of which sessions are waiting for your input. Sessions may sit idle for minutes or hours while you focus on another tab, unaware that Claude finished its task or needs permission approval.

## Solution Overview

A layered approach combining Claude Code's built-in notification system, iTerm2's native features, and custom hooks. The layers work independently so any single failure does not leave you blind.

| Layer | What It Does | Effort |
|-------|-------------|--------|
| 1. Terminal bell | Audible beep when Claude needs attention | 1 minute |
| 2. Desktop notifications | macOS Notification Center alerts with session context | 5 minutes |
| 3. iTerm2 tab color change | Waiting tabs turn visually distinct | 15 minutes |
| 4. iTerm2 badge overlay | Large text watermark showing session state | 15 minutes |
| 5. iTerm2 triggers | Pattern-matched visual cues in terminal output | 10 minutes |
| 6. Tab title propagation | Session name appears in the tab title | 5 minutes |

---

## Layer 1: Terminal Bell (Quickest Win)

Claude Code has a built-in notification channel that emits a terminal bell character (BEL / `\a`) when it needs your input or finishes a task.

### Setup

```bash
claude config set --global preferredNotifChannel terminal_bell
```

### iTerm2 Configuration

1. Open iTerm2 Settings (Cmd+,)
2. Go to **Profiles > Terminal**
3. Under **Notifications**, enable **"Send Growl/Notification Center alerts"**
4. Optionally enable **"Flash visual bell"** for a screen flash instead of / in addition to sound
5. Under **Filter Alerts**, enable **"Send escape sequence-generated alerts"**

### What Happens

When Claude finishes work or needs permission, it emits `\a`. iTerm2 picks this up and:
- Plays the system bell sound (or flashes the screen)
- Shows a bell icon on the tab if the tab is not focused
- Sends a macOS notification if configured

### Verification

```bash
echo -e "\a"
```

You should hear a beep or see a flash.

---

## Layer 2: Desktop Notifications via Hooks

The terminal bell is a blunt instrument -- it tells you *something* happened but not *what* or *where*. Notification hooks give you rich, contextual macOS notifications that identify which project and what kind of attention is needed.

### Setup

Create `~/.claude/hooks/notify-desktop.sh`:

```bash
#!/bin/bash
INPUT=$(cat)

NOTIF_TYPE=$(echo "$INPUT" | jq -r '.notification_type // empty')
MESSAGE=$(echo "$INPUT" | jq -r '.message // "Claude needs your attention"')
CWD=$(echo "$INPUT" | jq -r '.cwd // empty')
PROJECT=$(basename "$CWD")

case "$NOTIF_TYPE" in
  "idle_prompt")
    osascript -e "display notification \"$MESSAGE\" with title \"Claude Code\" subtitle \"$PROJECT\" sound name \"Ping\""
    ;;
  "permission_prompt")
    osascript -e "display notification \"$MESSAGE\" with title \"Claude Code - Permission\" subtitle \"$PROJECT\" sound name \"Submarine\""
    ;;
esac

exit 0
```

Make it executable and add to `~/.claude/settings.json` (applies to all projects):

```bash
chmod +x ~/.claude/hooks/notify-desktop.sh
```

```json
{
  "hooks": {
    "Notification": [
      {
        "matcher": "idle_prompt|permission_prompt",
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/notify-desktop.sh",
            "timeout": 5
          }
        ]
      }
    ]
  }
}
```

### How It Works

- The `Notification` hook event fires whenever Claude Code sends a notification
- The `matcher` field filters by notification type: `idle_prompt` (Claude finished and is waiting) or `permission_prompt` (Claude needs you to approve a tool use)
- The hook reads JSON from stdin, extracts the message and project directory, then runs `osascript` to trigger a native macOS notification
- Different sounds distinguish idle vs. permission notifications
- Using a separate script file avoids fragile nested quoting in settings.json

### Alternative: terminal-notifier

For more control (clickable notifications, custom icons), install `terminal-notifier`:

```bash
brew install terminal-notifier
```

Then create `~/.claude/hooks/notify-terminal-notifier.sh`:

```bash
#!/bin/bash
INPUT=$(cat)

MESSAGE=$(echo "$INPUT" | jq -r '.message // "Needs attention"')
SESSION=$(echo "$INPUT" | jq -r '.session_id // "default"')

terminal-notifier -title 'Claude Code' -message "$MESSAGE" -sound default -group "claude-$SESSION"
exit 0
```

Reference it from settings.json instead of `notify-desktop.sh`.

---

## Layer 3: iTerm2 Tab Color Change

This is the highest-impact visual cue. When Claude is waiting for input, the tab turns a distinct color (e.g., orange or red). When you respond and Claude starts working again, the tab returns to its default color.

### Implementation

Create `~/.claude/hooks/iterm2-tab-color.sh`:

```bash
#!/bin/bash
# Read hook input from stdin
INPUT=$(cat)

EVENT=$(echo "$INPUT" | jq -r '.hook_event_name')
NOTIF_TYPE=$(echo "$INPUT" | jq -r '.notification_type // empty')

set_tab_color() {
  local r=$1 g=$2 b=$3
  # Write escape sequences to /dev/tty, NOT stdout.
  # Hook stdout is captured by Claude Code and never reaches iTerm2.
  printf "\033]6;1;bg;red;brightness;%d\a" "$r" > /dev/tty
  printf "\033]6;1;bg;green;brightness;%d\a" "$g" > /dev/tty
  printf "\033]6;1;bg;blue;brightness;%d\a" "$b" > /dev/tty
}

reset_tab_color() {
  printf "\033]6;1;bg;*;default\a" > /dev/tty
}

case "$EVENT" in
  "Notification")
    case "$NOTIF_TYPE" in
      "permission_prompt")
        # Red for permission needed
        set_tab_color 220 60 60
        ;;
      "idle_prompt")
        # Orange for waiting for input
        set_tab_color 230 160 30
        ;;
    esac
    ;;
  "UserPromptSubmit")
    # User responded -- reset to default
    reset_tab_color
    ;;
esac

exit 0
```

Make it executable:

```bash
chmod +x ~/.claude/hooks/iterm2-tab-color.sh
```

Add to `~/.claude/settings.json`:

```json
{
  "hooks": {
    "Notification": [
      {
        "matcher": "permission_prompt|idle_prompt",
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-tab-color.sh"
          }
        ]
      }
    ],
    "UserPromptSubmit": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-tab-color.sh"
          }
        ]
      }
    ]
  }
}
```

### What Happens

- When Claude emits an `idle_prompt` notification, the tab turns **orange**
- When Claude emits a `permission_prompt` notification, the tab turns **red**
- When you submit a new prompt (`UserPromptSubmit`), the tab resets to its default color

Glancing at your tab bar instantly tells you which sessions need attention.

---

## Layer 4: iTerm2 Badge Overlay

Badges are large semi-transparent text labels in the top-right of a terminal session. They provide at-a-glance context without switching tabs.

### Implementation

Create `~/.claude/hooks/iterm2-badge.sh`:

```bash
#!/bin/bash
INPUT=$(cat)

EVENT=$(echo "$INPUT" | jq -r '.hook_event_name')
NOTIF_TYPE=$(echo "$INPUT" | jq -r '.notification_type // empty')
CWD=$(echo "$INPUT" | jq -r '.cwd // empty')
PROJECT=$(basename "$CWD")

set_badge() {
  local text="$1"
  local encoded=$(echo -n "$text" | base64)
  # Write escape sequences to /dev/tty, NOT stdout.
  # Hook stdout is captured by Claude Code and never reaches iTerm2.
  printf "\033]1337;SetBadgeFormat=%s\a" "$encoded" > /dev/tty
}

case "$EVENT" in
  "Notification")
    case "$NOTIF_TYPE" in
      "permission_prompt")
        set_badge "PERMISSION NEEDED"
        ;;
      "idle_prompt")
        set_badge "WAITING - $PROJECT"
        ;;
    esac
    ;;
  "UserPromptSubmit")
    set_badge "$PROJECT"
    ;;
  "SessionStart")
    set_badge "$PROJECT"
    ;;
esac

exit 0
```

Make executable and add it alongside the tab color hook in the same `Notification`, `UserPromptSubmit`, and `SessionStart` hook arrays.

### What Happens

- Session starts: badge shows the project name
- Claude waits: badge changes to "WAITING - projectname"
- Permission needed: badge changes to "PERMISSION NEEDED"
- You respond: badge resets to the project name

---

## Layer 5: iTerm2 Triggers

Triggers are regex-matched rules that fire when specific text appears in terminal output. They work without any Claude Code configuration -- they are purely iTerm2-side.

### Setup

1. Open iTerm2 Settings (Cmd+,)
2. Go to **Profiles > Advanced > Triggers > Edit**
3. Add triggers based on patterns you observe in your Claude Code output.

**Important:** The regex patterns below are examples to adapt, not exact matches. Claude Code's terminal output varies by version and context. Run a session, observe what text appears when Claude waits for input or needs permission, and craft your regex from that.

| Example Regex (adapt to your output) | Action | Parameter | Instant |
|---------------------------------------|--------|-----------|---------|
| `Waiting for your response` | Post Notification | Claude waiting for input | Yes |
| `Allow\|Deny\|permission` | Highlight Line | Yellow background | Yes |
| `Do you want to proceed` | Post Notification | Claude needs permission | Yes |

### Notes

- **Instant** triggers fire as soon as the text is matched, without waiting for a newline. This is important for prompts that stay on the same line.
- Triggers are a defense-in-depth layer. They catch cases that hooks might miss (e.g., if a hook fails or times out).
- The "Post Notification" action sends a macOS Notification Center alert.
- To discover the right patterns: run `claude --debug` or enable verbose mode (`Ctrl+O`) and note the exact text that appears when Claude stops and waits.

---

## Layer 6: Tab Title Propagation

Claude Code may support propagating session names to terminal titles via escape sequences (OSC 2). The feature request (#18326) was marked completed in January 2026, but verify your Claude Code version includes it. When supported, running `/rename` on a session updates the tab title automatically.

### Setup

In iTerm2 Settings:
1. Go to **Profiles > General > Title**
2. Set the title components to include **"Session Name"** or **"Profile & Session Name"**

Now when you run `/rename my-feature` in Claude Code, the iTerm2 tab title updates to show "Claude: my-feature" (or similar).

### Complementary: SessionStart Hook for Automatic Titles

If you want tabs to automatically show the working directory without manual renaming:

```bash
#!/bin/bash
# ~/.claude/hooks/set-title.sh
INPUT=$(cat)
CWD=$(echo "$INPUT" | jq -r '.cwd // empty')
PROJECT=$(basename "$CWD")
# Write escape sequences to /dev/tty, NOT stdout.
# Hook stdout is captured by Claude Code and never reaches iTerm2.
printf "\033]0;Claude: %s\033\\" "$PROJECT" > /dev/tty
exit 0
```

Add as a `SessionStart` hook:

```json
{
  "hooks": {
    "SessionStart": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/set-title.sh"
          }
        ]
      }
    ]
  }
}
```

---

## Combined Configuration

Here is a complete `~/.claude/settings.json` that implements layers 1-4 and 6:

```json
{
  "preferredNotifChannel": "terminal_bell",
  "hooks": {
    "SessionStart": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-badge.sh",
            "timeout": 5
          },
          {
            "type": "command",
            "command": "~/.claude/hooks/set-title.sh",
            "timeout": 5
          }
        ]
      }
    ],
    "Notification": [
      {
        "matcher": "permission_prompt|idle_prompt",
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-tab-color.sh",
            "timeout": 5
          },
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-badge.sh",
            "timeout": 5
          },
          {
            "type": "command",
            "command": "~/.claude/hooks/notify-desktop.sh",
            "timeout": 5
          }
        ]
      }
    ],
    "UserPromptSubmit": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-tab-color.sh",
            "timeout": 5
          },
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-badge.sh",
            "timeout": 5
          }
        ]
      }
    ]
  }
}
```

Layer 5 (triggers) is configured entirely within iTerm2 Settings, not in Claude Code configuration.

---

## Recommendation Priority

For most users, the optimal approach is:

1. **Start with Layer 1** (terminal bell) -- takes 1 minute, provides immediate audio feedback
2. **Add Layer 3** (tab colors) -- the single most impactful visual cue for multi-session workflows
3. **Add Layer 2** (desktop notifications) -- essential if you switch away from iTerm2 entirely
4. **Add Layer 6** (tab titles) -- helps identify sessions by name rather than guessing
5. **Add Layers 4 and 5** (badges, triggers) -- polish for power users running 4+ concurrent sessions

The tab color change (Layer 3) is the most transformative for the specific problem of "forgetting sessions are waiting." A glance at the tab bar tells you everything: orange tabs need input, red tabs need permission, default-colored tabs are working.

### Alternative: Using the `Stop` Hook Instead of `Notification`

The examples above use the `Notification` hook with `idle_prompt` matcher. An alternative is to use the `Stop` hook, which fires every time the main Claude Code agent finishes responding (regardless of notification type). The `Stop` hook is arguably more reliable for the "Claude is waiting" use case because it fires unconditionally when Claude stops, while `Notification` depends on Claude Code emitting a specific notification type.

The community [smart-notify](https://gist.github.com/WynnD/56ee07b1a88af61aa2390f514ee073b6) project uses `Stop` as its primary hook for this reason. To use `Stop` instead, replace the `Notification` matcher group with a `Stop` entry and check `stop_hook_active` in your script to avoid infinite loops:

```json
{
  "hooks": {
    "Stop": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hooks/iterm2-tab-color.sh",
            "timeout": 5
          }
        ]
      }
    ]
  }
}
```

Note: `Stop` does not support matchers -- it always fires. Your script must handle all cases. Also, `Stop` does not distinguish between "waiting for input" and "waiting for permission," so you lose the red/orange color differentiation. A hybrid approach using both `Stop` and `Notification` gives you the best of both worlds.

---

## Limitations and Caveats

- **Escape sequences must target `/dev/tty`**: Hook stdout is captured by Claude Code (parsed for JSON output), so `printf` to stdout will never reach iTerm2. All escape sequence scripts in this proposal write to `/dev/tty` to bypass Claude Code's stdout capture. If you adapt these scripts, always use `> /dev/tty` for terminal escape sequences.
- **iTerm2-specific**: Tab colors, badges, and proprietary escape sequences only work in iTerm2. Users of Terminal.app, Kitty, WezTerm, or VS Code terminal need different approaches.
- **Hook timing**: Notification hooks fire after Claude sends the notification. There is a small delay (typically <1 second) between Claude finishing and the hook executing.
- **Escape sequence interference**: If your shell prompt or other tools also set tab colors, they may conflict with the hook-based approach. Test carefully.
- **Badge formatting**: Badges use base64-encoded interpolated strings. Special characters in project names may require escaping.
- **Trigger regex fragility**: iTerm2 triggers match against terminal output text, which may change between Claude Code versions. Triggers are best as a supplementary layer, not the primary mechanism.
- **Hook timeouts**: The default hook timeout is 600 seconds for commands. Keep notification hooks lightweight (the examples above complete in milliseconds). Set an explicit `"timeout": 5` for notification hooks to avoid edge cases.
- **Session title propagation**: The `/rename` to terminal title feature was recently completed (issue #18326 closed January 2026). Ensure you are on a recent Claude Code version.

---

## Sources

- [Claude Code: Optimize your terminal setup](https://code.claude.com/docs/en/terminal-config)
- [Claude Code: Hooks reference](https://code.claude.com/docs/en/hooks)
- [Claude Code: Status line configuration](https://code.claude.com/docs/en/statusline)
- [iTerm2: Triggers documentation](https://iterm2.com/documentation-triggers.html)
- [iTerm2: Badges documentation](https://iterm2.com/documentation-badges.html)
- [iTerm2: Proprietary escape codes](https://iterm2.com/documentation-escape-codes.html)
- [Get notified when Claude Code needs your input (Martin Hjartmyr)](https://martin.hjartmyr.se/articles/claude-code-terminal-notifications/)
- [Get Notified When Claude Code Finishes With Hooks (alexop.dev)](https://alexop.dev/posts/claude-code-notification-hooks/)
- [smart-notify: per-session pane detection](https://gist.github.com/WynnD/56ee07b1a88af61aa2390f514ee073b6)
- [Claude Code terminal title feature (#18326)](https://github.com/anthropics/claude-code/issues/18326)
- [Claude Code Notifications That Don't Suck (Boris Buliga)](https://www.d12frosted.io/posts/2026-01-05-claude-code-notifications)
