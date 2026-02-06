# Critique: iTerm2 + Claude Code Session Awareness Proposal

## Overall Assessment

The proposal presents a well-structured layered approach with good progressive complexity. However, it contains a critical technical flaw that undermines two of its six layers, an unverified feature claim, and several smaller issues that would cause problems during implementation.

---

## Critical Issue: Layers 3 and 4 Will Not Work As Written

**Severity: Blocking**

The tab color script (Layer 3) and badge script (Layer 4) use `printf` to write iTerm2 escape sequences to stdout:

```bash
printf "\033]6;1;bg;red;brightness;%d\a" "$r"
printf "\033]1337;SetBadgeFormat=%s\a" "$encoded"
```

**This will not work.** Claude Code hooks capture stdout for JSON parsing. Hook stdout is consumed by Claude Code -- it does not reach the terminal emulator. This is a documented architectural limitation:

- [Issue #15082](https://github.com/anthropics/claude-code/issues/15082): Terminal escape sequences are captured as text output rather than being interpreted by the terminal emulator. Claude Code executes commands in a subprocess and captures stdout.
- [Issue #11120](https://github.com/anthropics/claude-code/issues/11120): Hook stdout is captured internally, visible only to Claude in system messages, not displayed in the terminal. This was closed as "not planned."
- The [hooks reference](https://code.claude.com/docs/en/hooks) confirms: "For most events, stdout is only shown in verbose mode (Ctrl+O)."

**The fix:** Write escape sequences to `/dev/tty` instead of stdout. This bypasses Claude Code's stdout capture and sends sequences directly to the terminal emulator:

```bash
printf "\033]6;1;bg;red;brightness;%d\a" "$r" > /dev/tty
printf "\033]1337;SetBadgeFormat=%s\a" "$encoded" > /dev/tty
```

This is a known workaround -- writing to `/dev/tty` sends output directly to the controlling terminal device. However, this workaround has its own caveat: Claude Code hooks reportedly run as detached processes with `TTY_NR: 0` (no controlling terminal). Whether `/dev/tty` resolves correctly from a hook subprocess depends on the OS and how Claude Code spawns the process. **This needs testing before being recommended as a working solution.**

The proposal should either:
1. Add the `/dev/tty` redirect and include a verification step confirming it works
2. Acknowledge this as an unverified approach and provide a tested fallback

---

## Unverified Claim: Feature #18326 (Layer 6)

**Severity: Moderate**

The proposal states: "Claude Code now supports propagating session names to terminal titles via escape sequences (OSC 2)."

This is partially misleading. [Issue #18326](https://github.com/anthropics/claude-code/issues/18326) was closed by the issue author on January 15, 2026, who noted that "terminal sequence propagation appeared to already be implemented" but it "wasn't connected to the `/rename` command." The issue was not closed by the Claude Code team with a confirmed fix -- the author self-closed it with uncertain language.

Furthermore, the underlying technical issue remains: [Issue #15082](https://github.com/anthropics/claude-code/issues/15082) documents that OSC escape sequences from Claude Code's Bash tool are captured rather than passed through to the terminal. If Claude Code itself sends OSC 2 sequences natively (outside of the Bash tool), it may work. But this should be verified with a concrete test, not assumed from a self-closed GitHub issue.

The SessionStart hook for automatic titles (the `set-title.sh` script) has the same stdout capture problem as Layers 3 and 4 -- `printf "\033]0;Claude: %s\033\\" "$PROJECT"` writes to stdout, which Claude Code captures.

---

## Layer 1: Terminal Bell -- Contradictory iTerm2 Instructions

**Severity: Minor**

The proposal's iTerm2 setup instructions say to enable "Send Growl/Notification Center alerts" (step 3), but the [official Claude Code terminal setup docs](https://code.claude.com/docs/en/terminal-config) say to enable "Silence bell" for iTerm2 notification setup. These are contradictory -- "Silence bell" suppresses the audible beep while still allowing escape-sequence-generated alerts.

The proposal's step 3 references "Send Growl/Notification Center alerts" which is outdated naming. Modern iTerm2 uses "Send Notification Center alerts" (Growl support was deprecated). The setting name in the proposal should be updated.

Step 5 references "Send escape sequence-generated alerts" under "Filter Alerts" which aligns with the official docs but may not exist as a separate checkbox in all iTerm2 versions. The exact UI path should be verified.

---

## Layer 2: terminal-notifier Command Has a Bug

**Severity: Minor**

The terminal-notifier alternative command contains:

```json
"command": "terminal-notifier -title 'Claude Code' -message \"$(echo $0 | jq -r '.message // \"Needs attention\"')\" -sound default -group claude-$SESSION_ID"
```

`$0` is the script name, not stdin. Hook input arrives via stdin, not as a positional argument. The correct approach is to read from stdin:

```bash
INPUT=$(cat)
MESSAGE=$(echo "$INPUT" | jq -r '.message // "Needs attention"')
terminal-notifier -title 'Claude Code' -message "$MESSAGE" -sound default
```

Additionally, `$SESSION_ID` is not a Claude Code environment variable. The session ID is available in the JSON input as `.session_id`, not as a shell environment variable. The grouping feature will silently fail (using an empty string), which means notifications from different sessions won't be grouped/replaced correctly.

---

## Layer 2: osascript Commands Lack Project Context

**Severity: Minor**

The osascript notification hooks embed `$CLAUDE_PROJECT_DIR` in the subtitle. However, `$CLAUDE_PROJECT_DIR` is an environment variable that Claude Code sets for hooks pointing to the project root. When configured in `~/.claude/settings.json` (global scope) and used across multiple projects, this should work correctly. But the proposal doesn't mention that this variable is only available inside hook execution context -- it won't work if someone tries to test the command directly in a terminal.

A more robust approach would extract the project name from the JSON input's `cwd` field:

```bash
INPUT=$(cat)
CWD=$(echo "$INPUT" | jq -r '.cwd')
PROJECT=$(basename "$CWD")
osascript -e "display notification \"Claude is waiting\" with title \"Claude Code - $PROJECT\""
```

This is what the badge script in Layer 4 correctly does, but Layer 2 doesn't follow the same pattern.

---

## Layer 3: Script Uses Wrong JSON Field Name

**Severity: Minor (only if /dev/tty fix is applied)**

The tab color script reads `notification_type` from the JSON input:

```bash
NOTIF_TYPE=$(echo "$INPUT" | jq -r '.notification_type // empty')
```

This is correct per the [official hooks reference](https://code.claude.com/docs/en/hooks), which documents the Notification event input as including a `notification_type` field. However, the script also reads `hook_event_name` to distinguish between Notification and UserPromptSubmit events. Since the same script is used for both hook events, it needs to handle the case where `notification_type` is absent (for UserPromptSubmit). The `// empty` fallback handles this correctly.

One concern: the script calls `jq` three times, spawning three subprocesses. For a notification hook that should complete in milliseconds, a single `jq` call extracting all needed fields would be more efficient:

```bash
read -r EVENT NOTIF_TYPE <<< $(echo "$INPUT" | jq -r '[.hook_event_name, (.notification_type // "")] | @tsv')
```

---

## Layer 5: Trigger Regex Patterns Are Speculative

**Severity: Moderate**

The proposed trigger patterns (`Waiting for your response`, `Allow|Deny|permission`, `Do you want to proceed`) are described as "approximate" with a note to "tune them to match the exact output format." This is honest, but the patterns themselves appear to be guesses -- Claude Code's TUI renders prompts with specific formatting that may not match these strings.

Claude Code uses a custom TUI (terminal user interface) that renders permission prompts with styled text, not plain "Allow|Deny|permission" strings. The actual text rendered depends on the Claude Code version and may include ANSI formatting codes that break regex matching.

The proposal correctly positions triggers as "defense-in-depth" but should either provide verified regex patterns or explicitly state that the user must inspect their own terminal output (e.g., with `cat -v` or iTerm2's "Log Session Output") to determine correct patterns.

---

## Combined Configuration: Structural Concern

**Severity: Minor**

The combined `settings.json` places `preferredNotifChannel` at the top level alongside `hooks`. This appears correct based on community usage, but the [official settings documentation](https://code.claude.com/docs/en/terminal-config) shows this value being set via `claude config set --global`, not as a direct JSON key. It likely works both ways, but the proposal should note that `claude config set` is the officially documented method.

---

## Missing Alternative: Built-in Notification Center Alerts

**Severity: Moderate (Missing approach)**

iTerm2 has a built-in feature that requires zero Claude Code configuration:

> Settings > Profiles > Terminal > Notifications > "Send Notification Center alerts"

With "notify when a terminal becomes idle" enabled, iTerm2 will natively send a macOS notification when a terminal session stops producing output (i.e., when Claude Code finishes and is waiting). This requires no hooks, no scripts, and no escape sequences. The idle delay is configurable.

The proposal doesn't mention this zero-configuration baseline. For users who just want basic awareness without setting up hooks, this is the simplest starting point and should be listed before Layer 1.

---

## Missing Alternative: iTerm2 Shell Integration "Alert on Next Mark"

**Severity: Minor (Missing approach)**

If Claude Code works with iTerm2 Shell Integration (which tracks command prompts via marks), the "Alert on Next Mark" feature (Edit > Marks and Annotations > Alert on Next Mark, or Cmd-Opt-A) provides a per-command notification. This is a manual per-command approach but worth mentioning as an available tool.

---

## Hook Timeout Claim

**Severity: Minor (Incorrect)**

The proposal states: "The default hook timeout is 600 seconds for commands." This is correct per the [hooks reference](https://code.claude.com/docs/en/hooks) which states command hooks default to 600 seconds. The proposal then recommends setting `"timeout": 5` for notification hooks, which is reasonable. No issue here.

---

## Summary of Required Changes

### Must Fix (blocking):
1. **Layers 3, 4, and 6 scripts**: Redirect escape sequence output to `/dev/tty` instead of stdout, and verify this works from a hook subprocess
2. **Layer 2 terminal-notifier**: Fix `$0` to read from stdin, fix `$SESSION_ID` reference

### Should Fix:
3. **Layer 6 claim**: Clarify that #18326 was self-closed, not confirmed shipped; verify `/rename` actually propagates to terminal titles
4. **Layer 5 triggers**: Either provide verified patterns or clearly state patterns are examples only
5. **Add zero-config baseline**: Mention iTerm2's built-in idle notification before Layer 1

### Nice to Fix:
6. **Layer 1 iTerm2 instructions**: Update "Growl" naming, align with official docs
7. **Layer 2 osascript**: Use JSON input `cwd` field instead of `$CLAUDE_PROJECT_DIR` for consistency
8. **Combined config**: Note `claude config set` as the documented method

---

## Sources Consulted

- [Claude Code Hooks Reference](https://code.claude.com/docs/en/hooks) - official hook event schemas and input/output formats
- [Claude Code Hooks Guide](https://code.claude.com/docs/en/hooks-guide) - setup walkthrough and examples
- [Claude Code Terminal Setup](https://code.claude.com/docs/en/terminal-config) - official terminal configuration guidance
- [Issue #15082: Terminal title escape sequences not passed through](https://github.com/anthropics/claude-code/issues/15082)
- [Issue #11120: Display startup hook stdout output in terminal](https://github.com/anthropics/claude-code/issues/11120) - closed as "not planned"
- [Issue #18326: Propagate session name to terminal title](https://github.com/anthropics/claude-code/issues/18326) - self-closed by author
- [iTerm2 Triggers Documentation](https://iterm2.com/documentation-triggers.html)
- [iTerm2 Proprietary Escape Codes](https://iterm2.com/documentation-escape-codes.html)
- [iTerm2 Badge Documentation](https://iterm2.com/documentation-badges.html)
- [iTerm2 Terminal Profile Preferences](https://iterm2.com/documentation-preferences-profiles-terminal.html)
- [Martin Hjartmyr: Get notified when Claude Code needs your input](https://martin.hjartmyr.se/articles/claude-code-terminal-notifications/)
- [Boris Buliga: Claude Code Notifications That Don't Suck](https://www.d12frosted.io/posts/2026-01-05-claude-code-notifications)
- [Ryan Patterson: Useful notifications in iTerm2](https://cgamesplay.com/post/2020/12/09/iterm-notifications/)
