# ui/terminal/ - MIRA Terminal Client

Interactive CLI for chatting with MIRA. Supports local (localhost:1993) and online (miraos.org) modes.

## Files
- `__init__.py` - Package init, exports `main()`
- `__main__.py` - `python -m ui.terminal` entry point
- `app.py` - `MiraTerminal` class: chat loop, prompt_toolkit, bottom toolbar, server management, `main()` entry point with argparse
- `api.py` - `MiraClient` class: all HTTP calls to MIRA API. Also `check_for_update()` and `strip_emotion_tag()` standalone functions
- `config.py` - `TerminalConfig` (Pydantic BaseModel): JSON config at `~/.mira/terminal.json`. Handles mode-based URL/token resolution
- `renderer.py` - Output functions: `print_user_message`, `print_mira_response`, `print_thinking`, `print_error`, `print_tool_use`, `print_info`. Uses prompt_toolkit `print_formatted_text` (thread-safe with `patch_stdout`)
- `commands.py` - Slash command dispatch via `COMMANDS` dict: /help, /tier, /status, /collapse, /domaindoc, /settings, /clear
- `settings.py` - Interactive settings UI with numbered enter-to-toggle list

## Key Patterns
- **Prompt layout**: Separator line (`─`) + `❯ ` input in `_get_prompt_message()`, status bar (tier, domaindocs, loading, update) in `_get_bottom_toolbar()` below the input. Both re-evaluated on `app.invalidate()`. `erase_when_done=True` prevents prompt echo duplication
- **Background API calls**: `threading.Thread` + `patch_stdout()` — response prints above the active prompt. Animation thread cycles `_loading_frame` and invalidates
- **Enter suppression**: Custom key binding swallows Enter while `self.loading` is True
- **Dual mode**: Config `mode` field ("local"/"online") controls URL and token resolution. Local uses Vault, online uses stored API key
- **Domaindoc editing**: Opens content in `$EDITOR` (config.editor, default nano) via tempfile, saves back through API
- **Server auto-start**: Local mode only. Shows `[MIRA IS STARTING]`, overwrites line when ready

## Entry Points
- `python talkto_mira.py` → `ui.terminal.main()`
- `python -m ui.terminal` → same
- `--headless "msg"` for one-shot, `--show-key` for API key display
