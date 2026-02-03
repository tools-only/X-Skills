# Revisit RP context_builder + chat_send for reviews

## Background

PR #73 migrated from RP's `--response-type review` back to manual prompt building because the two-step approach was too brittle.

## Current approach (manual prompt building)
1. `setup-review` → context building
2. `prompt-get` → get builder's handoff
3. Build custom prompt with explicit verdict tag requirement
4. `chat-send` → single message with our prompt

## Alternative tested (2026-01-26)

Tested `context_builder` (no response_type) → `chat_send` with verdict in message:
- Discovery worked: builder selected relevant files
- Verdict came back when included in chat message
- BUT: Previously tried and found brittle

## Why it was brittle (from PR #73)
- RP's `--response-type review` returns its own verdict format (request-changes, approve)
- Required follow-up message to get `<verdict>...</verdict>` tags
- Timing/reliability issues

## Eric's suggestion (pvncher)
- "You can edit the review prompt" in RP settings
- "Prompt your way there" - add XML tags in instructions
- Problem: Can't require all users to change RP settings

## Tasks when revisited
1. Test `context_builder` (no response_type) → `chat_send --mode review` with verdict in message
2. Verify reliability over multiple runs
3. Compare token usage vs manual approach
4. Check if RP has programmatic way to set review prompt defaults
