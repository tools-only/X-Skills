 # Activity Tracking System for Skill Discovery

 How can we use Claude Hooks to efficiently log our activities to allow it to recommend new skills
 to help with the most time consuming tasks.

 !!! prompt
  How can I use claude hooks to track my activity within Claude to help it suggest new skills? 


 ## 1. Hook Configuration

Create .claude/settings.json in your project (or ~/.claude/settings.json for global tracking):

```json
{
  "hooks": {
    "PostToolUse": [
      {
        "hooks": [
          {
            "type": "command",
            "command": ".claude/hooks/track-activity.sh",
            "description": "Log all tool usage for pattern analysis"
          }
        ]
      }
    ],
    "UserPromptSubmit": [
      {
        "hooks": [
          {
            "type": "command",
            "command": ".claude/hooks/track-prompts.sh",
            "description": "Log user prompts for task analysis"
          }
        ]
      }
    ],
      "SessionStart": [
        {
          "hooks": [
            {
              "type": "command",
              "command": ".claude/hooks/session-start.sh",
              "description": "Initialize activity tracking for session"
            }
          ]
        }
      ]
    }
  }
  ```

 ## 2. Activity Tracking Scripts

  .claude/hooks/track-activity.sh (logs tool usage):

```sh
#!/bin/bash
HOOK_INPUT=$(cat)
TOOL_NAME=$(echo "$HOOK_INPUT" | jq -r '.tool_name')
SESSION_ID=$(echo "$HOOK_INPUT" | jq -r '.session_id')
CWD=$(echo "$HOOK_INPUT" | jq -r '.cwd')
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Create log directory
mkdir -p ~/.claude/activity-logs

# Log as JSONL for easy processing
LOG_FILE=~/.claude/activity-logs/tool-usage.jsonl
echo
"{\"timestamp\":\"$TIMESTAMP\",\"tool\":\"$TOOL_NAME\",\"session\":\"$SESSION_ID\",\"project\":\"$CWD\"}" >>
"$LOG_FILE"

exit 0
```

  .claude/hooks/track-prompts.sh (logs user requests):

```sh
#!/bin/bash

HOOK_INPUT=$(cat)
PROMPT=$(echo "$HOOK_INPUT" | jq -r '.prompt')
SESSION_ID=$(echo "$HOOK_INPUT" | jq -r '.session_id')
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Create log directory
mkdir -p ~/.claude/activity-logs

# Log prompt with metadata
LOG_FILE=~/.claude/activity-logs/prompts.jsonl
jq -n \
  --arg ts "$TIMESTAMP" \
  --arg sid "$SESSION_ID" \
  --arg p "$PROMPT" \
  '{timestamp: $ts, session: $sid, prompt: $p}' >> "$LOG_FILE"

# Allow the prompt to continue
exit 0
```

  .claude/hooks/session-start.sh (session initialization):

```sh
#!/bin/bash

HOOK_INPUT=$(cat)
SESSION_ID=$(echo "$HOOK_INPUT" | jq -r '.session_id')
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Log session start
mkdir -p ~/.claude/activity-logs
LOG_FILE=~/.claude/activity-logs/sessions.jsonl
echo "{\"timestamp\":\"$TIMESTAMP\",\"session\":\"$SESSION_ID\",\"event\":\"start\"}" >> "$LOG_FILE"

exit 0
```

### 3. Analysis Script

  .claude/scripts/analyze-activity.py (identifies patterns):

```python
#!/usr/bin/env python3
  """Analyze Claude Code activity logs to suggest new skills."""

  import json
  from collections import Counter, defaultdict
  from datetime import datetime
  from pathlib import Path

  LOG_DIR = Path.home() / ".claude" / "activity-logs"

  def load_jsonl(filepath):
      """Load JSONL file into list of dicts."""
      if not filepath.exists():
          return []
      with open(filepath) as f:
          return [json.loads(line) for line in f if line.strip()]

  def analyze_tool_usage():
      """Find most-used tool combinations (potential skill patterns)."""
      tool_log = load_jsonl(LOG_DIR / "tool-usage.jsonl")

      # Group tools by session
      sessions = defaultdict(list)
      for entry in tool_log:
          sessions[entry['session']].append(entry['tool'])

      # Find common sequences (2-tool and 3-tool patterns)
      sequences_2 = []
      sequences_3 = []

      for tools in sessions.values():
          for i in range(len(tools) - 1):
              sequences_2.append(f"{tools[i]} → {tools[i+1]}")
          for i in range(len(tools) - 2):
              sequences_3.append(f"{tools[i]} → {tools[i+1]} → {tools[i+2]}")

      print("## Most Common Tool Sequences\n")
      print("### 2-Tool Patterns:")
      for seq, count in Counter(sequences_2).most_common(10):
          print(f"  {count:3d}x  {seq}")

      print("\n### 3-Tool Patterns:")
      for seq, count in Counter(sequences_3).most_common(10):
          print(f"  {count:3d}x  {seq}")

      return Counter(sequences_2), Counter(sequences_3)

  def analyze_prompts():
      """Identify common task types from user prompts."""
      prompt_log = load_jsonl(LOG_DIR / "prompts.jsonl")

      # Simple keyword analysis
      keywords = []
      for entry in prompt_log:
          prompt_lower = entry['prompt'].lower()
          keywords.extend(prompt_lower.split())

      # Find action verbs (common commands)
      action_verbs = ['create', 'generate', 'build', 'update', 'analyze',
                      'fix', 'debug', 'refactor', 'test', 'deploy',
                      'write', 'read', 'search', 'find', 'review']

      verb_counts = Counter()
      for verb in action_verbs:
          verb_counts[verb] = sum(1 for p in prompt_log if verb in p['prompt'].lower())

      print("\n## Common Task Types (Action Verbs)\n")
      for verb, count in verb_counts.most_common(15):
          if count > 0:
              print(f"  {count:3d}x  {verb}")

      return verb_counts

  def suggest_skills(tool_patterns, prompt_verbs):
      """Suggest potential new skills based on activity patterns."""
      print("\n## Suggested New Skills\n")

      # Look for repetitive patterns
      suggestions = []

      # Pattern-based suggestions
      for pattern, count in tool_patterns.most_common(5):
          if count >= 3:  # Repeated at least 3 times
              suggestions.append({
                  'name': f"Automated: {pattern.replace(' → ', '-')}",
                  'reason': f"You've done this sequence {count} times",
                  'pattern': pattern
              })

      # Verb-based suggestions
      common_verbs = [v for v, c in prompt_verbs.items() if c >= 5]
      if 'generate' in common_verbs:
          suggestions.append({
              'name': 'Content generator for your domain',
              'reason': f"You request generation {prompt_verbs['generate']} times",
              'pattern': 'generate-*'
          })

      if 'analyze' in common_verbs:
          suggestions.append({
              'name': 'Analysis automation skill',
              'reason': f"You run analysis {prompt_verbs['analyze']} times",
              'pattern': 'analyze-*'
          })

      for i, suggestion in enumerate(suggestions, 1):
          print(f"{i}. **{suggestion['name']}**")
          print(f"   - Reason: {suggestion['reason']}")
          print(f"   - Pattern: `{suggestion['pattern']}`\n")

  def main():
      print("# Claude Code Activity Analysis\n")
      print(f"Log directory: {LOG_DIR}\n")

      # Analyze tool usage patterns
      seq2, seq3 = analyze_tool_usage()

      # Analyze prompt patterns
      verbs = analyze_prompts()

      # Generate skill suggestions
      suggest_skills(seq3 if seq3 else seq2, verbs)

  if __name__ == "__main__":
      main()
```

## 4. How to Use This System

  Step 1: Set up the hooks
  # Create hook directory
  mkdir -p .claude/hooks

  # Create the tracking scripts (I can create these for you)
  chmod +x .claude/hooks/*.sh

  Step 2: Work normally in Claude Code
  - All tool usage and prompts are automatically logged
  - No disruption to your workflow
  - Logs accumulate in ~/.claude/activity-logs/

  Step 3: Analyze your patterns weekly
  python .claude/scripts/analyze-activity.py

  Step 4: Review with Claude
  - Share the analysis output with me
  - I'll identify specific repetitive workflows
  - We'll design custom skills for your most common patterns

  5. Example Output

  After a week of use, you might see:
  ## Most Common Tool Sequences

  ### 3-Tool Patterns:
     12x  Grep → Read → Edit
      8x  Glob → Read → Write
      6x  Read → Bash → Read

  ## Common Task Types
     15x  fix
     12x  update
      8x  analyze

  ## Suggested New Skills

  1. **Code Fix Automation**
     - Reason: You do Grep→Read→Edit 12 times
     - Pattern: `grep-read-edit`

  2. **Analysis Report Generator**
     - Reason: You run analysis 8 times
     - Pattern: `analyze-*`

  Then I can help you create skills like:
  - bug-fix-assistant - automates your grep→read→edit pattern
  - code-analysis-reporter - generates analysis reports in your style

