# Install Skill Tracker (Global)

The install-skill-tracker skill automates the installation of a **global** skill tracking
system for Claude Code. The system uses hooks installed in `~/.claude/` to automatically
log all skill invocations across **all projects**, including duration, token usage, and
user prompts for later analysis.

## Key Capabilities

This skill installs global tracking for:

- **Skill Usage**: Which skills are invoked and how often across all projects
- **Execution Duration**: Time spent in each skill
- **Token Usage**: API costs and efficiency metrics
- **User Prompts**: What triggers skill invocations
- **Session Tracking**: Group activities by session
- **Project Tracking**: Identify which project triggered each skill

## When to Use

Use this skill when:

- Setting up skill usage tracking that works across all your Claude Code projects
- Wanting to analyze which skills are used most frequently in your entire workflow
- Monitoring API costs and token usage across all skill executions
- Identifying time-consuming skills needing optimization
- Discovering patterns that could become new skills
- Tracking productivity gains from skill automation
- Understanding prompt cache effectiveness and optimization opportunities

## What Gets Installed

The skill creates this **global** directory structure in your home directory:

```
~/.claude/
├── settings.json             # Hook configuration (global)
├── hooks/
│   ├── track-prompts.sh      # Logs user prompts with project info
│   ├── track-skill-start.sh  # Logs skill start events
│   └── track-skill-end.sh    # Logs skill completion and tokens
├── scripts/
│   ├── analyze-skills.py     # Analysis script
│   └── show-skill-tokens.sh  # Token usage display
└── activity-logs/
    ├── prompts.jsonl         # All prompts across all projects
    └── skill-usage.jsonl     # All skill events across all projects
```

## Hook Scripts

- **track-prompts.sh**: Logs user prompts with timestamps, session IDs, and project directory
- **track-skill-start.sh**: Logs when skills begin execution with project context
- **track-skill-end.sh**: Logs skill completion, duration, and token usage

## Analysis Capabilities

The analysis scripts provide:

- Most frequently used skills (across all projects)
- Average execution time per skill
- Token usage breakdown (input, output, cache read, cache creation)
- Common prompt patterns
- Cost estimation
- Filter by project using the `project` field in logs

## Workflow

1. Create global directory structure (`~/.claude/hooks`, `~/.claude/scripts`, `~/.claude/activity-logs`)
2. Install hook scripts with proper permissions
3. Install analysis scripts
4. Configure hooks in `~/.claude/settings.json`
5. Verify installation

## Integration

This skill is typically run once to enable global tracking across all projects.
The centralized analysis data helps identify opportunities for new skills and
optimization of existing workflows. Since all logs include the project directory,
you can filter analysis by project when needed:

```bash
# View skills used in a specific project
jq 'select(.project | contains("my-project"))' ~/.claude/activity-logs/skill-usage.jsonl
```

## Version

**v2.0** - Global-only installation (December 2025)
