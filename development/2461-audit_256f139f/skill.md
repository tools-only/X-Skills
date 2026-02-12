---
description: Audit Namshub installation — show hooks, skills, permissions, data storage, and security posture
---

# /audit

Security and configuration audit for the Namshub toolkit. Shows exactly what's installed, what permissions are granted, and where data is stored.

## Purpose

Before activating autonomous mode, users should understand what Namshub does. This command provides full transparency.

## Execution

Run a comprehensive audit of the current Namshub installation:

```bash
python3 -c "
import json, os, sys
from pathlib import Path

home = Path.home()
claude_dir = home / '.claude'

print('=' * 60)
print('NAMSHUB AUDIT REPORT')
print('=' * 60)

# 1. Hooks
print('\n## Installed Hooks\n')
settings_path = claude_dir / 'settings.json'
if settings_path.exists():
    settings = json.loads(settings_path.read_text())
    hooks = settings.get('hooks', {})
    for event, hook_list in hooks.items():
        if hook_list:
            print(f'  {event}:')
            for entry in hook_list:
                for h in entry.get('hooks', []):
                    cmd = h.get('command', h.get('prompt', ''))[:80]
                    htype = h.get('type', 'command')
                    timeout = h.get('timeout', 'none')
                    matcher = entry.get('matcher', '*')
                    print(f'    [{htype}] matcher={matcher} timeout={timeout}s')
                    print(f'           {cmd}')
else:
    print('  No settings.json found')

# 2. Skills
print('\n## Installed Skills\n')
skills_dir = claude_dir / 'skills'
if skills_dir.exists():
    for skill_dir in sorted(skills_dir.iterdir()):
        if skill_dir.is_dir():
            skill_md = skill_dir / 'SKILL.md'
            if skill_md.exists():
                # Read first line after frontmatter for description
                lines = skill_md.read_text().split('\n')
                name = skill_dir.name
                print(f'  /{name}')
    print(f'\n  Total: {sum(1 for d in skills_dir.iterdir() if d.is_dir() and (d / \"SKILL.md\").exists())} skills')
else:
    print('  No skills directory found')

# 3. Commands
print('\n## Installed Commands\n')
commands_dir = claude_dir / 'commands'
if commands_dir.exists():
    for cmd_file in sorted(commands_dir.glob('*.md')):
        print(f'  /{cmd_file.stem}')
    print(f'\n  Total: {sum(1 for _ in commands_dir.glob(\"*.md\"))} commands')
else:
    print('  No commands directory found')

# 4. Permissions
print('\n## Permissions & Environment\n')
if settings_path.exists():
    settings = json.loads(settings_path.read_text())
    perms = settings.get('permissions', {})
    env = settings.get('env', {})
    print(f'  Default permission mode: {perms.get(\"defaultMode\", \"default\")}')
    print(f'  Agent teams enabled: {env.get(\"CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS\", \"0\")}')
    print(f'  Max output tokens: {env.get(\"CLAUDE_CODE_MAX_OUTPUT_TOKENS\", \"default\")}')
    print(f'  Always thinking: {settings.get(\"alwaysThinkingEnabled\", False)}')
else:
    print('  No settings.json found')

# 5. Data Storage
print('\n## Data Storage Locations\n')
memory_dir = home / '.claude' / 'memory'
if memory_dir.exists():
    project_count = sum(1 for d in memory_dir.iterdir() if d.is_dir())
    total_events = 0
    for proj in memory_dir.iterdir():
        events_dir = proj / 'events'
        if events_dir.exists():
            total_events += sum(1 for f in events_dir.glob('evt_*.json'))
    print(f'  Memory store: {memory_dir}')
    print(f'    Projects tracked: {project_count}')
    print(f'    Total events: {total_events}')
else:
    print(f'  Memory store: not initialized')

if memory_dir.exists():
    total_assertions = 0
    for proj in memory_dir.iterdir():
        ap = proj / 'core-assertions.jsonl'
        if ap.exists():
            total_assertions += sum(1 for line in ap.read_text().strip().split(chr(10)) if line.strip())
    if total_assertions:
        print(f'  Core assertions: {total_assertions} entries')

state_files = [
    '.claude/autonomous-state.json',
    '.claude/completion-checkpoint.json',
    '.claude/tool-usage-log.json',
    '.claude/session-snapshot.json',
]
print(f'\n  Session state files (project-local):')
for sf in state_files:
    exists = Path(sf).exists()
    print(f'    {sf}: {\"exists\" if exists else \"not present\"}')

# 6. Security Posture
print('\n## Security Posture\n')
guard_path = claude_dir / 'hooks' / 'cloud-command-guard.sh'
print(f'  Cloud command guard: {\"active\" if guard_path.exists() else \"not installed\"}')

state_guard = claude_dir / 'hooks' / 'state-file-guard.py'
print(f'  State file guard: {\"active\" if state_guard.exists() else \"not installed\"}')

deploy_guard = claude_dir / 'hooks' / 'deploy-enforcer.py'
print(f'  Deploy enforcer: {\"active\" if deploy_guard.exists() else \"not installed\"}')

stop_val = claude_dir / 'hooks' / 'stop-validator.py'
print(f'  Stop validator: {\"active\" if stop_val.exists() else \"not installed\"}')

print('\n' + '=' * 60)
print('Audit complete. Review the above before enabling autonomous mode.')
print('=' * 60)
"
```

## Output Sections

| Section | What It Shows |
|---------|---------------|
| **Installed Hooks** | Every hook registered in settings.json with event, matcher, type, and timeout |
| **Installed Skills** | All skills available via `/skillname` |
| **Installed Commands** | All slash commands available |
| **Permissions** | Default permission mode, agent teams, output limits |
| **Data Storage** | Memory event counts, project tracking, session state files |
| **Security Posture** | Which safety guards are active (cloud guard, state guard, deploy enforcer, stop validator) |

## When to Use

- **After installation**: Verify everything installed correctly
- **Before autonomous mode**: Understand what hooks will fire and what data is stored
- **Periodic check**: Ensure no unexpected changes to configuration
- **Debugging**: Verify hooks are registered and files exist
