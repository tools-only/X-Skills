# Git Operations in Autonomous Mode

When a skill activates autonomous mode (`autonomous-state.json` exists), the skill's
git instructions override the default "ask before committing" behavior. The user's
invocation of the skill IS their explicit permission for autonomous git operations.

**When autonomous mode is active:**
- Commit after each coherent logical change (not after every file edit — after each iteration, batch, or fix cycle)
- Push after successful lint/test verification
- Use conventional commit prefixes matching the skill: `feat:`, `fix:`, `refactor:`, `improve():`, `burndown:`, `appfix:`

**When autonomous mode is NOT active (normal interactive session):**
- Follow the default: commit and push only when the user asks
