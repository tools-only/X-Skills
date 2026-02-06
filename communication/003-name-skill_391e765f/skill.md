---
name: independent-research
description: "Use when about to ask the user a factual question, propose a solution, diagnose an error, or choose between approaches. Triggers on: 'Do you have X installed?', 'What version?', 'Is X configured?', 'We should...', 'The fix is...', 'Options: 1...', 'Based on my understanding...', 'I believe X supports...'. Before deciding anything, spin up parallel subagents to WebSearch for current docs, community solutions, framework best practices, and GitHub issues. Your memory is stale — verify everything."
version: 1.2.0
---

# Independent Research

Research-driven investigation. Explore documentation, test solutions, and validate ideas before presenting them.

## Critical Rules

🚨 **NEVER ASK FACTUAL QUESTIONS YOU CAN ANSWER YOURSELF.** If a command, config file, or web search can answer it, use that. Only ask users about preferences and priorities.

🚨 **VALIDATE BEFORE PRESENTING.** Test commands, verify syntax, check documentation. Never present untested recommendations.

## Lazy Question Detection

If you catch yourself about to ask any of these, STOP. Answer it yourself.

| Lazy Question | What To Do Instead |
|---|---|
| "Do you have X installed?" | Run `which X` or `X --version` |
| "What version of X?" | Run `X --version` |
| "Is X configured?" | Check the config file |
| "What OS are you using?" | Run `uname -a` |
| "Which tool do you use for X?" | Check `which`/`where`, look at project config |
| "Was this working before?" | Check git log, recent changes |
| "Do you want me to..." | Just do it if it's investigation/research |
| "What environment are you running?" | Check env vars, config files, runtime versions |
| "Have you tried X?" | Try it yourself first |

## Anti-pattern

### ❌ Asking Instead of Investigating

```
⏺ NX project graph failed — Node.js v24 compatibility issue.

⏺ Do you have nvm/fnm/volta installed? Was this working
  on a different Node version?
```

The user doesn't need to answer this. You have bash.

```bash
which nvm; which fnm; which volta
node --version
cat .nvmrc 2>/dev/null
cat package.json | grep -A2 '"engines"'
git log --oneline -20
```

## Premature Decision Detection

If you catch yourself about to recommend, diagnose, or decide without having researched, STOP.

| Premature Decision Signal | What To Do Instead |
|---|---|
| "We should downgrade/upgrade X" | Spin up parallel subagents: one to WebSearch the error + framework version, one to fetch the framework's migration/compatibility docs, one to search GitHub issues for the exact error |
| "The fix is to..." | Search for the exact error message first. Find the GitHub issue. Read the fix. Don't guess from memory |
| "This is a known issue with X" | Prove it. Fetch the GitHub issue URL, the changelog entry, the docs page. If you can't link to it, you're guessing |
| "Options: 1. X, 2. Y" | Before presenting options, WebSearch for what the framework/community actually recommends. Check migration guides, official docs, recent blog posts |
| "I believe X supports Y" / "Based on my understanding" | That's stale memory. WebFetch the actual docs page. Check the release notes. Verify the version compatibility matrix |
| Choosing between tools/versions/approaches | Launch parallel subagents: one per approach, each researching current docs, community consensus, and known issues |
| Error you haven't seen before | WebSearch the exact error message in quotes. Check GitHub issues. Check Stack Overflow. Check the framework's Discord/discussions |
| Recommending a config change | Fetch the framework's current configuration docs. Don't rely on what you remember the API being |

🚨 **Be paranoid that you're missing a better solution.** Your training data is stale. Libraries release weekly. Frameworks ship breaking changes. The answer you "know" may be outdated. Launch subagents to verify in parallel — it costs seconds and saves hours of wrong-direction work.

🚨 **Concrete actions, not vague research:**
- **WebSearch** the exact error message, the framework + version, the specific compatibility question
- **WebFetch** the official docs page, the changelog, the migration guide
- **Launch parallel subagents** when multiple angles need investigating simultaneously
- **Check GitHub issues** for the exact error — someone has probably hit this before
- **Read the actual release notes** — don't guess what version supports what

## When to Ask vs Research

**Research yourself (facts):**
- What's installed, what version, what's configured
- Whether something is compatible
- What the error means
- What solutions exist

**Ask the user (preferences):**
- Which approach they prefer
- What their priorities are
- Design decisions with multiple valid options
- Business context you can't infer

## Mandatory Checklist

Before asking the user a question:

1. [ ] Verify this is a preference/priority question, NOT a factual question
2. [ ] Verify you cannot answer it with a command, config file read, or web search
3. [ ] Verify you have already tried to answer it yourself

Do not ask until all checks pass.

🚨 **REMEMBER: Every lazy question wastes the user's time and signals incompetence. If you can look it up, look it up.**
