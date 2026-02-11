# Automated Code Reviews with Hodor

Hodor is optimized for **automated, READ-ONLY code reviews** in CI/CD pipelines. This document explains how to set up and use Hodor for continuous code quality monitoring.

## Overview

Hodor provides AI-powered code reviews that:
- ✅ Run automatically on every PR/MR (or on-demand)
- ✅ Use efficient search tools (grep, glob, planning_file_editor)
- ✅ Post review comments directly to GitHub/GitLab
- ✅ Focus on bugs, security, and performance (not style)
- ✅ Support repo-specific guidelines via skills system
- ✅ Work with self-hosted GitLab instances

## Key Features for CI/CD

### 1. Enhanced Tooling for Fast Reviews

Hodor uses specialized tools for efficient code analysis:

| Tool | Purpose | Example Use |
|------|---------|-------------|
| **grep** | Pattern search across codebase | Find all `TODO`, `FIXME`, null checks, error patterns |
| **glob** | File pattern matching | Find all `*.test.js`, `**/*.py`, config files |
| **planning_file_editor** | Read-optimized file viewer | View code with line numbers, search within files |
| **terminal** | Git and CLI commands | Run `git diff`, `git log`, language-specific linters |

These tools are **much faster** than having the LLM script search operations, reducing review time and token usage.

### 2. READ-ONLY by Design

Hodor is designed for **automated reviews without human intervention**:

- Workspace is a fresh clone (agent can't push changes)
- No confirmation policy needed (agent only analyzes, doesn't modify)
- Environment is isolated (temporary directory, cleaned up after)
- Safe for untrusted PRs (no risk of malicious code execution)

### 3. Repository-Specific Guidelines

Use the **skills system** to encode your team's standards:

```bash
# In your repository root
.cursorrules                    # Simple, single-file project guidelines
agents.md                       # Alternative single-file location
.hodor/skills/*.md              # Modular skills organized by topic
```

See [SKILLS.md](./SKILLS.md) for detailed documentation.

## Setup for GitHub Actions

### Quick Start

1. **Add Secrets** to your GitHub repository:
   - `LLM_API_KEY` or `ANTHROPIC_API_KEY`: Your Claude API key
   - `GITHUB_TOKEN` is automatically provided by GitHub Actions

2. **Copy the workflow file**:
   ```bash
   cp .github/workflows/pr-review.yml your-repo/.github/workflows/
   ```

3. **Trigger Options**:
   - **On-demand**: Add label `hodor-review` to a PR
   - **Auto-reviewer**: Request `hodor-agent` as a reviewer
   - **Always-on**: Modify workflow to trigger on all PRs

### Example Workflow

```yaml
# .github/workflows/hodor-review.yml
name: AI Code Review

on:
  pull_request:
    types: [opened, synchronize]  # Run on every PR

jobs:
  review:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Install Hodor
        run: |
          pip install uv
          git clone https://github.com/mr-karan/hodor /tmp/hodor
          cd /tmp/hodor && uv sync

      - name: Run Review
        env:
          LLM_API_KEY: ${{ secrets.LLM_API_KEY }}
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
        run: |
          cd /tmp/hodor
          PR_URL="${{ github.event.pull_request.html_url }}"
          uv run hodor "$PR_URL" --post --verbose
```

## Setup for GitLab CI

### Quick Start

1. **Set CI/CD Variables** in GitLab Settings > CI/CD > Variables:
   - `LLM_API_KEY` or `ANTHROPIC_API_KEY`: Your Claude API key
   - `GITLAB_TOKEN`: GitLab access token with `api` scope
   - `GITLAB_HOST` (if self-hosted): e.g., `gitlab.company.com`

2. **Copy the CI configuration**:
   ```bash
   cp .gitlab-ci-example.yml your-repo/.gitlab-ci.yml
   ```

3. **Configure triggers** (all PRs, specific labels, etc.)

### Example CI Configuration

```yaml
# .gitlab-ci.yml
hodor-review:
  image: python:3.13-slim
  stage: test

  before_script:
    - apt-get update && apt-get install -y git curl
    - curl -LsSf https://astral.sh/uv/install.sh | sh
    - export PATH="$HOME/.local/bin:$PATH"

  script:
    - git clone https://github.com/mr-karan/hodor /tmp/hodor
    - cd /tmp/hodor && uv sync
    - export MR_URL="${CI_MERGE_REQUEST_PROJECT_URL}/-/merge_requests/${CI_MERGE_REQUEST_IID}"
    - uv run hodor "$MR_URL" --post

  rules:
    - if: $CI_PIPELINE_SOURCE == "merge_request_event"
  allow_failure: true
```

## Configuration Options

### Model Selection

Choose the right model for your needs:

```yaml
# Fast and cost-effective (default)
LLM_MODEL: "anthropic/claude-sonnet-4-5-20250929"

# More thorough reviews
LLM_MODEL: "anthropic/claude-opus-4"

# OpenAI alternative
LLM_MODEL: "openai/gpt-4"

# Reasoning models for complex analysis
LLM_MODEL: "openai/o3-mini"
--reasoning-effort high
```

### Verbose Logging

Enable detailed logging for debugging:

```bash
uv run hodor $PR_URL --post --verbose
```

This shows:
- Which tools are configured
- Real-time command execution (grep, git, file reads)
- Token usage statistics
- Detailed error messages

### Custom Prompts

Override the default review prompt:

```bash
# Inline prompt
uv run hodor $PR_URL --prompt "Focus on security issues only..."

# From file (create your own custom prompt file)
uv run hodor $PR_URL --prompt-file custom-prompt.txt
```

## Cost Optimization

Hodor is designed to be cost-effective:

### Token Usage

| Review Type | Typical Tokens | Cost (Claude Sonnet 4.5) |
|-------------|---------------|--------------------------|
| Small PR (<5 files) | 5K-15K | $0.03-$0.09 |
| Medium PR (5-15 files) | 15K-40K | $0.09-$0.24 |
| Large PR (15+ files) | 40K-100K | $0.24-$0.60 |

*Costs based on Claude Sonnet 4.5 pricing ($3/M input, $15/M output)*

### Optimization Tips

1. **Use Efficient Tools**:
   - grep/glob are much faster than LLM-scripted searches
   - planning_file_editor is optimized for reading code
   - These tools reduce token usage by 30-50%

2. **Focus Reviews**:
   - Review only changed files (default behavior)
   - Use custom prompts for specific concerns
   - Skip trivial changes (docs, formatting) with CI rules

3. **Smart Triggers**:
   - Run on critical PRs only (use labels)
   - Skip draft PRs
   - Run once per PR (not on every commit)

Example GitLab CI rule:
```yaml
rules:
  - if: $CI_MERGE_REQUEST_LABELS !~ /skip-review/
    when: on_success
  - if: $CI_MERGE_REQUEST_DRAFT == "true"
    when: never
```

## Security Considerations

### Safe for Untrusted Code

Hodor's design makes it safe for reviewing untrusted PRs:

- **Isolated Environment**: Each review runs in a fresh, temporary workspace
- **No Write Access**: Agent can read code but can't push changes
- **No Credentials**: Agent doesn't have access to push/deploy credentials
- **Automatic Cleanup**: Workspace is deleted after review

### Protecting Secrets

Best practices for CI/CD:

1. **Never log secrets**: Hodor doesn't log API keys or tokens
2. **Use CI/CD variables**: Store sensitive data in GitHub Secrets / GitLab CI Variables
3. **Limit token scope**: GITLAB_TOKEN only needs `api` scope (not `write_repository`)
4. **Rotate tokens**: Periodically rotate API keys

### No Confirmation Policy Needed

Unlike interactive development, automated reviews don't need confirmation:
- Agent is READ-ONLY (can't modify repository)
- Agent can't access external resources (network isolation in CI)
- All actions are logged (audit trail in CI logs)

## Troubleshooting

### Review Fails with "command too long"

**Cause**: Large environment variables (DIRENV_DIFF, LS_COLORS)
**Solution**: Hodor automatically uses subprocess terminal (fixed in v1.0+)

### "No LLM API key found"

**Cause**: Missing CI/CD variable
**Solution**: Set `LLM_API_KEY` or `ANTHROPIC_API_KEY` in CI/CD settings

### "Failed to clone repository"

**Cause**: Missing or invalid GITLAB_TOKEN
**Solution**: Generate token with `api` scope in GitLab Settings > Access Tokens

### "Review finds no issues on obviously buggy code"

**Cause**: Model may need more context or different prompt
**Solution**: Try:
1. Use `--verbose` to see what agent is checking
2. Add repo-specific guidelines in `.cursorrules`
3. Use more capable model (`claude-opus-4`)
4. Increase `--reasoning-effort high`

### "glab: command not found" (GitLab CI)

**Cause**: GitLab CLI not installed
**Solution**: Install in `before_script`:
```bash
curl -fsSL https://gitlab.com/gitlab-org/cli/-/releases/permalink/latest/downloads/glab_linux_amd64.deb -o /tmp/glab.deb
dpkg -i /tmp/glab.deb
```

## Examples

### Review with Security Focus

```bash
# Custom prompt focusing on security
hodor $PR_URL --prompt "Review for security vulnerabilities: SQL injection, XSS, auth bypasses, secrets in code. Be thorough." --post
```

### Review Specific Files

```bash
# Review only backend changes
hodor $PR_URL --prompt "Review only files in src/backend/. Focus on database queries and API security." --post
```

### High-Thoroughness Review

```bash
# Use reasoning model for complex PRs
hodor $PR_URL --model anthropic/claude-opus-4 --reasoning-effort high --post
```

## Monitoring and Metrics

### Track Review Quality

Monitor review effectiveness:
- **True Positive Rate**: Issues found that are real bugs
- **False Positive Rate**: Issues reported that aren't actually problems
- **Coverage**: Percentage of bugs caught before production

### Track Costs

Monitor API usage:
```bash
# Verbose mode shows token usage
hodor $PR_URL --verbose --post | grep "Total tokens"

# Example output:
# Total tokens used: 23,456
# Cost estimate: $0.14
```

### Improve Over Time

1. **Analyze false positives**: Update `.cursorrules` to reduce noise
2. **Track missed bugs**: Add patterns to skills system
3. **Optimize prompts**: Focus on high-value checks
4. **Adjust triggers**: Run reviews where they provide most value

## Best Practices

### DO:
- ✅ Run automated reviews on all non-trivial PRs
- ✅ Use `.cursorrules` for project-specific guidelines
- ✅ Post reviews automatically (use `--post` flag)
- ✅ Enable verbose logging initially (debug issues)
- ✅ Monitor costs and adjust model based on budget
- ✅ Treat reviews as suggestions (not blockers)

### DON'T:
- ❌ Don't block PRs on review results (use as advisory)
- ❌ Don't review trivial changes (docs, typos, formatting)
- ❌ Don't use expensive models for small PRs
- ❌ Don't ignore review feedback (defeats the purpose)
- ❌ Don't commit secrets in `.cursorrules`

## Roadmap

Upcoming features for automated reviews:

- [ ] **Incremental reviews**: Review only new commits in updated PRs
- [ ] **Batch reviews**: Review multiple PRs in parallel
- [ ] **Learning mode**: Improve from human feedback on reviews
- [ ] **Custom MCP tools**: Integrate with your issue tracker, docs, etc.
- [ ] **Trend analysis**: Track code quality over time
- [ ] **Team metrics**: Compare review quality across teams

## Support

- **Documentation**: See [README.md](../README.md) and [SKILLS.md](./SKILLS.md)
- **Issues**: Report bugs at https://github.com/mr-karan/hodor/issues
- **Questions**: Start a discussion at https://github.com/mr-karan/hodor/discussions

## Conclusion

Hodor provides production-ready, automated code reviews that:
- Catch bugs before they reach production
- Reduce human review burden
- Maintain code quality consistently
- Cost less than $0.50 per typical PR
- Integrate seamlessly with GitHub/GitLab CI

Start with on-demand reviews (label-triggered), measure effectiveness, then graduate to always-on automated reviews once you're confident in the results.

Happy reviewing! 🚪
