---
skill_name: ci-cd-patterns
skill_category: devops
description: CI/CD pipeline patterns, GitHub Actions, build automation best practices
allowed_tools: [Read, Edit, Glob, Grep, Write]
token_estimate: 800
version: 1.0
last_updated: 2026-01-13
owner: DevOps Team
status: active

tags: [ci-cd, github-actions, pipeline, anti-pattern, best-practice, automation]
related_skills: [kubernetes, docker-patterns, terraform-patterns]

trigger_files: [".github/workflows/*.yml", ".github/workflows/*.yaml", "Jenkinsfile", ".gitlab-ci.yml", "azure-pipelines.yml", ".circleci/config.yml"]
trigger_keywords: [ci/cd, pipeline, github actions, workflow, build, deploy, continuous integration, continuous deployment]

quality_keywords: [anti-pattern, pattern, validation, best-practice, security, secrets]
---

# CI/CD Patterns

Best practices and anti-patterns for CI/CD pipelines, focusing on GitHub Actions but applicable to other CI systems.

## Purpose

- Ensure pipelines are fast, reliable, and secure
- Prevent common mistakes that cause build failures or security issues
- Establish patterns for maintainable pipeline configurations

---

## Core Patterns

### Pattern 1: Parallel Job Execution

**When to use:** Multiple independent tasks in pipeline.

**Implementation:**
```yaml
jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: npm run lint

  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: npm test

  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: npm run build

  # Only deploy after all checks pass
  deploy:
    needs: [lint, test, build]
    runs-on: ubuntu-latest
    steps:
      - run: ./deploy.sh
```

**Benefits:**
- Faster pipeline execution
- Early failure detection
- Clear dependency visualization

### Pattern 2: Dependency Caching

**When to use:** Repeated dependency installation across runs.

**Implementation:**
```yaml
- name: Cache dependencies
  uses: actions/cache@v4
  with:
    path: ~/.npm
    key: ${{ runner.os }}-node-${{ hashFiles('**/package-lock.json') }}
    restore-keys: |
      ${{ runner.os }}-node-

- name: Install dependencies
  run: npm ci
```

**Benefits:**
- 60-80% faster builds
- Reduced network usage
- Consistent dependency versions

### Pattern 3: Matrix Builds

**When to use:** Testing across multiple versions/platforms.

**Implementation:**
```yaml
strategy:
  matrix:
    node: [18, 20, 22]
    os: [ubuntu-latest, macos-latest]
  fail-fast: false

runs-on: ${{ matrix.os }}
steps:
  - uses: actions/setup-node@v4
    with:
      node-version: ${{ matrix.node }}
```

---

## Anti-Patterns

### Anti-Pattern 1: Hardcoded Secrets

| Aspect | Description |
|--------|-------------|
| **WHY** | Secrets in code are exposed in logs, version history, and forks. Security breach risk. |
| **DETECTION** | Look for: API keys, passwords, tokens directly in workflow files. Strings like `password: "..."` or `token: abc123`. |
| **FIX** | Use GitHub Secrets or environment-specific secret stores. Reference via `${{ secrets.NAME }}`. |

**Bad Example:**
```yaml
- name: Deploy
  run: |
    curl -H "Authorization: Bearer ghp_xxxxxxxxxxxx" \
      https://api.github.com/repos/owner/repo/deployments
```

**Good Example:**
```yaml
- name: Deploy
  run: |
    curl -H "Authorization: Bearer ${{ secrets.GITHUB_TOKEN }}" \
      https://api.github.com/repos/owner/repo/deployments
```

### Anti-Pattern 2: No Build Caching

| Aspect | Description |
|--------|-------------|
| **WHY** | Repeated full installs waste 5-15 minutes per build. Increases costs and slows feedback. |
| **DETECTION** | Look for `npm install` or `pip install` without preceding cache step. Check build times >5 min. |
| **FIX** | Add dependency caching. Use `actions/cache@v4` for npm, pip, cargo, etc. |

**Bad Example:**
```yaml
steps:
  - uses: actions/checkout@v4
  - run: npm install
  - run: npm test
```

**Good Example:**
```yaml
steps:
  - uses: actions/checkout@v4
  - uses: actions/cache@v4
    with:
      path: ~/.npm
      key: ${{ runner.os }}-node-${{ hashFiles('**/package-lock.json') }}
  - run: npm ci
  - run: npm test
```

### Anti-Pattern 3: Missing Health Checks

| Aspect | Description |
|--------|-------------|
| **WHY** | Deployment without verification can ship broken code to production silently. |
| **DETECTION** | Deploy steps without subsequent health/smoke tests. No rollback on failure. |
| **FIX** | Add health check after deploy. Implement automatic rollback on failure. |

**Bad Example:**
```yaml
deploy:
  steps:
    - name: Deploy
      run: kubectl apply -f deployment.yaml
    # No verification!
```

**Good Example:**
```yaml
deploy:
  steps:
    - name: Deploy
      run: kubectl apply -f deployment.yaml

    - name: Wait for rollout
      run: kubectl rollout status deployment/app --timeout=5m

    - name: Health check
      run: |
        curl -f https://app.example.com/health || exit 1

    - name: Rollback on failure
      if: failure()
      run: kubectl rollout undo deployment/app
```

### Anti-Pattern 4: Sequential Independent Jobs

| Aspect | Description |
|--------|-------------|
| **WHY** | Running independent tasks sequentially wastes time. 30-min pipeline could be 10-min. |
| **DETECTION** | Jobs with `needs:` dependencies that don't actually depend on each other's output. |
| **FIX** | Run independent jobs in parallel. Only use `needs:` for true dependencies. |

**Bad Example:**
```yaml
jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - run: npm run lint

  test:
    needs: lint  # Unnecessary dependency!
    runs-on: ubuntu-latest
    steps:
      - run: npm test
```

**Good Example:**
```yaml
jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - run: npm run lint

  test:
    runs-on: ubuntu-latest  # Runs in parallel
    steps:
      - run: npm test
```

### Anti-Pattern 5: No Timeout Configuration

| Aspect | Description |
|--------|-------------|
| **WHY** | Hanging builds waste resources and block pipelines. Default 6-hour timeout is too long. |
| **DETECTION** | No `timeout-minutes` at job or step level. Jobs that occasionally hang. |
| **FIX** | Set explicit timeouts based on expected duration + buffer (typically 2-3x normal). |

**Bad Example:**
```yaml
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - run: npm test  # Could hang forever
```

**Good Example:**
```yaml
jobs:
  build:
    runs-on: ubuntu-latest
    timeout-minutes: 30
    steps:
      - name: Run tests
        timeout-minutes: 15
        run: npm test
```

---

## Validation Checklist

### Pre-Implementation
- [ ] Review existing pipeline structure
- [ ] Identify parallelization opportunities
- [ ] Check current build times for baseline

### Implementation
- [ ] All secrets use `${{ secrets.NAME }}` syntax
- [ ] Dependency caching configured
- [ ] Independent jobs run in parallel
- [ ] Timeouts set on jobs and long-running steps
- [ ] Health checks after deployments
- [ ] Rollback mechanism in place

### Post-Implementation
- [ ] Build times improved or maintained
- [ ] No secrets exposed in logs
- [ ] Pipeline passes consistently
- [ ] Deployment verification works

---

## Quick Reference

| Task | Pattern |
|------|---------|
| Speed up builds | Add caching, parallelize jobs |
| Secure secrets | Use `${{ secrets.NAME }}` |
| Test multiple versions | Matrix builds |
| Safe deployments | Health checks + rollback |
| Prevent hangs | Set timeout-minutes |
