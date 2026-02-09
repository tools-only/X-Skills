---
name: langfuse-ci-integration
description: |
  Configure Langfuse CI/CD integration with GitHub Actions and automated testing.
  Use when setting up automated testing, configuring CI pipelines,
  or integrating Langfuse tests into your build process.
  Trigger with phrases like "langfuse CI", "langfuse GitHub Actions",
  "langfuse automated tests", "CI langfuse", "langfuse pipeline".
allowed-tools: Read, Write, Edit, Bash(gh:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse CI Integration

## Overview
Set up CI/CD pipelines for Langfuse integrations with automated testing and validation.

## Prerequisites
- GitHub repository with Actions enabled
- Langfuse test project API keys
- npm/pnpm project configured

## Instructions

### Step 1: Create GitHub Actions Workflow

Create `.github/workflows/langfuse-integration.yml`:

```yaml
name: Langfuse Integration Tests

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]
  schedule:
    # Run daily to catch Langfuse API changes
    - cron: '0 6 * * *'

env:
  NODE_VERSION: '20'

jobs:
  test:
    runs-on: ubuntu-latest
    environment: test

    env:
      LANGFUSE_PUBLIC_KEY: ${{ secrets.LANGFUSE_PUBLIC_KEY_TEST }}
      LANGFUSE_SECRET_KEY: ${{ secrets.LANGFUSE_SECRET_KEY_TEST }}
      LANGFUSE_HOST: https://cloud.langfuse.com

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: 'npm'

      - name: Install dependencies
        run: npm ci

      - name: Run unit tests
        run: npm test -- --coverage

      - name: Run Langfuse integration tests
        run: npm run test:integration
        timeout-minutes: 10

      - name: Upload coverage
        uses: codecov/codecov-action@v4
        with:
          token: ${{ secrets.CODECOV_TOKEN }}

  validate-traces:
    runs-on: ubuntu-latest
    needs: test
    if: github.event_name == 'push' && github.ref == 'refs/heads/main'

    env:
      LANGFUSE_PUBLIC_KEY: ${{ secrets.LANGFUSE_PUBLIC_KEY_TEST }}
      LANGFUSE_SECRET_KEY: ${{ secrets.LANGFUSE_SECRET_KEY_TEST }}

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: 'npm'

      - name: Install dependencies
        run: npm ci

      - name: Validate Langfuse traces
        run: npm run validate:langfuse
```

### Step 2: Configure Secrets

```bash
# Set test environment secrets
gh secret set LANGFUSE_PUBLIC_KEY_TEST --body "pk-lf-test-..."
gh secret set LANGFUSE_SECRET_KEY_TEST --body "sk-lf-test-..."

# Set production secrets (for deploy workflow)
gh secret set LANGFUSE_PUBLIC_KEY_PROD --body "pk-lf-prod-..."
gh secret set LANGFUSE_SECRET_KEY_PROD --body "sk-lf-prod-..."
```

### Step 3: Create Integration Tests

```typescript
// tests/langfuse.integration.test.ts
import { describe, it, expect, beforeAll, afterAll } from "vitest";
import { Langfuse } from "langfuse";

const SKIP_INTEGRATION = !process.env.LANGFUSE_PUBLIC_KEY;

describe.skipIf(SKIP_INTEGRATION)("Langfuse Integration", () => {
  let langfuse: Langfuse;

  beforeAll(() => {
    langfuse = new Langfuse({
      publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
      secretKey: process.env.LANGFUSE_SECRET_KEY!,
      baseUrl: process.env.LANGFUSE_HOST,
    });
  });

  afterAll(async () => {
    await langfuse.shutdownAsync();
  });

  it("should connect to Langfuse", async () => {
    const trace = langfuse.trace({
      name: "ci-connection-test",
      metadata: {
        ci: true,
        commit: process.env.GITHUB_SHA,
        workflow: process.env.GITHUB_WORKFLOW,
      },
    });

    expect(trace.id).toBeDefined();
    await langfuse.flushAsync();
  });

  it("should create complete trace with generation", async () => {
    const trace = langfuse.trace({
      name: "ci-full-trace-test",
      input: { query: "test query" },
    });

    const span = trace.span({
      name: "process",
      input: { step: 1 },
    });

    const generation = span.generation({
      name: "llm-call",
      model: "gpt-4",
      input: [{ role: "user", content: "test" }],
    });

    generation.end({
      output: "test response",
      usage: { promptTokens: 5, completionTokens: 10 },
    });

    span.end({ output: { processed: true } });
    trace.update({ output: { success: true } });

    await langfuse.flushAsync();

    // Verify trace URL is valid
    const traceUrl = trace.getTraceUrl();
    expect(traceUrl).toContain("langfuse.com");
  });

  it("should handle errors gracefully", async () => {
    const trace = langfuse.trace({
      name: "ci-error-test",
    });

    trace.update({
      level: "ERROR",
      statusMessage: "Test error for CI validation",
    });

    await langfuse.flushAsync();
    expect(trace.id).toBeDefined();
  });
});
```

### Step 4: Add Validation Script

```typescript
// scripts/validate-langfuse.ts
import { Langfuse } from "langfuse";

async function validateLangfuseSetup() {
  console.log("=== Langfuse CI Validation ===\n");

  const langfuse = new Langfuse();

  const checks = [
    {
      name: "Create trace",
      test: async () => {
        const trace = langfuse.trace({
          name: "ci-validation",
          metadata: { validation: true },
        });
        return !!trace.id;
      },
    },
    {
      name: "Create span",
      test: async () => {
        const trace = langfuse.trace({ name: "span-test" });
        const span = trace.span({ name: "test-span" });
        span.end();
        return !!span.id;
      },
    },
    {
      name: "Create generation",
      test: async () => {
        const trace = langfuse.trace({ name: "gen-test" });
        const gen = trace.generation({
          name: "test-gen",
          model: "test-model",
        });
        gen.end({ output: "test" });
        return !!gen.id;
      },
    },
    {
      name: "Flush data",
      test: async () => {
        await langfuse.flushAsync();
        return true;
      },
    },
    {
      name: "Shutdown gracefully",
      test: async () => {
        await langfuse.shutdownAsync();
        return true;
      },
    },
  ];

  let passed = 0;
  let failed = 0;

  for (const { name, test } of checks) {
    try {
      const result = await test();
      if (result) {
        console.log(`✓ ${name}`);
        passed++;
      } else {
        console.log(`✗ ${name}`);
        failed++;
      }
    } catch (error) {
      console.log(`✗ ${name}: ${error}`);
      failed++;
    }
  }

  console.log(`\nResults: ${passed}/${checks.length} passed`);

  if (failed > 0) {
    process.exit(1);
  }
}

validateLangfuseSetup();
```

### Step 5: Add package.json Scripts

```json
{
  "scripts": {
    "test": "vitest run",
    "test:integration": "vitest run tests/*.integration.test.ts",
    "validate:langfuse": "tsx scripts/validate-langfuse.ts"
  }
}
```

## Output
- Automated test pipeline
- PR checks configured
- Coverage reports uploaded
- Daily validation runs
- Trace validation on merge

## PR Status Checks

```yaml
# .github/branch-protection.yml (reference)
required_status_checks:
  strict: true
  contexts:
    - "test"
    - "validate-traces"
```

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Secret not found | Missing configuration | Add secret via `gh secret set` |
| Tests timeout | Network issues | Increase timeout or add retry |
| Auth failures | Invalid key | Check secret value |
| Flaky tests | Race conditions | Add proper flush/wait |

## Examples

### Matrix Testing (Multiple Node Versions)
```yaml
jobs:
  test:
    strategy:
      matrix:
        node-version: [18, 20, 22]

    steps:
      - uses: actions/setup-node@v4
        with:
          node-version: ${{ matrix.node-version }}
      # ... rest of steps
```

### Release Workflow with Langfuse Validation
```yaml
name: Release

on:
  push:
    tags: ['v*']

jobs:
  release:
    runs-on: ubuntu-latest

    steps:
      - uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: '20'
          registry-url: 'https://registry.npmjs.org'

      - name: Install dependencies
        run: npm ci

      - name: Validate Langfuse integration
        env:
          LANGFUSE_PUBLIC_KEY: ${{ secrets.LANGFUSE_PUBLIC_KEY_PROD }}
          LANGFUSE_SECRET_KEY: ${{ secrets.LANGFUSE_SECRET_KEY_PROD }}
        run: npm run validate:langfuse

      - name: Build
        run: npm run build

      - name: Publish
        run: npm publish
        env:
          NODE_AUTH_TOKEN: ${{ secrets.NPM_TOKEN }}
```

## Resources
- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Vitest Documentation](https://vitest.dev)
- [Langfuse SDK Testing](https://langfuse.com/docs/sdk)

## Next Steps
For deployment patterns, see `langfuse-deploy-integration`.
