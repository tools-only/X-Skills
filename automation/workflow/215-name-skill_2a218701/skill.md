---
name: twinmind-ci-integration
description: |
  Integrate TwinMind into CI/CD pipelines for automated testing and deployment.
  Use when setting up GitHub Actions, GitLab CI, or other CI systems
  with TwinMind API testing and validation.
  Trigger with phrases like "twinmind ci", "twinmind github actions",
  "twinmind pipeline", "automate twinmind testing".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind CI Integration

## Overview
Integrate TwinMind testing and validation into CI/CD pipelines.

## Prerequisites
- TwinMind Pro/Enterprise API access
- CI/CD system (GitHub Actions, GitLab CI, etc.)
- Test audio samples
- Secrets management

## Instructions

### Step 1: GitHub Actions Workflow

```yaml
# .github/workflows/twinmind-tests.yml
name: TwinMind Integration Tests

on:
  push:
    branches: [main, develop]
    paths:
      - 'src/twinmind/**'
      - 'tests/twinmind/**'
  pull_request:
    branches: [main]
  schedule:
    # Run daily at 6 AM UTC
    - cron: '0 6 * * *'

env:
  NODE_VERSION: '20'

jobs:
  lint-and-typecheck:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: 'npm'

      - name: Install dependencies
        run: npm ci

      - name: Run linter
        run: npm run lint

      - name: Run typecheck
        run: npm run typecheck

  unit-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: 'npm'

      - name: Install dependencies
        run: npm ci

      - name: Run unit tests
        run: npm run test:unit
        env:
          TWINMIND_API_KEY: ${{ secrets.TWINMIND_API_KEY_TEST }}

  integration-tests:
    runs-on: ubuntu-latest
    needs: [lint-and-typecheck, unit-tests]
    steps:
      - uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: 'npm'

      - name: Install dependencies
        run: npm ci

      - name: Run integration tests
        run: npm run test:integration
        env:
          TWINMIND_API_KEY: ${{ secrets.TWINMIND_API_KEY_TEST }}
          TWINMIND_WEBHOOK_SECRET: ${{ secrets.TWINMIND_WEBHOOK_SECRET }}
        timeout-minutes: 10

      - name: Upload test results
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: test-results
          path: coverage/

  api-health-check:
    runs-on: ubuntu-latest
    steps:
      - name: Check TwinMind API Health
        run: |
          RESPONSE=$(curl -s -w "\n%{http_code}" \
            -H "Authorization: Bearer ${{ secrets.TWINMIND_API_KEY_TEST }}" \
            https://api.twinmind.com/v1/health)

          HTTP_CODE=$(echo "$RESPONSE" | tail -1)
          BODY=$(echo "$RESPONSE" | head -n -1)

          echo "HTTP Status: $HTTP_CODE"
          echo "Response: $BODY"

          if [ "$HTTP_CODE" != "200" ]; then
            echo "API health check failed"
            exit 1
          fi

  transcription-smoke-test:
    runs-on: ubuntu-latest
    needs: [integration-tests, api-health-check]
    if: github.event_name == 'schedule' || github.ref == 'refs/heads/main'
    steps:
      - uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: 'npm'

      - name: Install dependencies
        run: npm ci

      - name: Run transcription smoke test
        run: npm run test:smoke
        env:
          TWINMIND_API_KEY: ${{ secrets.TWINMIND_API_KEY_TEST }}
          TEST_AUDIO_URL: ${{ secrets.TEST_AUDIO_URL }}
        timeout-minutes: 5
```

### Step 2: Unit Test Setup

```typescript
// tests/unit/twinmind/client.test.ts
import { describe, it, expect, beforeEach, vi } from 'vitest';
import { TwinMindClient } from '../../../src/twinmind/client';

describe('TwinMindClient', () => {
  let client: TwinMindClient;

  beforeEach(() => {
    client = new TwinMindClient({
      apiKey: 'test_api_key',
      baseUrl: 'https://api.twinmind.com/v1',
    });
  });

  describe('transcribe', () => {
    it('should send correct request for transcription', async () => {
      const mockFetch = vi.fn().mockResolvedValue({
        ok: true,
        json: () => Promise.resolve({
          id: 'tr_123',
          status: 'processing',
        }),
      });
      global.fetch = mockFetch;

      const result = await client.transcribe('https://example.com/audio.mp3');

      expect(mockFetch).toHaveBeenCalledWith(
        expect.stringContaining('/transcribe'),
        expect.objectContaining({
          method: 'POST',
          headers: expect.objectContaining({
            'Authorization': 'Bearer test_api_key',
          }),
        })
      );
    });

    it('should handle rate limit errors', async () => {
      const mockFetch = vi.fn().mockResolvedValue({
        ok: false,
        status: 429,
        headers: new Headers({ 'Retry-After': '60' }),
      });
      global.fetch = mockFetch;

      await expect(client.transcribe('https://example.com/audio.mp3'))
        .rejects.toThrow(/rate limit/i);
    });
  });

  describe('summarize', () => {
    it('should generate summary from transcript', async () => {
      const mockFetch = vi.fn().mockResolvedValue({
        ok: true,
        json: () => Promise.resolve({
          summary: 'Meeting discussed project timeline.',
          action_items: [{ text: 'Review budget', assignee: 'John' }],
        }),
      });
      global.fetch = mockFetch;

      const result = await client.summarize('tr_123');

      expect(result.summary).toContain('project timeline');
      expect(result.action_items).toHaveLength(1);
    });
  });
});
```

### Step 3: Integration Test Setup

```typescript
// tests/integration/twinmind/transcription.test.ts
import { describe, it, expect, beforeAll } from 'vitest';
import { TwinMindClient } from '../../../src/twinmind/client';

describe('TwinMind Transcription Integration', () => {
  let client: TwinMindClient;

  beforeAll(() => {
    const apiKey = process.env.TWINMIND_API_KEY;
    if (!apiKey) {
      throw new Error('TWINMIND_API_KEY required for integration tests');
    }
    client = new TwinMindClient({ apiKey });
  });

  describe('API Health', () => {
    it('should return healthy status', async () => {
      const health = await client.healthCheck();
      expect(health).toBe(true);
    });
  });

  describe('Account', () => {
    it('should return account info', async () => {
      const account = await client.getAccount();
      expect(account.id).toBeDefined();
      expect(account.plan).toBeDefined();
    });
  });

  describe('Transcription', () => {
    it('should transcribe test audio', async () => {
      const audioUrl = process.env.TEST_AUDIO_URL;
      if (!audioUrl) {
        console.log('Skipping: TEST_AUDIO_URL not set');
        return;
      }

      const transcript = await client.transcribe(audioUrl, {
        language: 'en',
        diarization: false,
      });

      expect(transcript.id).toBeDefined();
      expect(transcript.text.length).toBeGreaterThan(0);
    }, 60000); // 60 second timeout
  });
});
```

### Step 4: Smoke Test Script

```typescript
// tests/smoke/twinmind-smoke.ts
import { TwinMindClient } from '../../src/twinmind/client';

interface SmokeTestResult {
  name: string;
  passed: boolean;
  duration: number;
  error?: string;
}

async function runSmokeTests(): Promise<SmokeTestResult[]> {
  const results: SmokeTestResult[] = [];
  const client = new TwinMindClient({
    apiKey: process.env.TWINMIND_API_KEY!,
  });

  // Test 1: Health Check
  const healthStart = Date.now();
  try {
    await client.healthCheck();
    results.push({
      name: 'API Health Check',
      passed: true,
      duration: Date.now() - healthStart,
    });
  } catch (error: any) {
    results.push({
      name: 'API Health Check',
      passed: false,
      duration: Date.now() - healthStart,
      error: error.message,
    });
  }

  // Test 2: Account Access
  const accountStart = Date.now();
  try {
    const account = await client.getAccount();
    results.push({
      name: 'Account Access',
      passed: !!account.id,
      duration: Date.now() - accountStart,
    });
  } catch (error: any) {
    results.push({
      name: 'Account Access',
      passed: false,
      duration: Date.now() - accountStart,
      error: error.message,
    });
  }

  // Test 3: Transcription (if test audio available)
  if (process.env.TEST_AUDIO_URL) {
    const transcribeStart = Date.now();
    try {
      const transcript = await client.transcribe(process.env.TEST_AUDIO_URL, {
        model: 'ear-3',
      });
      results.push({
        name: 'Transcription',
        passed: !!transcript.id && transcript.text.length > 0,
        duration: Date.now() - transcribeStart,
      });
    } catch (error: any) {
      results.push({
        name: 'Transcription',
        passed: false,
        duration: Date.now() - transcribeStart,
        error: error.message,
      });
    }
  }

  return results;
}

// Main execution
async function main() {
  console.log('Running TwinMind Smoke Tests...\n');

  const results = await runSmokeTests();

  console.log('Results:');
  console.log('========');

  let allPassed = true;
  for (const result of results) {
    const status = result.passed ? 'PASS' : 'FAIL';
    console.log(`[${status}] ${result.name} (${result.duration}ms)`);
    if (result.error) {
      console.log(`       Error: ${result.error}`);
    }
    if (!result.passed) allPassed = false;
  }

  console.log('\n' + (allPassed ? 'All tests passed!' : 'Some tests failed.'));
  process.exit(allPassed ? 0 : 1);
}

main().catch(console.error);
```

### Step 5: GitLab CI Configuration

```yaml
# .gitlab-ci.yml
stages:
  - lint
  - test
  - integration
  - deploy

variables:
  NODE_VERSION: "20"

.node-template: &node-template
  image: node:${NODE_VERSION}
  cache:
    key: ${CI_COMMIT_REF_SLUG}
    paths:
      - node_modules/

lint:
  <<: *node-template
  stage: lint
  script:
    - npm ci
    - npm run lint
    - npm run typecheck

unit-tests:
  <<: *node-template
  stage: test
  script:
    - npm ci
    - npm run test:unit
  coverage: '/All files[^|]*\|[^|]*\s+([\d\.]+)/'
  artifacts:
    reports:
      coverage_report:
        coverage_format: cobertura
        path: coverage/cobertura-coverage.xml

integration-tests:
  <<: *node-template
  stage: integration
  variables:
    TWINMIND_API_KEY: $TWINMIND_API_KEY_TEST
  script:
    - npm ci
    - npm run test:integration
  rules:
    - if: $CI_COMMIT_BRANCH == "main"
    - if: $CI_PIPELINE_SOURCE == "schedule"

api-health:
  stage: integration
  image: curlimages/curl:latest
  script:
    - |
      curl -f -H "Authorization: Bearer $TWINMIND_API_KEY_TEST" \
        https://api.twinmind.com/v1/health
  rules:
    - if: $CI_COMMIT_BRANCH == "main"
```

### Step 6: Test Configuration

```typescript
// vitest.config.ts
import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    coverage: {
      provider: 'v8',
      reporter: ['text', 'json', 'html', 'cobertura'],
      exclude: [
        'node_modules/',
        'tests/',
        'dist/',
      ],
    },
    testTimeout: 30000,
    hookTimeout: 30000,
  },
});
```

```json
// package.json scripts
{
  "scripts": {
    "test": "vitest",
    "test:unit": "vitest run tests/unit",
    "test:integration": "vitest run tests/integration",
    "test:smoke": "ts-node tests/smoke/twinmind-smoke.ts",
    "test:coverage": "vitest run --coverage"
  }
}
```

## Output
- GitHub Actions workflow
- Unit test suite
- Integration test suite
- Smoke test script
- GitLab CI configuration

## Error Handling

| Issue | Cause | Solution |
|-------|-------|----------|
| Secret not found | Not configured | Add to GitHub Secrets |
| Test timeout | Large audio file | Use short test samples |
| Rate limited | Too many CI runs | Use test account limits |
| API unavailable | TwinMind outage | Add retry logic |

## Best Practices

1. **Use test API keys** - Never use production keys in CI
2. **Short test audio** - Use 10-30 second samples
3. **Cache dependencies** - Speed up pipeline runs
4. **Timeout limits** - Set reasonable timeouts
5. **Parallel tests** - Run independent tests concurrently
6. **Artifact storage** - Save test results and coverage

## Resources
- [GitHub Actions Documentation](https://docs.github.com/actions)
- [GitLab CI Documentation](https://docs.gitlab.com/ee/ci/)
- [Vitest Documentation](https://vitest.dev/)

## Next Steps
For deployment integration, see `twinmind-deploy-integration`.
