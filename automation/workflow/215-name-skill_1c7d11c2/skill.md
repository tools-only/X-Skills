---
name: openevidence-ci-integration
description: |
  Integrate OpenEvidence testing into CI/CD pipelines.
  Use when setting up automated testing, configuring GitHub Actions,
  or implementing continuous integration for clinical AI applications.
  Trigger with phrases like "openevidence ci", "openevidence github actions",
  "openevidence pipeline", "test openevidence ci", "automate openevidence tests".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence CI Integration

## Overview
Integrate OpenEvidence testing and validation into CI/CD pipelines for healthcare applications.

## Prerequisites
- GitHub, GitLab, or other CI platform
- OpenEvidence sandbox API credentials
- Test suite configured
- Secret management in CI platform

## Instructions

### Step 1: GitHub Actions Workflow
```yaml
# .github/workflows/openevidence-ci.yml
name: OpenEvidence CI

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

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

      - name: Lint
        run: npm run lint

      - name: Type check
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
        run: npm run test:unit -- --coverage

      - name: Upload coverage
        uses: codecov/codecov-action@v4
        with:
          files: ./coverage/lcov.info
          fail_ci_if_error: false

  integration-tests:
    runs-on: ubuntu-latest
    # Only run on main branch or explicit trigger
    if: github.ref == 'refs/heads/main' || contains(github.event.pull_request.labels.*.name, 'run-integration')
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
        env:
          OPENEVIDENCE_API_KEY: ${{ secrets.OPENEVIDENCE_SANDBOX_API_KEY }}
          OPENEVIDENCE_ORG_ID: ${{ secrets.OPENEVIDENCE_SANDBOX_ORG_ID }}
          OPENEVIDENCE_BASE_URL: https://api.sandbox.openevidence.com
        run: npm run test:integration
        timeout-minutes: 10

  clinical-validation:
    runs-on: ubuntu-latest
    needs: [unit-tests]
    if: github.ref == 'refs/heads/main'
    steps:
      - uses: actions/checkout@v4

      - name: Setup Node.js
        uses: actions/setup-node@v4
        with:
          node-version: ${{ env.NODE_VERSION }}
          cache: 'npm'

      - name: Install dependencies
        run: npm ci

      - name: Run clinical validation tests
        env:
          OPENEVIDENCE_API_KEY: ${{ secrets.OPENEVIDENCE_SANDBOX_API_KEY }}
          OPENEVIDENCE_ORG_ID: ${{ secrets.OPENEVIDENCE_SANDBOX_ORG_ID }}
        run: npm run test:clinical-validation

      - name: Upload validation report
        uses: actions/upload-artifact@v4
        with:
          name: clinical-validation-report
          path: reports/clinical-validation.json
```

### Step 2: Test Configuration
```typescript
// vitest.config.ts
import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['tests/**/*.test.ts'],
    exclude: ['tests/integration/**', 'tests/clinical-validation/**'],
    coverage: {
      provider: 'v8',
      reporter: ['text', 'lcov'],
      exclude: ['tests/**', 'node_modules/**'],
    },
    testTimeout: 10000,
  },
});

// vitest.config.integration.ts
import { defineConfig, mergeConfig } from 'vitest/config';
import baseConfig from './vitest.config';

export default mergeConfig(baseConfig, defineConfig({
  test: {
    include: ['tests/integration/**/*.test.ts'],
    exclude: [],
    testTimeout: 60000, // 60s for API calls
    retry: 2, // Retry flaky integration tests
    maxConcurrency: 1, // Avoid rate limits
  },
}));
```

### Step 3: Unit Tests with Mocks
```typescript
// tests/unit/clinical-query.test.ts
import { describe, it, expect, vi, beforeEach } from 'vitest';
import { queryClinicalEvidence } from '../../src/services/clinical-query';

// Mock the OpenEvidence SDK
vi.mock('@openevidence/sdk', () => ({
  OpenEvidenceClient: vi.fn().mockImplementation(() => ({
    query: vi.fn().mockResolvedValue({
      id: 'test-query-123',
      answer: 'Mock clinical answer for testing purposes.',
      citations: [
        { source: 'Test Journal', title: 'Test Article', year: 2025 }
      ],
      confidence: 0.95,
      lastUpdated: '2025-01-01',
    }),
  })),
}));

describe('Clinical Query Service', () => {
  beforeEach(() => {
    vi.clearAllMocks();
  });

  it('should return formatted clinical response', async () => {
    const result = await queryClinicalEvidence('Test question');

    expect(result.answer).toBeDefined();
    expect(result.citations).toHaveLength(1);
    expect(result.confidence).toBeGreaterThan(0);
  });

  it('should handle query errors gracefully', async () => {
    const { OpenEvidenceClient } = await import('@openevidence/sdk');
    vi.mocked(OpenEvidenceClient).mockImplementation(() => ({
      query: vi.fn().mockRejectedValue(new Error('API Error')),
    } as any));

    await expect(queryClinicalEvidence('Test question'))
      .rejects.toThrow('API Error');
  });
});
```

### Step 4: Integration Tests
```typescript
// tests/integration/openevidence-api.test.ts
import { describe, it, expect, beforeAll } from 'vitest';
import { OpenEvidenceClient } from '@openevidence/sdk';

describe('OpenEvidence API Integration', () => {
  let client: OpenEvidenceClient;

  beforeAll(() => {
    // Ensure we're using sandbox credentials
    const baseUrl = process.env.OPENEVIDENCE_BASE_URL;
    if (!baseUrl?.includes('sandbox')) {
      throw new Error('Integration tests must use sandbox environment');
    }

    client = new OpenEvidenceClient({
      apiKey: process.env.OPENEVIDENCE_API_KEY!,
      orgId: process.env.OPENEVIDENCE_ORG_ID!,
      baseUrl,
    });
  });

  it('should successfully query clinical evidence', async () => {
    const response = await client.query({
      question: 'What is the half-life of aspirin?',
      context: {
        specialty: 'pharmacology',
        urgency: 'routine',
      },
    });

    expect(response.id).toBeDefined();
    expect(response.answer).toBeDefined();
    expect(response.answer.length).toBeGreaterThan(50);
    expect(response.citations.length).toBeGreaterThan(0);
    expect(response.confidence).toBeGreaterThan(0.5);
  }, 30000);

  it('should handle invalid queries appropriately', async () => {
    await expect(client.query({
      question: '',
      context: { specialty: 'internal-medicine', urgency: 'routine' },
    })).rejects.toThrow();
  });

  it('should return citations with required fields', async () => {
    const response = await client.query({
      question: 'What are the first-line treatments for hypertension?',
      context: { specialty: 'cardiology', urgency: 'routine' },
    });

    response.citations.forEach(citation => {
      expect(citation.source).toBeDefined();
      expect(citation.title).toBeDefined();
      expect(citation.year).toBeGreaterThan(2000);
    });
  }, 30000);
});
```

### Step 5: Clinical Validation Tests
```typescript
// tests/clinical-validation/known-answers.test.ts
import { describe, it, expect } from 'vitest';
import { OpenEvidenceClient } from '@openevidence/sdk';

interface ClinicalValidationCase {
  question: string;
  specialty: string;
  expectedKeywords: string[];
  mustNotContain?: string[];
}

const VALIDATION_CASES: ClinicalValidationCase[] = [
  {
    question: 'What is the first-line treatment for type 2 diabetes in adults?',
    specialty: 'endocrinology',
    expectedKeywords: ['metformin', 'lifestyle', 'HbA1c'],
    mustNotContain: ['insulin-first'],
  },
  {
    question: 'What are the contraindications for aspirin?',
    specialty: 'internal-medicine',
    expectedKeywords: ['bleeding', 'allergy', 'ulcer'],
  },
  {
    question: 'What is the target LDL for secondary prevention of cardiovascular disease?',
    specialty: 'cardiology',
    expectedKeywords: ['LDL', 'mg/dL', 'statin'],
  },
];

describe('Clinical Validation - Known Answers', () => {
  const client = new OpenEvidenceClient({
    apiKey: process.env.OPENEVIDENCE_API_KEY!,
    orgId: process.env.OPENEVIDENCE_ORG_ID!,
  });

  const results: any[] = [];

  VALIDATION_CASES.forEach((testCase, index) => {
    it(`should answer correctly: ${testCase.question.substring(0, 50)}...`, async () => {
      const response = await client.query({
        question: testCase.question,
        context: {
          specialty: testCase.specialty,
          urgency: 'routine',
        },
      });

      const answerLower = response.answer.toLowerCase();

      // Check expected keywords
      const foundKeywords = testCase.expectedKeywords.filter(
        kw => answerLower.includes(kw.toLowerCase())
      );

      // Check must-not-contain
      const forbiddenFound = testCase.mustNotContain?.filter(
        kw => answerLower.includes(kw.toLowerCase())
      ) || [];

      // Record results for report
      results.push({
        question: testCase.question,
        expectedKeywords: testCase.expectedKeywords,
        foundKeywords,
        forbiddenFound,
        confidence: response.confidence,
        citationCount: response.citations.length,
        passed: foundKeywords.length >= testCase.expectedKeywords.length / 2 &&
                forbiddenFound.length === 0,
      });

      // Assert
      expect(foundKeywords.length).toBeGreaterThanOrEqual(
        Math.ceil(testCase.expectedKeywords.length / 2),
        `Missing keywords: ${testCase.expectedKeywords.filter(k => !foundKeywords.includes(k)).join(', ')}`
      );
      expect(forbiddenFound).toHaveLength(0);
      expect(response.confidence).toBeGreaterThan(0.7);
    }, 60000);
  });

  afterAll(async () => {
    // Write validation report
    const fs = await import('fs/promises');
    await fs.mkdir('reports', { recursive: true });
    await fs.writeFile(
      'reports/clinical-validation.json',
      JSON.stringify({
        timestamp: new Date().toISOString(),
        totalTests: results.length,
        passed: results.filter(r => r.passed).length,
        failed: results.filter(r => !r.passed).length,
        results,
      }, null, 2)
    );
  });
});
```

### Step 6: Secret Management
```yaml
# Store these as GitHub Secrets:
# - OPENEVIDENCE_SANDBOX_API_KEY
# - OPENEVIDENCE_SANDBOX_ORG_ID
# - OPENEVIDENCE_PROD_API_KEY (for deployment only)
# - OPENEVIDENCE_PROD_ORG_ID (for deployment only)

# For local development, use .env.test:
# OPENEVIDENCE_API_KEY=oe_sandbox_***
# OPENEVIDENCE_ORG_ID=org_sandbox_***
# OPENEVIDENCE_BASE_URL=https://api.sandbox.openevidence.com
```

## Output
- CI workflow running on push and PR
- Unit tests with mocked SDK
- Integration tests against sandbox
- Clinical validation test suite
- Coverage reports

## Package.json Scripts
```json
{
  "scripts": {
    "test": "vitest",
    "test:unit": "vitest run -c vitest.config.ts",
    "test:integration": "vitest run -c vitest.config.integration.ts",
    "test:clinical-validation": "vitest run tests/clinical-validation",
    "test:ci": "vitest run --coverage",
    "typecheck": "tsc --noEmit",
    "lint": "eslint src tests"
  }
}
```

## Error Handling
| CI Issue | Cause | Solution |
|----------|-------|----------|
| Integration test timeout | API slow or rate limited | Increase timeout, add retry |
| Secret not found | Missing GitHub secret | Add secret in repo settings |
| Flaky tests | Network variability | Add retries, improve assertions |
| Coverage drop | New code untested | Add tests, adjust thresholds |

## Resources
- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Vitest Documentation](https://vitest.dev/)

## Next Steps
For deployment integration, see `openevidence-deploy-integration`.
