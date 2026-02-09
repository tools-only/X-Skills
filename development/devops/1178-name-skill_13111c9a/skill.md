---
name: databricks-ci-integration
description: |
  Configure Databricks CI/CD integration with GitHub Actions and Asset Bundles.
  Use when setting up automated testing, configuring CI pipelines,
  or integrating Databricks deployments into your build process.
  Trigger with phrases like "databricks CI", "databricks GitHub Actions",
  "databricks automated tests", "CI databricks", "databricks pipeline".
allowed-tools: Read, Write, Edit, Bash(gh:*), Bash(databricks:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Databricks CI Integration

## Overview
Set up CI/CD pipelines for Databricks using GitHub Actions and Asset Bundles.

## Prerequisites
- GitHub repository with Actions enabled
- Databricks workspace with service principal
- Asset Bundles project structure

## Instructions

### Step 1: Configure Service Principal

```bash
# Create service principal in Databricks
databricks service-principals create --json '{
  "display_name": "GitHub Actions CI",
  "active": true
}'

# Note the application_id returned

# Create OAuth secret
databricks service-principal-secrets create \
  --service-principal-id <application_id>

# Grant permissions to service principal
databricks permissions update workspace --json '{
  "access_control_list": [{
    "service_principal_name": "<application_id>",
    "permission_level": "CAN_MANAGE"
  }]
}'
```

### Step 2: Configure GitHub Secrets

```bash
# Set GitHub secrets
gh secret set DATABRICKS_HOST --body "https://adb-1234567890.1.azuredatabricks.net"
gh secret set DATABRICKS_CLIENT_ID --body "your-client-id"
gh secret set DATABRICKS_CLIENT_SECRET --body "your-client-secret"

# For staging/prod environments
gh secret set DATABRICKS_HOST_STAGING --body "https://staging.azuredatabricks.net"
gh secret set DATABRICKS_HOST_PROD --body "https://prod.azuredatabricks.net"
```

### Step 3: Create GitHub Actions Workflow

```yaml
# .github/workflows/databricks-ci.yml
name: Databricks CI/CD

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

env:
  DATABRICKS_HOST: ${{ secrets.DATABRICKS_HOST }}
  DATABRICKS_CLIENT_ID: ${{ secrets.DATABRICKS_CLIENT_ID }}
  DATABRICKS_CLIENT_SECRET: ${{ secrets.DATABRICKS_CLIENT_SECRET }}

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Setup Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'
          cache: 'pip'

      - name: Install dependencies
        run: |
          pip install databricks-cli databricks-sdk pytest

      - name: Validate Asset Bundle
        run: databricks bundle validate

      - name: Run unit tests
        run: pytest tests/unit/ -v --tb=short

  deploy-staging:
    needs: validate
    if: github.ref == 'refs/heads/develop'
    runs-on: ubuntu-latest
    environment: staging
    env:
      DATABRICKS_HOST: ${{ secrets.DATABRICKS_HOST_STAGING }}
    steps:
      - uses: actions/checkout@v4

      - name: Setup Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'

      - name: Install Databricks CLI
        run: pip install databricks-cli

      - name: Deploy to Staging
        run: |
          databricks bundle deploy -t staging

      - name: Run Integration Tests
        run: |
          # Trigger test job and wait for completion
          RUN_ID=$(databricks bundle run -t staging integration-tests | jq -r '.run_id')
          databricks runs get --run-id $RUN_ID --wait
          # Check result
          RESULT=$(databricks runs get --run-id $RUN_ID | jq -r '.state.result_state')
          if [ "$RESULT" != "SUCCESS" ]; then
            echo "Integration tests failed!"
            exit 1
          fi

  deploy-production:
    needs: [validate, deploy-staging]
    if: github.ref == 'refs/heads/main'
    runs-on: ubuntu-latest
    environment:
      name: production
      url: ${{ secrets.DATABRICKS_HOST_PROD }}
    env:
      DATABRICKS_HOST: ${{ secrets.DATABRICKS_HOST_PROD }}
    steps:
      - uses: actions/checkout@v4

      - name: Setup Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'

      - name: Install Databricks CLI
        run: pip install databricks-cli

      - name: Deploy to Production
        run: |
          databricks bundle deploy -t prod

      - name: Verify Deployment
        run: |
          databricks bundle summary -t prod
          # Trigger smoke test
          databricks bundle run -t prod smoke-test
```

### Step 4: PR Validation Workflow

```yaml
# .github/workflows/pr-validation.yml
name: PR Validation

on:
  pull_request:
    branches: [main, develop]

jobs:
  lint-and-test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Setup Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'

      - name: Install dependencies
        run: |
          pip install ruff mypy pytest pytest-cov databricks-sdk

      - name: Lint with ruff
        run: ruff check src/

      - name: Type check with mypy
        run: mypy src/ --ignore-missing-imports

      - name: Run tests with coverage
        run: pytest tests/unit/ --cov=src --cov-report=xml

      - name: Upload coverage
        uses: codecov/codecov-action@v4
        with:
          files: coverage.xml

  bundle-validation:
    runs-on: ubuntu-latest
    env:
      DATABRICKS_HOST: ${{ secrets.DATABRICKS_HOST }}
      DATABRICKS_CLIENT_ID: ${{ secrets.DATABRICKS_CLIENT_ID }}
      DATABRICKS_CLIENT_SECRET: ${{ secrets.DATABRICKS_CLIENT_SECRET }}
    steps:
      - uses: actions/checkout@v4

      - name: Install Databricks CLI
        run: pip install databricks-cli

      - name: Validate bundle for all targets
        run: |
          databricks bundle validate -t dev
          databricks bundle validate -t staging
          databricks bundle validate -t prod

      - name: Check for breaking changes
        run: |
          # Compare job configurations
          databricks bundle summary -t prod --output json > current.json
          # Add logic to detect breaking changes
```

### Step 5: Nightly Test Workflow

```yaml
# .github/workflows/nightly-tests.yml
name: Nightly Tests

on:
  schedule:
    - cron: '0 2 * * *'  # 2 AM UTC daily
  workflow_dispatch:

jobs:
  integration-tests:
    runs-on: ubuntu-latest
    env:
      DATABRICKS_HOST: ${{ secrets.DATABRICKS_HOST_STAGING }}
      DATABRICKS_CLIENT_ID: ${{ secrets.DATABRICKS_CLIENT_ID }}
      DATABRICKS_CLIENT_SECRET: ${{ secrets.DATABRICKS_CLIENT_SECRET }}
    steps:
      - uses: actions/checkout@v4

      - name: Install dependencies
        run: pip install databricks-cli

      - name: Run full integration test suite
        run: |
          databricks bundle deploy -t staging
          RUN_ID=$(databricks bundle run -t staging full-integration-tests | jq -r '.run_id')
          databricks runs get --run-id $RUN_ID --wait

      - name: Generate test report
        if: always()
        run: |
          # Download test results
          databricks fs cp dbfs:/test-results/latest/ ./test-results/ --recursive

      - name: Upload test artifacts
        uses: actions/upload-artifact@v4
        if: always()
        with:
          name: test-results
          path: test-results/

      - name: Notify on failure
        if: failure()
        uses: slackapi/slack-github-action@v1
        with:
          channel-id: 'data-engineering-alerts'
          slack-message: 'Nightly tests failed! Check ${{ github.server_url }}/${{ github.repository }}/actions/runs/${{ github.run_id }}'
        env:
          SLACK_BOT_TOKEN: ${{ secrets.SLACK_BOT_TOKEN }}
```

## Output
- Automated test pipeline
- PR checks configured
- Staging deployment on merge to develop
- Production deployment on merge to main

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Auth failed | Invalid credentials | Regenerate service principal secret |
| Bundle validation failed | Invalid YAML | Run `databricks bundle validate` locally |
| Deployment timeout | Slow cluster startup | Use warm pools or increase timeout |
| Tests failed | Code regression | Fix code and re-run |

## Examples

### Matrix Testing (Multiple DBR Versions)
```yaml
jobs:
  test-matrix:
    strategy:
      matrix:
        dbr_version: ['13.3', '14.3', '15.1']
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Test on DBR ${{ matrix.dbr_version }}
        run: |
          databricks bundle deploy -t test-${{ matrix.dbr_version }}
          databricks bundle run -t test-${{ matrix.dbr_version }} tests
```

### Branch Protection Rules
```yaml
# Set via GitHub API or UI
required_status_checks:
  - "lint-and-test"
  - "bundle-validation"
required_reviews: 1
dismiss_stale_reviews: true
```

## Resources
- [Databricks Asset Bundles](https://docs.databricks.com/dev-tools/bundles/index.html)
- [GitHub Actions](https://docs.github.com/en/actions)
- [Service Principal Auth](https://docs.databricks.com/dev-tools/auth.html#oauth-machine-to-machine-m2m)

## Next Steps
For deployment patterns, see `databricks-deploy-integration`.
