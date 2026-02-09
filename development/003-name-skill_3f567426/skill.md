---
name: guidewire-ci-integration
description: |
  Configure CI/CD pipelines for Guidewire InsuranceSuite development.
  Use when setting up automated builds, tests, and deployments for Guidewire projects.
  Trigger with phrases like "guidewire ci", "guidewire pipeline",
  "automated testing guidewire", "jenkins guidewire", "github actions guidewire".
allowed-tools: Read, Write, Edit, Bash(gradle:*), Bash(git:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Guidewire CI Integration

## Overview

Configure continuous integration and continuous deployment pipelines for Guidewire InsuranceSuite projects using GitHub Actions, Jenkins, or Azure DevOps.

## Prerequisites

- Git repository with Guidewire project
- CI/CD platform access (GitHub Actions, Jenkins, or Azure DevOps)
- Guidewire Cloud Console access for deployment credentials
- JDK 17 and Gradle available in CI environment

## Instructions

### Step 1: GitHub Actions Workflow

```yaml
# .github/workflows/guidewire-ci.yml
name: Guidewire CI/CD

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

env:
  JAVA_VERSION: '17'
  GRADLE_VERSION: '8.5'
  GW_TENANT_ID: ${{ secrets.GW_TENANT_ID }}

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Set up JDK
        uses: actions/setup-java@v4
        with:
          java-version: ${{ env.JAVA_VERSION }}
          distribution: 'temurin'
          cache: 'gradle'

      - name: Grant Gradle execute permission
        run: chmod +x gradlew

      - name: Build
        run: ./gradlew build --no-daemon

      - name: Run Gosu checks
        run: ./gradlew gosucheck --no-daemon

      - name: Upload build artifacts
        uses: actions/upload-artifact@v4
        with:
          name: build-artifacts
          path: build/libs/

  test:
    needs: build
    runs-on: ubuntu-latest
    services:
      postgres:
        image: postgres:15
        env:
          POSTGRES_DB: pc_test
          POSTGRES_USER: postgres
          POSTGRES_PASSWORD: postgres
        ports:
          - 5432:5432
        options: >-
          --health-cmd pg_isready
          --health-interval 10s
          --health-timeout 5s
          --health-retries 5

    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Set up JDK
        uses: actions/setup-java@v4
        with:
          java-version: ${{ env.JAVA_VERSION }}
          distribution: 'temurin'
          cache: 'gradle'

      - name: Download build artifacts
        uses: actions/download-artifact@v4
        with:
          name: build-artifacts
          path: build/libs/

      - name: Set up test database
        run: |
          ./gradlew dbupgrade -PdbHost=localhost -PdbName=pc_test \
            -PdbUser=postgres -PdbPassword=postgres --no-daemon

      - name: Run unit tests
        run: ./gradlew test --no-daemon

      - name: Run integration tests
        run: ./gradlew integrationTest --no-daemon
        env:
          DB_HOST: localhost
          DB_NAME: pc_test
          DB_USER: postgres
          DB_PASSWORD: postgres

      - name: Upload test results
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: test-results
          path: build/reports/tests/

      - name: Upload coverage report
        uses: actions/upload-artifact@v4
        with:
          name: coverage-report
          path: build/reports/jacoco/

  security-scan:
    needs: build
    runs-on: ubuntu-latest
    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Run SAST scan
        uses: github/codeql-action/analyze@v2
        with:
          languages: java

      - name: Dependency vulnerability check
        run: ./gradlew dependencyCheckAnalyze --no-daemon

      - name: Upload security report
        uses: actions/upload-artifact@v4
        with:
          name: security-report
          path: build/reports/dependency-check/

  deploy-sandbox:
    needs: [test, security-scan]
    runs-on: ubuntu-latest
    if: github.ref == 'refs/heads/develop'
    environment: sandbox
    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Download build artifacts
        uses: actions/download-artifact@v4
        with:
          name: build-artifacts
          path: build/libs/

      - name: Deploy to Guidewire Cloud Sandbox
        env:
          GW_CLIENT_ID: ${{ secrets.GW_SANDBOX_CLIENT_ID }}
          GW_CLIENT_SECRET: ${{ secrets.GW_SANDBOX_CLIENT_SECRET }}
        run: |
          ./scripts/deploy-to-cloud.sh sandbox

      - name: Run smoke tests
        run: ./gradlew smokeTest -Penv=sandbox --no-daemon

  deploy-production:
    needs: deploy-sandbox
    runs-on: ubuntu-latest
    if: github.ref == 'refs/heads/main'
    environment:
      name: production
      url: https://your-tenant.cloud.guidewire.com
    steps:
      - name: Checkout code
        uses: actions/checkout@v4

      - name: Download build artifacts
        uses: actions/download-artifact@v4
        with:
          name: build-artifacts
          path: build/libs/

      - name: Deploy to Guidewire Cloud Production
        env:
          GW_CLIENT_ID: ${{ secrets.GW_PROD_CLIENT_ID }}
          GW_CLIENT_SECRET: ${{ secrets.GW_PROD_CLIENT_SECRET }}
        run: |
          ./scripts/deploy-to-cloud.sh production

      - name: Run production verification
        run: ./gradlew productionVerify -Penv=production --no-daemon
```

### Step 2: Jenkins Pipeline

```groovy
// Jenkinsfile
pipeline {
    agent {
        docker {
            image 'eclipse-temurin:17-jdk'
            args '-v gradle-cache:/root/.gradle'
        }
    }

    environment {
        GRADLE_OPTS = '-Dorg.gradle.daemon=false'
        GW_TENANT_ID = credentials('gw-tenant-id')
    }

    options {
        timeout(time: 1, unit: 'HOURS')
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '10'))
    }

    stages {
        stage('Checkout') {
            steps {
                checkout scm
            }
        }

        stage('Build') {
            steps {
                sh 'chmod +x gradlew'
                sh './gradlew clean build'
            }
            post {
                success {
                    archiveArtifacts artifacts: 'build/libs/*.jar', fingerprint: true
                }
            }
        }

        stage('Code Quality') {
            parallel {
                stage('Gosu Check') {
                    steps {
                        sh './gradlew gosucheck'
                    }
                }
                stage('Static Analysis') {
                    steps {
                        sh './gradlew spotbugsMain'
                    }
                    post {
                        always {
                            recordIssues(tools: [spotBugs(pattern: '**/build/reports/spotbugs/*.xml')])
                        }
                    }
                }
            }
        }

        stage('Test') {
            steps {
                sh './gradlew test integrationTest'
            }
            post {
                always {
                    junit '**/build/test-results/**/*.xml'
                    publishHTML([
                        reportDir: 'build/reports/tests/test',
                        reportFiles: 'index.html',
                        reportName: 'Test Report'
                    ])
                }
            }
        }

        stage('Security Scan') {
            steps {
                sh './gradlew dependencyCheckAnalyze'
            }
            post {
                always {
                    dependencyCheckPublisher pattern: '**/dependency-check-report.xml'
                }
            }
        }

        stage('Deploy Sandbox') {
            when {
                branch 'develop'
            }
            environment {
                GW_CLIENT_ID = credentials('gw-sandbox-client-id')
                GW_CLIENT_SECRET = credentials('gw-sandbox-client-secret')
            }
            steps {
                sh './scripts/deploy-to-cloud.sh sandbox'
                sh './gradlew smokeTest -Penv=sandbox'
            }
        }

        stage('Deploy Production') {
            when {
                branch 'main'
            }
            environment {
                GW_CLIENT_ID = credentials('gw-prod-client-id')
                GW_CLIENT_SECRET = credentials('gw-prod-client-secret')
            }
            steps {
                input message: 'Deploy to production?', ok: 'Deploy'
                sh './scripts/deploy-to-cloud.sh production'
            }
        }
    }

    post {
        always {
            cleanWs()
        }
        failure {
            emailext(
                subject: "Build Failed: ${env.JOB_NAME} [${env.BUILD_NUMBER}]",
                body: "Check console output at ${env.BUILD_URL}",
                recipientProviders: [developers()]
            )
        }
    }
}
```

### Step 3: Gradle Build Configuration

```groovy
// build.gradle - CI-optimized configuration
plugins {
    id 'com.guidewire.gradle' version '10.12.0'
    id 'org.owasp.dependencycheck' version '9.0.0'
    id 'com.github.spotbugs' version '6.0.0'
    id 'jacoco'
}

java {
    toolchain {
        languageVersion = JavaLanguageVersion.of(17)
    }
}

// Gosu compilation checks
tasks.named('gosucheck') {
    reports {
        xml.required = true
        html.required = true
    }
}

// Test configuration
test {
    useJUnitPlatform()
    maxParallelForks = Runtime.runtime.availableProcessors().intdiv(2) ?: 1

    testLogging {
        events 'passed', 'skipped', 'failed'
        showStandardStreams = false
        exceptionFormat = 'full'
    }

    finalizedBy jacocoTestReport
}

// Integration tests
tasks.register('integrationTest', Test) {
    description = 'Runs integration tests'
    group = 'verification'

    testClassesDirs = sourceSets.integrationTest.output.classesDirs
    classpath = sourceSets.integrationTest.runtimeClasspath

    shouldRunAfter test
}

// Code coverage
jacocoTestReport {
    dependsOn test
    reports {
        xml.required = true
        html.required = true
    }
}

jacocoTestCoverageVerification {
    violationRules {
        rule {
            limit {
                minimum = 0.7 // 70% coverage minimum
            }
        }
    }
}

// Security scanning
dependencyCheck {
    formats = ['HTML', 'XML', 'JSON']
    failBuildOnCVSS = 7.0
    suppressionFile = 'config/dependency-check-suppressions.xml'
}

// SpotBugs configuration
spotbugs {
    effort = 'max'
    reportLevel = 'medium'
    excludeFilter = file('config/spotbugs-exclude.xml')
}

// CI-specific tasks
tasks.register('ci') {
    description = 'Runs all CI checks'
    group = 'verification'
    dependsOn 'build', 'gosucheck', 'test', 'integrationTest',
              'jacocoTestCoverageVerification', 'dependencyCheckAnalyze'
}
```

### Step 4: Test Utilities

```gosu
// Test utilities for CI
package gw.test.ci

uses gw.testharness.v3.PLTestCase
uses gw.api.database.Query

class CITestBase extends PLTestCase {

  // Skip slow tests in CI
  static property get SkipSlowTests() : boolean {
    return System.getenv("CI") == "true" &&
           System.getenv("RUN_SLOW_TESTS") != "true"
  }

  // Database health check
  static function verifyDatabaseConnection() : boolean {
    try {
      var count = Query.make(Account).select().Count
      return true
    } catch (e : Exception) {
      return false
    }
  }

  // Clean test data
  override function beforeClass() {
    super.beforeClass()
    cleanupTestData()
  }

  protected function cleanupTestData() {
    // Remove test accounts created by previous runs
    Query.make(Account)
      .compare(Account#AccountNumber, StartsWith, "TEST-")
      .select()
      .each(\account -> {
        gw.transaction.Transaction.runWithNewBundle(\bundle -> {
          bundle.delete(bundle.add(account))
        })
      })
  }
}
```

### Step 5: Deployment Script

```bash
#!/bin/bash
# scripts/deploy-to-cloud.sh

set -e

ENV=${1:-sandbox}

echo "=== Deploying to Guidewire Cloud ($ENV) ==="

# Validate environment
if [[ "$ENV" != "sandbox" && "$ENV" != "production" ]]; then
    echo "Invalid environment: $ENV"
    exit 1
fi

# Get access token
TOKEN=$(curl -s -X POST "${GW_HUB_URL}/oauth/token" \
    -H "Content-Type: application/x-www-form-urlencoded" \
    -d "grant_type=client_credentials&client_id=${GW_CLIENT_ID}&client_secret=${GW_CLIENT_SECRET}" \
    | jq -r '.access_token')

if [ -z "$TOKEN" ] || [ "$TOKEN" == "null" ]; then
    echo "Failed to obtain access token"
    exit 1
fi

# Deploy configuration package
echo "Deploying configuration package..."
DEPLOYMENT_RESPONSE=$(curl -s -X POST "${GW_API_URL}/deployment/v1/packages" \
    -H "Authorization: Bearer ${TOKEN}" \
    -H "Content-Type: application/zip" \
    --data-binary @build/libs/configuration-package.zip)

DEPLOYMENT_ID=$(echo $DEPLOYMENT_RESPONSE | jq -r '.deploymentId')
echo "Deployment ID: $DEPLOYMENT_ID"

# Wait for deployment to complete
echo "Waiting for deployment to complete..."
while true; do
    STATUS=$(curl -s "${GW_API_URL}/deployment/v1/packages/${DEPLOYMENT_ID}" \
        -H "Authorization: Bearer ${TOKEN}" \
        | jq -r '.status')

    echo "Status: $STATUS"

    if [ "$STATUS" == "COMPLETED" ]; then
        echo "Deployment completed successfully"
        break
    elif [ "$STATUS" == "FAILED" ]; then
        echo "Deployment failed"
        exit 1
    fi

    sleep 30
done

echo "=== Deployment Complete ==="
```

## Output

- GitHub Actions workflow file
- Jenkins pipeline configuration
- Gradle build with CI tasks
- Test utilities for CI environment
- Cloud deployment script

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| Build timeout | Long-running tests | Increase timeout or parallelize |
| Database connection failed | Service not ready | Add health check wait |
| Deployment failed | Invalid package | Check build artifacts |
| Test flakiness | Race conditions | Add proper test isolation |

## Resources

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Jenkins Pipeline Documentation](https://www.jenkins.io/doc/book/pipeline/)
- [Guidewire Gradle Plugin](https://docs.guidewire.com/tools/gradle-plugin/)

## Next Steps

For deployment strategies, see `guidewire-deploy-integration`.
