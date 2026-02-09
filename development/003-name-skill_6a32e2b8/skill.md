---
name: guidewire-local-dev-loop
description: |
  Configure local Guidewire development workflow with Guidewire Studio, Gosu debugging,
  and hot reload capabilities.
  Use when setting up development environment, configuring IDE, or optimizing dev workflow.
  Trigger with phrases like "guidewire local dev", "guidewire studio setup",
  "gosu development", "guidewire debugging", "configure guidewire ide".
allowed-tools: Read, Write, Edit, Bash(java:*), Bash(gradle:*), Bash(git:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Guidewire Local Dev Loop

## Overview

Set up an efficient local development workflow with Guidewire Studio, including hot reload, Gosu debugging, and continuous testing.

## Prerequisites

- JDK 17 installed and configured
- IntelliJ IDEA (Ultimate recommended) or Guidewire Studio
- Gradle 8.x
- Git for version control
- Access to Guidewire Cloud sandbox environment

## Instructions

### Step 1: Configure IDE Settings

```bash
# IntelliJ IDEA settings for Guidewire
# File > Settings > Build, Execution, Deployment > Build Tools > Gradle

# Set Gradle JVM to JDK 17
# Enable auto-import for Gradle projects
# Set test runner to Gradle
```

**IntelliJ Plugins to Install:**
- Gosu Language Support
- Guidewire Studio Plugin
- EditorConfig

### Step 2: Project Structure Setup

```bash
# Standard Guidewire project structure
project-root/
├── build.gradle                 # Main build configuration
├── settings.gradle              # Multi-project settings
├── gradle.properties            # Gradle properties
├── modules/
│   ├── configuration/           # Gosu configuration code
│   │   ├── gsrc/               # Gosu source files
│   │   └── config/             # XML configuration
│   └── integration/             # Integration code
├── database/
│   └── upgrade/                 # Database upgrade scripts
└── build/
    └── idea/                    # IDE-specific files
```

### Step 3: Configure Local Server

```groovy
// build.gradle - Local server configuration
plugins {
    id 'com.guidewire.gradle' version '10.12.0'
}

guidewire {
    server {
        port = 8080
        debugPort = 5005
        jvmArgs = [
            '-Xmx4g',
            '-XX:+UseG1GC',
            '-Dgw.server.mode=dev'
        ]
    }

    database {
        server = 'localhost'
        port = 5432
        name = 'pc_dev'
        username = 'postgres'
        password = System.getenv('DB_PASSWORD') ?: 'password'
    }
}

// Hot reload configuration
tasks.register('devServer') {
    dependsOn 'classes'
    doLast {
        exec {
            commandLine 'java', '-agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=*:5005',
                '-jar', 'build/libs/server.jar'
        }
    }
}
```

### Step 4: Database Setup

```bash
# Create local development database
createdb -U postgres pc_dev

# Run database upgrade
./gradlew dbUpgrade

# Reset database (caution: destroys data)
./gradlew dbReset

# Generate sample data
./gradlew loadSampleData
```

### Step 5: Start Development Server

```bash
# Start server with hot reload
./gradlew runServer

# Start with debugging enabled
./gradlew runServer --debug-jvm

# Start specific application
./gradlew :policycenter:runServer
```

### Step 6: Configure Hot Reload

```gosu
// Enable Gosu hot swap in development
// config/dev-config.xml
<config>
  <development>
    <hot-swap enabled="true"/>
    <gosu-reload enabled="true"/>
    <pcf-reload enabled="true"/>
  </development>
</config>
```

**IDE Configuration for Hot Reload:**
1. Run > Edit Configurations > Remote JVM Debug
2. Set port to 5005
3. Start debug session after server is running

## Gosu Development Workflow

### Create New Gosu Class

```gosu
// gsrc/gw/custom/MyService.gs
package gw.custom

uses gw.api.database.Query
uses gw.api.util.Logger

class MyService {
  private static final var LOG = Logger.forCategory("MyService")

  static function processPolicy(policyNumber : String) : Policy {
    LOG.info("Processing policy: ${policyNumber}")

    var policy = Query.make(Policy)
      .compare(Policy#PolicyNumber, Equals, policyNumber)
      .select()
      .AtMostOneRow

    if (policy == null) {
      throw new IllegalArgumentException("Policy not found: ${policyNumber}")
    }

    // Business logic here
    return policy
  }
}
```

### Unit Testing

```gosu
// test/gsrc/gw/custom/MyServiceTest.gs
package gw.custom

uses gw.testharness.v3.PLTestCase
uses gw.testharness.v3.PLAssert

class MyServiceTest extends PLTestCase {

  function testProcessPolicy() {
    // Setup test data
    var account = createTestAccount()
    var policy = createTestPolicy(account)

    // Execute
    var result = MyService.processPolicy(policy.PolicyNumber)

    // Assert
    PLAssert.assertNotNull(result)
    PLAssert.assertEquals(policy.PolicyNumber, result.PolicyNumber)
  }

  private function createTestAccount() : Account {
    var account = new Account()
    account.AccountNumber = "TEST-" + System.currentTimeMillis()
    account.Bundle.commit()
    return account
  }
}
```

### Run Tests

```bash
# Run all tests
./gradlew test

# Run specific test class
./gradlew test --tests "gw.custom.MyServiceTest"

# Run with coverage
./gradlew test jacocoTestReport

# Continuous testing (watch mode)
./gradlew test --continuous
```

## PCF Development

```xml
<!-- pcf/AccountDetailScreen.pcf -->
<?xml version="1.0"?>
<PCF
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
  xsi:noNamespaceSchemaLocation="pcf.xsd">
  <Screen
    id="AccountDetailScreen"
    editable="true">
    <Require var="account" type="Account"/>

    <Toolbar>
      <ToolbarButton
        action="saveAccount()"
        id="SaveButton"
        label="Save"/>
    </Toolbar>

    <DetailViewPanel>
      <InputColumn>
        <TextInput
          editable="true"
          id="AccountNumber"
          label="Account Number"
          value="account.AccountNumber"/>
        <TextInput
          editable="true"
          id="AccountName"
          label="Account Name"
          value="account.AccountHolderContact.DisplayName"/>
      </InputColumn>
    </DetailViewPanel>

    <Code>
      <![CDATA[
      function saveAccount() {
        account.Bundle.commit()
        util.LocationUtil.addRequestScopedInfoMessage("Account saved")
      }
      ]]>
    </Code>
  </Screen>
</PCF>
```

## Output

- Running local development server
- Hot reload enabled for Gosu and PCF changes
- Debug session attached to IDE
- Unit tests executing successfully

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| Port already in use | Server already running | Kill existing process or change port |
| Database connection failed | Wrong credentials/missing DB | Check postgres is running |
| Gosu compilation error | Syntax error | Check IDE error highlights |
| Hot reload not working | Debug not attached | Reconnect debugger |
| Out of memory | Insufficient heap | Increase -Xmx in jvmArgs |

## Development Commands Cheatsheet

```bash
# Build
./gradlew clean build

# Run server
./gradlew runServer

# Database operations
./gradlew dbUpgrade
./gradlew dbReset
./gradlew loadSampleData

# Testing
./gradlew test
./gradlew test --tests "ClassName.methodName"
./gradlew test --continuous

# Code quality
./gradlew gosucheck
./gradlew spotlessApply

# Generate API documentation
./gradlew apiDoc
```

## IDE Keyboard Shortcuts

| Action | IntelliJ Shortcut |
|--------|-------------------|
| Hot swap code | Ctrl+Shift+F9 |
| Run to cursor | Alt+F9 |
| Evaluate expression | Alt+F8 |
| Find usages | Alt+F7 |
| Go to declaration | Ctrl+B |
| Refactor rename | Shift+F6 |

## Resources

- [Guidewire Studio Documentation](https://docs.guidewire.com/cloud/gcc-guide/insurer-developer/)
- [Gosu Language Reference](https://gosu-lang.github.io/)
- [Gradle Plugin Documentation](https://docs.guidewire.com/tools/gradle-plugin/)

## Next Steps

For SDK and API patterns, see `guidewire-sdk-patterns`.
