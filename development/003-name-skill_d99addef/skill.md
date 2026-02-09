---
name: langfuse-upgrade-migration
description: |
  Upgrade Langfuse SDK versions and migrate between API changes.
  Use when upgrading Langfuse SDK, handling breaking changes,
  or migrating between Langfuse versions.
  Trigger with phrases like "upgrade langfuse", "langfuse migration",
  "update langfuse SDK", "langfuse breaking changes", "langfuse version".
allowed-tools: Read, Write, Edit, Bash(npm:*), Bash(pip:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Upgrade & Migration

## Overview
Guide for upgrading Langfuse SDK versions and handling breaking changes.

## Prerequisites
- Existing Langfuse integration
- Access to test environment
- Version control (git)

## Instructions

### Step 1: Check Current Version and Updates

```bash
# Check current version (Node.js)
npm list langfuse
npm outdated langfuse

# Check current version (Python)
pip show langfuse
pip index versions langfuse

# View changelog
open https://github.com/langfuse/langfuse-js/releases
open https://github.com/langfuse/langfuse-python/releases
```

### Step 2: Review Breaking Changes

```markdown
## Common Breaking Changes by Version

### v2.x -> v3.x (TypeScript)
- `Langfuse.trace()` returns `Trace` instead of `Promise<Trace>`
- `flushAsync()` replaces `flush()` (now async)
- `observeOpenAI()` moved to main package export
- Generation `completionTokens` -> `completionTokens` (was `completion_tokens`)

### v1.x -> v2.x (Python)
- `langfuse.trace()` now returns synchronously
- Decorator `@observe()` replaces `@langfuse.observe()`
- `flush()` is now synchronous, use `shutdown()` for cleanup
```

### Step 3: Update SDK with Testing

```bash
# Create upgrade branch
git checkout -b chore/upgrade-langfuse

# Update package (Node.js)
npm install langfuse@latest

# Update package (Python)
pip install --upgrade langfuse

# Run tests
npm test
pytest

# If tests fail, check specific version
npm install langfuse@3.0.0  # Specific version
```

### Step 4: Handle TypeScript API Changes

```typescript
// BEFORE (v2.x)
import Langfuse from "langfuse";

const langfuse = new Langfuse({
  publicKey: "...",
  secretKey: "...",
});

// Old async pattern
const trace = await langfuse.trace({ name: "test" });
await langfuse.flush();

// AFTER (v3.x)
import { Langfuse } from "langfuse";  // Named export

const langfuse = new Langfuse({
  publicKey: "...",
  secretKey: "...",
});

// New sync pattern (returns immediately, batches in background)
const trace = langfuse.trace({ name: "test" });
await langfuse.flushAsync();  // Renamed method
```

### Step 5: Handle Python API Changes

```python
# BEFORE (v1.x)
from langfuse import Langfuse

langfuse = Langfuse()

@langfuse.observe()  # Old decorator syntax
def my_function():
    pass

langfuse.flush()

# AFTER (v2.x)
from langfuse import Langfuse
from langfuse.decorators import observe, langfuse_context  # New import

langfuse = Langfuse()

@observe()  # New decorator (no langfuse prefix)
def my_function():
    # Access context via langfuse_context
    langfuse_context.update_current_observation(
        metadata={"key": "value"}
    )
    pass

langfuse.flush()
```

### Step 6: Migration Script for Codemod

```typescript
// scripts/migrate-langfuse.ts
import { Project } from "ts-morph";

const project = new Project({
  tsConfigFilePath: "./tsconfig.json",
});

const sourceFiles = project.getSourceFiles();

for (const sourceFile of sourceFiles) {
  let modified = false;

  // Update imports
  const importDeclarations = sourceFile.getImportDeclarations();
  for (const importDecl of importDeclarations) {
    if (importDecl.getModuleSpecifierValue() === "langfuse") {
      // Change default import to named import
      const defaultImport = importDecl.getDefaultImport();
      if (defaultImport?.getText() === "Langfuse") {
        importDecl.removeDefaultImport();
        importDecl.addNamedImport("Langfuse");
        modified = true;
      }
    }
  }

  // Update flush() to flushAsync()
  sourceFile.forEachDescendant((node) => {
    if (
      node.getKindName() === "CallExpression" &&
      node.getText().includes(".flush()")
    ) {
      node
        .asKind(ts.SyntaxKind.CallExpression)
        ?.getExpression()
        .replaceWithText(
          node
            .getText()
            .replace(".flush()", ".flushAsync()")
        );
      modified = true;
    }
  });

  if (modified) {
    console.log(`Updated: ${sourceFile.getFilePath()}`);
    sourceFile.saveSync();
  }
}
```

### Step 7: Test Migration

```typescript
// tests/langfuse-migration.test.ts
import { describe, it, expect, beforeAll, afterAll } from "vitest";
import { Langfuse } from "langfuse";

describe("Langfuse Migration Tests", () => {
  let langfuse: Langfuse;

  beforeAll(() => {
    langfuse = new Langfuse({
      publicKey: process.env.LANGFUSE_PUBLIC_KEY!,
      secretKey: process.env.LANGFUSE_SECRET_KEY!,
    });
  });

  afterAll(async () => {
    await langfuse.shutdownAsync();
  });

  it("should create trace with new API", () => {
    const trace = langfuse.trace({
      name: "migration-test",
      metadata: { version: "3.x" },
    });

    expect(trace).toBeDefined();
    expect(trace.id).toBeDefined();
  });

  it("should create generation with correct usage format", () => {
    const trace = langfuse.trace({ name: "generation-test" });

    const generation = trace.generation({
      name: "test-gen",
      model: "gpt-4",
      input: [{ role: "user", content: "test" }],
    });

    generation.end({
      output: "response",
      usage: {
        promptTokens: 10,      // New format (camelCase)
        completionTokens: 20,  // New format (camelCase)
      },
    });

    expect(generation.id).toBeDefined();
  });

  it("should flush with new async method", async () => {
    langfuse.trace({ name: "flush-test" });

    // New method name
    await expect(langfuse.flushAsync()).resolves.not.toThrow();
  });
});
```

## Output
- Updated SDK to latest version
- Migrated deprecated API calls
- All tests passing
- No breaking changes in functionality

## Version Compatibility Matrix

| Feature | v2.x | v3.x | Migration |
|---------|------|------|-----------|
| Default import | `import Langfuse` | `import { Langfuse }` | Update imports |
| `trace()` return | `Promise<Trace>` | `Trace` | Remove `await` |
| `flush()` | Sync | N/A | Use `flushAsync()` |
| Usage keys | `snake_case` | `camelCase` | Update all usage objects |

## Error Handling
| Error | Cause | Solution |
|-------|-------|----------|
| Import error | Changed export | Use named import |
| Type error on usage | Key name change | Use camelCase keys |
| flush() not found | Method renamed | Use `flushAsync()` |
| Decorator error | New import path | Import from `langfuse.decorators` |

## Rollback Plan

```bash
# If upgrade fails, rollback
git checkout main
npm install langfuse@2.0.0  # Previous version

# Or with lock file
git checkout HEAD -- package-lock.json
npm ci
```

## Resources
- [Langfuse JS Changelog](https://github.com/langfuse/langfuse-js/releases)
- [Langfuse Python Changelog](https://github.com/langfuse/langfuse-python/releases)
- [Langfuse Migration Guide](https://langfuse.com/docs/sdk)
- [Langfuse Discord](https://langfuse.com/discord)

## Next Steps
For CI/CD integration, see `langfuse-ci-integration`.
