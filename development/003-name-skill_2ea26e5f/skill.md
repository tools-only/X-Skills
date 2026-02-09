---
name: obsidian-sdk-patterns
description: |
  Apply production-ready Obsidian plugin patterns for TypeScript.
  Use when implementing complex features, refactoring plugins,
  or establishing coding standards for Obsidian development.
  Trigger with phrases like "obsidian patterns", "obsidian best practices",
  "obsidian code patterns", "idiomatic obsidian plugin".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Obsidian SDK Patterns

## Overview
Production-ready patterns for Obsidian plugin development in TypeScript.

## Prerequisites
- Completed `obsidian-install-auth` setup
- Familiarity with TypeScript and async/await
- Understanding of Obsidian's plugin lifecycle

## Instructions

### Step 1: Type-Safe Settings Pattern
```typescript
// src/settings.ts
import { App, PluginSettingTab, Setting } from 'obsidian';
import type MyPlugin from './main';

export interface MyPluginSettings {
  apiEndpoint: string;
  enableFeatureX: boolean;
  maxItems: number;
  excludedFolders: string[];
}

export const DEFAULT_SETTINGS: MyPluginSettings = {
  apiEndpoint: 'https://api.example.com',
  enableFeatureX: true,
  maxItems: 100,
  excludedFolders: ['templates', 'archive'],
};

export class MyPluginSettingTab extends PluginSettingTab {
  plugin: MyPlugin;

  constructor(app: App, plugin: MyPlugin) {
    super(app, plugin);
    this.plugin = plugin;
  }

  display(): void {
    const { containerEl } = this;
    containerEl.empty();

    new Setting(containerEl)
      .setName('API Endpoint')
      .setDesc('The API endpoint for fetching data')
      .addText(text => text
        .setPlaceholder('https://api.example.com')
        .setValue(this.plugin.settings.apiEndpoint)
        .onChange(async (value) => {
          this.plugin.settings.apiEndpoint = value;
          await this.plugin.saveSettings();
        }));

    new Setting(containerEl)
      .setName('Enable Feature X')
      .setDesc('Toggle experimental feature')
      .addToggle(toggle => toggle
        .setValue(this.plugin.settings.enableFeatureX)
        .onChange(async (value) => {
          this.plugin.settings.enableFeatureX = value;
          await this.plugin.saveSettings();
        }));

    new Setting(containerEl)
      .setName('Max Items')
      .setDesc('Maximum number of items to display')
      .addSlider(slider => slider
        .setLimits(10, 500, 10)
        .setValue(this.plugin.settings.maxItems)
        .setDynamicTooltip()
        .onChange(async (value) => {
          this.plugin.settings.maxItems = value;
          await this.plugin.saveSettings();
        }));
  }
}
```

### Step 2: Service Layer Pattern
```typescript
// src/services/vault-service.ts
import { App, TFile, TFolder, Vault } from 'obsidian';

export class VaultService {
  constructor(private app: App) {}

  async getMarkdownFiles(folder?: string): Promise<TFile[]> {
    const files = this.app.vault.getMarkdownFiles();
    if (!folder) return files;

    return files.filter(f => f.path.startsWith(folder));
  }

  async readFile(file: TFile): Promise<string> {
    return this.app.vault.read(file);
  }

  async writeFile(file: TFile, content: string): Promise<void> {
    await this.app.vault.modify(file, content);
  }

  async createFile(path: string, content: string): Promise<TFile> {
    return this.app.vault.create(path, content);
  }

  async ensureFolder(path: string): Promise<TFolder> {
    const folder = this.app.vault.getAbstractFileByPath(path);
    if (folder instanceof TFolder) return folder;

    await this.app.vault.createFolder(path);
    return this.app.vault.getAbstractFileByPath(path) as TFolder;
  }

  getFileByPath(path: string): TFile | null {
    const file = this.app.vault.getAbstractFileByPath(path);
    return file instanceof TFile ? file : null;
  }
}
```

### Step 3: Event Management Pattern
```typescript
// src/events.ts
import { Plugin, EventRef, Events } from 'obsidian';

export class EventManager {
  private eventRefs: EventRef[] = [];

  constructor(private plugin: Plugin) {}

  register(events: Events, name: string, callback: (...args: any[]) => any): void {
    const ref = events.on(name as any, callback);
    this.eventRefs.push(ref);
    this.plugin.registerEvent(ref);
  }

  registerWorkspaceEvent(name: string, callback: (...args: any[]) => any): void {
    this.register(this.plugin.app.workspace, name, callback);
  }

  registerVaultEvent(name: string, callback: (...args: any[]) => any): void {
    this.register(this.plugin.app.vault, name, callback);
  }

  cleanup(): void {
    // Events are automatically cleaned up by Obsidian
    // But keep track for manual cleanup if needed
    this.eventRefs = [];
  }
}

// Usage in main.ts:
const eventManager = new EventManager(this);
eventManager.registerWorkspaceEvent('file-open', (file) => {
  if (file) console.log('Opened:', file.path);
});
eventManager.registerVaultEvent('modify', (file) => {
  console.log('Modified:', file.path);
});
```

### Step 4: Command Builder Pattern
```typescript
// src/commands.ts
import { Command, Editor, MarkdownView, Plugin } from 'obsidian';

interface CommandConfig {
  id: string;
  name: string;
  icon?: string;
  hotkeys?: { modifiers: string[]; key: string }[];
}

export class CommandBuilder {
  private commands: Command[] = [];

  constructor(private plugin: Plugin, private prefix: string) {}

  addSimple(config: CommandConfig, callback: () => void): this {
    this.commands.push({
      id: `${this.prefix}-${config.id}`,
      name: config.name,
      icon: config.icon,
      callback,
    });
    return this;
  }

  addEditor(
    config: CommandConfig,
    callback: (editor: Editor, view: MarkdownView) => void
  ): this {
    this.commands.push({
      id: `${this.prefix}-${config.id}`,
      name: config.name,
      icon: config.icon,
      editorCallback: callback,
    });
    return this;
  }

  addCheck(
    config: CommandConfig,
    check: () => boolean,
    callback: () => void
  ): this {
    this.commands.push({
      id: `${this.prefix}-${config.id}`,
      name: config.name,
      icon: config.icon,
      checkCallback: (checking) => {
        if (checking) return check();
        callback();
        return true;
      },
    });
    return this;
  }

  register(): void {
    this.commands.forEach(cmd => this.plugin.addCommand(cmd));
  }
}

// Usage:
new CommandBuilder(this, 'my-plugin')
  .addSimple({ id: 'greet', name: 'Show Greeting' }, () => {
    new Notice('Hello!');
  })
  .addEditor({ id: 'insert', name: 'Insert Text' }, (editor) => {
    editor.replaceSelection('Inserted text');
  })
  .register();
```

### Step 5: Async Queue Pattern
```typescript
// src/utils/async-queue.ts
export class AsyncQueue {
  private queue: (() => Promise<void>)[] = [];
  private processing = false;

  async add(task: () => Promise<void>): Promise<void> {
    return new Promise((resolve, reject) => {
      this.queue.push(async () => {
        try {
          await task();
          resolve();
        } catch (e) {
          reject(e);
        }
      });
      this.process();
    });
  }

  private async process(): Promise<void> {
    if (this.processing) return;
    this.processing = true;

    while (this.queue.length > 0) {
      const task = this.queue.shift();
      if (task) await task();
    }

    this.processing = false;
  }
}

// Usage for rate-limited operations:
const writeQueue = new AsyncQueue();

async function safeWrite(file: TFile, content: string) {
  await writeQueue.add(async () => {
    await this.app.vault.modify(file, content);
  });
}
```

## Output
- Type-safe settings management
- Service layer for vault operations
- Event registration with automatic cleanup
- Fluent command builder
- Async queue for rate limiting

## Error Handling
| Pattern | Use Case | Benefit |
|---------|----------|---------|
| Settings validation | User input | Prevents invalid config |
| Service layer | Vault access | Centralizes file ops |
| Event manager | Lifecycle events | Prevents memory leaks |
| Command builder | Plugin commands | Cleaner registration |
| Async queue | Bulk operations | Prevents race conditions |

## Examples

### Debounce Pattern
```typescript
// src/utils/debounce.ts
export function debounce<T extends (...args: any[]) => any>(
  func: T,
  wait: number
): (...args: Parameters<T>) => void {
  let timeout: NodeJS.Timeout | null = null;

  return (...args: Parameters<T>) => {
    if (timeout) clearTimeout(timeout);
    timeout = setTimeout(() => func(...args), wait);
  };
}

// Usage for search:
const debouncedSearch = debounce(async (query: string) => {
  const results = await performSearch(query);
  updateUI(results);
}, 300);
```

### Singleton Service Pattern
```typescript
// src/main.ts
export default class MyPlugin extends Plugin {
  private static instance: MyPlugin;
  private vaultService: VaultService;

  static getInstance(): MyPlugin {
    return MyPlugin.instance;
  }

  async onload() {
    MyPlugin.instance = this;
    this.vaultService = new VaultService(this.app);
  }

  getVaultService(): VaultService {
    return this.vaultService;
  }
}
```

## Resources
- [Obsidian API Reference](https://docs.obsidian.md/Reference/TypeScript+API)
- [Obsidian Plugin Developer Docs](https://docs.obsidian.md/Plugins)
- [TypeScript Best Practices](https://www.typescriptlang.org/docs/handbook/)

## Next Steps
Apply patterns in `obsidian-core-workflow-a` for vault operations.
