---
name: obsidian-performance-tuning
description: |
  Optimize Obsidian plugin performance for smooth operation.
  Use when experiencing lag, memory issues, or slow startup,
  or when optimizing plugin code for large vaults.
  Trigger with phrases like "obsidian performance", "obsidian slow",
  "optimize obsidian plugin", "obsidian memory usage".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Obsidian Performance Tuning

## Overview
Optimize Obsidian plugin performance for smooth operation in large vaults and resource-constrained environments.

## Prerequisites
- Working Obsidian plugin
- Developer Tools access (Ctrl/Cmd+Shift+I)
- Understanding of async JavaScript

## Performance Benchmarks

### Target Metrics
| Metric | Good | Warning | Critical |
|--------|------|---------|----------|
| Plugin load time | < 100ms | 100-500ms | > 500ms |
| Command execution | < 50ms | 50-200ms | > 200ms |
| File operation | < 10ms | 10-50ms | > 50ms |
| Memory increase | < 10MB | 10-50MB | > 50MB |
| Event handler | < 5ms | 5-20ms | > 20ms |

## Instructions

### Step 1: Profile Plugin Performance
```typescript
// src/utils/profiler.ts
export class PerformanceProfiler {
  private marks: Map<string, number> = new Map();
  private enabled: boolean;

  constructor(enabled: boolean = true) {
    this.enabled = enabled;
  }

  start(label: string): void {
    if (!this.enabled) return;
    this.marks.set(label, performance.now());
  }

  end(label: string): number {
    if (!this.enabled) return 0;

    const start = this.marks.get(label);
    if (!start) return 0;

    const duration = performance.now() - start;
    this.marks.delete(label);

    if (duration > 50) {
      console.warn(`[Performance] ${label}: ${duration.toFixed(2)}ms (slow)`);
    } else {
      console.log(`[Performance] ${label}: ${duration.toFixed(2)}ms`);
    }

    return duration;
  }

  async measure<T>(label: string, fn: () => Promise<T>): Promise<T> {
    this.start(label);
    try {
      return await fn();
    } finally {
      this.end(label);
    }
  }

  measureSync<T>(label: string, fn: () => T): T {
    this.start(label);
    try {
      return fn();
    } finally {
      this.end(label);
    }
  }
}

// Usage:
const profiler = new PerformanceProfiler(process.env.NODE_ENV !== 'production');

await profiler.measure('loadSettings', async () => {
  return this.loadData();
});
```

### Step 2: Lazy Initialization
```typescript
// src/services/lazy-service.ts
export class LazyService<T> {
  private instance: T | null = null;
  private initializing: Promise<T> | null = null;
  private factory: () => Promise<T>;

  constructor(factory: () => Promise<T>) {
    this.factory = factory;
  }

  async get(): Promise<T> {
    if (this.instance) return this.instance;

    if (this.initializing) return this.initializing;

    this.initializing = this.factory().then(instance => {
      this.instance = instance;
      this.initializing = null;
      return instance;
    });

    return this.initializing;
  }

  isInitialized(): boolean {
    return this.instance !== null;
  }

  clear(): void {
    this.instance = null;
    this.initializing = null;
  }
}

// Usage - defer expensive initialization
export default class MyPlugin extends Plugin {
  private indexService = new LazyService(() => this.buildIndex());

  async onload() {
    // Don't build index immediately
    // Build on first use instead

    this.addCommand({
      id: 'search',
      name: 'Search',
      callback: async () => {
        const index = await this.indexService.get();
        // Use index...
      },
    });
  }

  private async buildIndex(): Promise<SearchIndex> {
    // Expensive operation - only runs when first needed
    const files = this.app.vault.getMarkdownFiles();
    return new SearchIndex(files);
  }
}
```

### Step 3: Efficient File Processing
```typescript
// src/utils/file-processor.ts
import { TFile, Vault } from 'obsidian';

export class EfficientFileProcessor {
  private vault: Vault;
  private cache: Map<string, { content: string; mtime: number }> = new Map();

  constructor(vault: Vault) {
    this.vault = vault;
  }

  // Use cached read when possible
  async readWithCache(file: TFile): Promise<string> {
    const cached = this.cache.get(file.path);
    if (cached && cached.mtime === file.stat.mtime) {
      return cached.content;
    }

    const content = await this.vault.cachedRead(file);
    this.cache.set(file.path, {
      content,
      mtime: file.stat.mtime,
    });

    return content;
  }

  // Process files in chunks with pauses
  async processFilesInChunks<T>(
    files: TFile[],
    processor: (file: TFile) => Promise<T>,
    options: {
      chunkSize?: number;
      pauseMs?: number;
      onProgress?: (processed: number, total: number) => void;
    } = {}
  ): Promise<T[]> {
    const { chunkSize = 50, pauseMs = 10, onProgress } = options;
    const results: T[] = [];

    for (let i = 0; i < files.length; i += chunkSize) {
      const chunk = files.slice(i, i + chunkSize);

      // Process chunk in parallel
      const chunkResults = await Promise.all(
        chunk.map(file => processor(file))
      );
      results.push(...chunkResults);

      // Report progress
      onProgress?.(Math.min(i + chunkSize, files.length), files.length);

      // Pause to allow UI updates
      if (i + chunkSize < files.length) {
        await new Promise(r => setTimeout(r, pauseMs));
      }
    }

    return results;
  }

  // Generator for memory-efficient iteration
  async *iterateFiles(
    files: TFile[]
  ): AsyncGenerator<{ file: TFile; content: string }> {
    for (const file of files) {
      const content = await this.vault.cachedRead(file);
      yield { file, content };

      // Allow event loop to process
      await new Promise(r => setTimeout(r, 0));
    }
  }

  clearCache(): void {
    this.cache.clear();
  }

  removeFromCache(path: string): void {
    this.cache.delete(path);
  }
}
```

### Step 4: Memory-Efficient Data Structures
```typescript
// src/utils/efficient-structures.ts

// Use WeakMap for cached data that can be garbage collected
export class WeakFileCache<T> {
  private cache = new WeakMap<object, T>();
  private keyMap = new Map<string, WeakRef<object>>();

  set(path: string, file: object, value: T): void {
    this.cache.set(file, value);
    this.keyMap.set(path, new WeakRef(file));
  }

  get(file: object): T | undefined {
    return this.cache.get(file);
  }

  getByPath(path: string): T | undefined {
    const ref = this.keyMap.get(path);
    if (!ref) return undefined;

    const file = ref.deref();
    if (!file) {
      this.keyMap.delete(path);
      return undefined;
    }

    return this.cache.get(file);
  }

  has(file: object): boolean {
    return this.cache.has(file);
  }
}

// Efficient string storage for repeated values
export class StringInterner {
  private strings = new Map<string, string>();

  intern(str: string): string {
    const existing = this.strings.get(str);
    if (existing) return existing;
    this.strings.set(str, str);
    return str;
  }

  clear(): void {
    this.strings.clear();
  }

  get size(): number {
    return this.strings.size;
  }
}

// LRU cache for bounded memory usage
export class LRUCache<K, V> {
  private cache = new Map<K, V>();
  private maxSize: number;

  constructor(maxSize: number) {
    this.maxSize = maxSize;
  }

  get(key: K): V | undefined {
    const value = this.cache.get(key);
    if (value !== undefined) {
      // Move to end (most recently used)
      this.cache.delete(key);
      this.cache.set(key, value);
    }
    return value;
  }

  set(key: K, value: V): void {
    if (this.cache.has(key)) {
      this.cache.delete(key);
    } else if (this.cache.size >= this.maxSize) {
      // Remove oldest entry
      const firstKey = this.cache.keys().next().value;
      this.cache.delete(firstKey);
    }
    this.cache.set(key, value);
  }

  clear(): void {
    this.cache.clear();
  }
}
```

### Step 5: Optimize Event Handlers
```typescript
// src/utils/event-optimizer.ts
import { debounce, throttle } from 'lodash-es';

export class OptimizedEventManager {
  private plugin: Plugin;

  constructor(plugin: Plugin) {
    this.plugin = plugin;
  }

  // Debounced handler for file modifications
  registerDebouncedModify(
    handler: (file: TFile) => void,
    wait: number = 500
  ): void {
    const debouncedHandler = debounce(handler, wait);

    this.plugin.registerEvent(
      this.plugin.app.vault.on('modify', (file) => {
        if (file instanceof TFile) {
          debouncedHandler(file);
        }
      })
    );
  }

  // Throttled handler for frequent events
  registerThrottledScroll(
    element: HTMLElement,
    handler: (event: Event) => void,
    wait: number = 100
  ): void {
    const throttledHandler = throttle(handler, wait);
    this.plugin.registerDomEvent(element, 'scroll', throttledHandler);
  }

  // Batch multiple rapid events
  registerBatchedEvents(
    handler: (files: TFile[]) => void,
    wait: number = 100
  ): void {
    let pendingFiles: TFile[] = [];
    let timeoutId: NodeJS.Timeout | null = null;

    this.plugin.registerEvent(
      this.plugin.app.vault.on('modify', (file) => {
        if (file instanceof TFile) {
          pendingFiles.push(file);

          if (timeoutId) clearTimeout(timeoutId);

          timeoutId = setTimeout(() => {
            const files = [...new Set(pendingFiles)];
            pendingFiles = [];
            timeoutId = null;
            handler(files);
          }, wait);
        }
      })
    );
  }
}
```

### Step 6: UI Rendering Optimization
```typescript
// src/utils/render-optimizer.ts
export class RenderOptimizer {
  // Use DocumentFragment for batch DOM updates
  static batchRender(
    container: HTMLElement,
    items: string[],
    renderer: (item: string) => HTMLElement
  ): void {
    const fragment = document.createDocumentFragment();

    for (const item of items) {
      fragment.appendChild(renderer(item));
    }

    container.empty();
    container.appendChild(fragment);
  }

  // Virtual scrolling for long lists
  static createVirtualList(
    container: HTMLElement,
    items: any[],
    itemHeight: number,
    renderItem: (item: any) => HTMLElement
  ): void {
    const visibleCount = Math.ceil(container.clientHeight / itemHeight) + 2;
    let startIndex = 0;

    const render = () => {
      const scrollTop = container.scrollTop;
      const newStartIndex = Math.floor(scrollTop / itemHeight);

      if (newStartIndex !== startIndex) {
        startIndex = newStartIndex;
        container.empty();

        const fragment = document.createDocumentFragment();
        const spacer = document.createElement('div');
        spacer.style.height = `${startIndex * itemHeight}px`;
        fragment.appendChild(spacer);

        for (let i = startIndex; i < Math.min(startIndex + visibleCount, items.length); i++) {
          fragment.appendChild(renderItem(items[i]));
        }

        const bottomSpacer = document.createElement('div');
        bottomSpacer.style.height = `${(items.length - startIndex - visibleCount) * itemHeight}px`;
        fragment.appendChild(bottomSpacer);

        container.appendChild(fragment);
      }
    };

    container.addEventListener('scroll', render);
    render();
  }

  // Request animation frame for smooth updates
  static smoothUpdate(fn: () => void): number {
    return requestAnimationFrame(fn);
  }

  // Batch style changes
  static batchStyles(element: HTMLElement, styles: Record<string, string>): void {
    requestAnimationFrame(() => {
      Object.assign(element.style, styles);
    });
  }
}
```

## Output
- Performance profiler for identifying bottlenecks
- Lazy initialization patterns
- Efficient file processing with chunking
- Memory-efficient data structures
- Optimized event handlers
- UI rendering optimizations

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Plugin slow to load | Heavy initialization | Use lazy loading |
| UI freezes | Blocking operations | Use async + chunking |
| Memory growth | Unbounded caching | Use LRU cache |
| Event lag | Unthrottled handlers | Debounce/throttle |

## Examples

### Memory Usage Monitor
```typescript
function logMemoryUsage(label: string): void {
  if (performance.memory) {
    const used = performance.memory.usedJSHeapSize / 1048576;
    console.log(`[Memory] ${label}: ${used.toFixed(2)} MB`);
  }
}

// Use during development
logMemoryUsage('Before index build');
await buildIndex();
logMemoryUsage('After index build');
```

### Performance Checklist
```markdown
## Pre-Release Performance Checklist

- [ ] Plugin loads in < 100ms
- [ ] No blocking operations in onload()
- [ ] File operations use cachedRead when possible
- [ ] Event handlers are debounced/throttled
- [ ] Large lists use virtual scrolling
- [ ] Caches have size limits (LRU)
- [ ] Memory doesn't grow unboundedly
- [ ] Works smoothly with 1000+ files
```

## Resources
- [Chrome DevTools Performance](https://developer.chrome.com/docs/devtools/performance/)
- [JavaScript Performance Patterns](https://developer.mozilla.org/en-US/docs/Web/Performance)
- [Obsidian Performance Tips](https://docs.obsidian.md/Plugins/Guides/Performance)

## Next Steps
For resource optimization, see `obsidian-cost-tuning`.
