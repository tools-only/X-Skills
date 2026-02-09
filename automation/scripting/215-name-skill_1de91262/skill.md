---
name: obsidian-rate-limits
description: |
  Handle Obsidian file system operations and throttling patterns.
  Use when processing many files, handling bulk operations,
  or preventing performance issues from excessive operations.
  Trigger with phrases like "obsidian rate limit", "obsidian bulk operations",
  "obsidian file throttling", "obsidian performance limits".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Obsidian Rate Limits

## Overview
Manage file system operations and implement throttling to prevent performance issues in Obsidian plugins.

## Prerequisites
- Understanding of async JavaScript
- Familiarity with Obsidian vault operations
- Knowledge of file system performance considerations

## Key Concepts

### Obsidian Operation Limits
| Operation | Recommended Limit | Risk if Exceeded |
|-----------|------------------|------------------|
| File reads | 100/second | UI freeze |
| File writes | 10/second | Data corruption risk |
| Metadata cache reads | 1000/second | Memory pressure |
| DOM updates | 60/second | Visual lag |
| Event emissions | 100/second | Event queue backup |

## Instructions

### Step 1: Implement Async Queue
```typescript
// src/utils/async-queue.ts
export class AsyncQueue {
  private queue: Array<() => Promise<void>> = [];
  private processing = false;
  private concurrency: number;
  private activeCount = 0;

  constructor(concurrency: number = 1) {
    this.concurrency = concurrency;
  }

  async add<T>(task: () => Promise<T>): Promise<T> {
    return new Promise((resolve, reject) => {
      this.queue.push(async () => {
        try {
          const result = await task();
          resolve(result);
        } catch (error) {
          reject(error);
        }
      });
      this.process();
    });
  }

  private async process(): Promise<void> {
    if (this.activeCount >= this.concurrency) return;

    const task = this.queue.shift();
    if (!task) return;

    this.activeCount++;
    try {
      await task();
    } finally {
      this.activeCount--;
      this.process();
    }
  }

  get pending(): number {
    return this.queue.length;
  }

  get active(): number {
    return this.activeCount;
  }
}
```

### Step 2: Implement Rate Limiter
```typescript
// src/utils/rate-limiter.ts
export class RateLimiter {
  private timestamps: number[] = [];
  private limit: number;
  private window: number; // milliseconds

  constructor(limit: number, windowMs: number) {
    this.limit = limit;
    this.window = windowMs;
  }

  async acquire(): Promise<void> {
    const now = Date.now();

    // Remove timestamps outside window
    this.timestamps = this.timestamps.filter(t => now - t < this.window);

    if (this.timestamps.length >= this.limit) {
      // Calculate wait time
      const oldestTimestamp = this.timestamps[0];
      const waitTime = this.window - (now - oldestTimestamp);
      await this.sleep(waitTime);
      return this.acquire();
    }

    this.timestamps.push(now);
  }

  private sleep(ms: number): Promise<void> {
    return new Promise(resolve => setTimeout(resolve, ms));
  }
}

// Usage:
const writeLimiter = new RateLimiter(10, 1000); // 10 per second

async function safeWriteFile(file: TFile, content: string) {
  await writeLimiter.acquire();
  await this.app.vault.modify(file, content);
}
```

### Step 3: Batch Processing Pattern
```typescript
// src/utils/batch-processor.ts
export interface BatchOptions {
  batchSize: number;
  delayBetweenBatches: number;
  onProgress?: (processed: number, total: number) => void;
}

export async function processBatches<T, R>(
  items: T[],
  processor: (item: T) => Promise<R>,
  options: BatchOptions
): Promise<R[]> {
  const results: R[] = [];
  const total = items.length;

  for (let i = 0; i < items.length; i += options.batchSize) {
    const batch = items.slice(i, i + options.batchSize);

    // Process batch in parallel
    const batchResults = await Promise.all(
      batch.map(item => processor(item))
    );
    results.push(...batchResults);

    // Report progress
    options.onProgress?.(Math.min(i + options.batchSize, total), total);

    // Delay between batches (except for last batch)
    if (i + options.batchSize < items.length) {
      await sleep(options.delayBetweenBatches);
    }
  }

  return results;
}

function sleep(ms: number): Promise<void> {
  return new Promise(resolve => setTimeout(resolve, ms));
}

// Usage:
const files = this.app.vault.getMarkdownFiles();

const results = await processBatches(
  files,
  async (file) => {
    const content = await this.app.vault.read(file);
    return { path: file.path, length: content.length };
  },
  {
    batchSize: 50,
    delayBetweenBatches: 100, // 100ms between batches
    onProgress: (processed, total) => {
      console.log(`Processed ${processed}/${total}`);
    },
  }
);
```

### Step 4: Debounced Operations
```typescript
// src/utils/debounce.ts
export function debounce<T extends (...args: any[]) => any>(
  func: T,
  wait: number,
  options: { leading?: boolean; trailing?: boolean; maxWait?: number } = {}
): (...args: Parameters<T>) => void {
  let timeout: NodeJS.Timeout | null = null;
  let lastCallTime: number | null = null;
  let lastInvokeTime = 0;
  let result: ReturnType<T>;

  const { leading = false, trailing = true, maxWait } = options;

  function invokeFunc(time: number, args: Parameters<T>) {
    lastInvokeTime = time;
    result = func(...args);
    return result;
  }

  return function (...args: Parameters<T>) {
    const time = Date.now();
    const isInvoking = shouldInvoke(time);

    lastCallTime = time;

    if (isInvoking) {
      if (!timeout && leading) {
        return invokeFunc(time, args);
      }
    }

    if (!timeout) {
      timeout = setTimeout(() => {
        timeout = null;
        if (trailing && lastCallTime) {
          invokeFunc(Date.now(), args);
        }
      }, wait);
    }

    // Handle maxWait
    if (maxWait !== undefined) {
      const timeSinceLastInvoke = time - lastInvokeTime;
      if (timeSinceLastInvoke >= maxWait) {
        invokeFunc(time, args);
      }
    }
  };

  function shouldInvoke(time: number): boolean {
    const timeSinceLastCall = lastCallTime ? time - lastCallTime : 0;
    return !lastCallTime || timeSinceLastCall >= wait;
  }
}

// Usage for search input:
const debouncedSearch = debounce(
  async (query: string) => {
    const results = await performSearch(query);
    updateUI(results);
  },
  300,
  { leading: false, trailing: true }
);

// Typing triggers search only after 300ms pause
inputEl.addEventListener('input', (e) => {
  debouncedSearch(e.target.value);
});
```

### Step 5: Throttled Event Handling
```typescript
// src/utils/throttle.ts
export function throttle<T extends (...args: any[]) => any>(
  func: T,
  limit: number
): (...args: Parameters<T>) => void {
  let inThrottle = false;
  let lastArgs: Parameters<T> | null = null;

  return function (...args: Parameters<T>) {
    if (!inThrottle) {
      func(...args);
      inThrottle = true;
      setTimeout(() => {
        inThrottle = false;
        if (lastArgs) {
          func(...lastArgs);
          lastArgs = null;
        }
      }, limit);
    } else {
      lastArgs = args;
    }
  };
}

// Usage for scroll handling:
const throttledOnScroll = throttle(() => {
  // Update something based on scroll position
  console.log('Scroll position:', window.scrollY);
}, 100);

window.addEventListener('scroll', throttledOnScroll);
```

## Output
- Async queue for sequential operations
- Rate limiter for controlled throughput
- Batch processor for bulk operations
- Debounce for user input
- Throttle for frequent events

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| UI freezes | Too many sync operations | Use async with batching |
| Data loss | Write conflicts | Use async queue |
| Memory pressure | Too many cached reads | Use generators |
| Event storms | No debouncing | Debounce user input |
| Missed updates | Over-throttling | Reduce throttle time |

## Examples

### Progress Modal for Long Operations
```typescript
async function processAllFiles(app: App, plugin: Plugin) {
  const files = app.vault.getMarkdownFiles();
  const progressModal = new ProgressModal(app);
  progressModal.open();

  try {
    await processBatches(
      files,
      async (file) => {
        // Process file
        return processFile(file);
      },
      {
        batchSize: 20,
        delayBetweenBatches: 50,
        onProgress: (processed, total) => {
          const percent = Math.round((processed / total) * 100);
          progressModal.setProgress(percent, `Processing ${processed}/${total} files`);
        },
      }
    );
    new Notice('Processing complete!');
  } finally {
    progressModal.close();
  }
}
```

### Generator for Large File Sets
```typescript
async function* iterateFilesWithPause(
  files: TFile[],
  pauseEvery: number = 100,
  pauseMs: number = 10
): AsyncGenerator<TFile> {
  for (let i = 0; i < files.length; i++) {
    yield files[i];

    // Pause periodically to allow UI updates
    if (i > 0 && i % pauseEvery === 0) {
      await new Promise(r => setTimeout(r, pauseMs));
    }
  }
}

// Usage:
for await (const file of iterateFilesWithPause(files)) {
  await processFile(file);
}
```

## Resources
- [Obsidian Performance Tips](https://docs.obsidian.md/Plugins/Guides/Performance)
- [JavaScript Event Loop](https://developer.mozilla.org/en-US/docs/Web/JavaScript/EventLoop)

## Next Steps
For security practices, see `obsidian-security-basics`.
