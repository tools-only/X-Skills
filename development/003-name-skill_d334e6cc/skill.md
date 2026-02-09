---
name: lokalise-reference-architecture
description: |
  Implement Lokalise reference architecture with best-practice project layout.
  Use when designing new Lokalise integrations, reviewing project structure,
  or establishing architecture standards for Lokalise applications.
  Trigger with phrases like "lokalise architecture", "lokalise best practices",
  "lokalise project structure", "how to organize lokalise", "lokalise layout".
allowed-tools: Read, Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Lokalise Reference Architecture

## Overview
Production-ready architecture patterns for Lokalise integrations.

## Prerequisites
- Understanding of layered architecture
- Lokalise SDK knowledge
- TypeScript project setup
- Testing framework configured

## Project Structure

```
my-app/
├── src/
│   ├── i18n/
│   │   ├── index.ts              # i18n library setup
│   │   ├── config.ts             # Lokalise configuration
│   │   ├── types.ts              # TypeScript types
│   │   └── loaders/
│   │       ├── static.ts         # Bundled translations
│   │       ├── dynamic.ts        # Runtime loading
│   │       └── ota.ts            # Over-the-air (mobile)
│   ├── services/
│   │   └── lokalise/
│   │       ├── index.ts          # Service facade
│   │       ├── client.ts         # Lokalise client wrapper
│   │       ├── cache.ts          # Caching layer
│   │       ├── sync.ts           # Translation sync
│   │       └── webhooks.ts       # Webhook handlers
│   ├── locales/
│   │   ├── en.json               # English (source)
│   │   ├── es.json               # Spanish
│   │   ├── fr.json               # French
│   │   └── index.ts              # Locale exports
│   └── api/
│       └── webhooks/
│           └── lokalise.ts       # Webhook endpoint
├── scripts/
│   ├── lokalise-pull.sh          # Download translations
│   ├── lokalise-push.sh          # Upload source strings
│   └── check-translations.ts     # Validation script
├── tests/
│   ├── unit/
│   │   └── i18n/
│   └── integration/
│       └── lokalise/
├── config/
│   ├── lokalise.development.json
│   ├── lokalise.staging.json
│   └── lokalise.production.json
├── .env.example
├── lokalise.json                 # CLI configuration
└── package.json
```

## Layer Architecture

```
┌─────────────────────────────────────────┐
│           Application Layer              │
│   (React/Vue/Angular Components)         │
├─────────────────────────────────────────┤
│           i18n Library Layer             │
│   (i18next, react-intl, vue-i18n)        │
├─────────────────────────────────────────┤
│         Translation Service Layer        │
│  (Loading, Caching, Fallback Logic)      │
├─────────────────────────────────────────┤
│           Lokalise Layer                 │
│   (SDK Client, Sync, Webhooks)           │
├─────────────────────────────────────────┤
│         Infrastructure Layer             │
│    (Cache, Queue, Monitoring)            │
└─────────────────────────────────────────┘
```

## Key Components

### Step 1: Lokalise Client Wrapper
```typescript
// src/services/lokalise/client.ts
import { LokaliseApi } from "@lokalise/node-api";

let instance: LokaliseApi | null = null;

export interface LokaliseConfig {
  apiKey: string;
  projectId: string;
  enableCompression?: boolean;
}

export function getLokaliseClient(config?: Partial<LokaliseConfig>): LokaliseApi {
  if (!instance) {
    instance = new LokaliseApi({
      apiKey: config?.apiKey || process.env.LOKALISE_API_TOKEN!,
      enableCompression: config?.enableCompression ?? true,
    });
  }
  return instance;
}

export function getProjectId(): string {
  return process.env.LOKALISE_PROJECT_ID!;
}

// Reset for testing
export function resetClient(): void {
  instance = null;
}
```

### Step 2: Translation Service Facade
```typescript
// src/services/lokalise/index.ts
import { getLokaliseClient, getProjectId } from "./client";
import { TranslationCache } from "./cache";
import { syncTranslations } from "./sync";

export interface TranslationService {
  getTranslations(locale: string): Promise<Record<string, string>>;
  getKey(locale: string, key: string): Promise<string | null>;
  syncFromLokalise(): Promise<void>;
  invalidateCache(locale?: string): void;
}

const cache = new TranslationCache();

export const translationService: TranslationService = {
  async getTranslations(locale) {
    // Try cache first
    const cached = cache.get(locale);
    if (cached) return cached;

    // Load from bundled files or Lokalise
    const translations = await loadTranslations(locale);
    cache.set(locale, translations);
    return translations;
  },

  async getKey(locale, key) {
    const translations = await this.getTranslations(locale);
    return translations[key] ?? null;
  },

  async syncFromLokalise() {
    await syncTranslations(getProjectId());
  },

  invalidateCache(locale) {
    if (locale) {
      cache.delete(locale);
    } else {
      cache.clear();
    }
  },
};
```

### Step 3: Translation Loader Pattern
```typescript
// src/i18n/loaders/static.ts
// For bundled translations (build-time)
export async function loadStaticTranslations(
  locale: string
): Promise<Record<string, string>> {
  try {
    const translations = await import(`../../locales/${locale}.json`);
    return translations.default;
  } catch {
    console.warn(`Locale ${locale} not found, falling back to English`);
    const fallback = await import("../../locales/en.json");
    return fallback.default;
  }
}

// src/i18n/loaders/dynamic.ts
// For runtime loading (CDN/API)
export async function loadDynamicTranslations(
  locale: string,
  baseUrl = "/locales"
): Promise<Record<string, string>> {
  const response = await fetch(`${baseUrl}/${locale}.json`);

  if (!response.ok) {
    console.warn(`Failed to load ${locale}, falling back to English`);
    return loadDynamicTranslations("en", baseUrl);
  }

  return response.json();
}
```

### Step 4: i18n Library Integration
```typescript
// src/i18n/index.ts
import i18n from "i18next";
import { initReactI18next } from "react-i18next";
import { loadStaticTranslations } from "./loaders/static";

export const SUPPORTED_LOCALES = ["en", "es", "fr", "de", "ja"];
export const DEFAULT_LOCALE = "en";

export async function initI18n(locale = DEFAULT_LOCALE) {
  const resources: Record<string, { translation: any }> = {};

  // Load initial locale
  resources[locale] = {
    translation: await loadStaticTranslations(locale),
  };

  // Load English fallback if different
  if (locale !== "en") {
    resources.en = {
      translation: await loadStaticTranslations("en"),
    };
  }

  await i18n.use(initReactI18next).init({
    resources,
    lng: locale,
    fallbackLng: "en",
    supportedLngs: SUPPORTED_LOCALES,
    interpolation: {
      escapeValue: false,
    },
  });

  return i18n;
}

// Lazy load additional locales
export async function loadLocale(locale: string) {
  if (i18n.hasResourceBundle(locale, "translation")) {
    return;
  }

  const translations = await loadStaticTranslations(locale);
  i18n.addResourceBundle(locale, "translation", translations);
}
```

## Data Flow Diagram

```
Developer adds string
        │
        ▼
┌───────────────┐
│  en.json      │──────────────────┐
│  (source)     │                  │
└───────┬───────┘                  │
        │                          ▼
        │ git push           ┌───────────────┐
        │                    │   CI/CD       │
        ▼                    │   Pipeline    │
┌───────────────┐            └───────┬───────┘
│   Lokalise    │◀───────────────────┘
│   Project     │      lokalise push
└───────┬───────┘
        │
        │ Translators work
        ▼
┌───────────────┐
│  Translations │
│  Complete     │
└───────┬───────┘
        │
        │ Webhook / CI sync
        ▼
┌───────────────┐
│   App Build   │────▶ Production
│   with i18n   │
└───────────────┘
```

## Configuration Management

```typescript
// config/lokalise.ts
import devConfig from "./lokalise.development.json";
import stagingConfig from "./lokalise.staging.json";
import prodConfig from "./lokalise.production.json";

type Environment = "development" | "staging" | "production";

interface LokaliseEnvConfig {
  projectId: string;
  enableWebhooks: boolean;
  cacheEnabled: boolean;
  cacheTtlSeconds: number;
}

const configs: Record<Environment, LokaliseEnvConfig> = {
  development: devConfig,
  staging: stagingConfig,
  production: prodConfig,
};

export function getLokaliseEnvConfig(): LokaliseEnvConfig {
  const env = (process.env.NODE_ENV || "development") as Environment;
  return configs[env] || configs.development;
}
```

## Instructions

### Step 1: Create Directory Structure
Set up the project layout following the reference structure.

### Step 2: Implement Client Wrapper
Create the singleton client with caching support.

### Step 3: Build Translation Service
Implement the facade pattern for translation operations.

### Step 4: Integrate i18n Library
Connect Lokalise translations to your UI framework.

## Output
- Structured project layout
- Client wrapper with caching
- Translation service facade
- i18n library integration

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Circular imports | Wrong layering | Separate by layer |
| Missing locale | Not bundled | Add fallback logic |
| Stale translations | Cache not invalidated | Use webhooks to invalidate |
| Type errors | Missing types | Generate from source locale |

## Examples

### Quick Setup Script
```bash
#!/bin/bash
# Create reference structure

mkdir -p src/i18n/loaders
mkdir -p src/services/lokalise
mkdir -p src/locales
mkdir -p scripts
mkdir -p config
mkdir -p tests/{unit,integration}/lokalise

touch src/i18n/{index,config,types}.ts
touch src/i18n/loaders/{static,dynamic,ota}.ts
touch src/services/lokalise/{index,client,cache,sync,webhooks}.ts
touch scripts/{lokalise-pull.sh,lokalise-push.sh,check-translations.ts}
touch config/lokalise.{development,staging,production}.json
```

## Resources
- [Lokalise SDK Documentation](https://lokalise.github.io/node-lokalise-api/)
- [i18next Documentation](https://www.i18next.com/)
- [React i18next](https://react.i18next.com/)

## Flagship Skills
For multi-environment setup, see `lokalise-multi-env-setup`.
