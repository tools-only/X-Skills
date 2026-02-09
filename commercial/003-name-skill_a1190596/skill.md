---
name: obsidian-security-basics
description: |
  Implement secure Obsidian plugin development practices.
  Use when handling user data, implementing authentication,
  or ensuring plugin security best practices.
  Trigger with phrases like "obsidian security", "secure obsidian plugin",
  "obsidian data protection", "obsidian privacy".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Obsidian Security Basics

## Overview
Implement secure coding practices for Obsidian plugin development to protect user data and vault contents.

## Prerequisites
- Understanding of web security concepts
- Familiarity with Obsidian plugin architecture
- Knowledge of TypeScript

## Security Principles for Obsidian Plugins

### Core Security Rules
1. **Never store secrets in code** - Use settings or environment
2. **Validate all user input** - Sanitize paths, content, and settings
3. **Minimize permissions** - Request only what you need
4. **Protect vault data** - Don't leak content externally without consent
5. **Handle errors gracefully** - Don't expose stack traces to users

## Instructions

### Step 1: Secure Settings Storage
```typescript
// src/settings.ts
import { Plugin, PluginSettingTab, Setting } from 'obsidian';

interface SecureSettings {
  apiEndpoint: string;
  // Never store API keys in plain settings
  // Use Obsidian's requestUrl or secure methods
}

export class SecureSettingsTab extends PluginSettingTab {
  plugin: MyPlugin;

  display(): void {
    const { containerEl } = this;
    containerEl.empty();

    // API Key - use password field, don't log
    new Setting(containerEl)
      .setName('API Key')
      .setDesc('Your API key (stored locally, never sent in logs)')
      .addText(text => {
        text.inputEl.type = 'password';
        text.inputEl.autocomplete = 'off';
        text
          .setPlaceholder('Enter API key')
          .setValue(this.plugin.settings.apiKey || '')
          .onChange(async (value) => {
            // Don't log the value!
            this.plugin.settings.apiKey = value;
            await this.plugin.saveSettings();
          });
      });

    // Add show/hide toggle
    new Setting(containerEl)
      .setName('Show API Key')
      .addToggle(toggle => toggle
        .setValue(false)
        .onChange(value => {
          const input = containerEl.querySelector('input[type="password"], input[type="text"]');
          if (input) {
            (input as HTMLInputElement).type = value ? 'text' : 'password';
          }
        }));
  }
}
```

### Step 2: Input Validation and Sanitization
```typescript
// src/security/validation.ts
export class InputValidator {
  // Validate file path - prevent directory traversal
  static validatePath(path: string): { valid: boolean; error?: string } {
    // Check for directory traversal
    if (path.includes('..')) {
      return { valid: false, error: 'Path cannot contain ".."' };
    }

    // Check for absolute paths (platform-specific)
    if (path.startsWith('/') || /^[A-Za-z]:/.test(path)) {
      return { valid: false, error: 'Absolute paths not allowed' };
    }

    // Check for null bytes
    if (path.includes('\0')) {
      return { valid: false, error: 'Invalid characters in path' };
    }

    // Check allowed extensions
    const allowedExtensions = ['.md', '.txt', '.json', '.yaml', '.yml'];
    const ext = path.substring(path.lastIndexOf('.'));
    if (!allowedExtensions.includes(ext.toLowerCase())) {
      return { valid: false, error: 'File type not allowed' };
    }

    return { valid: true };
  }

  // Sanitize HTML content (for rendering in views)
  static sanitizeHtml(html: string): string {
    // Use DOMPurify or similar in production
    // This is a basic example
    const div = document.createElement('div');
    div.textContent = html;
    return div.innerHTML;
  }

  // Validate URL
  static validateUrl(url: string): { valid: boolean; error?: string } {
    try {
      const parsed = new URL(url);

      // Allow only HTTPS
      if (parsed.protocol !== 'https:') {
        return { valid: false, error: 'Only HTTPS URLs allowed' };
      }

      // Block localhost/internal IPs
      const hostname = parsed.hostname.toLowerCase();
      if (
        hostname === 'localhost' ||
        hostname === '127.0.0.1' ||
        hostname.startsWith('192.168.') ||
        hostname.startsWith('10.') ||
        hostname === '0.0.0.0'
      ) {
        return { valid: false, error: 'Internal URLs not allowed' };
      }

      return { valid: true };
    } catch {
      return { valid: false, error: 'Invalid URL format' };
    }
  }
}
```

### Step 3: Secure HTTP Requests
```typescript
// src/security/http.ts
import { requestUrl, RequestUrlParam } from 'obsidian';

export class SecureHttpClient {
  private apiKey: string;
  private baseUrl: string;

  constructor(apiKey: string, baseUrl: string) {
    const urlValidation = InputValidator.validateUrl(baseUrl);
    if (!urlValidation.valid) {
      throw new Error(`Invalid base URL: ${urlValidation.error}`);
    }
    this.apiKey = apiKey;
    this.baseUrl = baseUrl;
  }

  async request<T>(
    endpoint: string,
    options: Partial<RequestUrlParam> = {}
  ): Promise<T> {
    // Validate endpoint
    if (endpoint.includes('..') || endpoint.includes('//')) {
      throw new Error('Invalid endpoint');
    }

    const response = await requestUrl({
      url: `${this.baseUrl}${endpoint}`,
      method: options.method || 'GET',
      headers: {
        'Authorization': `Bearer ${this.apiKey}`,
        'Content-Type': 'application/json',
        ...options.headers,
      },
      body: options.body,
      throw: false, // Don't throw, handle errors manually
    });

    // Don't log response body (may contain sensitive data)
    if (response.status >= 400) {
      throw new Error(`HTTP ${response.status}: Request failed`);
    }

    return response.json as T;
  }
}
```

### Step 4: Data Protection
```typescript
// src/security/data-protection.ts
export class DataProtection {
  // Check if file path is in allowed directories
  static isPathAllowed(
    path: string,
    allowedFolders: string[]
  ): boolean {
    return allowedFolders.some(folder =>
      path.startsWith(folder + '/') || path === folder
    );
  }

  // Redact sensitive content before logging
  static redactForLogging(data: any): any {
    const sensitiveKeys = [
      'apiKey', 'api_key', 'token', 'password', 'secret',
      'authorization', 'auth', 'key', 'credential'
    ];

    if (typeof data !== 'object' || data === null) {
      return data;
    }

    const redacted = { ...data };
    for (const key of Object.keys(redacted)) {
      if (sensitiveKeys.some(sk =>
        key.toLowerCase().includes(sk.toLowerCase())
      )) {
        redacted[key] = '[REDACTED]';
      } else if (typeof redacted[key] === 'object') {
        redacted[key] = this.redactForLogging(redacted[key]);
      }
    }
    return redacted;
  }

  // Hash content for comparison without storing original
  static async hashContent(content: string): Promise<string> {
    const encoder = new TextEncoder();
    const data = encoder.encode(content);
    const hashBuffer = await crypto.subtle.digest('SHA-256', data);
    const hashArray = Array.from(new Uint8Array(hashBuffer));
    return hashArray.map(b => b.toString(16).padStart(2, '0')).join('');
  }
}
```

### Step 5: Permission Checks
```typescript
// src/security/permissions.ts
export class PermissionManager {
  private app: App;

  constructor(app: App) {
    this.app = app;
  }

  // Check if user has enabled the required setting
  async requestPermission(
    action: string,
    description: string
  ): Promise<boolean> {
    return new Promise((resolve) => {
      const modal = new ConfirmModal(
        this.app,
        `Allow "${action}"?\n\n${description}`,
        (confirmed) => resolve(confirmed)
      );
      modal.open();
    });
  }

  // Track external data access
  logExternalAccess(
    service: string,
    action: string,
    dataType: string
  ): void {
    // Log for audit purposes (without sensitive data)
    console.log(`[Audit] External access: ${service} - ${action} - ${dataType}`);
  }

  // Check before sending data externally
  async confirmExternalDataShare(
    service: string,
    dataDescription: string
  ): Promise<boolean> {
    return this.requestPermission(
      `Send data to ${service}`,
      `This will send the following data externally:\n${dataDescription}\n\nDo you want to proceed?`
    );
  }
}
```

## Output
- Secure settings storage with masked API keys
- Input validation for paths and URLs
- Safe HTTP client with request validation
- Data redaction for logging
- Permission prompts for sensitive operations

## Error Handling
| Risk | Mitigation |
|------|------------|
| API key exposure | Store in settings, mask in UI, never log |
| Path traversal | Validate all paths, block `..` |
| XSS in views | Sanitize HTML content |
| SSRF | Validate URLs, block internal addresses |
| Data leakage | Confirm before external transmission |

## Examples

### Content Security Policy for Custom Views
```typescript
// When creating custom views with external content
const iframe = document.createElement('iframe');
iframe.sandbox.add('allow-scripts');
iframe.sandbox.add('allow-same-origin');
// Don't add 'allow-top-navigation' or 'allow-forms' unless needed
```

### Secure Error Handling
```typescript
try {
  await riskyOperation();
} catch (error) {
  // Don't expose internal details
  console.error('Operation failed:', DataProtection.redactForLogging(error));

  // Show user-friendly message
  new Notice('Operation failed. Please check your settings.');

  // Don't re-throw with sensitive info
}
```

### Checklist Before Release
- [ ] No hardcoded secrets in code
- [ ] All user inputs validated
- [ ] API keys stored securely in settings
- [ ] External requests use HTTPS only
- [ ] Logging doesn't include sensitive data
- [ ] User consent for external data sharing
- [ ] Error messages don't leak internals

## Resources
- [OWASP Secure Coding Practices](https://owasp.org/www-project-secure-coding-practices-quick-reference-guide/)
- [Obsidian Plugin Guidelines](https://docs.obsidian.md/Plugins/Releasing/Plugin+guidelines)

## Next Steps
For pre-release checklist, see `obsidian-prod-checklist`.
