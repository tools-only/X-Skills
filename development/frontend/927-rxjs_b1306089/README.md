# RxJS Patterns

| Property | Value |
|----------|-------|
| **Name** | RxJS Patterns |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/angular-architect/references/rxjs.md) (⭐ 216) |
| **Original Path** | `skills/angular-architect/references/rxjs.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `b1306089e570e555...` |

## Description

typescript
import { Component, inject, signal } from '@angular/core';
import {
  map, filter, switchMap, catchError,
  debounceTime, distinctUntilChanged,
  tap, shareReplay, takeUntil
} from 'rxjs/operators';
import { Subject, of, EMPTY } from 'rxjs';

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/angular-architect/references/rxjs.md)*
