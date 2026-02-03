# Authentication & Guards

| Property | Value |
|----------|-------|
| **Name** | Authentication & Guards |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/nestjs-expert/references/authentication.md) (⭐ 216) |
| **Original Path** | `skills/nestjs-expert/references/authentication.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `7bdd11ae675c17f2...` |

## Description

typescript
// jwt.strategy.ts
import { Injectable } from '@nestjs/common';
import { PassportStrategy } from '@nestjs/passport';
import { ExtractJwt, Strategy } from 'passportjwt';
import { ConfigService } from '@nestjs/config';

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/nestjs-expert/references/authentication.md)*
