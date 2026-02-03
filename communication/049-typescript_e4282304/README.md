# Dtos Validation

| Property | Value |
|----------|-------|
| **Name** | Dtos Validation |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/nestjs-expert/references/dtos-validation.md) (⭐ 216) |
| **Original Path** | `skills/nestjs-expert/references/dtos-validation.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `e42823044324f592...` |

## Description

export class CreateUserDto {
  @ApiProperty({ example: 'user@example.com' })
  @IsEmail()
  email: string;

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/nestjs-expert/references/dtos-validation.md)*
