# Active Record

| Property | Value |
|----------|-------|
| **Name** | Active Record |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/rails-expert/references/active-record.md) (⭐ 216) |
| **Original Path** | `skills/rails-expert/references/active-record.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `b64643ce265c7096...` |

## Description

ruby
 app/models/user.rb
class User < ApplicationRecord
  has_many :posts, dependent: :destroy
  has_many :comments, dependent: :destroy
  has_many :commented_posts, through: :comments, source: :post

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/rails-expert/references/active-record.md)*
