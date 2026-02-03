# Docker Patterns

| Property | Value |
|----------|-------|
| **Name** | Docker Patterns |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/docker-patterns.md) (⭐ 216) |
| **Original Path** | `skills/devops-engineer/references/docker-patterns.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2025-12-14 |
| **Updated** | 2026-01-29 |
| **File Hash** | `b6ccc98c44aa64e9...` |

## Description

dockerfile
 Build stage
FROM node:20alpine AS builder
WORKDIR /app
COPY package.json ./
RUN npm ci only=production && npm cache clean force
COPY . .
RUN npm run build

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/devops-engineer/references/docker-patterns.md)*
