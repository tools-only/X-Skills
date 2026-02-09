# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/deepgram-pack/skills/deepgram-upgrade-migration/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/deepgram-pack/skills/deepgram-upgrade-migration/SKILL.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `f2544e25771d955e...` |

## Description

// v3 (new)
import { createClient } from '@deepgram/sdk';
const deepgram = createClient(apiKey);
const { result, error } = await deepgram.listen.prerecorded.transcribeUrl(
  { url: audioUrl },
  { punctuate: true }
);

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/deepgram-pack/skills/deepgram-upgrade-migration/SKILL.md)*
