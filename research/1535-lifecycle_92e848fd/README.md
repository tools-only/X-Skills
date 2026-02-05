# Lifecycle

| Property | Value |
|----------|-------|
| **Name** | Lifecycle |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/plugins/beagle-elixir/skills/liveview-code-review/references/lifecycle.md) (⭐ 17) |
| **Original Path** | `plugins/beagle-elixir/skills/liveview-code-review/references/lifecycle.md` |
| **Category** | research |
| **Subcategory** | data-gathering |
| **Tags** | research |
| **Created** | 2026-02-05 |
| **Updated** | 2026-02-05 |
| **File Hash** | `92e848fdb254f05d...` |

## Description

elixir
def mount(_params, _session, socket) do
   Only subscribe when actually connected (not during static render)
  if connected?(socket) do
    Phoenix.PubSub.subscribe(MyApp.PubSub, "updates")
  end

**Tags:** `research`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/plugins/beagle-elixir/skills/liveview-code-review/references/lifecycle.md)*
