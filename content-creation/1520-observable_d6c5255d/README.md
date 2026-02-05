# Observable

| Property | Value |
|----------|-------|
| **Name** | Observable |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/swift-code-review/references/observable.md) (⭐ 17) |
| **Original Path** | `skills/swift-code-review/references/observable.md` |
| **Category** | content-creation |
| **Subcategory** | editing |
| **Tags** | content creation |
| **Created** | 2026-01-11 |
| **Updated** | 2026-01-11 |
| **File Hash** | `d6c5255dea7e0c2c...` |

## Description

swift
// BAD  model recreated on view reconstruction
struct ContentView: View {
    @State private var viewModel = ExpensiveViewModel()  // init() runs repeatedly
}

**Tags:** `content creation`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/swift-code-review/references/observable.md)*
