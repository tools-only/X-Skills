# 3B Rust Compiler Analyzer

| Property | Value |
|----------|-------|
| **Name** | 3B Rust Compiler Analyzer |
| **Repository** | [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/3b-rust-compiler-analyzer.md) (⭐ 2.9k) |
| **Original Path** | `plugins/zeroize-audit/agents/3b-rust-compiler-analyzer.md` |
| **Category** | data-analysis |
| **Subcategory** | processing |
| **Tags** | data analysis |
| **Created** | 2026-02-26 |
| **Updated** | 2026-02-26 |
| **File Hash** | `f2a7fd92f233d4d7...` |

## Description

Performs crate-level MIR and LLVM IR analysis for Rust in zeroize-audit. A single instance runs per crate (unlike 3-tu-compiler-analyzer which runs one per C/C++ TU). Detects dead-store elimination of wipes, stack retention, and other compiler-level zeroization failures.

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [trailofbits/skills](https://raw.githubusercontent.com/trailofbits/skills/main/plugins/zeroize-audit/agents/3b-rust-compiler-analyzer.md)*
