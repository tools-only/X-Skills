# Provenance Report: “Methodology Comparison Template” Generation (2026-01-27)

This document records the **actual process executed** to produce:

- Output template file: `methodology_development/methodology_comparison_template.md`
- Supporting research summary: delivered in-chat as “Stage 1: Research & Information Gathering”

Goal: produce a reusable, universal comparison template grounded in **real-world comparison methodologies** (consumer, technical, business, academic), with an **AI-oriented pre-comparison reflection** and **bounded self-extension**.

---

## 1) Local repo context gathered (inputs read)

To match this repository’s conventions and constraints, I read:

- `sessions/CLAUDE.sessions.md` (cc-sessions workflow expectations; DAIC notes)
- `methodology_development/README.md` (how SAM docs are organized)
- `methodology_development/stateless-agent-methodology.md` (existing comparison style: matrices, stages, verification language)
- `methodology_development/sam-vs-get-shit-done.md` (comparison document structure and tone)

These reads informed the template’s emphasis on:

- gating steps (reflection before action),
- explicit failure modes / backpressure,
- structured outputs and checklists.

---

## 2) Web research execution log (what I actually did)

### 2.1 Initial broad sampling (search)

I ran broad web searches spanning:

- consumer review methodologies (Wirecutter, Consumer Reports, RTINGS, PCMag, Good Housekeeping),
- decision matrix frameworks (ASQ, Pugh),
- business methodology comparisons (Agile vs Waterfall, Lean vs Six Sigma),
- academic rubrics (AAC&U VALUE, NIH peer review),
- technical framework comparisons (React/Vue, Vue official comparisons).

Purpose: ensure the resulting domains weren’t biased toward software-only comparisons.

### 2.2 Primary-source retrieval (fetch)

I attempted direct retrieval via standard page fetch for the candidate primary sources.

Observed constraints:

- Some Wirecutter/NYT pages were blocked to direct fetch (403), so direct fetch was insufficient.
- Some RTINGS pages returned partial content with standard fetch (likely JS/content gating).

### 2.3 MCP tool selection decision (why Firecrawl MCP was used)

You requested MCP tools + Chrome browser when scraping fails.

I examined the available MCP servers in this workspace:

- `cursor-ide-browser` and `cursor-browser-extension` are present as MCP servers, but **do not provide local tool schema JSON files** in `./mcps/.../tools/` for me to read.
- This repo enforces “read schema first before calling MCP tools”; without tool schemas, I could not safely invoke those MCP tools via `CallMcpTool`.

Therefore, for blocked/JS-heavy sources, I switched to **Firecrawl MCP**, which _does_ provide tool schema JSON descriptors locally and can be invoked compliantly.

Tool schemas read (for compliance):

- `mcps/user-firecrawl/tools/firecrawl_scrape.json`
- `mcps/user-firecrawl/tools/firecrawl_search.json`
- `mcps/user-firecrawl/tools/firecrawl_extract.json`
- `mcps/user-firecrawl/tools/firecrawl_map.json`

### 2.4 MCP retrieval results (what succeeded / failed)

Succeeded via Firecrawl scrape:

- Wirecutter “Anatomy of a guide”: `https://www.nytimes.com/wirecutter/blog/anatomy-of-a-guide/`
- Wirecutter “Comparison tables”: `https://www.nytimes.com/wirecutter/blog/comparison-tables/`
- RTINGS “Test Benches and Scoring System”: `https://www.rtings.com/company/test-benches-and-scoring-system`
- RTINGS “Scoring Function Overhaul”: `https://www.rtings.com/company/scoring-function-overhaul`

Failed via Firecrawl (provider blocklist):

- Wirecutter washer/dryer/detergent testing article: `https://www.nytimes.com/wirecutter/reviews/testing-washers-dryers-and-detergents/`
  - Failure reason: Firecrawl reported the site/page as **blocklisted** (no content retrieved).
  - Mitigation: used other consumer-testing methodology sources (Consumer Reports, PCMag, Good Housekeeping) to cover the same “testing protocol + scoring” pattern without fabricating details.

Notes:

- Some attempted WebFetch requests were rejected due to the user choosing to skip those fetches (recorded in tool output at the time). Those sources were therefore **not used** as primary evidence.

---

## 3) Synthesis decisions (how research shaped the template)

The final template structure in `methodology_development/methodology_comparison_template.md` was directly driven by recurring patterns observed in the retrieved sources:

- **Structured scope framing** (“who it’s for”, “why trust us”) from Wirecutter.
- **Comparison tables / weighted attributes** from Wirecutter.
- **Category-dependent weighting + disqualifiers** (e.g., CR’s “cannot recommend if…”) from Consumer Reports.
- **Comparability constraints** (only compare within same script/protocol) from PCMag.
- **Weighted ratings as usage profiles** + explicit score interpretation bands from RTINGS.
- **Criteria selection + weighting + baseline scoring** + “keep scales aligned” from ASQ decision matrix guidance.
- **Context section separated from evaluation** and explicit rubric levels from Baldrige.
- **Mixed scoring types** (numeric + sufficiency/binary) and “not all criteria must be strong” from NIH peer review framework.
- **Bias acknowledgement** in stakeholder-authored comparisons from Vue’s comparison guide.

I then encoded these into:

- a mandatory **Pre-Comparison Reflection** gate,
- essential vs optional domains,
- bounded “domain-specific extension” with a hard cap,
- explicit stopping criteria,
- optional decision matrix with sensitivity check.

---

## 4) Outputs produced (filesystem provenance)

- Created: `methodology_development/methodology_comparison_template.md`

  - Contents include:
    - pre-comparison reflection gate,
    - universal domain worksheets,
    - optional decision matrix + sensitivity check,
    - evidence log for traceability.

- Created/updated: this file
  - `.claude/skills/research-and-compare/process-and-comparison-methodology-with-references.md`
  - This replaces a generic “how to do comparisons” spec with the concrete provenance record.

---

## 5) Reference set actually used (primary sources; access date 2026-01-27)

Consumer testing / review methodology:

1. Wirecutter, “The Anatomy of a Wirecutter Guide”. `https://www.nytimes.com/wirecutter/blog/anatomy-of-a-guide/` (accessed 2026-01-27)
2. Wirecutter, “Wirecutter’s Secret to Making Great Picks: Obsessive Spreadsheeting”. `https://www.nytimes.com/wirecutter/blog/comparison-tables/` (accessed 2026-01-27)
3. Consumer Reports Data Intelligence, “Rating Methods”. `https://data.consumerreports.org/rating-methods/` (accessed 2026-01-27)
4. PCMag, “How We Test Everything We Review”. `https://www.pcmag.com/about/how-we-test-everything-we-review` (accessed 2026-01-27)
5. Good Housekeeping (UK), “How the Good Housekeeping Institute Approved seal really works”. `https://www.goodhousekeeping.com/uk/consumer-advice/a563956/how-good-housekeeping-institute-approved-really-works/` (accessed 2026-01-27)
6. RTINGS, “Test Benches and Scoring System”. `https://www.rtings.com/company/test-benches-and-scoring-system` (accessed 2026-01-27)
7. RTINGS, “Scoring Function Overhaul”. `https://www.rtings.com/company/scoring-function-overhaul` (accessed 2026-01-27)

Decision matrices / selection frameworks: 8. ASQ, “What is a Decision Matrix?”. `https://asq.org/quality-resources/decision-matrix` (accessed 2026-01-27)

Business methodology comparisons / selection framing: 9. Atlassian, “Project management intro: Agile vs. waterfall methodologies”. `https://www.atlassian.com/agile/project-management/project-management-intro` (accessed 2026-01-27) 10. Lucidchart, “Lean vs. Six Sigma: Determining the Right Method for Your Business”. `https://lucidchart.com/blog/lean-vs-six-sigma` (accessed 2026-01-27)

Formal rubrics / peer review: 11. NIST Baldrige, “Award Criteria” (2025 edition page). `https://www.nist.gov/baldrige/baldrige-award/award-criteria` (accessed 2026-01-27) 12. NIH Grants, “Simplified Peer Review Framework”. `https://grants.nih.gov/policy-and-compliance/policy-topics/peer-review/simplifying-review/framework` (accessed 2026-01-27)

Technical framework self-comparison (bias-disclosure example): 13. Vue.js (Vue 2 docs), “Comparison with Other Frameworks”. `https://vuejs.org/v2/guide/comparison.html` (accessed 2026-01-27)

---

## 6) Known gaps / blocked evidence

- Wirecutter washer/dryer/detergent testing story was **not retrievable** via Firecrawl due to provider blocklisting:
  - `https://www.nytimes.com/wirecutter/reviews/testing-washers-dryers-and-detergents/`
  - Status: not used as evidence; replaced by other consumer-testing methodology sources listed above.
