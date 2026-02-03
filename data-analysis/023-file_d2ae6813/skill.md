# Book Metrics

This file contains overall metrics for the intelligent textbook.

| Metric Name | Value | Link | Notes |
|-------------|-------|------|-------|
| Chapters | 0 | [Chapters](../chapters/index.md) | Number of chapter directories |
| Concepts | 0 | [Concept List](./concept-list.md) | Concepts from learning graph |
| Glossary Terms | 0 | [Glossary](../glossary.md) | Defined terms |
| FAQs | 0 | [FAQ](../faq.md) | Frequently asked questions |
| Quiz Questions | 0 | - | Questions across all chapters |
| Diagrams | 78 | - | Level 4 headers starting with '#### Diagram:' |
| Equations | 268 | - | LaTeX expressions (inline and display) |
| MicroSims | 0 | [Simulations](../sims/index.md) | Interactive MicroSims |
| Total Words | 357,571 | - | Words in all markdown files |
| Links | 1268 | - | Hyperlinks in markdown format |
| Equivalent Pages | 1449 | - | Estimated pages (250 words/page + visuals) |

## Metrics Explanation

- **Chapters**: Count of chapter directories containing index.md files
- **Concepts**: Number of rows in learning-graph.csv
- **Glossary Terms**: H4 headers in glossary.md
- **FAQs**: H3 headers in faq.md
- **Quiz Questions**: H4 headers with numbered questions (e.g., '#### 1.') or H2 headers in quiz.md files
- **Diagrams**: H4 headers starting with '#### Diagram:'
- **Equations**: LaTeX expressions using $ and $$ delimiters
- **MicroSims**: Directories in docs/sims/ with index.md files
- **Total Words**: All words in markdown files (excluding code blocks and URLs)
- **Links**: Markdown-formatted links `[text](url)`
- **Equivalent Pages**: Based on 250 words/page + 0.25 page/diagram + 0.5 page/MicroSim
