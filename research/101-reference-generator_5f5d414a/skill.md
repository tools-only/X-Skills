# Reference Generator

## Overview

The reference-generator skill generates curated, verified reference lists for educational textbooks with level-appropriate resources. It creates 10-40 references depending on target audience (junior-high to graduate level), with links, publication details, and relevance descriptions.

## Purpose

This skill automates the creation of high-quality, academically appropriate reference lists that enhance textbook credibility and provide students with pathways for deeper learning at their comprehension level.

## Key Features

- **Level-Appropriate Quantity**: 10 (junior-high), 20 (senior-high), 30 (college), 40 (graduate)
- **Verified URLs**: Every link tested with WebFetch before inclusion
- **Publication Details**: Dates, sources, and relevance descriptions
- **Quality Filtering**: Age-appropriate content and academic rigor
- **Two Modes**: Book-level or chapter-level references
- **ISO Format Dates**: Consistent YYYY-MM-DD formatting

## When to Use

Use this skill when:
- Creating a new intelligent textbook that needs comprehensive references
- Adding references to existing textbook
- Updating or expanding reference sections
- A user explicitly requests reference generation

## Reference Quantity by Level

### Junior-High (Middle School) - 10 References
- Educational websites with interactive content
- Videos from reputable educational channels
- Visual resources, infographics, and animations
- Age-appropriate articles from educational publishers
- Museums, science centers, and educational organizations

### Senior-High (High School) - 20 References
- Mix of educational websites and academic sources
- Reputable news organizations and science journalism
- Educational videos and documentaries
- Introduction to academic journals (accessible papers)
- Government and NGO educational resources

### College (Undergraduate) - 30 References
- Peer-reviewed journal articles (50%+ of references)
- Academic textbooks and monographs
- University course materials and lectures
- Research institution publications
- Industry white papers and technical reports

### Graduate (Masters/PhD) - 40 References
- Heavily weighted toward peer-reviewed journals (70%+ of references)
- Seminal papers in the field
- Recent research (last 5 years) showing current state
- Meta-analyses and systematic reviews
- Academic books from university presses

## Workflow Steps

### Step 1: Analyze Course Description
Read `/docs/course-description.md` to determine:
- Grade level or target audience
- Prerequisites (indicates reader sophistication)
- Subject matter (determines reference topics)
- Learning objectives (guides reference selection)

### Step 2: Check for Chapter-Level Content
Search for chapter content:
```bash
find /docs/chapters
find /docs -name "chapter*.md" -o -name "*-chapter-*.md"
```

If chapters exist, ask user: "Book-level or chapter-level references?"

### Step 3: Generate References with Verification
For each reference:
1. **Search** for authoritative sources using WebSearch tool
2. **Verify** each URL using WebFetch to ensure accessibility
3. **Format** according to standard template

### Step 4: Format Each Reference

Standard format:
```markdown
1. [Link Title](URL) - YYYY-MM-DD - Publication Name - Brief description of resource and specific relevance to the textbook topic.
```

**Format Specifications:**
- **Link Title**: Exact title of article/paper/video/resource
- **URL**: Verified, working link
- **Date**: Publication date in YYYY-MM-DD format (or YYYY-MM / YYYY if unavailable)
- **Publication Name**: Journal, website, organization, or publisher
- **Description**: 1-2 sentences explaining content and relevance

### Step 5: Write References to File

**Book-level references:**
Create `/docs/references.md`:

```markdown
# References

This textbook draws upon the following high-quality resources:

[Generated numbered list of references]

---
*References last updated: [Current Date]*
```

**Chapter-level references:**
Append to each chapter file:

```markdown

## References

[Generated numbered list of references for this chapter]
```

### Step 6: Validation and Reporting
1. Count references to ensure correct quantity
2. Verify all URLs were checked with WebFetch
3. Report summary to user with any failed verifications

## Reference Format Examples

```markdown
1. [How Neural Networks Really Work](https://distill.pub/2020/circuits/zoom-in/) - 2020-03-10 - Distill - Interactive visualization explaining the inner workings of neural networks through explorable explanations, perfect for visual learners beginning their ML journey.

2. [Attention Is All You Need](https://arxiv.org/abs/1706.03762) - 2017-06-12 - arXiv - Seminal paper introducing the Transformer architecture that revolutionized natural language processing and forms the foundation for modern LLMs like GPT and BERT.

3. [Khan Academy: Introduction to Algorithms](https://www.khanacademy.org/computing/computer-science/algorithms) - 2024-01-15 - Khan Academy - Free, interactive course covering fundamental algorithms including sorting and searching, with visualizations and practice exercises suitable for high school students.
```

## URL Verification Process

**Critical**: Every URL must be verified before inclusion.

```python
# Use WebFetch for each URL
WebFetch(url=reference_url, prompt="Is this page accessible? Provide the title and a brief description of the content.")
```

If a URL returns an error or redirect:
- Try to find an updated or archived version
- Use Internet Archive / Wayback Machine if appropriate
- Skip the reference if no valid URL exists
- Note in the report any references that couldn't be verified
- For academic papers behind paywalls, reference the citation page
- For academic textbooks, prefer highly-cited works

## Quality Checklist

Before finalizing references, ensure:
- [ ] Correct quantity for target level (10/20/30/40)
- [ ] All URLs verified and accessible
- [ ] Publication dates included
- [ ] Mix of resource types (articles, videos, papers)
- [ ] Descriptions explain relevance to textbook
- [ ] Academic rigor matches target audience
- [ ] No duplicate sources
- [ ] Proper formatting throughout

## Example Usage Scenarios

### Scenario 1: New Textbook
```
User: "Generate references for my textbook"
→ Read /docs/course-description.md
→ Identify level (e.g., college)
→ Check for chapters (none found)
→ Generate 30 verified references
→ Write to /docs/references.md
```

### Scenario 2: Existing Textbook with Chapters
```
User: "Add references to my course"
→ Read /docs/course-description.md
→ Find chapter files exist
→ Ask: "Book-level or chapter-level references?"
→ User selects chapter-level
→ Generate references for each chapter
→ Append to each chapter file
```

## Best Practices

### Source Selection
1. **Verify Authority**: Use established publishers, recognized experts
2. **Check Recency**: Prefer recent sources for rapidly evolving fields
3. **Balance Types**: Mix videos, articles, papers, books
4. **Cross-Reference**: Include multiple perspectives on key topics
5. **Accessibility**: Consider open-access resources when available

### Academic Papers
- For college/graduate levels, prefer Google Scholar citations
- Include seminal papers (highly cited, foundational work)
- Balance classic papers with recent research
- Check if full text is available or just abstract

### Educational Websites
- Verify reputation (edu domains, established organizations)
- Check for regular updates and maintenance
- Ensure mobile-friendly, accessible design
- Prefer interactive or multimedia content

### URLs and Link Rot
- Test all links before inclusion
- Note if archived version used
- Include DOIs for academic papers when available
- Consider adding archive.org links as backup

## Output Summary

After generation, report:
- Number of references generated
- Target level identified
- File location (book or chapter-level)
- Any URLs that failed verification
- Suggestion to use citation graph skill for academic papers

## Integration with Other Skills

- **course-description-analyzer**: Determines appropriate reference level and topics
- **chapter-content-generator**: References support chapter content
- **glossary-generator**: Reference definitions align with glossary
- **learning-graph-generator**: References support concept dependencies

## Tools Used

- **WebSearch**: Find authoritative sources on topics
- **WebFetch**: Verify URLs are accessible and extract metadata
- **AskUserQuestion**: Clarify book-level vs chapter-level preference

## Advanced Features

### Citation Graph Analysis
After generating references, suggest:
- Use citation graph skill to find highly-cited papers
- Identify influential works in the field
- Discover seminal papers that shape the domain

### Multiple Formats
Generate references in various formats:
- Markdown (default)
- BibTeX for LaTeX integration
- RIS for reference managers
- JSON for programmatic access

## References

- **WebSearch and WebFetch**: Built-in Claude Code tools
- **Dublin Core**: Metadata standards for resources
- **Academic Citation Standards**: MLA, APA, Chicago styles
