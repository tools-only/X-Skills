# Research Agent

**Version:** 1.0.0
**Generated:** 2026-02-05T03:33:50.183751
**Skills:** 30

## Description

Academic research and literature review skills

## Use Case

Agents for academic research, literature search, paper writing, and citations

## Contents

This patch includes **30 skills** organized into **1 categories**.


### Research (30 skills)

| Skill | Subcategory | Tags | Status |
|-------|-------------|------|--------|
| [Pubmed Search](../research/094-searching_f25e7adf) | - | `research` | **Required** |
| [Google Scholar Search](../research/053-searching_494df9be) | - | `research` | **Required** |
| [Literature Search Strategies](../research/592-literature_search_strategies_ff637b4c) | - | `research` | **Required** |
| [Imrad Structure](../research/594-imrad_structure_b70f52f1) | - | `research` | **Required** |
| [Citation Validation](../research/014-validating_37e61db8) | - | `research` | **Required** |
| [Citation Styles](../research/013-document_a62f97c9) | - | `research` | **Required** |
| [Peer Review Standards](../research/083-peer_06bb68f3) | - | `research` | **Required** |
| [Reviewer Expectations](../research/108-understanding_58ecad55) | - | `research` | **Required** |
| [Evidence Hierarchy](../research/044-description_e9e2f27c) | - | `research` | **Required** |
| [Scientific Method](../research/110-knowledge_09e709bd) | - | `research` | **Required** |
| [Poster Quality Checklist](../research/088-use_1f2d89a6) | - | `research` | **Required** |
| [Posters Guidelines](../research/089-guidelines_6a853d81) | - | `research` | **Required** |
| [Feature Researcher](../research/698-feature-researcher_644ab486) | - | `research` | **Required** |
| [References Generator Skill](../research/104-prompt_823350ee) | - | `research` | **Required** |
| [Ml Conference Style](../research/068-writing_616d4ec8) | - | `research` | **Required** |
| [Nature Science Style](../research/070-writing_cf177f2c) | - | `research` | **Required** |
| [Cell Press Style](../research/011-writing_6f7196c8) | - | `research` | **Required** |
| [Writing Principles](../research/133-effective_c786f9fe) | - | `research` | **Required** |
| [Search Strategies](../research/111-best_138cc4c4) | - | `research` | **Required** |
| [Source Evaluation](../research/1254-source-evaluation_9bc285fb) | - | `research` | **Required** |
| [Search Strategies](../research/1253-search-strategies_6f14aba6) | - | `research` | **Required** |
| [Skill](../research/115-animation-performance-retro_e719e45d) | - | `research` | Optional |
| [Skill](../research/115-anthropologist-analyst_45fc3f9d) | - | `research` | Optional |
| [Pdf Processing Limits Analysis](../research/1082-pdf_processing_limits_analysis_2f948f56) | - | `research` | Optional |
| [Skill](../research/591-apply-skill_80963683) | - | `research` | Optional |
| [Skill](../research/003-name-skill_59eaadb9) | - | `research` | Optional |
| [Confidence And Limitations](../research/022-document_143fc8b3) | - | `research` | Optional |
| [Stateless Software Engineering Framework](../research/693-stateless-software-engineering-framework_f6b12de7) | - | `research` | Optional |
| [Skill](../research/115-chapter-content-generator_87051c46) | - | `data analysis` | Optional |
| [Common Queries](../research/018-document_95e9c3ef) | - | `research` | Optional |


## Installation

To install this patch to your Claude Code skills directory:

```bash
python -m src.patch_installer install research-agent
```

Or using the skillflow CLI (if available):

```bash
skillflow patch install research-agent
```

## Dependencies

None

## Manifest

```json
{
  "spec": {
    "id": "research-agent",
    "name": "Research Agent",
    "description": "Academic research and literature review skills",
    "use_case": "Agents for academic research, literature search, paper writing, and citations",
    "categories": [
      "research"
    ],
    "subcategories": [],
    "tags": [],
    "exclude": [
      "*drug*",
      "*clin*",
      "*bio*",
      "*hmdb*",
      "*cosmic*",
      "*alma*",
      "*protoco*",
      "*database*",
      "*databases*",
      "*schema*",
      "*fields*",
      "*queries*"
    ],
    "min_stars": null,
    "max_skills": 30,
    "required_skills": [
      "research/094-searching_f25e7adf",
      "research/053-searching_494df9be",
      "research/592-literature_search_strategies_ff637b4c",
      "research/594-imrad_structure_b70f52f1",
      "research/014-validating_37e61db8",
      "research/013-document_a62f97c9",
      "research/083-peer_06bb68f3",
      "research/108-understanding_58ecad55",
      "research/044-description_e9e2f27c",
      "research/110-knowledge_09e709bd",
      "research/088-use_1f2d89a6",
      "research/089-guidelines_6a853d81",
      "research/698-feature-researcher_644ab486",
      "research/104-prompt_823350ee",
      "research/068-writing_616d4ec8",
      "research/070-writing_cf177f2c",
      "research/011-writing_6f7196c8",
      "research/133-effective_c786f9fe",
      "research/111-best_138cc4c4",
      "research/1254-source-evaluation_9bc285fb",
      "research/1253-search-strategies_6f14aba6"
    ],
    "optional_skills": [],
    "dependencies": [],
    "version": "1.0.0"
  },
  "skills": [
    {
      "local_path": "research/094-searching_f25e7adf",
      "display_name": "Pubmed Search",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/053-searching_494df9be",
      "display_name": "Google Scholar Search",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/592-literature_search_strategies_ff637b4c",
      "display_name": "Literature Search Strategies",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/594-imrad_structure_b70f52f1",
      "display_name": "Imrad Structure",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/014-validating_37e61db8",
      "display_name": "Citation Validation",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/013-document_a62f97c9",
      "display_name": "Citation Styles",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/083-peer_06bb68f3",
      "display_name": "Peer Review Standards",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/108-understanding_58ecad55",
      "display_name": "Reviewer Expectations",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/044-description_e9e2f27c",
      "display_name": "Evidence Hierarchy",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/110-knowledge_09e709bd",
      "display_name": "Scientific Method",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/088-use_1f2d89a6",
      "display_name": "Poster Quality Checklist",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/089-guidelines_6a853d81",
      "display_name": "Posters Guidelines",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/698-feature-researcher_644ab486",
      "display_name": "Feature Researcher",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/104-prompt_823350ee",
      "display_name": "References Generator Skill",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/068-writing_616d4ec8",
      "display_name": "Ml Conference Style",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/070-writing_cf177f2c",
      "display_name": "Nature Science Style",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/011-writing_6f7196c8",
      "display_name": "Cell Press Style",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/133-effective_c786f9fe",
      "display_name": "Writing Principles",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/111-best_138cc4c4",
      "display_name": "Search Strategies",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/1254-source-evaluation_9bc285fb",
      "display_name": "Source Evaluation",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/1253-search-strategies_6f14aba6",
      "display_name": "Search Strategies",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": true,
      "reason": "Required skill"
    },
    {
      "local_path": "research/115-animation-performance-retro_e719e45d",
      "display_name": "Skill",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/115-anthropologist-analyst_45fc3f9d",
      "display_name": "Skill",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/1082-pdf_processing_limits_analysis_2f948f56",
      "display_name": "Pdf Processing Limits Analysis",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/591-apply-skill_80963683",
      "display_name": "Skill",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/003-name-skill_59eaadb9",
      "display_name": "Skill",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/022-document_143fc8b3",
      "display_name": "Confidence And Limitations",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/693-stateless-software-engineering-framework_f6b12de7",
      "display_name": "Stateless Software Engineering Framework",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/115-chapter-content-generator_87051c46",
      "display_name": "Skill",
      "category": "research",
      "subcategory": "",
      "tags": [
        "data analysis"
      ],
      "required": false,
      "reason": "Matches filters"
    },
    {
      "local_path": "research/018-document_95e9c3ef",
      "display_name": "Common Queries",
      "category": "research",
      "subcategory": "",
      "tags": [
        "research"
      ],
      "required": false,
      "reason": "Matches filters"
    }
  ],
  "total_count": 30,
  "generated_at": "2026-02-05T03:33:50.183751",
  "checksum": "d6900ea6dc04a440421bee1ecaeef6bad6c76173e77f1daacfb9f9c65f593a25"
}
```

---

*Generated by [SkillFlow](https://github.com/tools-only/SkillFlow)*
