---
argument-hint: "provide the full text to be processed, or a list of files to be processed."
description: "Text to LLM Optimized Text workflow instructions"
disable-model-invocation: false
---

**Situation**

The assistant is operating in a text transformation workflow where source material must be converted into machine-optimized instruction sets. The target consumer is an LLM or AI agent with expert-level domain knowledge that requires concise, deterministic rules for execution, and reference materials, tables, instructions, documentation for verifying details against.

**Task**

The assistant must transform the provided text into one or more structured markdown files containing technical instructions, reference material, and rules optimized for LLM consumption. Each file must follow a strict format with YAML frontmatter and ASCII-only content structured as decision triggers and pattern-matching rules.

**Objective**

Produce deterministic, high-density instruction sets that enable an expert-level AI to efficiently parse, apply, and execute rules without ambiguity. The output must maximize information density while maintaining clarity through structured XML-tagged sections and explicit rule hierarchies.

**Knowledge**

The source text or list of files to process will be provided between XML `<text>` tokens:

<text>
$ARGUMENTS
</text>

The assistant must assume the consuming AI possesses comprehensive technical knowledge of all domain concepts, internals, and industry-specific terminology. Content should be written as lookup references and decision triggers rather than explanatory material.

XML tags provide critical parsing benefits:

- `<instructions>`, `<example>`, `<formatting>` tags prevent the LLM from conflating different prompt components
- Tags enable modular prompt construction and easier extraction of specific response elements
- Structured separation improves accuracy and reduces interpretation errors

File naming convention: all-lower-case-with-hyphens.md

YAML frontmatter structure: <frontmatter_example>

---

name: [skill name or file identifier] description: The model must use this [SKILL/RULE/INSTRUCTION] when: ACTION when TRIGGER to OUTCOME. [Additional trigger examples within 600 characters] version: "1.0.0" last_updated: "<ISO 8601 date and time>" [optional_metadata_key]: [value for disambiguation or context]

---

</frontmatter_example>

Standard glob patterns must be unquoted (e.g., _.js, src/\*\*/_.{ts,js})

**Instructions**

The assistant must apply the following transformation rules in order:

1. The assistant must begin content with a directive instructing how to read and apply all subsequent rules.

2. The assistant must maximize information density by employing technical jargon, domain-specific terminology, dense vocabulary, equations, and industry-standard nomenclature appropriate to the source material's subject matter.

3. The assistant must rephrase all content for maximum accuracy and specificity, eliminating vague guidance and ambiguous phrasing.

4. The assistant must write using declarative phrasing with "The model must" construction for all imperative instructions.

5. The assistant must produce flat ASCII text without markdown formatting for bold, italic, or stylistic emphasis. Headings for hierarchy and lists are permitted.

6. The assistant must structure content into explicit sections: Table of Contents (when applicable), and References (when applicable), identity, intent, task rules, issue handling, triggers.

7. The assistant must preserve or expand any structured examples present in source text, wrapping them in appropriate XML tags (`<example>`, `<positive_example>`, `<negative_example>`).

8. The assistant must use XML tags (`<instructions>`, `<context>`, `<constraints>`, `<formatting>`, `<example>`) to separate distinct prompt components and improve parseability.

9. The assistant must format rules as ACTION TRIGGER OUTCOME patterns optimized for decision-making rather than education.

10. The assistant must establish clear priority levels between conflicting rules and instructions.

11. The assistant must target under 500 lines per file, splitting significant concepts into multiple composable rules when necessary.

12. The assistant must preserve any output structure specifications described in the original text.

13. The assistant must omit greetings, preambles, unnecessary explanatory text, and all markdown styling.

14. The assistant must use only visible ASCII characters in output.

15. The assistant must optimize for AI context window efficiency by removing redundant information and non-essential content.

16. The assistant must populate frontmatter descriptions with comprehensive TRIGGERS that specify when rules should be applied, maintaining clear intent for rule selection.

17. The assistant must limit examples to essential patterns only, providing both positive and negative examples when helpful for rule application.

18. The assistant must output content in the following format:

## <file_content>

name: [name] description: The model must use this [type] when: [trigger conditions with specific examples, max 600 chars] version: "1.0.0" last_updated: "[ISO 8601 timestamp]" [optional_metadata]: [value]

---

[Transformed content following all rules above] </file_content>

19. The assistant must create additional files when source material contains multiple distinct skill domains or rule categories that exceed 500 lines when combined.

20. The assistant must ensure all instructions are focused, imperative, actionable, and properly scoped to their domain.
