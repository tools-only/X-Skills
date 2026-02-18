# Integrating Agent Skills into Your Agent

> Source: <https://agentskills.io/integrate-skills.md>

**Use this guide when** building or modifying an AI agent or development tool to support the Agent Skills format — discovering skills, loading metadata, injecting into context, and handling script execution safely.

This guide explains how to add skills support to an AI agent or development tool.

---

## Contents

- [Integration approaches](#integration-approaches)
- [Overview](#overview)
- [Skill discovery](#skill-discovery)
- [Loading metadata](#loading-metadata)
- [Injecting into context](#injecting-into-context)
- [Security considerations](#security-considerations)
- [Reference implementation](#reference-implementation)

---

## Integration Approaches

**Filesystem-based agents** operate within a computer environment (bash/unix). Skills are activated when models issue shell commands like `cat /path/to/my-skill/SKILL.md`. Bundled resources are accessed through shell commands. This is the most capable option.

**Tool-based agents** function without a dedicated computer environment. They implement tools allowing models to trigger skills and access bundled assets. The specific tool implementation is up to the developer.

---

## Overview

A skills-compatible agent needs to:

1. **Discover** skills in configured directories
2. **Load metadata** (name and description) at startup
3. **Match** user tasks to relevant skills
4. **Activate** skills by loading full instructions
5. **Execute** scripts and access resources as needed

---

## Skill Discovery

Skills are folders containing a `SKILL.md` file. The agent should scan configured directories for valid skills.

---

## Loading Metadata

At startup, parse only the frontmatter of each `SKILL.md` file. This keeps initial context usage low.

### Parsing Frontmatter

```
function parseMetadata(skillPath):
    content = readFile(skillPath + "/SKILL.md")
    frontmatter = extractYAMLFrontmatter(content)

    return {
        name: frontmatter.name,
        description: frontmatter.description,
        path: skillPath
    }
```

---

## Injecting into Context

Include skill metadata in the system prompt so the model knows what skills are available.

For Claude models, the recommended format uses XML:

```xml
<available_skills>
  <skill>
    <name>pdf-processing</name>
    <description>Extracts text and tables from PDF files, fills forms, merges documents.</description>
    <location>/path/to/skills/pdf-processing/SKILL.md</location>
  </skill>
  <skill>
    <name>data-analysis</name>
    <description>Analyzes datasets, generates charts, and creates summary reports.</description>
    <location>/path/to/skills/data-analysis/SKILL.md</location>
  </skill>
</available_skills>
```

**For filesystem-based agents:** Include the `location` field with the absolute path to the SKILL.md file.

**For tool-based agents:** The location can be omitted.

Keep metadata concise. Each skill should add roughly 50-100 tokens to the context. **Reason:** This keeps the system prompt within typical limits for skill discovery; excessive metadata can cause truncation and prevent skills from being auto-invoked.

---

## Security Considerations

Script execution requires safeguards to remain safe:

- **Sandboxing** — Run scripts in isolated environments
- **Allowlisting** — Execute scripts only from trusted skills
- **Confirmation** — Ask users before running potentially dangerous operations
- **Logging** — Record all script executions for auditing

---

## Reference Implementation

The [skills-ref](https://github.com/agentskills/agentskills/tree/main/skills-ref) library provides Python utilities and a CLI for working with skills.

**Validate a skill directory:**

```bash
skills-ref validate <path>
```

**Generate `<available_skills>` XML for agent prompts:**

```bash
skills-ref to-prompt <path>...
```

Use the library source code as a reference implementation for building your own skills integration.

---

## Adoption

Agent Skills are supported by 25+ products including:

Claude Code, Cursor, VS Code, Gemini CLI, OpenAI Codex, GitHub, Roo Code, Amp, Goose, Factory, Databricks, Spring AI, Letta, Firebender, OpenCode, Autohand, Mux, Piebald, TRAE, Qodo, Agentman, Mistral Vibe, Command Code, Ona, VT Code.
