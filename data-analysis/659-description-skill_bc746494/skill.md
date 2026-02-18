---
description: Summarize files by reading content, extracting key passages, and applying type-specific strategies. Activates on summarize this file, what's in this file, describe this codebase, file summary, analyze this file, tl;dr this file, what does this code do, explain this config, break down this script. Routes to strategies for code, config, data, documentation, markup, and binary files based on extension and word count.
---
# File Summarization

Apply this methodology when summarizing files of any type. This skill provides the routing logic and type-specific strategies for faithful file summarization.

## Pre-Summarization Assessment

Before summarizing any file, the model MUST:

1. **Read the file** - Use the Read tool to access the actual content. Never guess from the filename.

2. **Assess size** - Run `$CLAUDE_PLUGIN_ROOT/scripts/file_metrics.py` to determine word count and file type. If the script is unavailable, use the Read tool and manually estimate word count from line count.

3. **Select strategy** - Based on size thresholds from the table below.

4. **Verify file type** - Use file extension and content inspection to determine which type-specific strategy to apply.

## Size-Based Strategy Selection

| File Size | Strategy | Approach |
|-----------|----------|----------|
| Small (< 2,000 words) | Full read with extractive summarization | Read entire file, extract key passages, summarize from extracts |
| Medium (2,000-10,000 words) | Section-based extraction | Read full file, identify sections/modules, extract from each section, synthesize |
| Large (> 10,000 words) | Chunk and map-reduce | Split into chunks, summarize each chunk, synthesize chunk summaries |

SOURCE: Size thresholds adapted from Anthropic knowledge-synthesis skill (knowledge-work-plugins repository, accessed 2026-02-06). Strategy patterns informed by Map-Reduce Summarization methodology.

## File Type Strategies

### Code Files

File extensions: `.py`, `.js`, `.ts`, `.jsx`, `.tsx`, `.rs`, `.go`, `.java`, `.c`, `.cpp`, `.h`, `.rb`, `.php`, `.swift`, `.kt`, `.scala`, `.sh`, `.bash`, `.zsh`

The model MUST extract:

1. **Imports/dependencies** - List external modules and standard library imports
2. **Structure** - Classes, functions, methods with signatures
3. **Purpose** - Inferred from docstrings, comments, function names
4. **Key logic** - Core algorithms, state machines, data transformations
5. **Entry points** - `main()`, CLI argument parsing, exported functions
6. **Configuration** - Environment variables, config file references

**Extraction method**: Read sequentially. Capture top-level definitions with their line numbers. Extract docstrings verbatim. Quote complex logic rather than paraphrasing.

**Example summary structure**:

```text
## Summary

Python module for HTTP client authentication. Implements JWT token refresh flow with retry logic. Exports `AuthClient` class and `refresh_token()` function.

## What Was Found

- Class `AuthClient` (lines 15-87): JWT-based HTTP client with automatic token refresh
- Function `refresh_token()` (lines 92-105): Retries up to 3 times on 401 errors
- Dependencies: `httpx`, `jwt`, `tenacity` (lines 1-3)
- Environment variables: `AUTH_BASE_URL`, `AUTH_CLIENT_ID` (lines 10-11)

## What Was NOT Found

- No test coverage information in this file
- No error handling for network failures
- Configuration schema not documented
```

### Configuration Files

File extensions: `.json`, `.yaml`, `.yml`, `.toml`, `.ini`, `.env`, `.conf`, `.cfg`, `.properties`

The model MUST extract:

1. **Top-level keys** - All root keys with their value types
2. **Nested structure** - Hierarchy depth and organization
3. **Settings categories** - Group keys by purpose if clear
4. **Notable values** - Endpoints, file paths, feature flags, credentials (note presence, do not expose values)
5. **Validation constraints** - Type requirements, enums, ranges if documented

**Extraction method**: Parse structure. For small files, include all keys. For large files, sample representative sections and note structure patterns.

**Example summary structure**:

```text
## Summary

Application configuration in YAML format. Defines database connection, API endpoints, feature flags, and logging settings. 47 configuration keys across 5 top-level sections.

## What Was Found

- `database.host`, `database.port`, `database.name` (lines 2-4): PostgreSQL connection settings
- `api.base_url`, `api.timeout` (lines 7-8): External API configuration
- `features.experimental_mode: false` (line 12): Feature flag for beta features
- `logging.level: INFO`, `logging.format` (lines 15-16): Logging configuration

## What Was NOT Found

- No schema validation rules present
- No environment-specific overrides documented
- API authentication credentials not in this file
```

### Data Files

File extensions: `.csv`, `.tsv`, `.parquet`, `.json` (when data-structured), `.jsonl`, `.ndjson`

The model MUST extract:

1. **Row count** - Exact number of records
2. **Column names** - All column headers
3. **Data types** - Inferred from first N rows
4. **Sample values** - Representative examples from each column
5. **Missing data** - Columns with null/empty values
6. **Unique identifiers** - Primary key columns if evident

**Extraction method**: For CSV/TSV, read header row and first 10 data rows. For Parquet, note that binary inspection is limited. For JSON, inspect array structure.

**Example summary structure**:

```text
## Summary

CSV file containing user activity logs. 1,247 rows with 8 columns. Timestamps range from 2025-01-01 to 2026-02-06. No missing values detected.

## What Was Found

- Column `user_id` (integer): User identifiers, range 1001-5432
- Column `timestamp` (ISO 8601): Activity timestamps
- Column `action` (string): Values include "login", "logout", "view_page", "click_button"
- Column `duration_ms` (integer): Range 0-45000
- 1,247 total records (line count: 1,248 including header)

## What Was NOT Found

- No schema documentation in file
- Column `referrer` is present but not documented
- No indication of data collection methodology
```

### Documentation Files

File extensions: `.md`, `.rst`, `.txt`, `.adoc`, `.org`

The model MUST extract:

1. **Topic hierarchy** - Top-level headings and structure
2. **Key sections** - Main topics covered
3. **Commands/examples** - Code blocks, shell commands, API calls
4. **Links** - External references and internal cross-references
5. **Definitions** - Technical terms defined in the text

**Extraction method**: Read sequentially. Extract headings to build table of contents. Quote key passages that define core concepts. Note code examples.

**Example summary structure**:

```text
## Summary

User guide for deploying containerized applications. Covers Docker setup, image building, registry configuration, and troubleshooting. 5 main sections with 23 subsections. Includes 12 shell command examples.

## What Was Found

- Section "Getting Started" (lines 10-45): Docker installation on Linux and macOS
- Section "Building Images" (lines 47-89): Dockerfile syntax and multi-stage builds
- Section "Troubleshooting" (lines 200-245): Common errors with solutions
- 12 shell command examples throughout document

## What Was NOT Found

- No Windows deployment instructions
- Security best practices not covered
- Performance tuning section mentioned but not written (line 15: "TODO")
```

### Binary and Unknown Files

File extensions: `.pdf`, `.zip`, `.tar`, `.gz`, `.bin`, `.exe`, `.so`, `.dylib`, `.dll`, or unrecognized extensions

The model MUST:

1. **Attempt to read** - Use the Read tool. If the tool returns binary content or an error, note this.

2. **State limitation** - Do NOT guess contents. State: "Binary file, cannot extract text content."

3. **Provide file metadata** - File size, extension, location.

4. **For PDFs**: Use the Read tool with `pages` parameter to extract text from specific page ranges. Summarize text content if extraction succeeds.

**Example for unreadable binary**:

```text
## Summary

Binary file, cannot extract text content.

## What Was Found

- File path: ./build/output.bin
- File size: 2.3 MB
- Extension: .bin

## What Was NOT Found

Unable to determine contents without binary inspection tools.

## Uncertain

File may be compiled binary, compressed archive, or proprietary format.
```

## Quote-Grounding Technique

For all text-based files, the model MUST apply the quote-grounding technique:

1. **First pass** - Read file, identify key passages
2. **Extract** - Copy exact quotes with line numbers
3. **Organize extracts** - Group by theme or importance
4. **Summarize from extracts** - Write summary grounded in the extracted quotes
5. **Verify** - Ensure every claim in summary traces to an extract

SOURCE: Technique adapted from [Fidelity Rules](../summarizer/references/fidelity-rules.md) Rule 2 (lines 27-41).

## Output Format

All file summaries MUST use the structured output format defined in [Structured Summary](../summarizer/templates/structured.md).

Required sections:

1. **YAML frontmatter** - Include `source_type: file`, `source_path`, `method`, `confidence`, word counts
2. **Summary** - Condensed content (BLUF style)
3. **What Was Found** - Items discovered with line number references
4. **What Was NOT Found** - Expected items that were absent
5. **Uncertain** - Ambiguous items requiring interpretation
6. **Sources** - Full file path, access date

## Fidelity Rules

The model MUST follow all fidelity rules defined in [Fidelity Rules](../summarizer/references/fidelity-rules.md).

Critical rules for file summarization:

- **Rule 1**: Read the file before summarizing. Never guess from filename.
- **Rule 2**: Extract before abstracting. Identify key passages first.
- **Rule 3**: Preserve counts and specifics. "7 functions" not "several functions."
- **Rule 4**: Distinguish absence from nonexistence. "Not in file" not "doesn't exist."
- **Rule 6**: State confidence explicitly. Full read of small file = high confidence. Truncated large file = medium/low confidence.

## Multi-File Summarization

When the user requests summarization of multiple files:

1. Summarize each file individually using this methodology
2. Write each summary to a separate output file or section
3. Do NOT merge file summaries into a single combined summary without explicit user request
4. If synthesis across files is requested, load the multi-source-synthesis skill after completing individual summaries

SOURCE: Multi-source synthesis approach from [Summarizer](../summarizer/SKILL.md) lines 33-37.

## Error Handling

If a file cannot be read:

1. Attempt to read with the Read tool
2. If read fails, report the error: "Unable to read [file path]: [error message]"
3. Do NOT speculate about file contents
4. Do NOT proceed with summarization
5. Ask user if they want to try alternative access methods

## Output Rendering

1. **Read template** - Load the template file at `../summarizer/templates/{format_id}.md` (default: `structured`). The template defines the schema, required sections, and fidelity constraints for the selected format.
2. **Render** - Produce output following the template's Schema section. Use the template's Example as a reference for structure and style.
3. **Verify fidelity** - Confirm the output satisfies the template's Fidelity Constraints and all applicable [Fidelity Rules](../summarizer/references/fidelity-rules.md).

## Anti-Patterns

The model MUST NOT:

- Summarize a file based on its name without reading it
- Guess file contents from directory structure or naming conventions
- Assume file type from extension without verifying contents
- Summarize from partial reads (head/tail/grep) without disclosing the limitation
- Upgrade "not found in file" to "file doesn't contain" in a way that implies certainty about what the file should contain
- Present interpretation as observation
- Skip the "What Was NOT Found" section
- Omit line number references for key findings
