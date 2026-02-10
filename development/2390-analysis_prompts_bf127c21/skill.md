# AI Analysis Prompts for Git Changes

This file contains the prompts used to analyze git changes and categorize them into different domains. These prompts are adapted from the JIRA release notes analysis system but work purely with git data (commits and diffs).

## Primary Analysis Prompt

Use this prompt to analyze git changes and understand intent, context, and implementation:

```
You are analyzing git changes to generate a comprehensive merge request description.

<input_data>
<commits>{commit_details}</commits>
<diff>{changes_diff}</diff>
<file_list>{changed_files}</file_list>
<stats>
Commits: {commit_count}
Files Changed: {files_changed}
Lines Added: {lines_added}
Lines Deleted: {lines_deleted}
</stats>
</input_data>

For each logical change unit (commit or related file group):

1. **Understand the Intent**
   - Extract purpose from commit messages
   - Identify change type from diff patterns
   - Determine affected components/modules
   - Look for conventional commit prefixes (fix:, feat:, refactor:, etc.)

2. **Categorize Changes**
   Primary categories:

   - **Bug Fix**: Corrects defective behavior
     * Indicators: "fix", "bug", "issue", "correct", defensive code additions, error handling improvements
     * Diff patterns: Conditional checks added, error handling blocks, validation logic

   - **Enhancement**: Adds/improves functionality
     * Indicators: "feat", "add", "improve", "enhance", new functions/classes, expanded APIs
     * Diff patterns: New public methods, additional parameters, new modules/files

   - **Technical Debt**: Internal improvements
     * Indicators: "refactor", "cleanup", "simplify", renamed files, reduced duplication
     * Diff patterns: Code reorganization, extraction of functions, dependency updates

   - **Documentation**: Updates to docs/comments
     * Indicators: .md files, .rst files, docstrings, README, inline comments
     * Diff patterns: Only documentation files changed, or only comment lines in code

   - **Testing**: Test additions/modifications
     * Indicators: test files, test utilities, spec files, coverage improvements
     * Diff patterns: Files in test/ or spec/ directories, new test cases

   - **Non-Functional**: Formatting, config, minor cleanup
     * Indicators: whitespace changes, import reordering, formatting, linting fixes
     * Diff patterns: Only whitespace/formatting diffs, config file updates (.yml/.toml/.json)

   - **Build/CI**: Build system and CI/CD changes
     * Indicators: .github/workflows, .gitlab-ci.yml, Dockerfile, package.json scripts
     * Diff patterns: CI configuration files, build scripts, deployment configs

3. **Analyze Implementation**
   For each change:
   - **Logic changes**: Algorithms, control flow, business logic modifications
   - **API changes**: New methods, modified signatures, deprecated functions
   - **Configuration**: Environment variables, settings files, feature flags
   - **Database/Schema**: Migrations, model changes, query updates
   - **Dependencies**: Packages added, removed, upgraded (package.json, requirements.txt, etc.)
   - **Breaking changes**: Incompatibilities, removed features, changed interfaces

4. **Identify Impact**
   - **User-facing**: Changes visible to end users
   - **Developer-facing**: API changes, developer experience improvements
   - **Performance**: Speed improvements, optimization, resource usage
   - **Security**: Security fixes, vulnerability patches, auth improvements
   - **Migration required**: Steps needed to adopt changes

<output_format>
Return structured JSON:
```json
{
  "summary": "High-level overview of all changes in 2-3 sentences",
  "change_categories": {
    "bug_fixes": [
      {
        "title": "Clear, specific description of bug fixed",
        "description": "What was broken and how it's now fixed",
        "files_affected": ["path/to/file.py", "another/file.js"],
        "commits": ["abc123f", "def456g"],
        "impact": "Who/what this affects and why it matters",
        "technical_details": ["Implementation specifics if relevant"]
      }
    ],
    "enhancements": [
      {
        "title": "Feature added or capability improved",
        "description": "What new functionality was added",
        "benefits": "Value provided to users/developers",
        "usage": "How to use the new feature (if applicable)",
        "files_affected": ["path/to/file.py"],
        "commits": ["abc123f"],
        "technical_details": ["API additions", "Configuration options"]
      }
    ],
    "tech_debt": [
      {
        "title": "Refactoring or cleanup performed",
        "purpose": "Why this internal update was needed",
        "changes": "What was modified or reorganized",
        "impact": "How this improves code quality/maintainability",
        "files_affected": ["path/to/file.py"],
        "commits": ["abc123f"]
      }
    ],
    "documentation": [
      {
        "title": "Documentation added or updated",
        "updates": "What content was added/changed",
        "location": "Where to find the documentation",
        "importance": "Why this helps developers/users",
        "files_affected": ["README.md", "docs/api.md"],
        "commits": ["abc123f"]
      }
    ],
    "testing": [
      {
        "title": "Tests added or improved",
        "coverage": "What's now tested that wasn't before",
        "test_type": "unit|integration|e2e",
        "files_affected": ["tests/test_foo.py"],
        "commits": ["abc123f"]
      }
    ],
    "build_ci": [
      {
        "title": "Build or CI/CD change",
        "changes": "What was modified in build/deployment",
        "impact": "How this improves the development workflow",
        "files_affected": [".github/workflows/test.yml"],
        "commits": ["abc123f"]
      }
    ],
    "non_functional": [
      {
        "type": "formatting|linting|config|deps",
        "description": "Brief description of minor changes",
        "files_affected": ["path/to/file.py"],
        "commits": ["abc123f"]
      }
    ]
  },
  "components_affected": ["module1", "module2", "api-layer"],
  "breaking_changes": [
    {
      "change": "Description of breaking change",
      "migration": "Steps to migrate or adapt to the change",
      "files_affected": ["path/to/file.py"],
      "commits": ["abc123f"]
    }
  ],
  "statistics": {
    "commits": 15,
    "files_changed": 42,
    "lines_added": 1523,
    "lines_deleted": 847
  }
}
```
</output_format>

<guidelines>
- Focus on "why" behind changes, not just "what"
- Use functional/business terms where applicable
- Group related changes logically
- Highlight breaking changes prominently
- Provide migration guidance for API changes
- Quantify improvements when visible in diffs (performance, coverage, etc.)
- Connect technical changes to user/developer value
</guidelines>
```

## Domain-Specific Formatting Prompts

After categorization, use these prompts to format each category into polished MR description text:

### Bug Fixes Formatting

```
Format the following bug fix data into clear, professional MR description entries.

<analyzed_bugs>{bug_fixes_json}</analyzed_bugs>

For each bug fix, create an entry following this structure:

- **{title}**
  - **Description:** {what_was_broken_and_how_fixed}
  - **Impact:** {who_affected_and_why_it_matters}
  - **Files:** {affected_files}
  - **Technical Notes:** {optional_implementation_details}

Requirements:
1. Title: Use cause→effect pattern, avoid generic "Fixed X"
2. Description: Explain the bug behavior and resolution clearly
3. Impact: Focus on user/system consequences
4. Files: List key files changed (not exhaustive if many)
5. Technical Notes: Include only if adds value

Return formatted markdown text ready for MR description.
```

### Enhancements Formatting

```
Format the following enhancement data into compelling MR description entries.

<analyzed_enhancements>{enhancements_json}</analyzed_enhancements>

For each enhancement, create an entry following this structure:

- **{clear_feature_description}**
  - **Feature:** {what_was_added_or_changed}
  - **Benefits:** {how_this_improves_product_or_performance}
  - **Usage:** {how_to_use_if_applicable}
  - **Files:** {affected_files}
  - **Technical Details:** {optional_configuration_or_API_info}

Requirements:
1. Title: Clear, action-oriented description of capability
2. Feature: Specific functionality added or improved
3. Benefits: User/developer value proposition
4. Usage: Practical guidance if user-facing, omit if internal
5. Technical Details: API additions, configuration options if relevant

Return formatted markdown text ready for MR description.
```

### Technical Debt Formatting

```
Format the following technical debt items into MR description entries that highlight improvement value.

<analyzed_tech_debt>{tech_debt_json}</analyzed_tech_debt>

For each item, create an entry following this structure:

- **{summary_of_refactor}**
  - **Purpose:** {why_this_internal_update_was_needed}
  - **Changes:** {what_was_modified_or_removed}
  - **Impact:** {how_this_improves_system}
  - **Files:** {affected_files}

Requirements:
1. Title: Concise summary of refactoring (improvement-focused)
2. Purpose: Business or technical driver
3. Changes: High-level overview of modifications
4. Impact: Tangible improvements (maintainability, performance, etc.)

Return formatted markdown text ready for MR description.
```

### Documentation Formatting

```
Format the following documentation updates into MR description entries.

<analyzed_documentation>{documentation_json}</analyzed_documentation>

For each documentation update, create an entry:

- **{summary_of_doc_update}**
  - **Updates:** {what_content_was_updated}
  - **Location:** {where_to_find_it}
  - **Importance:** {why_it_matters}

Return formatted markdown text ready for MR description.
```

## Title Generation Prompt

Use this to generate a clear, descriptive MR title:

```
Generate a clear, concise merge request title based on the analysis below.

<change_summary>{json_analysis_summary}</change_summary>

Requirements:
- Keep under 80 characters
- Use imperative mood ("Add feature" not "Added feature" or "Adds feature")
- Be specific about what changed
- If multiple categories, lead with the most significant
- Use conventional commit prefixes if appropriate (feat:, fix:, refactor:, etc.)

Examples:
- "feat: Add user authentication with JWT and refresh tokens"
- "fix: Resolve race condition in payment processing"
- "refactor: Simplify database query layer and improve performance"
- "docs: Add comprehensive API documentation with examples"

Return only the title text, no additional explanation.
```
