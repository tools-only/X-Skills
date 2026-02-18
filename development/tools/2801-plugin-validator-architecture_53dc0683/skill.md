# Plugin Validator Architecture

**Analysis Date:** 2026-02-13
**Target:** `plugins/plugin-creator/scripts/plugin_validator.py`
**Total Lines:** 3045
**Purpose:** Comprehensive validation framework for Claude Code plugins with token-based complexity measurement

---

## 1. Validator Protocol

**Location:** Lines 249-286

**Interface Contract:**

```python
class Validator(Protocol):
    def validate(self, path: Path) -> ValidationResult: ...
    def can_fix(self) -> bool: ...
    def fix(self, path: Path) -> list[str]: ...
```

**Protocol Requirements:**

- `validate()` — Runs validation check on file/directory, returns ValidationResult
- `can_fix()` — Returns True if validator supports auto-fixing, False otherwise
- `fix()` — Applies automatic fixes, returns list of human-readable fix descriptions

**Design Principle:** All validators implement this protocol for polymorphic usage in validation pipelines.

---

## 2. Data Models

### 2.1 Core Models

**ValidationIssue** (lines 168-196):

```python
@dataclass(frozen=True)
class ValidationIssue:
    field: str                      # Which field has the issue
    severity: Literal["error", "warning", "info"]
    message: str                    # Human-readable description
    code: str                       # Error code (e.g., "FM001", "SK006")
    line: int | None = None         # Line number if applicable
    suggestion: str | None = None   # How to fix the issue
    docs_url: str | None = None     # Link to error code documentation
```

**ValidationResult** (lines 199-206):

```python
@dataclass(frozen=True)
class ValidationResult:
    passed: bool                    # True if no errors, False otherwise
    errors: list[ValidationIssue]   # Blocking issues (fail validation)
    warnings: list[ValidationIssue] # Non-blocking issues (pass with warnings)
    info: list[ValidationIssue]     # Informational messages
```

**ComplexityMetrics** (lines 209-242):

```python
@dataclass(frozen=True)
class ComplexityMetrics:
    total_tokens: int               # Total tokens in file
    frontmatter_tokens: int         # Tokens in frontmatter section
    body_tokens: int                # Tokens in body content
    encoding: str = "cl100k_base"   # Tokenizer encoding used

    @property
    def status(self) -> Literal["ok", "warning", "error"]:
        # error: body_tokens > 10000 (TOKEN_ERROR_THRESHOLD)
        # warning: body_tokens > 4000 (TOKEN_WARNING_THRESHOLD)
        # ok: body_tokens <= 4000
```

**Why Token-Based Metrics:** Token count directly measures what Claude processes, unlike line counts which vary by verbosity.

**Thresholds** (lines 44-45):
- `TOKEN_WARNING_THRESHOLD = 4000` (~500 lines equivalent)
- `TOKEN_ERROR_THRESHOLD = 10000` (~1250 lines equivalent)

### 2.2 FileType Enum

**Location:** Lines 138-166

**Values:**
- `SKILL` — Detected by filename `SKILL.md`
- `AGENT` — Detected by `agents/` in path
- `COMMAND` — Detected by `commands/` in path
- `PLUGIN` — Detected by `plugin.json` or `.claude-plugin/plugin.json`
- `UNKNOWN` — Cannot determine type

**Detection Logic** (`detect_file_type`, lines 147-165):

```python
if path.name == "SKILL.md":
    return FileType.SKILL
if path.name == "plugin.json" or (path / ".claude-plugin/plugin.json").exists():
    return FileType.PLUGIN
if "agents" in path.parts:
    return FileType.AGENT
if "commands" in path.parts:
    return FileType.COMMAND
return FileType.UNKNOWN
```

---

## 3. Pydantic Schema Models

### 3.1 SkillFrontmatter (lines 1030-1101)

**Purpose:** Validates YAML frontmatter for SKILL.md files

**Fields:**

| Field | Type | Required | Validation |
|-------|------|----------|------------|
| `name` | `str \| None` | No | Max 64 chars, pattern `^[a-z][a-z0-9-]*$` |
| `description` | `str \| None` | No | No colons except in URLs |
| `argument-hint` | `str \| None` | No | None |
| `allowed-tools` | `str \| None` | No | Comma-separated string |
| `model` | `str \| None` | No | None |
| `skills` | `str \| None` | No | Comma-separated string |
| `context` | `Literal["fork"] \| None` | No | None |
| `agent` | `str \| None` | No | None |
| `user-invocable` | `bool \| None` | No | None |
| `disable-model-invocation` | `bool \| None` | No | None |
| `hooks` | `dict[str, Any] \| None` | No | None |

**Validators:**

1. **normalize_comma_separated** (lines 1052-1062): Converts YAML arrays to comma-separated strings
   - Applied to: `skills`, `allowed_tools`
   - Example: `["Read", "Write"]` → `"Read, Write"`

2. **normalize_single_line** (lines 1064-1074): Collapses multiline descriptions
   - Applied to: `description`
   - Example: `"Line 1\nLine 2"` → `"Line 1 Line 2"`

3. **validate_no_colons** (lines 1076-1100): Rejects colons outside URLs
   - Applied to: `description`
   - Rationale: Colons trigger YAML quoting issues

### 3.2 CommandFrontmatter (lines 1103-1164)

**Purpose:** Validates frontmatter for command .md files

**Fields:**

| Field | Type | Required | Validation |
|-------|------|----------|------------|
| `description` | `str` | **Yes** | No colons except in URLs |
| `argument-hint` | `str \| None` | No | None |
| `allowed-tools` | `str \| None` | No | Comma-separated string |
| `model` | `str \| None` | No | None |
| `context` | `Literal["fork"] \| None` | No | None |
| `agent` | `str \| None` | No | None |
| `hooks` | `dict[str, Any] \| None` | No | None |

**Validators:** Same as SkillFrontmatter (comma-separated, single-line, no-colons)

### 3.3 AgentFrontmatter (lines 1167-1235)

**Purpose:** Validates frontmatter for agent .md files

**Fields:**

| Field | Type | Required | Validation |
|-------|------|----------|------------|
| `name` | `str` | **Yes** | Max 64 chars, pattern `^[a-z][a-z0-9-]*$` |
| `description` | `str` | **Yes** | No colons except in URLs |
| `tools` | `str \| None` | No | Comma-separated string |
| `disallowedTools` | `str \| None` | No | Comma-separated string |
| `model` | `Literal["sonnet", "opus", "haiku", "inherit"] \| None` | No | Enum validation |
| `permissionMode` | `Literal[...] \| None` | No | Enum: default, acceptEdits, dontAsk, bypassPermissions, plan |
| `skills` | `str \| None` | No | Comma-separated string |
| `hooks` | `dict[str, Any] \| None` | No | None |
| `color` | `str \| None` | No | None |

**Note:** Field names use camelCase to match official agent schema (lines 1172-1173)

**Validators:** Same as SkillFrontmatter (comma-separated, single-line, no-colons)

---

## 4. Validator Implementations

### 4.1 FrontmatterValidator (lines 1243-1600)

**Purpose:** Validates and auto-fixes YAML frontmatter in all capability files

**Validation Steps:**

1. **Extract frontmatter** (lines 1286): Uses `extract_frontmatter()` utility
2. **Check delimiter presence** (lines 1287-1300): FM003 if missing `---` markers
3. **Detect forbidden multiline** (lines 1308-1318): FM004 if uses `|`, `>`, `|-`, `>-`
4. **Parse YAML** (lines 1321-1335): FM002 if YAML syntax invalid
5. **Select Pydantic model** (lines 1351-1356): SkillFrontmatter, AgentFrontmatter, or CommandFrontmatter
6. **Validate with Pydantic** (lines 1358-1377): Convert Pydantic errors to ValidationIssue

**Auto-Fix Capabilities** (`can_fix() = True`, lines 1430-1590):

- Converts YAML arrays to comma-separated strings
- Normalizes multiline descriptions to single-line
- Quotes descriptions containing colons
- **Bug workaround:** Removes `name` field from plugin skills (lines 1541-1560)

**Bug Workaround Detail:** Plugin skills with `name:` field don't appear as slash commands in Claude Code v2.1.23. Validator auto-removes this field with explanation (line 1558).

### 4.2 NameFormatValidator (lines 1603-1787)

**Purpose:** Validates name field format (lowercase, hyphens only)

**Validation Rules:**

- SK001: Uppercase characters not allowed
- SK002: Underscores not allowed (use hyphens)
- SK003: Leading/trailing/consecutive hyphens not allowed
- FM010: Must match pattern `^[a-z0-9][a-z0-9-]*[a-z0-9]$` or `^[a-z0-9]$`

**Auto-Fix:** No (`can_fix() = False`) — name changes require semantic decisions

**Location References:**
- Pattern constant: Line 52
- Error codes: Lines 83-86

### 4.3 DescriptionValidator (lines 1790-1936)

**Purpose:** Validates description field length and trigger phrases

**Validation Rules:**

- SK004: Description too short (< 20 characters minimum, line 48)
- SK005: Missing trigger phrases (required phrases, lines 56-64)
- Warning: Description exceeds 1024 characters (recommended max, line 49)

**Trigger Phrases Required** (lines 56-64):
```python
REQUIRED_TRIGGER_PHRASES = [
    "use when", "use this", "used when", "used by",
    "when ", "trigger", "activate"
]
```

**Rationale:** Descriptions must include activation context for effective skill selection.

**Auto-Fix:** No (`can_fix() = False`) — trigger phrases require semantic context

### 4.4 ComplexityValidator (lines 1939-2122)

**Purpose:** Token-based complexity measurement for skills

**Validation Flow:**

1. **Count tokens** using `tiktoken` library (cl100k_base encoding)
2. **Split frontmatter vs body** tokens for separate measurement
3. **Compare to thresholds:**
   - SK006 (warning): body_tokens > 4000
   - SK007 (error): body_tokens > 10000

**Why Token-Based:** Line count varies by verbosity; token count measures actual Claude processing cost.

**Auto-Fix:** No (`can_fix() = False`) — splitting skills requires semantic analysis

**Location References:**
- Thresholds: Lines 44-45
- Error codes: Lines 89-90

### 4.5 InternalLinkValidator (lines 448-631)

**Purpose:** Validates markdown links in SKILL.md files

**Validation Rules:**

- LK001 (error): Broken link (target file doesn't exist)
- LK002 (warning): Link missing `./` prefix

**Link Extraction** (lines 457-490):
- Regex pattern: `\[([^\]]+)\]\(([^)]+)\)` (line 458)
- **Critical:** Strips code blocks and inline code before scanning (lines 468-490)
- **Rationale:** Prevents false positives from code examples

**Ignored Links:**
- External URLs (`http://`, `https://`, `ftp://`)
- Anchor links (`#heading`)
- Absolute paths (`/absolute/path`)

**Auto-Fix:** No (`can_fix() = False`) — broken links require file creation or manual path correction

### 4.6 ProgressiveDisclosureValidator (lines 334-441)

**Purpose:** Checks for progressive disclosure directory structure

**Checked Directories** (line 345):
- `references/`
- `examples/`
- `scripts/`

**Issue Severity:** INFO only (lines 376-386) — directories are optional organizational aids

**Error Codes:**
- PD001: No `references/` directory
- PD002: No `examples/` directory
- PD003: No `scripts/` directory

**Auto-Fix:** No (`can_fix() = False`) — creating directories requires content creation decisions

### 4.7 NamespaceReferenceValidator (lines 638-1023)

**Purpose:** Validates namespace-qualified references in plugin files

**Detected Patterns:**

| Pattern | Regex | Example | Target Type |
|---------|-------|---------|-------------|
| Skill command | `Skill\(command:\s*"([^"]+):([^"]+)"` | `Skill(command: "plugin:skill")` | skill |
| Skill skill | `Skill\(skill="([^"]+):([^"]+)"` | `Skill(skill="plugin:skill")` | skill |
| Task agent | `Task\(agent[=:]\s*"([^"]+):([^"]+)"` | `Task(agent="plugin:agent")` | agent |
| At-agent | `@([a-z0-9-]+):([a-z0-9-]+)` | `@plugin:agent` | agent |
| Slash command | `(?<!\w)/([a-z0-9-]+):([a-z0-9-]+)` | `/plugin:skill` | command |

**Validation Flow:**

1. **Parse references** from body text (lines 993-1022)
2. **Resolve plugin directory** from `~/.claude/settings.json` or `.claude/settings.json`
3. **Check target file exists:**
   - Skills: `{plugin}/skills/{name}/SKILL.md`
   - Agents: `{plugin}/agents/{name}.md`
   - Commands: `{plugin}/commands/{name}.md`

**Builtin Agent Filtering** (lines 661-697): Ignores references to built-in agents like `Explore`, `general-purpose`, `context-gathering`

**Error Codes:**
- NR001: Namespace reference target does not exist
- NR002: Namespace reference points outside plugin directory

**Auto-Fix:** No (`can_fix() = False`) — requires creating missing files or fixing references

### 4.8 PluginStructureValidator (lines 2125-2414)

**Purpose:** Validates plugin.json schema and component references

**Validation Steps:**

1. **Check plugin.json exists** (PL001)
2. **Parse JSON syntax** (PL002)
3. **Validate required fields** (PL003 — `name` field)
4. **Validate component paths** (PL004 — must start with `./`)
5. **Check referenced files exist** (PL005)

**Integration:** Optionally runs `claude plugin validate` CLI if available (lines 2421-2500)

**Auto-Fix:** No (`can_fix() = False`) — structural issues require manual intervention

---

## 5. Validator Registration

**Location:** `validate_path()` function, lines 2865-2951

**Registration Logic:**

```python
def _validate_single_path(path: Path, *, check: bool, fix: bool, verbose: bool):
    file_type = FileType.detect_file_type(path)

    validators: list[Validator] = []

    if file_type in {FileType.SKILL, FileType.AGENT, FileType.COMMAND}:
        # All capability files
        validators.extend([
            FrontmatterValidator(),
            NameFormatValidator(),
            DescriptionValidator(),
            NamespaceReferenceValidator(),
        ])

        # Skill-specific validators
        if file_type == FileType.SKILL:
            validators.extend([
                ComplexityValidator(),
                InternalLinkValidator(),
                ProgressiveDisclosureValidator(),
            ])

    elif file_type == FileType.PLUGIN:
        validators.append(PluginStructureValidator())

    # Run validation
    for validator in validators:
        result = validator.validate(path)
        results.append((path, result))

    # Apply fixes if requested
    if fix:
        for validator in validators:
            if validator.can_fix():
                fixes = validator.fix(path)
                # Re-validate after fixes
```

**Selection Strategy:**

| File Type | Validators Applied |
|-----------|-------------------|
| SKILL.md | Frontmatter, NameFormat, Description, NamespaceReference, Complexity, InternalLink, ProgressiveDisclosure |
| Agent .md | Frontmatter, NameFormat, Description, NamespaceReference |
| Command .md | Frontmatter, NameFormat, Description, NamespaceReference |
| Plugin directory | PluginStructure |

---

## 6. CLI Interface

**Location:** Lines 2954-3045

**Typer App Structure:**

```python
app = typer.Typer(
    name="plugin-validator",
    help="Validate Claude Code plugins with comprehensive checks",
    no_args_is_help=True,
)

@app.command()
def main(
    paths: list[Path],              # Multiple paths supported
    check: bool = False,            # Validate only, no auto-fix
    fix: bool = False,              # Auto-fix issues
    verbose: bool = False,          # Show info messages
    no_color: bool = False,         # CI mode (no Rich formatting)
):
    """Validate Claude Code plugins, skills, agents, and commands."""
```

**Flags:**

| Flag | Short | Purpose | Default |
|------|-------|---------|---------|
| `--check` | - | Validate only, don't auto-fix | False |
| `--fix` | - | Auto-fix issues where possible | False |
| `--verbose` | `-v` | Show detailed output including info messages | False |
| `--no-color` | - | Disable color output for CI environments | False |

**Exit Codes:**

| Code | Meaning |
|------|---------|
| 0 | Success (all checks passed) |
| 1 | Validation errors found |
| 2 | Usage error (invalid arguments) |
| 130 | Interrupted by user (Ctrl+C) |

**Entry Point:** `if __name__ == "__main__": app()`

---

## 7. Report Architecture

### 7.1 Reporter Protocol (lines 2561-2592)

```python
class Reporter(Protocol):
    def report(self, results: list[tuple[Path, ValidationResult]], verbose: bool = False) -> None: ...
    def summarize(self, total_files: int, passed: int, failed: int, warnings: int) -> None: ...
```

### 7.2 ConsoleReporter (lines 2599-2721)

**Purpose:** Rich-formatted terminal output with colors and panels

**Features:**
- Uses `rich` library for formatting
- Displays issues with icons (✓, ✗, ⚠, ℹ)
- Shows summary in colored panel
- Respects `NO_COLOR` environment variable

**Known Bug** (lines 2681-2709): Counting logic counts files, not issues. Summary shows "Passed: X files" but counts individual ValidationResult objects, which can be multiple per file.

**Impact:** When validating a single file with multiple validators, each validator produces a ValidationResult. The reporter counts these as separate "files" in the summary.

### 7.3 CIReporter (lines 2728-2807)

**Purpose:** Plain text output for CI environments

**Format:** `file:line [CODE] field: message`

**Example:**
```
✓ path/to/SKILL.md - PASSED

path/to/agent.md
  ✗ ERROR [FM001] name: Missing required field: name
    → Add 'name' field to frontmatter
```

### 7.4 SummaryReporter (lines 2814-2858)

**Purpose:** Single-line status for quick checks

**Output:** `✓ 5/5 files passed (2 with warnings)`

---

## 8. Test Patterns

**Test Directory:** `plugins/plugin-creator/tests/`

### 8.1 Fixtures (conftest.py)

**Location:** Lines 1-342

**Key Fixtures:**

| Fixture | Purpose | Returns |
|---------|---------|---------|
| `cli_runner` | CliRunner for testing Typer app | `CliRunner(mix_stderr=False)` |
| `sample_skill_dir` | Valid skill directory | `Path` to skill with SKILL.md |
| `sample_agent_dir` | Valid agent directory | `Path` to agents/ with agent.md |
| `sample_plugin_dir` | Valid plugin directory | `Path` to plugin with plugin.json |
| `invalid_frontmatter_samples` | Invalid frontmatter examples | `dict[str, str]` mapping scenario to content |
| `broken_link_samples` | Broken link test cases | `dict[str, tuple[Path, str]]` |
| `no_color_env` | Sets NO_COLOR=1 | Generator for env var context |

**Example Fixture Usage:**

```python
@pytest.fixture
def sample_skill_dir(tmp_path: Path) -> Path:
    skill_dir = tmp_path / "test-skill"
    skill_dir.mkdir()
    (skill_dir / "SKILL.md").write_text("""---
description: Use this skill when testing validation
---
# Test Skill
""")
    return skill_dir
```

### 8.2 Test Structure (test_frontmatter_validator.py)

**Location:** Lines 1-150

**Test Class Pattern:**

```python
class TestFrontmatterValidatorBasic:
    """Test basic FrontmatterValidator functionality."""

    def test_validator_instantiation(self) -> None:
        """Test FrontmatterValidator can be instantiated."""
        validator = FrontmatterValidator()
        assert validator is not None
        assert validator.can_fix() is True

    def test_valid_skill_frontmatter(self, tmp_path: Path) -> None:
        """Test validation passes for valid skill frontmatter.

        Tests: Valid skill SKILL.md with minimal frontmatter
        How: Create file with description only, validate
        Why: Ensure validator accepts valid minimal skill frontmatter
        """
        skill_md = tmp_path / "SKILL.md"
        skill_md.write_text("""---
description: Use this skill when testing validation
---
# Test Skill Content
""")

        validator = FrontmatterValidator()
        result = validator.validate(skill_md)

        assert result.passed is True
        assert len(result.errors) == 0
```

**Test Organization:**

| Test File | Coverage |
|-----------|----------|
| `test_frontmatter_validator.py` | FrontmatterValidator validation and auto-fix |
| `test_plugin_structure_validator.py` | PluginStructureValidator |
| `test_progressive_disclosure_validator.py` | ProgressiveDisclosureValidator |
| `test_internal_link_validator.py` | InternalLinkValidator |
| `test_complexity_validator.py` | ComplexityValidator token-based metrics |
| `test_description_validator.py` | DescriptionValidator trigger phrases |
| `test_name_format_validator.py` | NameFormatValidator naming rules |
| `test_external_tools.py` | Claude CLI integration tests |
| `test_cli.py` | CLI argument parsing and exit codes |
| `test_token_counting.py` | Token counting accuracy |

**What's Covered:**
- Valid input acceptance
- Invalid input rejection with correct error codes
- Auto-fix correctness
- Boundary conditions
- Integration with external tools (claude CLI)

**What's Not Covered:**
- NamespaceReferenceValidator (no dedicated test file)
- Multi-file validation workflows
- Pre-commit hook integration

---

## 9. Extension Points

### 9.1 Adding a New FileType

**Steps:**

1. **Add enum value** to `FileType` (line 138):
   ```python
   class FileType(StrEnum):
       HOOK = "hook"
   ```

2. **Update detection logic** in `FileType.detect_file_type()` (lines 147-165):
   ```python
   if "hooks" in path.parts and path.suffix == ".js":
       return FileType.HOOK
   ```

3. **Create Pydantic model** (if frontmatter validation needed):
   ```python
   class HookFrontmatter(BaseModel):
       name: str
       type: Literal["SubagentStart", "SubagentStop", ...]
   ```

4. **Update FrontmatterValidator._get_model_class()** (line 1352):
   ```python
   if file_type == FileType.HOOK:
       return HookFrontmatter
   ```

5. **Add validator registration** in `_validate_single_path()` (lines 2896-2913):
   ```python
   if file_type == FileType.HOOK:
       validators.extend([
           FrontmatterValidator(),
           # Hook-specific validators
       ])
   ```

### 9.2 Adding a New Validator

**Steps:**

1. **Implement Validator protocol:**
   ```python
   class MyValidator:
       def validate(self, path: Path) -> ValidationResult:
           errors: list[ValidationIssue] = []
           # Validation logic
           return ValidationResult(passed=len(errors) == 0, errors=errors, warnings=[], info=[])

       def can_fix(self) -> bool:
           return False

       def fix(self, path: Path) -> list[str]:
           raise NotImplementedError("...")
   ```

2. **Define error codes** (lines 67-110):
   ```python
   MV001 = "MV001"  # My validator error 001
   ```

3. **Register in `_validate_single_path()`** (lines 2896-2913):
   ```python
   if file_type == FileType.SKILL:
       validators.append(MyValidator())
   ```

4. **Add tests:**
   ```python
   # tests/test_my_validator.py
   class TestMyValidator:
       def test_valid_input(self, tmp_path: Path) -> None: ...
       def test_invalid_input(self, tmp_path: Path) -> None: ...
   ```

### 9.3 Adding a New Frontmatter Field

**Steps:**

1. **Update Pydantic model** (e.g., SkillFrontmatter, lines 1030-1101):
   ```python
   class SkillFrontmatter(BaseModel):
       new_field: str | None = None
   ```

2. **Add field validator** if needed:
   ```python
   @field_validator("new_field", mode="before")
   @classmethod
   def validate_new_field(cls, v: object) -> str | None:
       # Custom validation logic
       return v
   ```

3. **Update documentation** (ERROR_CODES.md)

4. **Add tests** for new field validation

---

## 10. Known Issues

### 10.1 Reporter Counting Bug

**Location:** ConsoleReporter.summarize(), lines 2681-2709

**Issue:** Counts ValidationResult objects instead of unique files

**Impact:** When a single file has multiple validators, summary shows inflated "total files" count

**Example:**
```
Input: 1 file validated by 7 validators
Output: "Total files: 7" (should be "Total files: 1")
```

**Root Cause:** Each validator produces a ValidationResult, and the reporter counts these as separate files.

**Fix Approach:**
1. Change `results` parameter to include unique file paths
2. Count `set(path for path, _ in results)` instead of `len(results)`
3. Update test expectations

### 10.2 Name Field Bug Workaround

**Location:** FrontmatterValidator.fix(), lines 1541-1560

**Issue:** Claude Code v2.1.23 bug — plugin skills with `name:` field don't appear as slash commands

**Workaround:** Validator auto-removes `name:` field from plugin skills

**Detection:** Checks if skill is in a plugin directory structure

**Status:** Temporary workaround until Claude Code bug fixed

---

## Summary

**Architecture Pattern:** Protocol-based validation pipeline with Pydantic schema models and token-based complexity measurement.

**Key Design Decisions:**

1. **Protocol over inheritance** — Validator protocol enables polymorphic composition without tight coupling
2. **Token-based metrics** — More accurate than line counts for measuring Claude processing cost
3. **Pydantic validation** — Declarative schema validation with automatic error conversion
4. **Auto-fix separation** — `can_fix()` and `fix()` clearly separate read-only validators from mutating ones
5. **Reporter protocol** — Pluggable output formatting (Rich terminal, CI, summary)
6. **FileType-based dispatch** — Validates different component types with appropriate validators

**Extension Pattern:** Add new validators by implementing Validator protocol and registering in `_validate_single_path()`. Add new file types by extending FileType enum and detection logic.

**Test Coverage:** Comprehensive unit tests for individual validators, CLI integration tests, fixture-based test data generation.

**Dependencies:**
- `typer` — CLI framework
- `rich` — Terminal formatting
- `tiktoken` — Token counting (cl100k_base encoding)
- `pyyaml` — YAML parsing
- `pydantic` — Schema validation

---

**Analysis Complete:** All 8 requested architecture components documented with line number citations.
