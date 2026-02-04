---
name: am-i-complete
description: Self-verification protocol that drives completion of testing steps to properly finish the task
---

# IS IT WORKING? - Task Verification Protocol

**Workflow Reference**: See [Master Workflow](./../knowledge/workflow-diagrams/master-workflow.md) for how this command fits into the verification stage of the agentic workflow.

You're in the middle of a task. Before claiming completion, execute this verification protocol.

## IMMEDIATE ACTION REQUIRED

<task_identification> First, understand your project's test structure:

- Check for existing test directories: tests/, test/, spec/, **tests**/
- Identify test categories: unit/, integration/, e2e/, functional/, regression/
- Note the testing framework: pytest, jest, mocha, go test, etc.
- Find naming conventions: test\__.py,_.test.js, _\_test.go,_.spec.ts

Your task type:

- FIX: Resolving a bug → Execute FIX_VERIFICATION
- FEATURE: Adding functionality → Execute FEATURE_VERIFICATION
- REFACTOR: Restructuring code → Execute REFACTOR_VERIFICATION
- OTHER: Any other task → Execute GENERAL_VERIFICATION </task_identification>

---

## FIX_VERIFICATION

<verification_protocol type="fix">

<step_1_reproduce> **REPRODUCE THE BUG NOW**

Find or create the appropriate test file:

- If fixing src/module/component.py → Create/update tests/unit/test_component.py
- If fixing a integration issue → Create/update tests/integration/test\_[issue_description].py
- Add a test method that REPRODUCES the bug:

<example>
def test_task_completed_should_print_progress(capsys: pytest.CaptureFixture[str]) -> None:
    """Regression test for issue where task_completed() prints nothing"""
    # This test should FAIL before the fix
    task_completed("Building")
    captured = capsys.readouterr()
    assert "Building" in captured.out  # This will fail, proving the bug exists
</example>

Run this test NOW and confirm it fails:

<commands>
# -x: stop at first failure (fast feedback for TDD)
# -v: verbose output with test names
# -s: show print statements and stdout (disable output capture)
pytest tests/unit/test_component.py::test_task_completed_should_print_progress -xvs
</commands>

Save the failure output for reference. </step_1_reproduce>

<step_2_fix> **APPLY YOUR FIX**

Make your code changes to fix the bug. The test you just wrote should guide your fix. </step_2_fix>

<step_3_verify> **VERIFY THE FIX WORKS**

Run the SAME test - it must pass now: <commands>

# The regression test must now pass

pytest tests/unit/test_component.py::test_task_completed_should_print_progress -xvs

# Run all related tests to ensure you didn't break anything

pytest tests/unit/test_component.py -xvs

# If this was an integration issue, run integration tests

pytest tests/integration/ -xvs -k "relevant_keyword" </commands>

Show actual output demonstrating the test passes. </step_3_verify>

<step_4_regression> **CHECK FOR REGRESSIONS**

Run the full test suite relevant to your changes: <commands>

# Run unit tests for the module you changed

pytest tests/unit/test\_[affected_module].py -xvs

# Run integration tests if you changed interfaces

pytest tests/integration/ -xvs

# Run the full suite if changes are significant

pytest tests/ -xvs

# Or use the project's test command

make test npm test go test ./... </commands>

Document any failures and fix them. </step_4_regression>

<step_5_commit> **COMMIT YOUR FIX WITH THE TEST**

Your commit should include:

- The fix itself
- The regression test that prevents this bug from recurring

<commands>
git add src/module/component.py tests/unit/test_component.py
git commit -m "fix: task_completed() now prints progress messages

- Added print statement to task_completed()
- Added regression test to prevent silent execution
- Fixes #123" </commands> </step_5_commit>

</verification_protocol>

**STATUS CHECK**: Is your regression test committed and passing? If not, continue working.

---

## FEATURE_VERIFICATION

<verification_protocol type="feature">

<step_1_criteria> **DEFINE TEST REQUIREMENTS**

Based on your feature, determine where tests belong:

- Unit tests: tests/unit/test\_[new_module].py for new functions/classes
- Integration tests: tests/integration/test\_[feature_name].py for component interactions
- E2E tests: tests/e2e/test\_[user_workflow].py for user-facing features

Write the test signatures FIRST (they will fail):

<example>
# tests/unit/test_json_exporter.py

def test_export_to_json_creates_valid_file() -> None: """Test that export creates a valid JSON file""" pass # TODO: implement

def test_export_handles_special_characters() -> None: """Test that special characters are properly escaped""" pass # TODO: implement

def test_export_raises_on_invalid_input() -> None: """Test that appropriate errors are raised""" pass # TODO: implement </example> </step_1_criteria>

<step_2_implement> **IMPLEMENT FEATURE WITH TESTS**

For each test you defined:

1. Implement the test fully
2. Run it - confirm it fails
3. Implement the feature code to make it pass
4. Confirm the test passes

<commands>
# Run each test as you implement
pytest tests/unit/test_json_exporter.py::test_export_to_json_creates_valid_file -xvs

# Once all unit tests pass, run them together

pytest tests/unit/test_json_exporter.py -xvs </commands>

Don't move forward until all tests pass. </step_2_implement>

<step_3_integration> **ADD INTEGRATION TESTS**

Create tests/integration/test\_[feature]\_integration.py:

<example>
def test_json_export_full_workflow() -> None:
    """Test the complete export workflow as a user would use it"""
    # Setup
    data = create_realistic_test_data()

    # Execute
    result = export_workflow(data, format='json')

    # Verify
    assert result.success
    assert valid_json(result.output)
    assert all_data_preserved(data, result.output)

</example>

Run integration tests: <commands> pytest tests/integration/test\_\*\_integration.py -xvs </commands> </step_3_integration>

<step_4_e2e> **ADD E2E TEST IF APPLICABLE**

For user-facing features, add tests/e2e/test\_[feature]\_workflow.py:

- Test the complete user journey
- Include error scenarios
- Test edge cases

<commands>
pytest tests/e2e/ -xvs
</commands>
</step_4_e2e>

<step_5_commit> **COMMIT FEATURE WITH ALL TESTS**

<commands>
git add src/[feature_files] tests/
git commit -m "feat: add JSON export functionality

- Implemented JsonExporter class with validation
- Added comprehensive unit tests
- Added integration tests for full workflow
- Added E2E tests for user scenarios
- Closes #456" </commands> </step_5_commit>

</verification_protocol>

**STATUS CHECK**: Are all your tests committed and passing in the proper test directories? If not, continue working.

---

## REFACTOR_VERIFICATION

<verification_protocol type="refactor">

<step_1_baseline> **ENSURE EXISTING TESTS PASS**

Before refactoring, verify the current state: <commands>

# Run tests for the code you're about to refactor

pytest tests/unit/test\_[module_to_refactor].py -xvs > /tmp/baseline_tests.txt pytest tests/integration/ -xvs -k [module] >> /tmp/baseline_tests.txt

# Ensure all pass before starting

</commands>

If tests are missing for the code you're refactoring, ADD THEM FIRST. </step_1_baseline>

<step_2_refactor> **PERFORM REFACTOR**

Make your refactoring changes. After each significant change, run the tests: <commands> pytest tests/unit/test\_[refactored_module].py -xvs </commands>

Tests should continue passing throughout the refactor. </step_2_refactor>

<step_3_verify> **VERIFY BEHAVIOR UNCHANGED**

<commands>
# Run the same tests as baseline
pytest tests/unit/test_[module_to_refactor].py -xvs > /tmp/after_tests.txt
pytest tests/integration/ -xvs -k [module] >> /tmp/after_tests.txt

# Compare test results

diff /tmp/baseline_tests.txt /tmp/after_tests.txt

# Run performance tests if applicable

pytest tests/performance/ -xvs </commands>

All tests should pass with similar or better performance. </step_3_verify>

<step_4_commit> **COMMIT REFACTOR**

<commands>
git add -A
git commit -m "refactor: [description of refactoring]

- [What was changed]
- [Why it's better]
- All existing tests still pass
- No functional changes" </commands> </step_4_commit>

</verification_protocol>

**STATUS CHECK**: Do all existing tests still pass after refactoring? If not, continue working.

---

## GENERAL_VERIFICATION

<verification_protocol type="general">

<thinking>
What type of task am I doing?
Does it need tests?
Where do those tests belong in the existing structure?
</thinking>

<checklist>
Determine what validation is needed:
- [ ] If code changes: are there tests?
- [ ] If bug fix: is there a regression test?
- [ ] If feature: are there unit + integration tests?
- [ ] If refactor: do existing tests still pass?
- [ ] If config/docs: have you validated the changes work?
</checklist>

<execute>
Based on your task:
1. Find or create appropriate test files in the existing test structure
2. Write tests that will live in the codebase (not throwaway tests)
3. Ensure all tests pass
4. Commit tests with your changes
</execute>

</verification_protocol>

---

## CRITICAL REMINDERS

<test_location_rules> **Where Tests Belong**

- Unit tests: tests/unit/test\_[module_name].py - test individual functions/classes
- Integration tests: tests/integration/test\_[feature]\_integration.py - test component interactions
- E2E tests: tests/e2e/test\_[workflow].py - test complete user workflows
- Regression tests: Add to existing test files as new test methods
- Performance tests: tests/performance/test\_[module]\_perf.py if applicable

NEVER create random test files in the project root. ALWAYS follow the project's existing test structure. </test_location_rules>

<permanent_tests> **Tests Are Permanent**

- Every bug fix MUST include a regression test
- Every feature MUST include unit and integration tests
- Tests are part of the deliverable, not verification artifacts
- Tests prevent future breaks and document behavior
- Tests must be committed with the code changes </permanent_tests>

<temp_artifacts> **Temporary Artifacts Only** Use /tmp/ or .verification/ (gitignored) ONLY for:

- Output captures for comparison
- Log files for debugging
- Temporary state files

Never commit temporary verification artifacts. </temp_artifacts>

---

## YOUR IMMEDIATE NEXT STEP

<action>
Answer these questions:
1. What is your specific task?
2. What test files already exist for the code you're modifying?
3. What new test files need to be created in the test structure?
4. Have you written the failing test FIRST?

Then execute: <commands>

# Find relevant existing tests

find tests -name "_[relevant_module]_" -type f

# Run existing tests to establish baseline

pytest tests/[relevant_path] -xvs

# Create or update test file in proper location

# Write the test that validates your change

# Run it and show it fails (for fixes) or passes (for features)

</commands>
</action>

**This command ensures you CREATE PERMANENT TESTS, not temporary verification files. Keep working until tests are committed.**
