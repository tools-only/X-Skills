# MicroSim Utils - TODO

## Pending Enhancements

### Add iframe Example Path Validation to Standardization Check

**Priority:** Medium
**Added:** 2026-01-01

**Description:**
Add a new standardization check that validates sample iframe example URLs in `index.md` files to ensure they reference the correct GitHub Pages deployment path.

**Problem Discovered:**
During a routine check of `docs/sims/*/index.md` files, several MicroSims had incorrect iframe example paths:
- Wrong repository name (e.g., `microsims` instead of `intro-to-physics-course`)
- Typos in repository name (e.g., `intro-to-physics` missing `-course`)
- Wrong directory name (e.g., `basic-fft` instead of `fft-basic`)

**Proposed Implementation:**

1. Add to `references/standardization.md` checklist:
   ```markdown
   #### 13. Sample iframe URL Validation
   - Check that the sample iframe example URL matches the project's GitHub Pages URL
   - Expected format: `https://<username>.github.io/<repo-name>/sims/<microsim-name>/main.html`
   - The `<microsim-name>` in the URL must match the directory name
   - The `<repo-name>` should be extracted from `mkdocs.yml` or detected from git remote
   - If mismatch found: Add TODO to fix iframe example URL
   ```

2. Add scoring to the rubric table:
   ```markdown
   |iframe URL correct|The sample iframe URL matches the project's GitHub Pages path|3|
   ```

3. Implementation details:
   - Extract `site_url` from `mkdocs.yml` to determine the correct base URL
   - Parse the iframe src attribute from the code block
   - Verify the URL contains the correct repository name
   - Verify the microsim directory name in the URL matches the actual directory

**Detection Command:**
```bash
# Find incorrect paths (adjust exclusions as needed)
grep -r "github.io" docs/sims/*/index.md | grep -v "<expected-repo-name>" | grep -v "external-references"
```

**Reference:**
See `/logs/check-iframe-example-path.md` in the intro-to-physics-course project for the full analysis that identified this issue.
