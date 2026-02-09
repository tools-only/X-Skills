# MIRA Roadmap

This document outlines planned development for MIRA. If you're interested in contributing, look for issues tagged with `help wanted` or `good first issue`.

## Want to Pitch In? Here are tasks that could be done and PR'd:
- Import/Export API endpoint that collects messages and memories and puts them in a file that can be imported. Put them on the /actions endpoint. 
- Identify dead code and include a reference to a grep that it is truly orphaned in the PR
- OpenWebUI integration. the hosted version has the bespoke UI but some basic web UI would help lower the barrier to entry for new users
- ADVANCED: writing tests. lord, if you can work through the directories writing pytests with strong precise assertations I will be forever grateful. 


## Q1 2025

### January: Codebase Cleanup

Focus on reducing technical debt and improving code quality.

- [ ] Remove dead code and unused imports
- [ ] Standardize error handling patterns
- [ ] Consolidate duplicate logic across modules
- [ ] Improve type annotations coverage
- [ ] Update and clean up documentation
- [ ] Refactor inconsistent naming conventions

**Contribution opportunities**: Code cleanup tasks are excellent entry points for new contributors.

---

### February: LT_Memory Refinement

Enhance the long-term memory system and self-reflection capabilities.

- [ ] Refine LT_Memory extraction and linking logic
- [ ] Adjust end-of-segment self-reflection step
- [ ] Implement blind scratchpad for MIRA observations
  - MIRA generates unstructured observations before structured extraction
  - Captures intuitions that might be lost in direct-to-structure processing
- [ ] Improve memory refinement pipeline

---

### March: TUI Enhancement

Build out a more complete terminal user interface.

- [ ] Expand TUI feature set
- [ ] Improve navigation and usability
- [ ] Add configuration options accessible via TUI
- [ ] Enhanced conversation display and controls

---

## How to Contribute

1. Check the [Issues](../../issues) page for tasks tagged `help wanted`
2. Read [CONTRIBUTING.md](CONTRIBUTING.md) for development setup
3. For larger changes, open an issue first to discuss the approach

Questions? Open a Discussion or reach out in Issues.
