## [5.0.0] - 2026-02-10 - "Antigravity Workflows Foundation"

> First-class Workflows are now available to orchestrate multiple skills through guided execution playbooks.

### 🚀 New Skills

### 🧭 [antigravity-workflows](skills/antigravity-workflows/)

**Orchestrates multi-step outcomes using curated workflow playbooks.**
This new skill routes users from high-level goals to concrete execution steps across related skills and bundles.

- **Key Feature 1**: Workflow routing for SaaS MVP, Security Audit, AI Agent Systems, and Browser QA.
- **Key Feature 2**: Explicit step-by-step outputs with prerequisites, recommended skills, and validation checkpoints.

> **Try it:** `Use @antigravity-workflows to run ship-saas-mvp for my project.`

---

## 📦 Improvements

- **Workflow Registry**: Added `data/workflows.json` for machine-readable workflow metadata.
- **Workflow Docs**: Added `docs/WORKFLOWS.md` to distinguish Bundles vs Workflows and provide practical execution playbooks.
- **Trinity Sync**: Updated `README.md`, `docs/GETTING_STARTED.md`, and `docs/FAQ.md` for workflow onboarding.
- **Go QA Path**: Added optional `@go-playwright` wiring in QA/E2E workflow steps.
- **Registry Update**: Catalog regenerated; repository now tracks 714 skills.

## 👥 Credits

A huge shoutout to our community and maintainers:

- **@Walapalam** for the Workflows concept request ([Issue #72](https://github.com/sickn33/antigravity-awesome-skills/issues/72))
- **@sickn33** for workflow integration, release preparation, and maintenance updates

---

_Upgrade now: `git pull origin main` to fetch the latest skills._
