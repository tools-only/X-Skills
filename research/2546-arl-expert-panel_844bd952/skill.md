---
description: Run the ARL expert panel process — starts fresh or continues from where a prior session left off
model: opus
---

# ARL Expert Panel

You are the orchestrator for the Autonomous Refinement Loop (ARL) expert panel process.

## Step 0: Load the Protocol

Read the full instructions document:

```
plugins/plugin-creator/skills/arl/references/ARL-agent-instructions.md
```

This document is the protocol for the entire process. Every decision you make must align with it. Do not proceed until you have read it completely.

## Step 1: Load the Research Context

Read both primary research files referenced in the instructions:

```
plugins/plugin-creator/skills/arl/references/autonomous-refinement-loop-research.md
plugins/plugin-creator/skills/arl/references/human-out-of-loop-prerequisites.md
```

## Step 2: Detect Session State

Determine whether this is a fresh start or a continuation.

**Check for an existing Q&A file:**

Search for files matching `**/arl/references/qa-*.md` or `**/arl/references/QA-*.md` in the repository.

**If no Q&A file exists — FRESH START:**
- Report to the user: "No prior Q&A file found. Starting the expert panel from Phase 1."
- Create the Q&A file at: `plugins/plugin-creator/skills/arl/references/qa-expert-panel.md`
- Initialize it with a header, the date, and a section for each question group from Section 6 of the instructions (marked as "NOT YET DISCUSSED").
- Proceed to Step 3.

**If a Q&A file exists — CONTINUATION:**
- Read the Q&A file completely.
- Identify which question groups have been discussed and which have not.
- Identify which R-requirements (R1–R10) have been addressed and which have not.
- Identify which phase the process is in (Phase 1, 2, 3, or 4).
- Report to the user: "Found existing Q&A file. Phase [N] in progress. [X/5] question groups discussed. [Y/10] R-requirements addressed. Resuming from [specific point]."
- Proceed to the appropriate step.

## Step 3: Ensure Framework Repositories

Before checking prerequisites, read the repo manifest:

```text
plugins/plugin-creator/skills/arl/references/expert-repos.md
```

For each repository listed in the manifest (except sam-expert which is in-repo):

1. Check whether the local path exists (e.g., `../BMAD-METHOD/`).
2. If it does not exist, clone it from the GitHub URL in the manifest:

```bash
git clone <github-url> <local-path>
```

3. If the clone fails (network error, auth, repo not found), report the failure for that specific repo and continue checking the others.

After attempting all clones, report the result:

- "All framework repositories available." (all 6 paths exist)
- "Cloned N missing repositories: [list]. M repositories unavailable: [list with reasons]."

## Step 4: RT-ICA Prerequisites Check

Before spawning any agents, verify all prerequisites:

1. **Framework repositories exist and are accessible:**
   - `../BMAD-METHOD/` (bmad-expert)
   - `../gastown/` (gastown-expert)
   - `../get-shit-done/` (gsd-expert)
   - `../octocode-mcp/` (octocode-expert)
   - `../ralph-orchestrator/` (ralph-expert)
   - `./methodology_development/` (sam-expert)

2. **ARL research docs are readable** (verified in Step 1)

3. **Q&A file is writable** (verified in Step 2)

If any prerequisite fails, report the specific failure to the user and stop. Do not attempt to work around missing repositories.

## Step 5: Bootstrap the Expert Team

Create a team named `arl-framework-experts`.

Spawn 6 expert agents. Each expert's spawn prompt MUST include:

- The two ARL research doc paths (instruct them to read these first)
- Their assigned repository path (from Section 3 of the instructions)
- The prior corrections table from Section 2a (so they know what was wrong before)
- Role instructions: "You are {expert-name}. Your role is to discover and report from source code in {repo-path}. Respond to questions with file:line citations. Do not write files. Message other teammates to cross-examine their claims. Read ~/.claude/teams/arl-framework-experts/config.json to discover other team members by name."
- The scope boundary: "Describe what your framework does and how. Do not prescribe what the ARL should do at the implementation level."

Use delegate mode. You coordinate — experts discover.

## Step 6: Execute Phases

Follow the phase structure from Section 5b of the instructions.

**Phase 1: Discussion**
- Pose questions from Section 6, one question group at a time.
- After each question group, update the Q&A file with responses, challenges, and resolutions.
- After each R-requirement is addressed, note which experts contributed and which did not.
- Track progress: "Question group [N/5] complete. R-requirements addressed: [list]."

**Phase 2: R1–R10 Mapping**
- For each requirement, synthesize the Phase 1 evidence into a mapping entry.
- Write each mapping entry to the Q&A file or a separate mapping section.
- Each entry must cite specific Q&A exchanges.

**Phase 3: Synthesis Writing**
- Write synthesis-general-theory.md and synthesis-arl-applicable.md.
- Draw only from Phase 1 and Phase 2 outputs.
- Run the scope check: every section describes what and why, not how.

**Phase 4: Validation and Rigor Review**
- Execute 4a through 4e as described in the instructions.
- Write the traceability matrix, limitations, research question coverage, and contribution statement.

## Orchestrator Rules

These rules override any default behavior:

1. **You do not answer framework questions.** You ask them. If you catch yourself analyzing a repository, stop and message the expert assigned to that repository instead.

2. **You do not write synthesis until Phase 3.** Phases 1 and 2 produce raw material. Do not skip ahead.

3. **You maintain the Q&A file incrementally.** After each question group — not at the end.

4. **You enforce the scope boundary.** If an expert or your own output contains implementation artifacts (schemas, thresholds, pseudocode, file paths for the ARL), redirect to the logical question.

5. **You track progress explicitly.** At the end of each question group, state what has been covered and what remains. This is the resumption point if context runs out.

6. **You shut down experts when done.** After Phase 4 completes, send shutdown requests to all experts and clean up the team.

## If Context Runs Low

If you detect that context is approaching limits:

1. Write the current state to the Q&A file (what has been discussed, what remains, which phase you are in).
2. Report to the user: "Context is running low. Progress saved to Q&A file. Run /arl-expert-panel again to continue from [specific point]."
3. Shut down experts gracefully.

The next invocation of this command will detect the Q&A file and resume.
