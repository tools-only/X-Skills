
The workflow begins the moment you run `/zerg:init`. What happens next depends entirely on what ZERG finds.

ZERG's orchestrator spawns a lightweight reconnaissance subagent using Claude's Task tool. This subagent has read-only permissions. It scans the directory for project markers: package.json, pyproject.toml, Cargo.toml, go.mod, Makefile, Dockerfile, any signal that code already lives here. The scan takes seconds. The subagent returns a structured assessment: either "existing project with characteristics X, Y, Z" or "empty directory, nothing to analyze."

That fork determines everything that follows.

---

**The Empty Directory Path: Inception Mode**

When reconnaissance returns empty, ZERG shifts from discovery to conversation. There's nothing to analyze, so it asks questions instead.

ZERG opens a dialogue about what you're building. Not technical specifications. The problem itself. Who uses this system? What problem does it solve? What constraints exist? Regulatory requirements? Performance targets? Integration points with existing systems? Your answers get captured in `.gsd/PROJECT.md` as a project charter. This document becomes the source of truth for all downstream decisions.

From the charter, ZERG asks about technology. Do you have preferences? Organizational mandates? Existing infrastructure you must integrate with? If you have opinions, ZERG records them. If you're open to recommendations, ZERG makes them based on the problem domain. A REST API serving mobile clients might get Python FastAPI with PostgreSQL. A real-time collaboration tool might get TypeScript, Next.js, and WebSockets. A CLI utility might get Go for single-binary distribution. ZERG explains its reasoning so you can override with cause.

Once technology is locked, ZERG spawns a scaffolding subagent. This subagent generates the entire project skeleton in a single atomic operation. Directory structure following community conventions for the chosen stack. Configuration files for the language runtime, package manager, and build tools. A README with setup instructions, architecture overview, and contribution guidelines. A Dockerfile for the application itself, not just development. CI pipeline configuration for GitHub Actions or GitLab CI. License file. Gitignore tuned to the stack.

In parallel, ZERG generates the devcontainer. It selects a base image appropriate to the stack. Python projects get the Microsoft Python devcontainer image with the detected Python version. Node projects get the Node image. Polyglot projects get a composite. The devcontainer configuration includes VS Code extensions for the languages involved, debugger configurations, and port forwarding for common development servers.

The devcontainer also receives secure coding rules. ZERG pulls from TikiTribe/claude-secure-coding-rules public github repo based on the technology choices. A Python backend gets OWASP Top 10 rules, SQL injection prevention patterns, and authentication best practices. A React frontend gets XSS prevention rules, CSP configuration guidance, and secure state management patterns, etc. These rules get installed in the .claude/rules/ directory or subdirectory CLAUDE.md files to organize modular, relevant rules inside the devcontainer, so every worker that runs in that container inherits security constraints automatically. No explicit prompting required.

Code quality tooling installs in the same pass. Linters configured with strict settings: ESLint for TypeScript, Ruff for Python, clippy for Rust. Formatters with opinionated defaults: Prettier, Black, rustfmt. Static analyzers with security focus: Semgrep with OWASP patterns, Bandit for Python, npm audit for Node. These tools get cached in the devcontainer image layer so worker startup is fast.

ZERG doesn't stop at infrastructure. With an empty directory, there's no code to extend. So it continues into backlog generation automatically. The scaffolding subagent takes the charter and produces an initial product backlog. These aren't implementation tasks yet. They're user stories at the epic level. User authentication. Data persistence. Core domain logic. API surface. Administrative interface. Each story has a rough priority based on dependency order and risk. High-risk items surface early because you want to validate uncertain assumptions before building on them.

ZERG presents everything for approval. The charter summary. The technology choices with rationale. The project structure. The initial backlog. You review. You modify. You challenge the recommendations. ZERG adjusts. Only after explicit approval does ZERG commit the scaffolding. It initializes a git repository if one doesn't exist, creates the initial commit with the skeleton, and tags it as the inception point.

At this moment, the empty directory contains a runnable project. You can build it. You can run it. The tests execute but fail because there's no implementation. The structure exists. The tooling works. The security rules are active. Sprint Zero is complete.

---

**The Existing Directory Path: Discovery Mode**

When reconnaissance returns a populated assessment, ZERG has material to work with. It doesn't ask what you're building. It figures that out from what exists.

The reconnaissance subagent's report includes detected languages, frameworks, test runners, build systems, and deployment configurations. ZERG spawns additional discovery subagents in parallel, each focused on a different aspect. One maps the directory structure and identifies architectural patterns. One reads configuration files to understand dependencies and their versions. One scans for existing tests and measures coverage. One checks for security configurations, environment variable handling, and secret management.

These subagents run concurrently using Claude's Task parallelism. Each operates in its own context window with read-only tool access. They can't modify anything. They report findings to the orchestrator, which synthesizes a comprehensive project profile.

The profile feeds devcontainer generation. Unlike inception mode, ZERG doesn't choose technologies. It matches them. If the project uses Python 3.11, the devcontainer gets Python 3.11. If it uses Node 18, Node 18. ZERG detects the test framework and includes its runner. It finds the linter configuration and preserves those settings rather than imposing defaults.

Secure coding rules still get injected, but they're filtered to the detected stack. A Django project gets Django-specific rules in addition to Python core rules. A NestJS project gets NestJS patterns alongside TypeScript basics. The rules don't conflict with existing project conventions because they focus on security constraints, not style preferences.

Code quality gates configure themselves from existing tooling. If the project already has ESLint, ZERG uses that configuration. If it has pytest with coverage thresholds, ZERG adopts those thresholds. ZERG augments but doesn't replace. It adds Semgrep security scanning if absent. It adds pre-commit hooks if missing. It fills gaps without disrupting established workflows.

ZERG captures the infrastructure assessment in `.gsd/INFRASTRUCTURE.md`. This document records what exists, what ZERG added, and what assumptions it made. If you disagree with an assumption, you modify the document and ZERG adjusts.

Discovery mode ends with a working devcontainer that matches the existing project. You can open the project in the container and develop normally. All existing scripts work. All existing tests run. ZERG added security and quality tooling, but the core development experience is unchanged.

---

**Convergence: Planning Forward**

Both paths converge at `/zerg:plan`. Whether you started from nothing or from a mature codebase, planning works identically.

You describe what you want to build next. A new feature. A refactor. A bug fix. A migration. ZERG spawns a product owner subagent that engages in requirements elicitation. It asks clarifying questions. It identifies edge cases. It extracts acceptance criteria. It validates that the scope is achievable and testable.

The output is a feature specification in `.gsd/specs/{feature}/REQUIREMENTS.md`. This document contains the problem statement, user stories with acceptance criteria, constraints, and verification commands that prove the feature works. The spec is human-readable and version-controlled. It becomes the contract that implementation must satisfy.

For projects that started empty, the first planning session typically pulls from the initial backlog. You refine those rough stories into precise specifications. For existing projects, planning addresses whatever you need next, whether that's new functionality or technical debt.

---

**Design: Task Graph Generation**

`/zerg:design` transforms specifications into executable work.

ZERG spawns an architect subagent that reads all specs for the feature and produces the technical blueprint. This subagent makes structural decisions: where new files live, what interfaces look like, how components interact, what the data model requires. It produces architecture documentation in `.gsd/specs/{feature}/ARCHITECTURE.md`.

The architect also generates the task graph. Each task is atomic: one clear objective, one verification command, exclusive file ownership. No two tasks modify the same file. Dependencies between tasks determine the wave structure.

Level 0 tasks are foundations. Type definitions. Interface declarations. Configuration schemas. Database migrations. These have no dependencies on other tasks in the feature.

Level 1 tasks implement core logic against Level 0 interfaces. Service classes. Repository implementations. Business rule engines.

Level 2 tasks build on Level 1. API endpoints that call services. Event handlers that trigger workflows. Integration points with external systems.

Level 3 tasks are consumers. Frontend components that call APIs. CLI commands that orchestrate services. Background jobs that process queues.

Level 4 tasks verify the whole. Integration tests. End-to-end tests. Performance benchmarks. Documentation updates.

The architect subagent also generates stub files. Interfaces with method signatures but no implementation. Type definitions with structures but no logic. These stubs compile and type-check. They form the contract that workers implement against. Workers in later levels can import these stubs and write code that satisfies the interfaces without knowing how other workers will implement them.

For projects that started empty, the first design session produces substantial scaffolding because there's no existing code to extend. For existing projects, design produces only the delta needed for the feature.

---

**Execution: Parallel Workers in Isolated Containers**

`/zerg:rush` launches the workforce.

The orchestrator reads the task graph and identifies all Level 0 tasks. For each task, it provisions a worker. Provisioning means:

First, creating a git worktree. The orchestrator runs `git worktree add .zerg-worktrees/{feature}/worker-{N} -b zerg/{feature}/worker-{N}`. Each worker gets its own branch and its own directory. Changes in one worktree don't affect others until merge.

Second, launching a devcontainer instance. The orchestrator runs `docker-compose` with environment variables that identify the worker: `ZERG_WORKER_ID`, `ZERG_FEATURE`, `ZERG_TASK_ID`. The container mounts the worker's worktree as `/workspace`. It mounts `.gsd/specs/` as read-only so the worker can read specifications but not modify them. It connects to the `zerg-net` Docker network for inter-worker communication if needed.

Third, spawning a Claude Code subagent inside the container. The orchestrator uses Claude's Task tool to create a specialized subagent. The subagent's system prompt comes from `.claude/agents/worker.md`, which includes the secure coding rules, the quality standards, and the execution protocol. The subagent receives the task specification as its initial prompt.

Workers execute in parallel. Claude's Task tool supports up to 10 concurrent subagents. If you have more than 10 Level 0 tasks, they queue and execute in batches. Workers don't know about each other. They read specifications, write code, run verification commands, and report results.

Each worker follows the same protocol. Read the task specification. Understand the acceptance criteria. Read any files listed as dependencies. Write the implementation. Write tests that prove the acceptance criteria. Run the verification command. If it passes, commit to the worker branch. If it fails, analyze the failure, fix the code, and retry. After three failures, the worker marks the task as blocked and halts.

Workers can't spawn their own subagents. Claude's Task tool prohibits subagent nesting. This keeps execution predictable. The orchestrator maintains the single point of coordination.

---

**Wave Boundaries: Quality Gates and Merges**

When all Level 0 workers complete, the orchestrator runs quality gates before any Level 1 work begins.

The gate sequence is fixed. First, the orchestrator pulls all worker branches. It runs `git merge-base --is-ancestor` checks to verify branches haven't diverged unexpectedly.

Second, the orchestrator merges worker branches into a staging branch using sequential fast-forward merges. If exclusive file ownership held, merges succeed without conflict. If there's a conflict, something went wrong in task decomposition. The orchestrator halts and reports the conflict for human review.

Third, quality checks run against the merged staging branch. Linters execute with the project's configured rules. Type checkers verify that interfaces match implementations. Test suites run with coverage measurement. Security scanners check for vulnerabilities. Each check produces a structured report.

Fourth, the orchestrator evaluates gate results against thresholds. If coverage dropped below 80%, the gate fails. If security scanners found high-severity issues, the gate fails. If linters reported errors (not warnings), the gate fails.

Failed gates produce specific feedback. The orchestrator identifies which worker produced the failing code and which rule was violated. It spawns a remediation subagent for that worker with the failure context injected. The remediation subagent attempts to fix the issue. If it succeeds, the orchestrator reruns the affected checks. If it fails three times, the task marks as blocked.

Passed gates trigger the merge to main. The orchestrator runs `git checkout main && git merge --no-ff staging -m "Level 0: {feature}"`. The commit message lists all tasks completed in the level. The staging branch resets for the next wave.

Now Level 1 workers can start. They pull from main, which contains all Level 0 artifacts. They can import the types and interfaces their predecessors created. The dependency ordering is enforced by wave boundaries, not by complex scheduling logic.

This pattern repeats for each level. Parallel execution within the wave. Quality gates at wave end. Merge before proceeding. Each wave builds on the verified foundation of previous waves.

---

**Coordination: Claude Tasks as Shared State**

Throughout execution, Claude's native Tasks feature provides the coordination layer.

When ZERG generates the task graph, it creates a Claude Task for each ZERG task. The Task includes structured metadata: dependencies, assigned worker, current status, verification results, retry count. All workers connect to the same Task list via the `CLAUDE_CODE_TASK_LIST_ID` environment variable, which is set to the feature name.

Workers query the Task list to find their assignment. They update Task status as they progress: `pending` to `in_progress` when starting, `in_progress` to `completed` when verification passes, `in_progress` to `blocked` when retries exhaust.

The orchestrator watches the Task list for state changes. When all tasks in a level reach `completed` or `blocked`, the orchestrator triggers the wave boundary logic. If any task is `blocked`, the orchestrator reports which tasks and why before deciding whether to proceed with partial completion or halt entirely.

This shared state mechanism means workers are stateless and restartable. If a worker crashes, the orchestrator sees the Task stuck in `in_progress`, spawns a replacement worker, and the Task resumes from the last known state. Workers don't carry conversation history. They read from specifications and write to Tasks. The spec files and Task list are the memory.

---

**Completion: Deliverable Software**

The final wave is always integration. These tasks have dependencies on all previous levels because they verify the whole system works together.

Integration workers run end-to-end tests that exercise complete user flows. They generate coverage reports that aggregate unit and integration coverage. They update documentation with API references and usage examples. They verify the CI pipeline passes. They produce a deployment artifact: a Docker image, a compiled binary, a deployable package.

When the integration wave passes quality gates and merges, the feature is complete. The orchestrator generates a summary report. Tasks completed with time elapsed. Test coverage achieved. Security findings addressed. Token usage across all workers. The report goes to `.gsd/specs/{feature}/COMPLETION.md`.

You have working software. Not prototype software. Working, tested, secure, documented software. The git history shows clean wave-by-wave commits, each verified before merge. You can trace any line of code back to the task that created it, the specification that required it, and the acceptance criteria that validated it.

---

**Iteration: Sprint Over Sprint**

For the next sprint, infrastructure persists. The devcontainer exists. The quality gates are configured. The secure coding rules are active. You don't rebuild the foundation.

You run `/zerg:plan` with new stories. Maybe they come from the initial backlog if you started empty. Maybe they come from user feedback if you shipped the previous sprint. The product owner subagent refines them into specifications.

You run `/zerg:design`. The architect subagent reads existing code alongside new specs. It extends the architecture rather than replacing it. It generates tasks that modify existing files where necessary, always ensuring exclusive ownership within each wave.

You run `/zerg:rush`. Workers execute against the current codebase. They build on what exists. Each sprint adds a vertical slice of functionality. Each slice is tested, secure, and deployable.

The cycle repeats. Plan, design, execute, deliver. Parallel execution within sprints compresses timelines. Wave boundaries enforce dependency ordering. Quality gates catch regressions before they merge. Secure coding rules prevent vulnerabilities at generation time rather than detection time.

Whether you started from an empty directory or a mature codebase, the steady state is identical. ZERG adapts to what exists, provisions what's missing, and executes in parallel against a verified foundation. The inception cost differs. The ongoing workflow converges.

</worflow narrative>

The output will be a .md file that is optimized for consumption by claude code so that there is no mistake in how zerg should operate.