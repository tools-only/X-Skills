# Troubleshooting

Common issues encountered when running ZERG, organized by symptom. Each entry follows a Problem / Cause / Solution format.

For deeper investigation of complex failures, see [[Debug Guide]].

---

## Table of Contents

- [Workers Not Starting](#workers-not-starting)
- [Docker and Container Issues](#docker-and-container-issues)
- [API Key and Authentication](#api-key-and-authentication)
- [Task Failures](#task-failures)
- [Port Conflicts](#port-conflicts)
- [State and Coordination Issues](#state-and-coordination-issues)
- [Git and Worktree Issues](#git-and-worktree-issues)
- [Performance Issues](#performance-issues)

---

## Workers Not Starting

### Workers launch but immediately exit

**Problem:** `/zerg:rush` reports that workers were spawned, but they terminate within seconds. No task progress is recorded.

**Cause:** The worker cannot find the spec files or task graph it needs to execute. Workers are stateless and rely entirely on spec files for their instructions.

**Solution:**

1. Verify the spec directory exists and is populated:
   ```bash
   ls .gsd/specs/<feature>/
   ```
2. Confirm that `task-graph.json` is present and valid:
   ```bash
   python -m json.tool .gsd/specs/<feature>/task-graph.json
   ```
3. Check the worker logs for the specific error:
   ```bash
   cat .zerg/logs/worker-*.stderr.log | tail -20
   ```
4. Re-run `/zerg:design` if the spec files are missing or corrupt.

### Workers hang on "in_progress" without producing output

**Problem:** `/zerg:status` shows workers as `in_progress`, but no files are being created or modified. The workers appear stalled.

**Cause:** Workers may be waiting on a dependency that was not properly resolved, or Claude Code itself may have hit a context limit.

**Solution:**

1. Check if there is a dependency deadlock in the task graph:
   ```bash
   python -c "import json; g=json.load(open('.gsd/specs/<feature>/task-graph.json')); [print(t['id'], '->', t.get('dependencies',[])) for t in g['tasks']]"
   ```
2. Review the worker logs for context-related errors:
   ```bash
   cat .zerg/logs/worker-*.stderr.log
   ```
3. Stop the stalled workers with `/zerg:stop` and retry with `/zerg:rush --resume`.

### Workers fail to claim tasks

**Problem:** Workers start but report "no tasks available" despite pending tasks existing.

**Cause:** The Task ecosystem is out of sync. Workers use `TaskUpdate` to claim tasks, and if the `CLAUDE_CODE_TASK_LIST_ID` is mismatched, workers cannot see the correct task list.

**Solution:**

1. Verify the task list ID matches across all workers. All workers for a feature must share the same `CLAUDE_CODE_TASK_LIST_ID`.
2. Run `/zerg:status` to check for mismatches between the Task system and `.zerg/state/` files.
3. If state is corrupted, stop all workers and re-launch with `/zerg:rush --resume`, which calls `TaskList` first to reconcile state.

---

## Docker and Container Issues

### Docker daemon not running

**Problem:** `/zerg:rush --mode container` fails with "Cannot connect to the Docker daemon."

**Cause:** The Docker Desktop application or Docker daemon service is not running.

**Solution:**

1. Start Docker Desktop (macOS/Windows) or the Docker service (Linux):
   ```bash
   # macOS
   open -a Docker

   # Linux
   sudo systemctl start docker
   ```
2. Verify Docker is accessible:
   ```bash
   docker info
   ```
3. Retry the rush command.

### Container authentication failure

**Problem:** Container workers start but fail with authentication errors. Claude Code inside the container cannot connect to the API.

**Cause:** Container workers authenticate via one of two methods -- OAuth (mounting `~/.claude`) or API key (environment variable). Neither is configured correctly.

**Solution:**

For OAuth (Claude Pro/Team accounts):
1. Verify `~/.claude` exists and contains valid session data.
2. Confirm the volume mount is working:
   ```bash
   docker run --rm -v ~/.claude:/home/user/.claude alpine ls /home/user/.claude
   ```

For API key:
1. Verify the `ANTHROPIC_API_KEY` environment variable is set:
   ```bash
   echo $ANTHROPIC_API_KEY
   ```
2. Ensure the key is valid and has not expired.
3. Check that the key is being passed into the container environment.

### Container image build fails

**Problem:** ZERG cannot build the worker container image. Errors occur during the Docker build step.

**Cause:** Network issues preventing package downloads, or missing Dockerfile context.

**Solution:**

1. Try building the image manually to see full error output:
   ```bash
   docker build -t zerg-worker .
   ```
2. If packages fail to download, check your network and any proxy settings.
3. Clear the Docker build cache and retry:
   ```bash
   docker builder prune
   ```

### Container mode substituted with subprocess

**Problem:** You specified `--mode container` but workers are running as subprocesses instead of Docker containers.

**Cause:** ZERG may fall back to subprocess mode if Docker is unavailable. Container mode is a first-class execution path and should not be silently substituted.

**Solution:**

1. Verify Docker is running (see above).
2. Check `.zerg/config.yaml` for the execution mode setting.
3. Re-run with explicit container mode: `/zerg:rush --mode container --workers=5`.
4. If the fallback persists, run `/zerg:debug --deep` to diagnose infrastructure issues.

---

## API Key and Authentication

### ANTHROPIC_API_KEY not set

**Problem:** Workers fail immediately with an authentication error mentioning a missing API key.

**Cause:** The `ANTHROPIC_API_KEY` environment variable is not set in the shell where ZERG is running.

**Solution:**

1. Set the key in your shell:
   ```bash
   export ANTHROPIC_API_KEY="sk-ant-..."
   ```
2. To persist across sessions, add it to your shell profile (`~/.zshrc`, `~/.bashrc`).
3. Verify it is set:
   ```bash
   echo $ANTHROPIC_API_KEY | head -c 10
   ```

### API rate limiting

**Problem:** Multiple workers start failing simultaneously with HTTP 429 errors.

**Cause:** Running many parallel workers can exceed your API rate limits. Each worker is an independent Claude Code session making API calls.

**Solution:**

1. Reduce the number of parallel workers:
   ```
   /zerg:rush --workers=3
   ```
2. Check your API plan's rate limits and adjust the worker count accordingly.
3. ZERG has built-in retry with backoff for transient rate limit errors, but sustained over-limit usage requires fewer workers.

---

## Task Failures

### Task verification command fails

**Problem:** A task's code was generated but the verification command (`verify` field in `task-graph.json`) returns a non-zero exit code.

**Cause:** The generated code has errors (syntax, import, logic) or the verification command itself is misconfigured.

**Solution:**

1. Check which task failed and its verification command:
   ```bash
   /zerg:status
   ```
2. Run the verification command manually to see the full error:
   ```bash
   # Example: the verify command from the task definition
   python -m py_compile path/to/file.py
   ```
3. If the code has errors, use `/zerg:retry <task-id>` to re-attempt the task.
4. If the verification command is wrong, fix it in `task-graph.json` and retry.

### Tasks stuck in "pending" after level completion

**Problem:** All tasks at the current level are complete, but the next level's tasks remain pending.

**Cause:** The merge step between levels may have failed, blocking level advancement. ZERG requires all tasks at level N to complete and merge before level N+1 begins.

**Solution:**

1. Check if the merge is pending or failed:
   ```
   /zerg:status
   ```
2. If the merge failed, check the quality gate output:
   ```bash
   cat .zerg/logs/merge-*.log
   ```
3. Manually trigger the merge: `/zerg:merge`.
4. Fix any quality gate failures (lint, typecheck, test errors) and retry.

### Dependency errors in generated code

**Problem:** Workers generate code that imports modules not present in the project.

**Cause:** Workers operate from spec files and may not have full visibility into the project's installed dependencies.

**Solution:**

1. Verify the project's dependencies are listed in `requirements.txt` or `package.json`.
2. Check if the design spec (`design.md`) correctly identifies required packages.
3. Install missing dependencies and retry the failed task.

---

## Port Conflicts

### Address already in use

**Problem:** Workers or services fail to start because a required port is already in use.

**Cause:** Another process (or a previous ZERG run that did not clean up) is occupying the port.

**Solution:**

1. Identify what is using the port:
   ```bash
   lsof -i :<port-number>
   ```
2. Stop the conflicting process, or configure ZERG to use a different port range in `.zerg/config.yaml`.
3. If leftover from a previous run, clean up with `/zerg:cleanup`.

### Port allocation exhaustion

**Problem:** ZERG reports it cannot allocate ports for new workers.

**Cause:** The configured port range has been exhausted by active or orphaned workers.

**Solution:**

1. Stop all workers: `/zerg:stop`.
2. Verify no orphaned processes remain:
   ```bash
   ps aux | grep claude
   ```
3. Expand the port range in `.zerg/config.yaml` if running many workers.
4. Run `/zerg:cleanup` to release any held resources.

---

## State and Coordination Issues

### State JSON and Task system disagree

**Problem:** `/zerg:status` reports mismatches between `.zerg/state/<feature>.json` and the Claude Code Task system.

**Cause:** A worker crashed or was killed before it could update both state stores. The Task system is the source of truth; state JSON files are supplementary.

**Solution:**

1. Trust the Task system output from `/zerg:status`.
2. Run `/zerg:rush --resume` to reconcile. The resume flag calls `TaskList` first and only creates tasks that do not already exist.
3. If corruption is severe, stop all workers, back up the state directory, and re-launch.

### Orphaned worktrees

**Problem:** `git worktree list` shows worktrees that are no longer associated with active workers.

**Cause:** Workers were terminated without cleaning up their git worktrees.

**Solution:**

1. List all worktrees:
   ```bash
   git worktree list
   ```
2. Remove orphaned worktrees:
   ```bash
   git worktree remove <path>
   ```
3. If the worktree is locked:
   ```bash
   git worktree unlock <path>
   git worktree remove <path>
   ```
4. Run `/zerg:cleanup` to automate this process.

---

## Git and Worktree Issues

### Merge conflicts during level advancement

**Problem:** The merge step between levels fails due to git merge conflicts.

**Cause:** Despite file ownership rules, merge conflicts can arise from configuration files, lock files, or overlapping generated content.

**Solution:**

1. Review the conflict:
   ```bash
   git diff --name-only --diff-filter=U
   ```
2. Resolve conflicts manually or use `/zerg:merge` with manual resolution.
3. If the task graph has overlapping file assignments, fix the design with `/zerg:design`.

### Worktree creation fails

**Problem:** Workers cannot create git worktrees. Errors mention the branch already exists or the path is already a worktree.

**Cause:** Previous run left behind branches or worktrees that were not cleaned up.

**Solution:**

1. Clean up stale worktrees:
   ```bash
   git worktree prune
   ```
2. Delete leftover branches:
   ```bash
   git branch -D zerg-worker-<id>
   ```
3. Run `/zerg:cleanup` before starting a new rush.

---

## Performance Issues

### Workers running slower than expected

**Problem:** Individual workers are taking much longer than estimated.

**Cause:** Large context windows, complex tasks, or system resource contention.

**Solution:**

1. Check system resources (CPU, memory, disk):
   ```bash
   top -l 1 | head -10
   df -h .
   ```
2. Reduce the number of parallel workers to decrease resource contention.
3. Enable context engineering in `.zerg/config.yaml` to reduce token usage:
   ```yaml
   plugins:
     context_engineering:
       enabled: true
       command_splitting: true
       task_context_budget_tokens: 4000
   ```

### Disk space exhaustion

**Problem:** Workers or git operations fail with "no space left on device."

**Cause:** Each worker creates a git worktree (a full copy of the working directory). Many workers on a large repo can exhaust disk space.

**Solution:**

1. Check available disk space:
   ```bash
   df -h .
   ```
2. Clean up completed worktrees: `/zerg:cleanup`.
3. Reduce the number of simultaneous workers.
4. Consider using container mode, which isolates filesystem usage.

---

## Still Stuck?

If none of the above solutions resolve your issue:

1. Run a full diagnostic: `/zerg:debug --deep --env`
2. Check the generated report in `claudedocs/debug-<timestamp>.md`
3. See [[Debug Guide]] for advanced investigation techniques
4. Review the [[Command Reference]] for correct syntax and flags
