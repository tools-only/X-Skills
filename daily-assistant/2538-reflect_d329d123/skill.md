# Reflect Command

Review and process unreviewed corrections captured by the correction detection system.

## Overview

The `/reflect` command implements the second stage of the two-stage correction workflow:
1. **Auto-Capture** (Stage 1): Corrections are automatically detected in user messages
2. **Manual Review** (Stage 2): User confirms, rejects, or refines corrections via this command

## Arguments

- No arguments: Show all pending corrections
- `--agent <id>`: Filter by agent (me, ta, qa, doc, etc.)
- `--target <type>`: Filter by target (skill, agent, memory, preference)
- `--all`: Include already-reviewed corrections
- `--dedupe`: Consolidate similar corrections before display

## Step 1: Retrieve Corrections

Call these tools to gather correction data:

1. **Get pending corrections:**
   ```
   correction_list({ status: 'pending', limit: 20 })
   ```

2. **Get correction statistics:**
   ```
   correction_stats()
   ```

3. **If --agent flag provided:**
   ```
   correction_list({ status: 'pending', agentId: '<agent-id>' })
   ```

4. **If --target flag provided:**
   ```
   correction_list({ status: 'pending', target: '<target-type>' })
   ```

5. **If --all flag provided:**
   ```
   correction_list({ includeExpired: true, limit: 50 })
   ```

## Step 2: Display Corrections Summary

Format the output as a review dashboard:

```
## Correction Review

**Statistics:**
| Pending | Approved | Rejected | Applied | Expired |
|---------|----------|----------|---------|---------|
| 5       | 12       | 3        | 10      | 2       |

### Pending Corrections (5)

---
**#1** | ID: `abc123` | Confidence: 0.87 | Target: skill
**Pattern:** explicit_correction
**Original:** "use var for declarations"
**Corrected:** "use const for declarations"
**Context:** "Correction: use const instead of var for declarations"
**Agent:** @agent-me | **Created:** 2025-01-15

**Actions:** [Approve] [Reject] [Modify]

---
**#2** | ID: `def456` | Confidence: 0.72 | Target: preference
**Pattern:** preference
**Original:** (not specified)
**Corrected:** "tabs over spaces"
**Context:** "I prefer tabs over spaces"
**Agent:** (global) | **Created:** 2025-01-14

**Actions:** [Approve] [Reject] [Modify]

---
```

## Step 3: Interactive Review

For each correction, prompt the user:

```
Review correction #1 (abc123)?
1. Approve - Store as confirmed correction
2. Reject - Mark as false positive
3. Modify - Edit before approving
4. Skip - Review later
5. Apply - Approve and immediately apply to target

Choice [1-5 or Enter to skip]:
```

### Handle User Response

**Approve (1):**
```
correction_update({ correctionId: 'abc123', status: 'approved' })
```
Output: "Correction #1 approved. Will be applied during next skill/agent update."

**Reject (2):**
Ask for rejection reason, then:
```
correction_update({
  correctionId: 'abc123',
  status: 'rejected',
  rejectionReason: '<user-provided-reason>'
})
```
Output: "Correction #1 rejected: <reason>"

**Modify (3):**
Prompt user for modified correction content, then create new correction with modified content and approve it.

**Skip (4):**
Move to next correction without changing status.

**Apply (5):**
```
correction_update({ correctionId: 'abc123', status: 'applied' })
```
Then delegate to appropriate agent to apply the correction:
- If target is 'skill': Route to @agent-me to update skill file
- If target is 'agent': Route to @agent-me to update agent file
- If target is 'memory': Store as lesson via memory_store
- If target is 'preference': Store as preference via memory_store

## Step 4: Handle Deduplication (--dedupe)

If `--dedupe` flag is present:

1. Group corrections by target + correctedContent similarity
2. Show grouped view:

```
### Similar Corrections (grouped)

**Group 1:** "use const instead of var" (3 occurrences)
- abc123: Confidence 0.87, 2025-01-15
- def789: Confidence 0.82, 2025-01-14
- ghi012: Confidence 0.75, 2025-01-13

**Actions:** [Approve All] [Reject All] [Review Individually]
```

3. If "Approve All", update all corrections in group to approved

## Step 5: Summary

After review session:

```
## Review Session Complete

**Processed:** 5 corrections
- Approved: 3
- Rejected: 1
- Applied: 1
- Skipped: 0

**Next Steps:**
- Run `/reflect --all` to see full correction history
- Approved corrections will be applied during next agent invocation
- Use `correction_stats()` to track correction patterns over time
```

## Edge Cases

### No Pending Corrections

```
## Correction Review

No pending corrections to review.

**Statistics:**
- Total captured: 25
- Applied: 20
- Rejected: 5

Use `/reflect --all` to see correction history.
```

### All Corrections Expired

```
## Correction Review

No active corrections. All pending corrections have expired.

**Tip:** Corrections expire if not reviewed promptly.
Consider reviewing corrections more frequently.
```

### High Confidence Auto-Apply

If a correction has confidence >= 0.95 and autoStore was true:

```
**Note:** This correction was auto-stored due to high confidence (0.95).
It will be applied automatically unless you reject it.
```

## Correction Targets Explained

| Target | Where Applied | Agent Responsible |
|--------|---------------|-------------------|
| skill | `.claude/skills/*.md` files | @agent-me |
| agent | `.claude/agents/*.md` files | @agent-me |
| memory | Stored as lesson/decision | Memory Copilot |
| preference | Stored as user preference | Memory Copilot |

## Display Notes

- Show corrections in chronological order (newest first)
- Truncate long content to 100 characters with "..."
- Color-code confidence: >= 0.85 (green), 0.6-0.84 (yellow), < 0.6 (red)
- Show pattern type that triggered detection
- Include raw context (user message) for verification
- Highlight corrections that are about to expire (< 24h remaining)

## End

Present the review dashboard and begin interactive review process.
