# Pattern: Tool/Extension Evaluation Framework

**Pattern ID:** `tool-evaluation-framework-001`
**Category:** Decision Making
**Confidence:** High
**Last Updated:** 2025-11-01
**Reuse Count:** 0 (newly created)

---

## Pattern Summary

Systematic framework for deciding whether to install, keep, or remove development tools, MCP servers, Claude Extensions, and Skills.

**Use When:**
- Evaluating new tools/servers/extensions
- Optimizing installed components
- Reducing system overhead
- Choosing between alternatives

---

## Evaluation Criteria (5 Dimensions)

### 1. Practical Value
**Question:** Does it solve a real, current problem?

**Scale:**
- 🟢 **High:** Solves frequent, concrete problems (daily/weekly use)
- 🟡 **Medium:** Solves occasional problems (monthly use)
- 🔴 **Low:** Theoretical benefit, no observed use

**Examples:**
- 🟢 `code-runner` → Execute code directly (high value)
- 🟡 `doc-generator` → Generate docs when needed (medium value)
- 🔴 `memory-server` → Theoretical learning, no observable benefit (low value)

---

### 2. Transparency
**Question:** Can I see what it's doing? Is output inspectable?

**Scale:**
- 🟢 **Clear:** Actions visible, logs readable, output inspectable
- 🟡 **Partial:** Some visibility, some black box
- 🔴 **Opaque:** Black box, no insight into operations

**Examples:**
- 🟢 `browser` → Clear actions (navigate, click, scrape)
- 🟡 `sequential-thinking` → Partial visibility (thinking logs)
- 🔴 `memory-server` → Opaque (hidden knowledge graph)

---

### 3. Redundancy
**Question:** Does existing system already provide this capability?

**Scale:**
- 🟢 **Unique:** No overlap with existing tools
- 🟡 **Partial:** Some overlap but adds unique value
- 🔴 **Redundant:** Existing system handles it

**Examples:**
- 🟢 `notion` → Unique Notion integration
- 🟡 `everything` → Overlaps with native filesystem, but adds convenience
- 🔴 `memory-server` → Redundant with `.claude/projects/`, `AGENT_MEMORY.md`

---

### 4. Stability
**Question:** Does it work reliably without errors?

**Scale:**
- 🟢 **Stable:** No errors, fast startup
- 🟡 **Occasional Issues:** Rare errors, acceptable
- 🔴 **Problematic:** Frequent crashes, errors, warnings

**Check:**
```bash
# Log size (larger = more issues)
ls -lh logs/mcp-server-[name].log

# Error patterns
grep -i "error\|warning\|crash\|fail" logs/mcp-server-[name].log
```

**Examples:**
- 🟢 `browser` → Stable, minimal logs
- 🟡 `notion` → Occasional connection timeouts (acceptable)
- 🔴 `aws-api` → 256KB error logs, constant crashes

---

### 5. Overhead
**Question:** What's the cost (startup time, memory, complexity)?

**Scale:**
- 🟢 **Light:** Fast startup, low memory
- 🟡 **Medium:** Noticeable but acceptable
- 🔴 **Heavy:** Slow startup, high memory, complex config

**Examples:**
- 🟢 `code-runner` → Lightweight Node.js
- 🟡 `browser` → Heavier (Puppeteer) but justified
- 🔴 `aws-api` → Heavy Python + botocore

---

## Decision Matrix

| Score | Decision | Action |
|-------|----------|--------|
| 4-5 Green | **Keep** | Essential tool |
| 2-3 Green + rest Yellow/Green | **Keep** | Useful tool |
| Mixed with 1-2 Red | **Evaluate** | Consider alternatives |
| 3+ Red | **Remove** | Not worth the cost |

---

## Example Evaluations

### Example 1: AWS API MCP Server
| Criterion | Score | Reason |
|-----------|-------|--------|
| Practical Value | 🟡 Medium | Would be useful IF stable |
| Transparency | 🟢 Clear | AWS CLI is well-documented |
| Redundancy | 🟡 Partial | AWS CLI exists, but integration adds value |
| Stability | 🔴 Problematic | 256KB error logs, DeprecationWarnings |
| Overhead | 🔴 Heavy | Large Python dependencies |

**Verdict:** ❌ **Remove** (2 Red flags: Stability + Overhead)

---

### Example 2: Memory Server
| Criterion | Score | Reason |
|-----------|-------|--------|
| Practical Value | 🔴 Low | No observable benefit |
| Transparency | 🔴 Opaque | Hidden knowledge graph |
| Redundancy | 🔴 Redundant | `.claude/projects/`, `AGENT_MEMORY.md` |
| Stability | 🟢 Stable | No errors |
| Overhead | 🟡 Medium | Node.js, acceptable size |

**Verdict:** ❌ **Remove** (3 Red flags: Value + Transparency + Redundancy)

---

### Example 3: Continual Learning Skill
| Criterion | Score | Reason |
|-----------|-------|--------|
| Practical Value | 🟢 High | Captures reusable patterns |
| Transparency | 🟢 Clear | Markdown files, readable |
| Redundancy | 🟢 Unique | No existing pattern capture system |
| Stability | 🟢 Stable | Local files only |
| Overhead | 🟢 Light | Just markdown files |

**Verdict:** ✅ **Keep** (5 Green: Perfect fit)

---

## Quick Evaluation Template

```markdown
### [Tool Name] Evaluation

| Criterion | Score | Reason |
|-----------|-------|--------|
| Practical Value | 🟢/🟡/🔴 | [Why?] |
| Transparency | 🟢/🟡/🔴 | [Can I see what it does?] |
| Redundancy | 🟢/🟡/🔴 | [Unique or redundant?] |
| Stability | 🟢/🟡/🔴 | [Errors in logs?] |
| Overhead | 🟢/🟡/🔴 | [Cost?] |

**Verdict:** ✅ Keep / ❌ Remove / 🤔 Evaluate
**Reason:** [Primary factors]
```

---

## Special Cases

### Case 1: Experimental Tools
**Approach:** Install → Use for 1-2 weeks → Re-evaluate
**Threshold:** If no value observed after 2 weeks → Remove

### Case 2: High-Value but Unstable
**Options:**
1. Check for updates (may fix stability)
2. Configure more conservatively
3. Report issues to maintainer
4. Temporarily disable, check back later

### Case 3: Low-Value but Stable
**Decision:** Remove anyway
**Reason:** Every installed component adds:
- Startup time
- Complexity
- Maintenance burden
- Cognitive load

---

## Optimization Guidelines

1. **Target: 3-5 Essential Tools**
   - More than 10 → too much overhead
   - Less than 3 → missing key capabilities

2. **Prefer Single-Purpose over Multi-Purpose**
   - Easier to understand
   - Easier to debug
   - Easier to replace

3. **Transparent over "Smart"**
   - Readable logs/output
   - Predictable behavior
   - User-controllable

4. **Active Use over "Just in Case"**
   - If not used in 2 weeks → probably not needed
   - Exception: Disaster recovery tools

---

## Related Patterns

- `claude-desktop-troubleshooting` (debugging bad tools)
- `skill-selection-framework` (choosing Claude Code skills)
- `system-config-optimization` (general config cleanup)

---

## Version History

- **v1.0 (2025-11-01):** Initial framework from MCP server optimization session
