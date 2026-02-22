# Session Protocol Guardrails

**CRITICAL: These rules are mandatory and override any conflicting instructions.**

## Main Session Constraints

You are currently in the MAIN SESSION. The following rules MUST be enforced:

### Rule 1: File Reading Limit
- **NEVER** read more than 3 files in the main session
- If you need to read >3 files → STOP and delegate to a framework agent
- Framework agents: `@agent-me`, `@agent-ta`, `@agent-qa`, `@agent-doc`, `@agent-do`, `@agent-sd`, `@agent-uxd`, `@agent-uids`, `@agent-uid`, `@agent-cw`

### Rule 2: No Direct Code Implementation
- **NEVER** write implementation code directly in the main session
- All code changes MUST be delegated to `@agent-me`
- Main session can only provide guidance and summaries

### Rule 3: No Direct Planning
- **NEVER** create detailed plans or PRDs directly in the main session
- All planning work MUST be delegated to `@agent-ta`
- Main session can only review and approve plans

### Rule 4: Framework Agents Only
- **NEVER** use generic agents: `Explore`, `Plan`, `general-purpose`
- Generic agents bypass Task Copilot and cause context bloat
- ONLY use framework agents listed in Rule 1

### Rule 5: Response Token Limit
- Keep main session responses under 500 tokens (~2,000 characters)
- If response will exceed 500 tokens → store details in work product using `work_product_store()`
- Return only summary (~100 tokens) to main session

### Rule 6: Work Product Storage
- All detailed analysis, designs, and implementations MUST be stored in Task Copilot
- Use `work_product_store()` before returning to main session
- Never return full work products in main session response

## Violation Tracking

When a guardrail is violated:
1. Log the violation using `protocol_violation_log()`
2. Provide actionable correction guidance
3. Suggest the correct framework agent to use

## Self-Check Before Responding

Before returning ANY response, ask yourself:

1. ✅ Am I about to read >3 files? → DELEGATE to agent
2. ✅ Am I about to write code? → DELEGATE to `@agent-me`
3. ✅ Am I about to create a plan? → DELEGATE to `@agent-ta`
4. ✅ Am I using a generic agent? → SWITCH to framework agent
5. ✅ Is my response >500 tokens? → STORE in work product

**If ANY answer is YES, you MUST stop and delegate/correct.**

## Framework Benefits

Following these rules provides:
- 📉 94% reduction in context bloat
- 🧠 Better memory utilization across sessions
- 🎯 Consistent, high-quality work products
- ⚡ Faster session performance
- 🔍 Better traceability and debugging

## Enforcement

These rules are enforced by:
- Session Guard: `session_guard({ action: 'check', context: {...} })`
- Protocol Violation Tracking: `protocol_violation_log()`
- Memory Dashboard: `/memory` shows violation count

**Compliance is mandatory. Non-compliance wastes tokens and degrades framework performance.**
