# fn-8 Explore OpenProse for Ralph orchestration

## Overview

Evaluate OpenProse (prose.md) as a potential replacement or enhancement for Ralph's bash loop orchestration. OpenProse is a DSL where the LLM interprets "prose programs" as a virtual machine, enabling smarter agent coordination without external frameworks.

## Context

From Twitter discussion with @irl_danB and @tiagoefreitas:
- OpenProse could be "the missing piece" for bitter-lesson aligned architecture
- Allows complex loops where edge cases break dumb bash
- Coordinator agent pattern: one agent manages VM context, sub-agents execute
- Heavy use of sub-agents for clean context separation

## Scope

Research and prototype:
- Understand OpenProse syntax and capabilities
- Evaluate fit for Ralph use cases
- Prototype a simple flow (e.g., single task with review)
- Compare with current bash approach

## Approach

```prose
# Example: Ralph-style task loop in OpenProse
agent implementer:
  model: sonnet
  skills: ["bash", "edit", "read"]

agent reviewer:
  model: opus
  skills: ["read"]

loop until **all tasks done or blocked** (max: 10):
  task = session: coordinator
    prompt: "Pick next ready task via flowctl"

  session: implementer
    prompt: "Implement {task}"
    context: { task }

  review = session: reviewer
    prompt: "Review the implementation"
    context: { task }

  if **review passed**:
    session: coordinator
      prompt: "Mark task done"
```

## Considerations

**Pros:**
- Smarter loop conditions (natural language)
- Built-in sub-agent spawning
- Context management via coordinator pattern
- No bash edge cases

**Cons:**
- Beta software with telemetry
- Requires Opus (expensive for long runs)
- Another abstraction layer
- May be overkill for simple loops

## Quick commands

- Visit https://prose.md for cloud playground
- Clone https://github.com/openprose/prose for local exploration

## Acceptance

- [ ] Understand OpenProse syntax and VM model
- [ ] Document pros/cons vs current Ralph approach
- [ ] Prototype single-task flow in OpenProse
- [ ] Decision: adopt, defer, or reject

## References

- [OpenProse GitHub](https://github.com/openprose/prose)
- [prose.md Cloud](https://www.prose.md/)
- [Dan's Twitter thread](https://x.com/irl_danB/status/2009871120892342707)
- Twitter discussion with @tiagoefreitas about coordinator agent pattern
