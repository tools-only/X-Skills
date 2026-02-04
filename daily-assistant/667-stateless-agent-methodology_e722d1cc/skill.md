# Stateless Agent Methodology (SAM)

A constraint-driven development framework that compensates for LLM limitations through architectural structure rather than behavioral instructions.

---

This methodology has companion documentation at:

- **Canonical specification**: [Stateless Software Engineering Framework](./stateless-software-engineering-framework.md)
- **Prompt/automation helper**: [SAM framework generator](./sam-framework-generator.md)

---

## The Core Insight

**Claude is not a knowledge worker - Claude is a stateless computation engine** that happens to have noisy, stale priors baked in.

Treat Claude like a pure function:

- **Input**: Complete context (task file with all answers)
- **Output**: Verified result
- **No side effects**: Fresh context each time
- **No memory**: Everything externalized to artifact files

---

## Part 1: The Problem

This section is maintained in the canonical framework document:

- [Stateless Software Engineering Framework: Part 1 (The Problem)](./stateless-software-engineering-framework.md#part-1-the-problem)

---

## Part 2: The Architectural Solution

This section is maintained in the canonical framework document:

- [Stateless Software Engineering Framework: Part 2 (Architectural Solution)](./stateless-software-engineering-framework.md#part-2-the-architectural-solution)

---

## Part 3: The Pipeline Architecture

This section is maintained in the canonical framework document:

- [Stateless Software Engineering Framework: Pipeline architecture + artifact flow](./stateless-software-engineering-framework.md#part-2-the-architectural-solution)
- [Stateless Software Engineering Framework: Agent specifications (artifact templates)](./stateless-software-engineering-framework.md#part-3-agent-specifications)

---

## Part 4: Input Constraints

This section is maintained in the canonical framework document:

- [Stateless Software Engineering Framework: Appendix H (SAM Input Constraints)](./stateless-software-engineering-framework.md#appendix-h-sam-input-constraints-canonicalized)

---

## Part 5: Stage Specifications

This section is maintained in the canonical framework document:

- [Stateless Software Engineering Framework: Agent specifications (stage roles + templates)](./stateless-software-engineering-framework.md#part-3-agent-specifications)

---

## Further Reading (canonical)

- [Framework: Why this works](./stateless-software-engineering-framework.md#part-5-why-this-works)
- [Framework: Theoretical foundations](./stateless-software-engineering-framework.md#part-4-theoretical-foundations)
- [Framework: Anti-patterns](./stateless-software-engineering-framework.md#appendix-b-anti-patterns-to-avoid)
- [Framework: Success metrics](./stateless-software-engineering-framework.md#appendix-c-success-metrics)

---

## Summary

SAM treats the model as a stateless computation engine and enforces correctness via externalized artifacts + verification gates. The detailed normative specification lives in the framework document linked above.
