# Auditor's Guide: Inspecting Lár Agents

> **Audience**: Compliance Officers, Quality Assurance (QA), and external Auditors.
> **Scope**: This guide explains how to verify the behavior of a "High-Risk AI System" built with the Lár engine, in accordance with **EU AI Act Article 13 (Transparency)**.

---

## 1. The "Glass Box" Concept

Most AI systems are "Black Boxes"—a prompt goes in, and an answer comes out, with no visibility into the intermediate steps.

**Lár is different.** It is a "Glass Box." 
Every Lár agent is a **Graph** (a flowchart) of explicit steps. Nothing happens "magically."

*   **Nodes**: The boxes in the flowchart (Steps).
*   **Edges**: The arrows connecting them (Logic).
*   **State**: The memory passed between them.

As an auditor, you have the right to inspect all three.

---

## 2. The Artifacts

For every execution of a Lár Agent, the system produces two artifacts. You should request these from the engineering team for any incident you are investigating.

### A. The Graph Definition (`architecture.mermaid`)

This is the map of *what could happen*. It is usually a Mermaid diagram or a Python file defining the nodes.

**What to look for:**

* **Loops**: Are there infinite loops? (Lár automatically caps loops at 100 steps to prevent this).
* **Gates**: Is there a "Human-in-the-Loop" before sensitive actions (e.g., "Refund Tool")?
* **Separation**: Is the "Reasoning" (LLM) separated from the "Action" (Tool)?

### B. The Flight Recorder (`state.json` / `flight_recorder.json`)

This is the map of *what actually happened*.

Lár generates a **State-Diff Ledger**. It does not just dump logs; it records exactly what changed at every step.

**Example Trace:**

```json
[
  {
    "step": 1,
    "node": "TriageNode",
    "diff": { "classification": "high_risk" }
  },
  {
    "step": 2,
    "node": "SupervisorNode",
    "diff": { "next_action": "escalate_to_human" }
  }
]
```

**Verification Steps:**

1. **Traceability**: Can you follow the `step` numbers 1, 2, 3...?
2. **Causality**: Did Step 2 happen *because* Step 1 output "high_risk"?
3. **Completeness**: Are there any gaps in the sequence?

---

## 3. Proving Determinism

The "Golden Rule" of Lár is **Determinism**.
*   **Rule**: If you provide the same Input and the same Random Seed (for the LLM), the agent MUST produce the exact same Output.

**The Test:**

1. Ask the engineer for the **Input State** of the incident.
2. Ask them to re-run the agent with that Input.
3. Compare the new `flight_recorder.json` with the original.
4. **They must match.**

If they do not match, the system is not compliant with "Repeatability" standards. Lár is designed to pass this test.

---

## 4. Common Failure Modes

As an auditor, look for these common "Anti-Patterns":

*   **"The Hallucinating Jury"**: In a Juried Layer (Proposer -> Jury), did the Jury ignore a policy violation? Check the `JuryNode` log in the flight recorder.
*   **"The Magic Loop"**: Did the agent spin in a loop (Planner -> Executor -> Planner) without making progress? Check if the `step` count hit the limit (100).
*   **"The Silent Fail"**: Did a tool fail (e.g., API Error) but the LLM ignored it? Check the `ToolNode` output in the JSON.
