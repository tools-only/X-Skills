# Lár v1.1 Release Notes: "The Metacognition Update"

We are proud to announce **Lár v1.1**, the "Metacognition Update". This release introduces Level 4 Agency capabilities, allowing agents to modify their own execution logic at runtime while maintaining strict "Glass Box" compliance.

## Major Features

### 1. Dynamic Graphs (Metacognition)
Agents can now inspect and modify their own topology during execution.
- **`DynamicNode`**: A new primitive that allows a node to return a *new subgraph* or *topology modification* instead of just a state update.
- **`TopologyValidator`**: A deterministic validation engine that ensures any self-assembled graph is valid, safe, and terminates ensuring compliance with safety constraints.
- **Audit-Logged Self-Modification**: Every change to the graph structure is recorded in the `State-Diff Ledger`, preserving the "Glass Box" guarantee even for self-modifying code.

### 2. New Examples & "Metacognition" Labs
We have added 5 advanced examples demonstrating dynamic graph capabilities:
- **`25_dynamic_depth.py`**: An agent that decides *how deep* to think based on the question's complexity.
- **`26_tool_inventor.py`**: An agent that writes its own Python tools on the fly to solve novel problems.
- **`27_self_healing.py`**: A pipeline that detects failures (e.g., auth errors), dynamically inserts a "fix" node (e.g., credential rotation), and resumes execution.
- **`28_adaptive_deep_dive.py`**: An agent that spawns parallel sub-agents to research topics it identifies as important.
- **`29_expert_summoner.py`**: An agent that dynamically instantiates specialized "expert" nodes based on the user's query topic.

### 3. Compliance & Safety Updates
- **EU AI Act & FDA Compatibility**: Updated compliance documentation to address "Self-Modifying AI". By using `TopologyValidator` and the immutable `State-Diff Ledger`, Lár v1.1 ensures that even dynamic behavior is fully traceable and reproducible.
- **Safety First**: New documentation alerts and guidelines for deploying dynamic agents in production environments.

## Improvements

- **Example Reorganization**: All examples have been moved to `lar/examples/` and organized into clear categories:
  - `basic/`: Fundamental primitives (Linear, Branching).
  - `patterns/`: Common agentic patterns (RAG, ReAct, Map-Reduce).
  - `compliance/`: Safety and auditing examples (Human Jury, Audit Logs).
  - `scale/`: High-performance and distributed examples.
  - `metacognition/`: The new v1.1 dynamic graph examples.
- **Website Updates**:
  - **Red Teaming Hook**: Added a prominent hook to the Compliance section of the homepage.
  - **Doc Improvements**: Fixed documentation rendering for alerts/admonitions and audited all links.
  - **Metacognition Section**: Added a dedicated section for dynamic graphs.
- **Version Bump**: Framework version updated to `1.1.0`.

## Upgrading

To upgrade to Lár v1.1, simply run:
```bash
pip install lar-engine==1.1.0
```

## Documentation
- [New Metacognition Guide](https://docs.snath.ai/core-concepts/9-metacognition)
- [Compliance & Safety](https://docs.snath.ai/compliance)
