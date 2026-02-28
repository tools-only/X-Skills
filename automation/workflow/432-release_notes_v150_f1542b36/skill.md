# Lár v1.5.0 Release Notes
*(Major Feature Release - February 19, 2026)*

**"The Fractal Update"** 

v1.5.0 introduces **Fractal Agency** and **True Parallelism**, transforming Lár from a linear workflow engine into a system capable of massive, concurrent, self-organizing complexity.

##  Major Features

### 1. Fractal Agency (Recursive Dynamic Nodes)
Agents can now design **Sub-Agents** at runtime. A `DynamicNode` can instantiate another `DynamicNode` (or a whole graph) which executes and returns results up the chain.
- **Use Case:** A "Manager" agent spawns a "Coder" agent and a "Reviewer" agent on the fly to solve a problem it wasn't pre-programmed for.
- **Implementation:** `DynamicNode` now supports recursive instantiation and execution.

### 2. True Parallelism (`BatchNode`)
Lár now supports **multithreaded execution**. The `BatchNode` spins up independent threads for each branch of your graph.
- **Performance:** Run 5+ sub-agents simultaneously (e.g., searching web, writing code, checking logs) instead of sequentially.
- **Safety:** Each thread gets a deep-copied `State`, preventing race conditions. Merging is handled automatically.

##  Examples
- **NEW:** `examples/advanced/fractal_polymath.py`
  - A comprehensive showcase where a Manager orchestrates a "Cipher Specialist" and a "Poet" in parallel to solve a complex task.

##  Internal Changes
- **Tests Reorganized:** The `tests/` directory is now structured into `unit`, `integration`, and `scenarios` for better maintainability.
- **Core Patch:** `node.py` and `dynamic.py` updated to support `_pending_concurrent_ids` and recursive execution loops.

---
*Lár v1.5.0 represents a paradigm shift. Your agents are no longer just chains; they are living organisms that can grow.*
