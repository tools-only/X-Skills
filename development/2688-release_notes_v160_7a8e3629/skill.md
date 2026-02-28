# Lár v1.6.0: Memory Compression & Economic Constraints 

We are excited to announce **Lár v1.6.0**, introducing powerful new defensive mechanisms for mission-critical Agents.

## The Problem
As Agents become more autonomous (like the `DynamicNode` introduced in v1.5), they run the risk of two things:
1. **Context Bloat:** Parallel workers returning huge documents that overwhelm the downstream context window (The "Black Hole").
2. **Infinite Loops:** An agent getting stuck in a self-correction loop and burning through API credits.

## The Solution

### 1. Memory Compression (`ReduceNode`)
The new `ReduceNode` natively supports the Map-Reduce pattern to prevent context window bloat.
It reads multiple heavy documents from the `GraphState`, generates a concise executive summary, and then **explicitly deletes the raw data keys from memory.**
This ensures the state "baton" remains light and prevents downstream models from choking or hallucinating on massive context sizes.

### 2. Economic Constraints (Token Budgets)
You can now set a strict, mathematical ceiling on an agent's execution cost. 
Simply add a `token_budget` integer to the `initial_state`. Every time an `LLMNode` runs, it deducts its usage from this budget. If the budget hits zero, the node refuses to execute. You can now mathematically guarantee an agent will never exceed a specific dollar amount, regardless of how complex its dynamic graph becomes.

### 3. Node Fatigue (Infinite Loop Protection)
The `GraphExecutor` now tracks precisely how many times a node is visited during a single execution trace. If an agent loops over the same node (like a "Fix Code" node) more times than the `max_node_fatigue` limit (default: 20), the executor trips a circuit breaker and safely terminates the run, logging a fatigue error.

## New Showcase Example
- Added `examples/advanced/11_map_reduce_budget.py` to demonstrate `ReduceNode`, map-reduce parallel fan-out, and token budget exhaustion.

## Minor Fixes
- Added `ReduceNode` to `lar/__init__.py` for top-level import.
