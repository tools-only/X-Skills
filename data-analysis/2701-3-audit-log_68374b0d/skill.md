

# The Event-Sourced Audit Log

**The `GraphExecutor`'s "killer feature" is the structured history log it generates. This is not just a `print` statement; it's a forensic trail.**

For complex agent workflows, you need detailed provenance.

**Lár delivers:**
1.  **Immutability**: Every step is a discrete event.
2.  **Causality**: We log the *exact state change* (`state_diff`) that caused the new reality.
3.  **Reproducibility**: Replaying these events with the same seed guarantees the same outcome.

## The "Event Sourcing" Model

To be infinitely scalable, `lar` (as of v2.0) does not log the entire state at every step. That is slow and expensive.

Instead, it logs **"events" (state diffs)**.

The `GraphExecutor` `yields` a `step_log` object for each step. This log is a small, efficient JSON object:
```json
{
  "step": 0,
  "node": "LLMNode",
  "outcome": "success",
  
  "state_before": {
    "task": "What is the capital of France?"
  },
  
  "state_diff": {
    "added": {
      "plan": "TEXT"
    },
    "removed": {},
    "modified": {}
  },

  "run_metadata": {
    "prompt_tokens": 42,
    "output_tokens": 1,
    "total_tokens": 43
  }
}
```

**What This Log Tells You**

- `step` **&** `node`: Where you are in the graph.

- `outcome`: What happened ("success" or "error").

- `state_before`: A "snapshot" of the agent's entire memory before this node ran.

- `state_diff`: The "**killer feature.**" This is the exact change the node made to the state. You can see it `added` the `"plan": "TEXT" key`.

- `run_metadata`: The **cost audit**. You can see this single `LLMNode` step cost 43 tokens.

When your agent fails, the "glass box" gives you a perfect record:
```json
{
  "step": 3,
  "node": "LLMNode",
  "outcome": "error",
  "error": "429 LLMNode failed after 3 retries.",
  "state_before": { ... },
  "state_diff": { ... },
  "run_metadata": {}
}
```

This is the `lar` difference. You don't have to guess why it failed. Your log tells you the exact node, the exact error, and the exact state it failed with.

---

## New in v1.3.0: Modular Observability

**Lár now separates logging and tracking into dedicated, injectable components:**

### AuditLogger

Manages the audit trail with GxP-compliant JSON persistence.

**Default usage:**
```python
executor = GraphExecutor(log_dir="my_logs")
# AuditLogger created automatically
history = executor.logger.get_history()
```

**Advanced (custom injection):**
```python
from lar import AuditLogger

custom_logger = AuditLogger(log_dir="compliance_logs")
executor = GraphExecutor(logger=custom_logger)

# Direct access to audit history
for step in custom_logger.get_history():
    print(f"Step {step['step']}: {step['node']} - {step['outcome']}")
```

### TokenTracker

Aggregates token usage across models and workflows.

**Default usage:**
```python
executor = GraphExecutor()
summary = executor.tracker.get_summary()
print(f"Total tokens: {summary['total_tokens']}")
print(f"By model: {summary['tokens_by_model']}")
```

**Advanced (shared tracker for cost aggregation):**
```python
from lar import TokenTracker

# Share one tracker across multiple executors
shared_tracker = TokenTracker()

executor1 = GraphExecutor(log_dir="workflow1", tracker=shared_tracker)
executor2 = GraphExecutor(log_dir="workflow2", tracker=shared_tracker)

# Run both workflows
executor1.run(agent1, state1)
executor2.run(agent2, state2)

# Get aggregated cost
total = shared_tracker.get_summary()
print(f"Total cost across both workflows: {total['total_tokens']} tokens")
```

**See full example:** `examples/patterns/16_custom_logger_tracker.py`