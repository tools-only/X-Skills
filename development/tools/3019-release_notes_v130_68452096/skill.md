# Lár v1.3.0 Release Notes

**"The Compliance & Architecture Update"**

This release focuses on hardening the framework for enterprise and regulatory compliance (EU AI Act), while refactoring the core execution engine for better separation of concerns.

## Key Features

### 1. "Human-in-the-Loop" Primitive (`HumanJuryNode`)
A new node type that pauses execution to request explicit human feedback via the CLI.
- **Why**: Directly satisfies EU AI Act Article 14 ("Human Oversight") requirements.
- **How**:
  ```python
  jury = HumanJuryNode(
      prompt="Approve sensitive action?",
      choices=["approve", "reject"],
      output_key="approval_status"
  )
  ```

### 2. Static Safety Analysis (`TopologyValidator`)
We've added a `TopologyValidator` (backed by NetworkX) that runs comprehensive checks on **Dynamic Graphs** *before* they execute.
- **Cycle Detection**: Catches infinite loops in generated subgraphs.
- **Structural Integrity**: Validates `next_node` references.
- **Tool Allowlisting**: Enforces strict boundaries on what tools a dynamic agent can access.

### 3. Core Refactor: Modular Observability
The `GraphExecutor` has been refactored to delegate responsibilities to dedicated components, keeping the engine lightweight.
- **`AuditLogger`**: Centralizes audit trail logging and file persistence.
- **`TokenTracker`**: Aggregates token usage across multiple providers and models with precision.

## Usage Updates

### Breaking Changes
- `GraphExecutor` constructor now accepts optional `logger` and `tracker` instances for dependency injection.
- **Compliance**: The "Glass Box" is now even more transparent with improved metadata fidelity in logs.

### How to Use (New in v1.3.0)

#### Option 1: Automatic (Default Behavior)
```python
from lar import GraphExecutor, LLMNode

# Logger and Tracker are created automatically
executor = GraphExecutor(log_dir="my_logs")

node = LLMNode(model_name="ollama/phi4", prompt_template="test", output_key="result")
result = executor.run(node, {})

# Access automatically created instances
print(executor.logger.get_history())  # Audit trail
print(executor.tracker.get_summary())  # Token usage
```

#### Option 2: Custom Injection (Advanced)
```python
from lar import GraphExecutor, AuditLogger, TokenTracker

# Create custom instances
custom_logger = AuditLogger(log_dir="advanced_logs")
custom_tracker = TokenTracker()

# Inject into executor
executor = GraphExecutor(
    logger=custom_logger,
    tracker=custom_tracker
)

# Share tracker across multiple executors for aggregated cost tracking
executor2 = GraphExecutor(
    logger=AuditLogger(log_dir="other_logs"),
    tracker=custom_tracker  # Same tracker = aggregated tokens
)
```

**Why Custom Injection?**
- Centralized audit trail management
- Cost aggregation across workflows
- Custom log formatting/persistence
- Integration with existing monitoring systems

**See:** `examples/patterns/16_custom_logger_tracker.py` for full demo


## Chains to State Machines
With `v1.3`, the debate is settled. Lár's state machine architecture now offers native compliance features (Breakpoints, Static Analysis) that Chain-based frameworks struggle to implement.

## Changelog
- **[NEW]** `HumanJuryNode` in `src/lar/node.py`
- **[NEW]** `TopologyValidator` in `src/lar/dynamic.py`
- **[NEW]** `AuditLogger` in `src/lar/logger.py`
- **[NEW]** `TokenTracker` in `src/lar/tracker.py`
- **[REFACTOR]** `GraphExecutor` to use new helper classes.
