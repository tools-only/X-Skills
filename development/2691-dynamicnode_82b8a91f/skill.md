# DynamicNode API Reference

## Overview

`DynamicNode` enables **metacognitive agents** that can modify their own behavior at runtime. It asks an LLM to generate a JSON graph specification, validates it for safety, and then executes the generated subgraph.

> [!CAUTION]
> **Self-Modifying Code is Risky**: `DynamicNode` introduces inherent security risks. Always use with `TopologyValidator` to enforce safety policies.

## Class Signature

```python
class DynamicNode(BaseNode):
    def __init__(
        self,
        llm_model: str,
        prompt_template: str,
        validator: TopologyValidator,
        next_node: BaseNode = None,
        context_keys: List[str] = [],
        system_instruction: str = None
    )
```

## Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `llm_model` | `str` | Yes | Model to generate the graph JSON (e.g., "gpt-4", "gemini-pro") |
| `prompt_template` | `str` | Yes | Prompt asking the LLM to design a graph. Must include schema instructions. |
| `validator` | `TopologyValidator` | Yes | Safety validator to enforce security policies (cycle detection, tool allowlisting) |
| `next_node` | `BaseNode` | No | Node to execute after the dynamic subgraph completes |
| `context_keys` | `List[str]` | No | State keys to include in the LLM's context when designing the graph |
| `system_instruction` | `str` | No | System prompt for the LLM |

## Behavior

### Execution Flow

1. **Generate Graph JSON**: Call LLM with `prompt_template` and `context_keys`
2. **Parse JSON**: Extract graph specification from LLM response
3. **Validate**: Run `TopologyValidator` to check for cycles, unauthorized tools, etc.
4. **Build**: Instantiate nodes from the JSON spec
5. **Execute**: Run the generated subgraph
6. **Continue**: Proceed to `next_node`

### JSON Graph Schema

The LLM must generate JSON matching this format:

```json
{
  "nodes": {
    "node_1": {
      "type": "LLMNode",
      "config": {
        "model_name": "gpt-4",
        "prompt_template": "Analyze: {input}",
        "output_key": "analysis"
      },
      "next_node": "node_2"
    },
    "node_2": {
      "type": "ToolNode",
      "config": {
        "tool_name": "approved_tool",  // Must be in validator allowlist
        "input_keys": ["analysis"],
        "output_key": "result"
      },
      "next_node": null
    }
  },
  "start_node": "node_1"
}
```

## Safety & Validation

### Why Validation is Critical

Self-modifying code can:
- Create infinite loops
- Call blacklisted APIs
- Exfiltrate data
- Execute arbitrary code

### TopologyValidator Protections

1. **Cycle Detection**: Ensures graph is a DAG (no infinite loops)
2. **Tool Allowlisting**: Only permits pre-approved ToolNode functions
3. **Structural Integrity**: Validates all `next_node` references exist

See [`TopologyValidator` API](topologyvalidator.md) for details.

## Example Usage

### 1. Adaptive Depth (Variable Complexity)

```python
from lar import DynamicNode, TopologyValidator, GraphExecutor

# Define allowed tools
def simple_search(query): 
    return f"Quick result for {query}"

def deep_research(query):
    return f"Detailed analysis of {query}"

# Create validator
validator = TopologyValidator(allowed_tools=[simple_search, deep_research])

# Create metacognitive node
adapter = DynamicNode(
    llm_model="gpt-4",
    prompt_template="""
    Based on the query complexity, design a graph:
    - Simple query: Use simple_search
    - Complex query: Use deep_research

    Query: {user_query}

    Output JSON graph spec with 'nodes' and 'start_node'.
    """,
    validator=validator,
    context_keys=["user_query"]
)

executor = GraphExecutor()
result = list(executor.run_step_by_step(adapter, {"user_query": "What is 2+2?"}))
```

**What Happens**:
- For "What is 2+2?" → LLM generates 1-node graph with `simple_search`
- For "Explain quantum entanglement" → LLM generates 3-node graph with `deep_research` + synthesis

### 2. Tool Inventor (Runtime Code Generation)

```python
def execute_generated_function(code: str):
    """Executes LLM-generated Python code (DANGEROUS - use with caution)"""
    exec(code)  # Extreme trust required

validator = TopologyValidator(allowed_tools=[execute_generated_function])

inventor = DynamicNode(
    llm_model="gpt-4",
    prompt_template="""
    The user needs a tool to: {task_description}

    Generate Python code for a function that does this, then create a ToolNode graph to execute it.

    Return JSON with:
    - 'code': The Python function
    - 'nodes': Graph spec using execute_generated_function
    """,
    validator=validator,
    context_keys=["task_description"]
)
```

> [!WARNING]
> **Code Execution Risk**: This pattern executes arbitrary code. Only use in sandboxed environments or with strict human-in-the-loop oversight.

### 3. Self-Healing (Error Recovery)

```python
from lar import DynamicNode, TopologyValidator

# Allow error correction tools
validator = TopologyValidator(allowed_tools=[retry_with_backoff, use_fallback_api])

healer = DynamicNode(
    llm_model="gpt-4",
    prompt_template="""
    An error occurred: {last_error}

    Design a recovery subgraph using allowed tools:
    - retry_with_backoff
    - use_fallback_api

    Return JSON graph to fix the error.
    """,
    validator=validator,
    context_keys=["last_error"]
)
```

## Audit Trail

`DynamicNode` logs **the exact generated JSON** in the audit trail:

```json
{
  "step": 5,
  "node": "DynamicNode",
  "state_diff": {
    "added": {
      "_generated_graph_spec": {...},  // Full JSON logged
      "_dynamic_result": "Success"
    }
  }
}
```

This ensures **full compliance auditability** - you can always see what the agent decided to do.

## Best Practices

### ✅ Do

1. **Always use TopologyValidator** - No exceptions
2. **Limit tool allowlist** - Principle of least privilege
3. **Log generated graphs** - For debugging and compliance
4. **Use system instructions** - Guide the LLM's design choices
5. **Add human oversight** - For high-risk decisions

### ❌ Don't

1. **Allow unrestricted `exec()`** - Extreme security risk
2. **Skip validation** - Enables injection attacks
3. **Use in production without testing** - Test generated graphs first
4. **Allow network tools in validator** - Data exfiltration risk

## Real-World Examples

See the [`examples/metacognition/`](https://github.com/snath-ai/lar/blob/main/examples/metacognition/) directory:

| Example | Capability |
|---------|------------|
| [`1_dynamic_depth.py`](https://github.com/snath-ai/lar/blob/main/examples/metacognition/1_dynamic_depth.py) | Adaptive complexity (1 node vs N nodes) |
| [`2_tool_inventor.py`](https://github.com/snath-ai/lar/blob/main/examples/metacognition/2_tool_inventor.py) | Self-coding (writes tools at runtime) |
| [`3_self_healing.py`](https://github.com/snath-ai/lar/blob/main/examples/metacognition/3_self_healing.py) | Error recovery (Injects fix subgraphs) |
| [`4_adaptive_deep_dive.py`](https://github.com/snath-ai/lar/blob/main/examples/metacognition/4_adaptive_deep_dive.py) | Recursive research (spawns sub-agents) |
| [`5_expert_summoner.py`](https://github.com/snath-ai/lar/blob/main/examples/metacognition/5_expert_summoner.py) | Dynamic persona instantiation |

## Compliance

### EU AI Act Article 13 (Transparency)

`DynamicNode` satisfies transparency requirements by:
- Logging the exact generated graph JSON
- Recording which tools were invoked
- Providing deterministic validation results

### Security Auditing

Every `DynamicNode` execution creates an audit entry with:
- `input`: Context keys and prompt
- `output`: Generated JSON graph
- `validation_result`: Pass/fail and reason

## See Also

- [TopologyValidator API](topologyvalidator.md) - Safety enforcement
- [Metacognition Guide](../core-concepts/9-metacognition.md) - Deep dive into self-modifying agents
- [Red Teaming](../case-studies/red-teaming.md) - Security testing dynamic graphs
