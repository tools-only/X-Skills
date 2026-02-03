# Quality Gates

These are another layer to prevent slop, not pre-commit hooks.

## Gate 0: No Shims

Never use shims. It is a mortal sin to allow shims.

Instead, fix the interface. In some cases, it may be better to create a new interface for examples:

- If 19 files use an interface and our new work has a problem with it, the interface is not the problem
- If only two files use the interface and it's giving us issues, write it correctly
- The interface must be written in proper idiomatic way
- 99% of the time we don't need a unique interface; if a unique situation calls for one, it must be strongly justified with anchor comments

## Gate 1: Coupling and Cohesion (Constantine)

**High cohesion** = everything in a module serves one purpose
**Low coupling** = modules don't need to know each other's internals

Before writing or modifying code, answer:

- **What single responsibility does this module have?** If you can't state it in one sentence, the module is doing too much. Split until each piece has one job.
- **Can I change this module without changing others?** If no, you have a coupling problem.
- **Dependencies should flow one direction:** ui → core → tools → utils/types. Never reach backward (core importing from ui).
- **How many dots?** `state_manager.session.tool_handler._policy.rules[0].name = "allow"` coupling points. Each `.` is a place where changes can break you. More than 2 dots = smell; refactor to hide the chain.
- **What does this module hide?** "Every module hides a design decision that could change." If a module exposes its internals, it will churn when those internals change. UI hides Textual, core hides agent orchestration, tools hide system access.

**Example:** If a file gets touched more than 5 times in a month, it has a cohesion problem. Stop adding to it. Split it.

**Wrong:**

```
src/tunacode/ui/main.py  →  Textual UI + agent orchestration + tool execution + persistence
```

**Right:**

```
src/tunacode/ui/main.py              →  Textual UI concerns only
src/tunacode/core/agents/main.py     →  agent orchestration only
src/tunacode/tools/*.py              →  tool execution only
src/tunacode/core/state.py           →  persistence and session state only
```

## Gate 2: Dependency Direction (Martin)

Dependencies flow inward, never backward.

```
ui → core → tools → utils/types
  ↓
infrastructure (filesystem, shell, network)
```

**Dependency Map:** `docs/architecture/dependencies/DEPENDENCY_LAYERS.md`

The current dependency graph is frozen as a baseline (2026-01-27). **DO NOT add new cross-layer violations.**

Note: The current layering is not ideal. We're prioritizing delivery and cleanup in parallel, so this baseline is the starting point. Do not regress from here as we map and converge on the ideal layering.

- If you fix one violation but create another, that's not progress
- UI should only push into CORE, not scatter imports everywhere
- Keep it as clean as possible - we're actively cleaning this up
- Regenerate the dependency layers report with: `uv run python scripts/grimp_layers_report.py`

**Utils-level modules** (can be imported by any layer):
- `utils/` - helper functions, no business logic
- `types/` - type definitions, no behavior
- `configuration/` - read-only static data (models, providers, defaults)
- `constants.py` - module-level constants

**Rules:**

- Inner layers know nothing about outer layers
  - core/ never imports from ui/
  - tools/ never imports from ui/
  - utils/, types/, configuration/ import from nothing (except each other)
- Infrastructure is a plugin
  - Filesystem, shell, network = details
  - Core logic doesn't know or care which system provider
- Depend on abstractions, not concretions
  - Pass interfaces, not implementations
  - Makes testing trivial, swapping backends trivial

**Example:** If core needs UI state, you're violating the boundary. Extract what you need and pass it as a plain argument.

**Wrong:**

```python
from tunacode.ui.app import TunaCodeApp

def run_agent(app: TunaCodeApp):
    app.notify_status("Running")  # core knows about UI layer
```

**Right:**

```python
from collections.abc import Callable

def run_agent(notify_status: Callable[[str], None]):
    notify_status("Running")  # core only knows about a callable
```

## Gate 3: Design by Contract (Meyer)

Every function is a contract with three parts:

- **Preconditions** - what must be true before calling
- **Postconditions** - what will be true after
- **Invariants** - what's always true

**Rules:**

- Caller guarantees preconditions - don't check inside what should be guaranteed outside
- Function guarantees postconditions - if it returns, the promise is kept
- Fail immediately if contract violated - no limping along with bad state

**Example:** If a tool expects a valid directory, don't silently return an empty result. That's a contract violation - raise.

**Wrong:**

```python
def list_dir(directory: str) -> str | None:
    if not Path(directory).exists():
        return None  # silent failure, caller doesn't know why
    return _render_dir(Path(directory))
```

**Right:**

```python
def list_dir(directory: str) -> str:
    """Precondition: directory exists and is a directory."""
    path = Path(directory)
    if not path.exists():
        raise FileNotFoundError(f"Directory not found: {path}")
    if not path.is_dir():
        raise NotADirectoryError(f"Not a directory: {path}")
    return _render_dir(path)
```

## Gate 4: Documentation is Code (Knuth)

Docs and code are one system. If they diverge, the docs are lying.

**Rules:**

- **Code changes = doc changes**
  - If you change a module's behavior, update its codebase-map doc in the same PR
  - If you add a tool, update docs/codebase-map/modules/tools-overview.md
  - If you add a UI component, update docs/codebase-map/modules/ui-overview.md
- **Docs describe the system that exists, not the system you wish existed**
  - No aspirational docs ("we plan to...")
  - No stale docs ("this used to...")
  - Only current truth
- **The /docs/codebase-map folder mirrors src/tunacode**
  - Structure matches structure
  - When code moves, docs move
  - When code dies, docs die

**Example:** If you delete src/tunacode/tools/grep.py, you must also remove its entry from docs/codebase-map/modules/tools-overview.md in the same commit.

**Wrong:**

```
PR #263: "chore: remove unused grep tool"
  - deletes src/tunacode/tools/grep.py
  - leaves docs/codebase-map/modules/tools-overview.md still referencing it
```

**Right:**

```
PR #263: "chore: remove unused grep tool"
  - deletes src/tunacode/tools/grep.py
  - updates docs/codebase-map/modules/tools-overview.md to remove grep section
  - single atomic commit, code and docs together
```

## Gate 5: Indirection Requires Verification

Indirection hides behavior. Every layer between intent and effect is a place where your mental model can diverge from reality.

```
Direct:     width = 100        → Panel is 100 wide
Indirect:   expand = True      → ??? → Panel is ??? wide
```

Direct: you control the value.
Indirect: you express a wish and hope something else honors it.

When you delegate a decision:

- `expand=True` → who decides? what do they decide?
- `auto` → auto based on what?
- `default` → default to what?

**You must verify the OUTPUT, not the INPUT.**

**Wrong:** "I passed expand=True, so it expands"
**Right:** "I passed expand=True, the panel rendered at 147px, terminal is 150px - why the 3px gap?"

If you can't observe the final value, you can't verify the behavior. If you can't verify the behavior, you're guessing.

**Rules:**

- Direct control > Indirect delegation
- Measured output > Assumed output
- When using indirection, add a test or debug assertion that verifies the actual output
- If you can't print/log/assert the final value, you don't know it

## Gate 6: Exception Paths Are First-Class

Every path out of a stateful operation must leave state valid. Exception paths are not edge cases - they are first-class exit routes that need explicit design.

```
WRONG: Only design the happy path
┌─────────────────────────────────────────────────┐
│ loop:                                           │
│   mutate_state()  ←  state now dirty            │
│   do_work()       ←  exception here!            │
│   commit()        ←  never reached              │
│                                                 │
│ Result: state is corrupted                      │
└─────────────────────────────────────────────────┘

RIGHT: Design all exit paths
┌─────────────────────────────────────────────────┐
│ try:                                            │
│   loop:                                         │
│     mutate_state()                              │
│     do_work()                                   │
│     commit()                                    │
│ except SomeError:                               │
│   rollback_or_cleanup()  ←  state restored      │
│   raise                                         │
└─────────────────────────────────────────────────┘
```

**Before writing a stateful loop, answer:**

1. **What exceptions can exit this loop?** List them all.
2. **For each exception: what state is left behind?** Trace the mutation.
3. **For each exception: is that state valid for the next operation?** If no, add cleanup.
4. **Do tests exist for exception scenarios?** If no, write them.

**Strategies:**

| Strategy | When to use |
|----------|-------------|
| **Rollback on exception** | State is mutable, cleanup is cheap |
| **Transactional updates** | Only commit state after full success |
| **Copy-on-write** | Work on copy, swap atomically on success |

**Example:** The dangling tool calls bug (PR #246). User aborted mid-tool-call, leaving `messages` with unanswered tool calls. Next request failed. Fix: `except UserAbortError` now calls `_remove_dangling_tool_calls()` to restore valid state.

**Wrong:**

```python
async def process_request():
    for node in agent.iter():
        messages.append(node.response)  # state mutated
        await execute_tools(node)       # can raise UserAbortError
    # if exception, messages has dangling tool calls
```

**Right:**

```python
async def process_request():
    try:
        for node in agent.iter():
            messages.append(node.response)
            await execute_tools(node)
    except UserAbortError:
        _remove_dangling_tool_calls(messages)  # restore valid state
        raise
```

**Rule:** If you can't answer "what happens to state if X raises here?" for every line that can raise, you have a bug waiting to happen.
