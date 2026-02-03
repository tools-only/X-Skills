---
title: Core Module
path: tunacode/core
type: directory
depth: 1
description: Agent orchestration, state management, and prompting system
seams: [agents, state, prompting, setup]
---

# Core Module (`src/tunacode/core`)

## Where
`src/tunacode/core/` - Central business logic layer for the TunaCode agent system.

## What
Implements the **agent orchestration engine** that powers TunaCode's AI capabilities:
- **Agent System**: Multi-agent architecture with specialized roles (`agents/`)
- **State Management**: Centralized application state (`state/`)
- **Prompting**: Dynamic prompt construction and template system (`prompting/`)
- **Setup Flow**: Onboarding and configuration wizards (`setup/`)
- **Background Tasks**: Async job execution and monitoring (`background/`)
- **Token Usage**: Cost tracking and usage analytics (`token_usage/`)
- **Logging**: Structured logging system (`logging/`)

## Directory Structure
```
core/
├── agents/
│   ├── main.py                    # Primary agent orchestration
│   ├── research_agent.py          # Codebase research specialist
│   ├── delegation_tools.py        # Agent tool delegation logic
│   └── agent_components/          # Agent building blocks
│       ├── (state, tools)
├── prompting/
│   ├── (prompt builders)
├── state/
│   ├── (application state models)
├── setup/
│   ├── (setup wizard logic)
├── background/
│   ├── (async task execution)
├── token_usage/
│   ├── (cost calculation)
└── logging/
    ├── (structured logging)
```

## How
The core module implements a **multi-agent system with layered orchestration**:

### Agent Architecture
The agent system follows a **specialized agent pattern**:

1. **Main Agent** (`agents/main.py`)
   - Primary conversational interface
   - Routes tasks to specialized sub-agents
   - Manages tool execution and response synthesis
   - Handles user interaction loop

2. **Research Agent** (`agents/research_agent.py`)
   - Specialized in codebase analysis
   - Deep exploration of code structure
   - Generates research documentation
   - Used for `/research` slash commands

3. **Agent Components** (`agents/agent_components/`)
   - Modular agent capabilities
   - State management components
   - Tool orchestration logic

### State Management
- **Centralized State**: Single source of truth for application state
- **Immutable Updates**: State changes create new state objects
- **Persistence**: State snapshots saved to `.tunacode/sessions/`
- **Recovery**: State restoration after crashes/interruptions

### Prompting System
- **Dynamic Construction**: Prompts built from templates + context
- **Context Injection**: File contents, search results, tool outputs
- **Template Library**: Reusable prompt patterns in `prompts/`
- **Token Budgeting**: Smart context window management

### Background Tasks
- **Async Execution**: Non-blocking tool execution
- **Progress Tracking**: Real-time status updates
- **Cancellation**: Clean shutdown of background jobs
- **Error Isolation**: Background failures don't crash main app

### Token Usage Tracking
- **Real-Time Monitoring**: Token counting during streaming
- **Cost Calculation**: Per-model pricing from `configuration/pricing.py`
- **Session Totals**: Cumulative usage tracking
- **Budget Enforcement**: Optional spending limits

## Why
**Separation of Business Logic**: The core module contains pure business logic:
- **UI Agnostic**: No dependencies on Textual or UI components
- **Testability**: Isolated logic enables comprehensive unit testing
- **Reusability**: Core can be imported by other frontends (API, web UI)
- **Maintainability**: Clear boundaries prevent feature creep

**Multi-Agent Design**:
- **Specialization**: Different agents optimize for different tasks
- **Parallel Execution**: Independent agents can run concurrently
- **Resilience**: Agent failure doesn't crash the system
- **Extensibility**: New agent types can be added without refactoring

**Design Principles**:
1. **Agent Isolation**: Each agent is self-contained
2. **Tool Abstraction**: Agents interact with tools via unified interface
3. **State Immutability**: Predictable state transitions
4. **Streaming First**: All model interactions support streaming responses
5. **Error Recovery**: Graceful degradation when components fail

## Integration Points
- **Tools**: `tunacode.tools` provides callable tool implementations
- **Configuration**: `tunacode.configuration` provides model settings and API keys
- **UI**: `tunacode.ui` consumes agent responses and state updates
- **Services**: `tunacode.services` provides external integrations (Git, LSP)
- **Indexing**: `tunacode.indexing` provides code search capabilities

## Code Quality Notes
- **Type Safety**: Extensive use of `TypedDict` and dataclasses for state
- **Async/Await**: Proper async patterns throughout
- **Error Handling**: Comprehensive exception handling with recovery
- **Documentation**: Inline comments explain complex agent reasoning
- **Testing**: Agent logic designed for deterministic testing

## Key Patterns
- **Dependency Injection**: Components receive dependencies via constructors
- **Strategy Pattern**: Different agent types implement common interface
- **Observer Pattern**: State changes trigger UI updates via events
- **Builder Pattern**: Prompt construction via fluent builders
- **Middleware Pattern**: Token usage logging wraps model calls
