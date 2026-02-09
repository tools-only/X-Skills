# MIRA Architecture Overview Skill

**Quick Reference**: MIRA is a multi-tenant conversational AI system that maintains persistent working and long-term memory across interactions. This skill provides instant understanding of system architecture, data flows, and core subsystems.

## What is MIRA?

**Core Concept**: Memory-Integrated Reasoning Assistant - an event-driven conversational AI that learns from every interaction and maintains persistent context across conversations.

**Key Differentiator**: Unlike traditional chatbots that treat each message as isolated, MIRA builds a continuous understanding through:
- Autonomous memory extraction and consolidation
- Proactive context surfacing based on conversation relevance
- Event-driven architecture coordinating all system components
- Complete user data isolation via Row Level Security

## System Architecture (30,000ft View)

### Technology Stack

**Core Framework**:
- Python 3.11+ with FastAPI (async web framework)
- Hypercorn ASGI server with HTTP/2 support
- Pydantic for data validation
- 269 Python files total

**Databases & Caching**:
- PostgreSQL 14+ with pgvector extension (vector search)
- Valkey (Redis-compatible) for caching and sessions
- Single unified database: `mira_service`

**AI & ML**:
- Anthropic Claude Sonnet 4.5 (reasoning model)
- Groq/OpenRouter (execution models for simple tasks)
- mdbr-leaf-ir-asym (768-dim asymmetric embeddings for retrieval)
- BGE reranker (semantic result ranking)
- SpaCy en_core_web_lg (entity extraction)
- LlamaIndex (temporal RAG)
- Letta (domain knowledge blocks)

**Security**:
- HashiCorp Vault (credential storage)
- WebAuthn (passkey authentication)
- Row Level Security (RLS) enforced at database layer
- Magic link authentication via direct SMTP

### Project Structure

```
mira/
├── main.py                  # Application entry point
├── config/                  # Configuration and prompts
│   ├── config.py           # Pydantic config schemas
│   ├── config_manager.py   # Config loader
│   ├── system_prompt.txt   # Base system prompt
│   └── prompts/            # LLM prompt templates
├── cns/                     # Central Nervous System
│   ├── api/                # REST & WebSocket endpoints
│   ├── core/               # Domain models, events, continuum
│   ├── infrastructure/     # Repositories, caching, pools
│   ├── integration/        # Event bus, factory patterns
│   └── services/           # Orchestrator, services
├── lt_memory/              # Long-Term Memory system
│   ├── extraction.py       # LLM-based memory extraction
│   ├── refinement.py       # Memory consolidation
│   ├── linking.py          # Relationship classification
│   ├── proactive.py        # Context-aware surfacing
│   ├── batching.py         # Anthropic Batch API orchestration
│   └── factory.py          # Dependency injection
├── working_memory/         # Working Memory system
│   ├── core.py            # Event-driven orchestrator
│   ├── composer.py        # System prompt composition
│   └── trinkets/          # Context providers
├── tools/                  # Tool framework
│   ├── repo.py            # Tool base class & repository
│   ├── registry.py        # Tool config registry
│   └── implementations/   # Concrete tools
├── auth/                   # Authentication system
│   ├── service.py         # Auth orchestrator
│   ├── session.py         # Session management
│   └── user_credentials.py # Encrypted credentials
├── clients/                # External service clients
│   ├── llm_provider.py    # Anthropic/OpenAI/Groq
│   ├── hybrid_embeddings_provider.py # mdbr-leaf-ir-asym + BGE
│   ├── postgres_client.py # PostgreSQL with RLS
│   ├── valkey_client.py   # Valkey cache
│   └── vault_client.py    # HashiCorp Vault
├── utils/                  # Shared utilities
└── web/                    # Frontend assets
```

## Core Subsystems

### 1. CNS (Central Nervous System)

**Purpose**: Event-driven conversation orchestration following Domain-Driven Design principles.

**Key Components**:
- **Continuum Aggregate** (`cns/core/continuum.py`): Immutable conversation state, one per user (replaces discrete "conversations")
- **Orchestrator** (`cns/services/orchestrator.py`): Coordinates entire message processing flow
- **Event Bus** (`cns/integration/event_bus.py`): Synchronous pub/sub for component coordination
- **Continuum Pool** (`cns/infrastructure/continuum_pool.py`): In-memory pool with Valkey caching
- **Unit of Work Pattern**: Batches database operations (metadata + messages) for transactional consistency

**Domain Events**:
- `MessageReceivedEvent` - User message added
- `TurnCompletedEvent` - Assistant response generated
- `ComposeSystemPromptEvent` - Request system prompt composition
- `SystemPromptComposedEvent` - Structured prompt sections ready
- `UpdateTrinketEvent` - Refresh specific trinket content

**Request → Response Flow**:
```
1. User input → WebSocket/HTTP API endpoint
2. Get/create continuum from pool (Valkey → DB if miss)
3. Add user message to continuum cache
4. Generate fingerprint (retrieval-optimized query expansion)
5. Generate 768-dim embedding of fingerprint
6. Surface relevant memories via hybrid search (BM25 + vector)
7. Log retrieval results (JSONL for quality evaluation)
8. Compose system prompt via event bus + trinkets
9. Get enabled tool definitions from tool repository
10. Call LLM with streaming events API
11. Parse tags from response (<mira:memory_ref>, <mira:my_emotion>)
12. Add assistant message to cache
13. Batch persist via Unit of Work (user msg + assistant msg + metadata)
14. Publish TurnCompletedEvent (triggers auto-unload, Letta buffering, etc.)
15. Return response to client
```

**Segment Collapse System**:
- Segments collapse after context-aware timeout (30min morning, 60min normal, 120min late-night)
- Segment summaries generated via LLM
- Summaries embedded (768-dim) for semantic segment search
- Manifest display shows recent segments (configurable depth)
- Memory extraction triggered per-segment

### 2. Long-Term Memory (LT_Memory)

**Purpose**: Autonomous memory lifecycle from extraction through deletion, following "Earning Your Keep" model.

**Key Components**:
- **ExtractionService** (`lt_memory/extraction.py`): LLM-based memory extraction with deduplication
- **RefinementService** (`lt_memory/refinement.py`): Consolidation and verbose trimming
- **LinkingService** (`lt_memory/linking.py`): Relationship classification (supports, conflicts, supersedes, related)
- **ProactiveService** (`lt_memory/proactive.py`): Context-aware memory surfacing
- **BatchingService** (`lt_memory/batching.py`): Anthropic Batch API orchestration (50% cost savings)
- **EntityGC** (`lt_memory/entity_gc.py`): Dormant entity detection and merge candidate scoring

**Memory Lifecycle**:
```
1. Segment Collapse → Extraction Event
2. Chunk messages (max 100 per chunk)
3. Submit to Anthropic Batch API (async)
4. Poll batch status (scheduled job)
5. Parse extracted memories, validate structure
6. Deduplicate via fuzzy + vector similarity
7. Generate embeddings (768-dim mdbr-leaf-ir-asym)
8. Store in PostgreSQL with RLS
9. Submit relationship classification batch
10. Create memory links (inbound_links, outbound_links JSONB arrays)
11. Link entities (SpaCy NER extraction)
12. Consolidation review (identify duplicate clusters)
13. Refinement (trim verbose memories after 7 days + 3 accesses)
14. Proactive surfacing during turns (vector + graph traversal + scoring)
15. Access tracking updates importance_score
16. Decay/deletion via activity-based scoring
```

**Scoring Formula** (SQL function in `lt_memory/scoring_formula.sql`):
- Base importance score (0.0-1.0)
- Recency decay (activity-day-based, vacation-proof)
- Access pattern boost (logarithmic scaling)
- Hub score (memory with many quality links)
- Result: Dynamic importance for retrieval ranking

**Memory Link Types**:
- `supports`: Memory reinforces/adds evidence to another
- `conflicts`: Direct contradiction requiring resolution
- `supersedes`: Newer information replaces older
- `related`: Thematic connection without explicit relationship

**Entity Linking**:
- SpaCy NER extracts entities (PERSON, ORG, GPE, PRODUCT, etc.)
- Entities stored with 300-dim word vectors (en_core_web_lg)
- Hub topology: memories link to entities, not to each other
- Garbage collection merges dormant duplicates

### 3. Working Memory

**Purpose**: Event-driven system prompt composition via specialized "trinkets".

**Key Components**:
- **WorkingMemory** (`working_memory/core.py`): Event-driven orchestrator
- **SystemPromptComposer** (`working_memory/composer.py`): Assembles sections with cache breakpoints
- **Trinkets** (`working_memory/trinkets/`): Specialized context providers

**Available Trinkets**:
- `TimeManager`: Current date, time, timezone
- `ReminderManager`: Active reminders and deadlines
- `ProactiveMemoryTrinket`: Surfaced memories formatted for LLM
- `ToolGuidanceTrinket`: Tool usage hints and patterns
- `ToolLoaderTrinket`: Available tool hints for invokeother_tool
- `ManifestTrinket`: Recent segment summaries
- `DomainKnowledgeTrinket`: Letta agent memory blocks
- `PunchclockTrinket`: Active time tracking status

**System Prompt Structure**:
```
Block 1 (CACHED): Base prompt + stable trinket content
    └─ cache_control: ephemeral (Anthropic prompt caching)

Block 2 (NON-CACHED): Dynamic trinket content + temporal context
    └─ No cache control (changes every turn)
```

**Event Flow**:
```
1. Orchestrator publishes ComposeSystemPromptEvent
2. WorkingMemory routes UpdateTrinketEvent to each trinket
3. Trinkets publish TrinketContentEvent with content + cache_policy
4. Composer assembles sections (base → cached → non-cached)
5. WorkingMemory publishes SystemPromptComposedEvent
6. Orchestrator receives structured prompt sections
```

### 4. Tool System

**Purpose**: Extensible tool framework with dynamic loading via `invokeother_tool` pattern.

**Architecture**:
- **Tool Base Class** (`tools/repo.py`): Abstract base with user context injection
- **ToolRepository** (`tools/repo.py`): Manages tool lifecycle and definitions
- **InvokeOtherTool** (`tools/implementations/invokeother_tool.py`): Dynamic loader meta-tool

**Tool Isolation**:
- `self.user_id` - Automatic user context via context variable
- `self.db` - User-scoped database access (RLS enforced)
- `self.user_data_path` - Per-user file storage directory

**InvokeOther Pattern**:
```
Essential tools (always loaded):
- web_tool
- reminder_tool
- invokeother_tool

Secondary tools (loaded on demand):
- contacts_tool, punchclock_tool, continuum_tool
- weather_tool, maps_tool, calendar_tool
- email_tool, kasa_tool, pager_tool
- square_tool, customerdatabase_tool

Auto-unload:
- Tools unused for 5 turns automatically unloaded
- ToolLoaderTrinket shows hints in working memory
- LLM can see all available tools without context overhead
```

**Available Tools**:
- `web_tool`: Web search, fetch, and HTTP requests (Playwright JS rendering)
- `reminder_tool`: Natural language reminder management
- `contacts_tool`: Contact storage with geocoding
- `punchclock_tool`: Time tracking for activities/projects
- `continuum_tool`: Semantic search across conversation history
- `invokeother_tool`: Dynamic tool loader

### 5. Authentication & Security

**Purpose**: Multi-tenant security with complete user data isolation.

**Authentication Flow** (Magic Link):
```
1. User requests magic link via email
2. Rate limiting check (5 requests/minute per email)
3. Generate cryptographically secure token
4. Send via direct SMTP (no SendGrid dependency)
5. User clicks link with token
6. Verify token, check expiry, mark as used
7. Create session in Valkey (with TTL)
8. Return session cookie (httponly, secure, samesite=strict)
```

**Session Management** (`auth/session.py`):
- Valkey-backed sessions with TTL
- Idle timeout: 7 days (configurable)
- Max lifetime: 30 days (configurable)
- CSRF protection via session-derived tokens
- API tokens: Long-lived sessions (90 days) with kind="api"

**User Isolation** (Row Level Security):
```sql
-- Example RLS policy (enforced at PostgreSQL level)
CREATE POLICY memories_user_policy ON memories
    FOR ALL TO PUBLIC
    USING (user_id = current_setting('app.current_user_id')::uuid);
```

**Context Variable Pattern**:
```python
# Set once at request start (API endpoints, WebSocket)
from utils.user_context import set_current_user_id
set_current_user_id(user_id)

# PostgreSQL session variable automatically set
# All queries filtered by RLS
# Tools automatically scoped via self.user_id property
```

**Credential Storage**:
- System credentials: HashiCorp Vault (`clients/vault_client.py`)
- User credentials: Encrypted in PostgreSQL (`auth/user_credentials.py`)
- Per-user encryption keys derived from user_id
- Automatic Vault token renewal (prevents expiration)

### 6. Database Schema

**Key Tables**:

**Users & Auth**:
- `users` - User accounts with activity tracking
- `sessions` - Active sessions (also in Valkey)
- `magic_links` - One-time authentication tokens
- `user_credentials` - Encrypted per-user credentials
- `user_activity_days` - Granular activity tracking (vacation-proof scoring)

**Continuum & Messages**:
- `continuums` - One per user, stores metadata (linked_days, etc.)
- `messages` - All messages with role (user/assistant/tool), content, metadata
- `continuum_segments` - Time-bounded segments with summaries and embeddings

**Long-Term Memory**:
- `memories` - Core memory storage with 768-dim embeddings
- `entities` - SpaCy-extracted entities with 300-dim word vectors
- `extraction_batches` - Batch API job tracking
- `post_processing_batches` - Relationship classification batches

**Domain Knowledge**:
- `domain_knowledge_blocks` - Letta agent configurations
- `domain_knowledge_block_content` - Block content with sync tracking

**Vector Indexes**:
```sql
-- Memory similarity search (mdbr-leaf-ir-asym 768-dim)
CREATE INDEX idx_memories_embedding ON memories
    USING ivfflat(embedding vector_cosine_ops);

-- Segment summary search (768-dim)
CREATE INDEX idx_continuum_segments_embedding ON continuum_segments
    USING ivfflat(summary_embedding vector_cosine_ops);

-- Entity similarity (SpaCy 300-dim)
CREATE INDEX idx_entities_embedding ON entities
    USING ivfflat(embedding vector_cosine_ops);
```

## Data Flow Patterns

### Message Processing (Detailed)

```
┌─────────────────────────────────────────────────────────────┐
│ 1. API Layer (WebSocket/HTTP)                               │
│    - Authenticate user (Bearer token or session cookie)     │
│    - Set user context (RLS enforcement)                     │
│    - Validate input (sanitize, rate limit)                  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. Continuum Pool                                           │
│    - get_or_create(user_id)                                 │
│    - Check Valkey cache                                     │
│    - On miss: load from PostgreSQL                          │
│    - Load recent messages into cache                        │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. Orchestrator.process_message()                           │
│    - Add user message to continuum cache                    │
│    - Generate fingerprint (query expansion via Groq)        │
│    - Generate 768-dim embedding of fingerprint              │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. Memory Surfacing                                         │
│    - Vector similarity search (embedding)                   │
│    - Graph traversal (follow memory links)                  │
│    - Scoring (importance + recency + access)                │
│    - Return top N memories                                  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 5. System Prompt Composition (Event-Driven)                 │
│    - Publish ComposeSystemPromptEvent                       │
│    - Trinkets update content                                │
│    - Composer assembles sections                            │
│    - Receive SystemPromptComposedEvent                      │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 6. Tool Selection                                           │
│    - Get enabled tools from repository                      │
│    - Format tool definitions for LLM                        │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 7. LLM Call (Streaming)                                     │
│    - Build structured messages                              │
│    - Stream events API (TextEvent, ThinkingEvent, etc.)     │
│    - Tool execution loop (if tools called)                  │
│    - Circuit breaker (max 10 iterations)                    │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 8. Response Processing                                      │
│    - Extract final text content                             │
│    - Parse tags (<mira:memory_ref>, <mira:my_emotion>)      │
│    - Clean response text                                    │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 9. Persistence (Unit of Work)                               │
│    - Batch user message + assistant message                 │
│    - Update continuum metadata                              │
│    - Commit transaction                                     │
│    - Update Valkey cache                                    │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 10. Event Publishing                                        │
│    - TurnCompletedEvent (triggers downstream)               │
│    - Tool auto-unload check                                 │
│    - Letta message buffering                                │
│    - Segment timeout check                                  │
└─────────────────────────────────────────────────────────────┘
```

### Memory Extraction (Batch Processing)

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Segment Collapse (Timeout Detection)                     │
│    - Scheduled job checks for inactive segments             │
│    - Context-aware timeouts (time of day)                   │
│    - Generate segment summary via LLM                       │
│    - Embed summary (768-dim)                                │
│    - Mark segment as collapsed                              │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 2. Batching Orchestration                                   │
│    - Check min message threshold (20 messages)              │
│    - Load segment messages                                  │
│    - Chunk messages (max 100 per chunk)                     │
│    - Capture memory context (surfaced/referenced IDs)       │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 3. Extraction Payload Generation                            │
│    - Format messages (Human:/Assistant:)                    │
│    - Load memory context texts                              │
│    - Shorten UUIDs (8 chars) for prompt efficiency          │
│    - Build extraction prompt                                │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 4. Batch API Submission                                     │
│    - Create Anthropic Message Batch                         │
│    - Store batch metadata in extraction_batches             │
│    - 24-hour expiry window                                  │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 5. Polling (Scheduled Job)                                  │
│    - Check batch status every 5 minutes                     │
│    - Retrieve completed results                             │
│    - Parse JSON responses                                   │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 6. Memory Processing                                        │
│    - Remap short UUIDs to full UUIDs                        │
│    - Validate memory structure                              │
│    - Deduplicate (fuzzy + vector similarity)                │
│    - Generate embeddings (768-dim)                          │
│    - Store in memories table                                │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────┐
│ 7. Post-Processing Batches                                  │
│    - Relationship classification (batch)                    │
│    - Consolidation review (batch)                           │
│    - Create memory links (JSONB arrays)                     │
│    - Entity extraction (SpaCy)                              │
└─────────────────────────────────────────────────────────────┘
```

## Key Design Patterns

### 1. Event-Driven Architecture

All major state changes publish domain events:
```python
# Event definition
class TurnCompletedEvent(ContinuumEvent):
    turn_number: int
    continuum: Continuum

# Publisher (orchestrator)
self.event_bus.publish(TurnCompletedEvent.create(...))

# Subscriber (tool repository)
self.event_bus.subscribe('TurnCompletedEvent', self._handle_turn_completed)
```

### 2. User Context Injection

User context flows automatically through the entire stack:
```python
# API layer sets context
from utils.user_context import set_current_user_id
set_current_user_id(user_id)

# PostgreSQL session variable set
postgres_client.execute_with_context(user_id, query)

# Tools get automatic scoping
class MyTool(Tool):
    def run(self, **params):
        user_id = self.user_id  # Automatic via context variable
        db = self.db  # RLS enforced automatically
```

### 3. Singleton Factory Pattern

Expensive resources initialized once:
```python
# LT_Memory factory
from lt_memory.factory import get_lt_memory_factory
factory = get_lt_memory_factory()  # Singleton
extraction_service = factory.extraction  # Shared instance

# Embeddings provider
from clients.hybrid_embeddings_provider import get_hybrid_embeddings_provider
embeddings = get_hybrid_embeddings_provider()  # Singleton
```

### 4. Unit of Work Pattern

Batch database operations for transactional consistency:
```python
# Create UoW
uow = continuum_pool.begin_work(continuum)

# Add operations
uow.add_messages(user_msg, assistant_msg)
uow.mark_metadata_updated()

# Commit atomically
uow.commit()
```

### 5. Repository Pattern

Clean separation of domain logic from persistence:
```python
# Domain object (Continuum) has no database dependencies
continuum = Continuum.create_new(user_id)
continuum.add_user_message("Hello")

# Repository handles persistence
continuum_repo.save(continuum)
```

## Configuration

**Main Config** (`config/config.py`):
- `ApiConfig` - LLM provider settings, model selection
- `ApiServerConfig` - FastAPI server settings
- `ToolConfig` - Tool behavior, invokeother settings
- `EmbeddingsConfig` - Fast/deep model configuration
- `SystemConfig` - Global settings, segment timeouts
- `LTMemoryConfig` - Memory extraction, linking, refinement

**Environment Variables**:
- `MIRA_DEV` - Development mode (hot reload)
- All sensitive values stored in HashiCorp Vault

**Prompts** (`config/prompts/`):
- `memory_extraction_system.txt` - Memory extraction guidance
- `memory_extraction_user.txt` - Template with context
- `memory_refinement_system.txt` - Consolidation rules
- `segment_summary_system.txt` - Segment summary generation

## API Endpoints

**Authentication** (`/v0/auth`):
- `POST /v0/auth/request-magic-link` - Request magic link
- `POST /v0/auth/verify-magic-link` - Verify token, create session
- `POST /v0/auth/logout` - Revoke session
- `GET /v0/auth/me` - Get current user profile

**Chat** (`/v0/api`):
- `POST /v0/api/chat` - Simple HTTP request/response
- `WS /v0/ws/chat` - WebSocket streaming (preferred)
- `GET /v0/api/health` - Health check

**Data** (`/v0/api`):
- `GET /v0/api/data/manifest` - Get segment manifest
- `GET /v0/api/data/domain-knowledge-blocks` - List Letta blocks
- `POST /v0/api/data/domain-knowledge-blocks/{id}/toggle` - Enable/disable

**Actions** (`/v0/api`):
- `POST /v0/api/actions/create-reminder` - Quick reminder creation
- Additional action endpoints for common workflows

## Performance Optimizations

### 1. Prompt Caching (Anthropic)
- Base system prompt cached (~90% cost reduction)
- Stable trinket content cached
- Dynamic content uncached

### 2. Model Routing
- Simple tools (reminder, punchclock, weather) → Groq (fast)
- Complex reasoning → Claude Sonnet 4.5
- Emergency fallback → OpenAI/Ollama

### 3. Embedding Efficiency
- Generate once per turn (orchestrator)
- mdbr-leaf-ir-asym for all operations (768-dim asymmetric)
- Query encoding (realtime) for fingerprints
- Document encoding (deep) for memory storage

### 4. Connection Pooling
- PostgreSQL connection pools (max 100)
- Valkey connection reuse
- Thread pool for synchronous endpoints

### 5. Batch Processing
- Anthropic Message Batches API (50% cost savings)
- Background processing for extraction
- Scheduled consolidation jobs

### 6. Continuum Pool
- In-memory LRU cache (100 conversations)
- Valkey backing for persistence
- Lazy loading from PostgreSQL

## Deployment

**Requirements**:
- Python 3.11+
- PostgreSQL 14+ with pgvector
- Valkey/Redis
- HashiCorp Vault

**Deploy Script** (`deploy.sh`):
1. Create virtual environment
2. Install dependencies
3. Download SpaCy models
4. Set up databases
5. Configure Vault
6. Create initial user
7. Optional systemd service

**Running MIRA**:
```bash
# Development mode (hot reload)
MIRA_DEV=true python main.py

# Production mode
python main.py

# With firehose logging (debug LLM calls)
python main.py --firehose
```

## Testing

**Test Framework**: pytest with real components (minimal mocking)

**Test Structure**:
- `tests/api/` - API endpoint tests
- `tests/auth/` - Authentication tests
- `tests/clients/` - External client tests
- `tests/cns/` - CNS component tests
- `tests/lt_memory/` - Memory system tests
- `tests/fixtures/` - Shared test fixtures

**Test Principles** (from `.claude/skills/pytest-real-testing/`):
- Use real components, avoid mocking
- Test user isolation via RLS
- Validate complete workflows
- Clean up test data

## Key Files Reference

**Entry Points**:
- `/Users/taylut/Programming/GitHub/botwithmemory/main.py` - Application bootstrap

**Core Orchestration**:
- `/Users/taylut/Programming/GitHub/botwithmemory/cns/services/orchestrator.py` - Message processing
- `/Users/taylut/Programming/GitHub/botwithmemory/cns/core/continuum.py` - Conversation aggregate
- `/Users/taylut/Programming/GitHub/botwithmemory/cns/integration/event_bus.py` - Event coordination

**Memory Systems**:
- `/Users/taylut/Programming/GitHub/botwithmemory/lt_memory/extraction.py` - Memory extraction
- `/Users/taylut/Programming/GitHub/botwithmemory/lt_memory/proactive.py` - Memory surfacing
- `/Users/taylut/Programming/GitHub/botwithmemory/working_memory/core.py` - Working memory orchestrator

**Database**:
- `/Users/taylut/Programming/GitHub/botwithmemory/deploy/mira_service_schema.sql` - Complete schema

**Configuration**:
- `/Users/taylut/Programming/GitHub/botwithmemory/config/config.py` - All config schemas
- `/Users/taylut/Programming/GitHub/botwithmemory/config/system_prompt.txt` - Base system prompt

## Common Development Tasks

### Adding a New Tool
1. Create tool class inheriting from `Tool` (`tools/repo.py`)
2. Implement `run()` method with type hints
3. Register config in `tools/registry.py`
4. Add to secondary tools directory if not essential
5. Update `ToolLoaderTrinket` hints if needed
6. Write tests in `tests/tools/`

### Adding a New Trinket
1. Create trinket class in `working_memory/trinkets/`
2. Implement `handle_update_request(event)` method
3. Publish `TrinketContentEvent` with content
4. Set appropriate cache_policy (True/False)
5. Register in factory (`cns/integration/factory.py`)
6. Test via event bus subscription

### Modifying Database Schema
1. Update schema in `deploy/mira_service_schema.sql`
2. Create migration SQL file if needed
3. Test migration on fresh database
4. Update Pydantic models if applicable
5. Verify RLS policies still apply

### Adding a New API Endpoint
1. Create endpoint function in appropriate router (`cns/api/`)
2. Use `BaseHandler` pattern for consistency
3. Add authentication dependency (`get_current_user`)
4. Set user context for RLS
5. Return standardized responses (`APIResponse`)
6. Add tests in `tests/api/`

## Debugging Tips

**Enable Firehose Logging**:
```bash
python main.py --firehose
# Logs all LLM API calls to firehose_output.json
```

**Check Event Bus Flow**:
```python
# Add debug logging to specific events
logger.debug(f"Event {event_type} published to {len(subscribers)} subscribers")
```

**Inspect Continuum State**:
```python
# In orchestrator or handlers
logger.info(f"Continuum state: {continuum.to_dict()}")
logger.info(f"Messages in cache: {len(continuum.messages)}")
```

**Verify RLS**:
```sql
-- Check RLS policies
SELECT schemaname, tablename, policyname, qual
FROM pg_policies
WHERE tablename = 'memories';

-- Test RLS enforcement
SET app.current_user_id = 'some-uuid';
SELECT * FROM memories;  -- Should only return that user's memories
```

**Monitor Batch Processing**:
```sql
-- Check extraction batch status
SELECT batch_id, status, created_at, completed_at
FROM extraction_batches
WHERE status != 'completed'
ORDER BY created_at DESC;
```

## Architecture Principles (from CLAUDE.md)

1. **Evidence-Based Decisions**: Form technical positions based on evidence, maintain them despite pushback
2. **Brutal Honesty**: Reject unsound ideas directly, propose superior alternatives
3. **No Backwards Compatibility**: Breaking changes preferred at this stage (greenfield)
4. **Thoughtful Component Design**: Build components that eliminate repetitive work
5. **Integrate Rather Than Invent**: Prefer established patterns over custom solutions
6. **Simple Solutions First**: Never sacrifice correctness for simplicity
7. **UTC Everywhere**: All timestamps in UTC via `utils/timezone_utils.py`
8. **Vault for Credentials**: All sensitive values in HashiCorp Vault
9. **Synchronous Over Async**: Only use async when genuine concurrency benefit exists

## Performance Characteristics

**Typical Message Processing**: 2-5 seconds
- LLM generation: 1-3 seconds
- Memory surfacing: 100-300ms
- System prompt composition: 50-100ms
- Persistence: 50-100ms

**Memory Extraction**: Async (background)
- Batch submission: <1 second
- Anthropic processing: 5-30 minutes
- Per-user: ~20-100 memories extracted per session

**Database Performance**:
- Vector similarity search: 50-200ms (ivfflat index)
- Message retrieval: 10-50ms (cached in pool)
- Segment collapse: 500ms-2s (includes LLM summary)

## Known Limitations & Future Enhancements

**Current Limitations**:
- Single-region deployment (no multi-region support)
- No conversation export/import
- Limited multimedia support (images only)
- Temporal RAG features experimental

**Planned Enhancements** (see `junk_drawer/FUTURE_ENHANCEMENTS.md`):
- Multi-modal memory extraction
- Voice input/output
- Conversation clustering
- Advanced memory consolidation
- Real-time collaboration features

## Developer Notes

**From the Developer** (Taylor Satula):
> "I strive to build software I would personally use. I loath garbo half-implemented buggy slopcode. All of the functionality described works properly and in concert with the other components."

**On Claude Sonnet 4.5**:
> "Claude Sonnet 4.5 has a gestalt quality to it that has enabled an incredibly accurate mimicry of a real person. It speaks to you with such realism and is so good at creative tool use combos. MIRA will work just fine with other models but they'll never have that OEM feel that Sonnet MIRA has."

**On the System Prompt** (`config/system_prompt.txt`):
> "Check out the system_prompt.txt for MIRA in the config/ folder. It is very different than many system prompts from other LLM based projects in its content and its brevity."

---

**Last Updated**: 2025-10-23
**Total Files**: 269 Python files
**Lines of Code**: ~50,000 (estimated)
**Test Coverage**: Comprehensive (API, auth, memory, tools)
