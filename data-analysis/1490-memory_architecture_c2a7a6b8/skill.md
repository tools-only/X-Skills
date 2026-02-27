# 记忆系统架构

> 最后更新: 2026-02-23 (v2 架构)

## 一、整体架构

```mermaid
flowchart TB
    subgraph storage [存储层]
        SQLITE[(SQLite v2 主存储)]
        FTS[FTS5 全文索引]
        MEMMD[MEMORY.md 精华]
        USERMD[USER.md 用户档案]
    end

    subgraph semantic [语义记忆]
        PREF[用户偏好 PREFERENCE]
        RULE[规则约束 RULE]
        SKILL[成功模式 SKILL]
        ERROR[错误教训 ERROR]
        FACT[重要事实 FACT]
    end

    subgraph episode [情节记忆]
        EP[Episode 历史操作记录]
        SP[Scratchpad 当前任务]
    end

    subgraph session [Session 级别]
        MSG[对话历史 messages]
        TASK[当前任务 current_task]
        VAR[会话变量 variables]
    end

    subgraph daily [每日归纳]
        TURNS[conversation_turns SQLite]
        QUEUE[extraction_queue 重试队列]
    end

    semantic --> SQLITE
    semantic --> FTS
    episode --> SQLITE
    SQLITE -->|每日刷新| MEMMD
    SQLITE -->|每日刷新| USERMD
    MSG -->|v2 LLM 提取| semantic
    MSG -->|会话结束| EP
    TURNS -->|凌晨批量提取| semantic
```

## 二、存储设计 (v2)

| 存储 | 用途 | 技术 | 更新频率 |
|------|------|------|----------|
| **memories 表** | 语义记忆主存储 | SQLite + FTS5 全文索引 | 实时（每轮提取） |
| **episodes 表** | 情节记忆（对话摘要） | SQLite | 会话结束时生成 |
| **scratchpad 表** | 工作记忆（当前任务） | SQLite | 每轮同步更新 |
| **conversation_turns 表** | 对话原文索引 | SQLite | 实时 |
| **extraction_queue 表** | 提取重试队列 | SQLite (UNIQUE去重) | 异步 |
| **MEMORY.md** | 核心记忆精华 | 文件 | 每日归纳 |
| **USER.md** | 用户档案 | 文件 | 每日归纳 |

## 三、数据流 — 渐进式披露

```mermaid
sequenceDiagram
    participant User
    participant Agent
    participant Compiler as PromptCompiler
    participant Builder as PromptBuilder
    participant MM as MemoryManager
    participant LLM

    User->>Agent: 发送消息
    Agent->>Compiler: 编译任务描述
    Compiler-->>Agent: task_description

    Note over Agent: 同步更新 Scratchpad 当前任务
    Agent->>MM: 更新 Scratchpad.current_focus

    Agent->>Builder: 构建 System Prompt
    Builder->>Builder: 注入记忆系统自描述
    Builder->>MM: 读取 Scratchpad (当前任务+近期完成)
    Builder->>MM: 读取 Core Memory (MEMORY.md)
    Builder-->>Agent: 极简 System Prompt

    Agent->>LLM: 请求 (不自动注入动态记忆)

    alt 第一级: LLM 觉得需要历史经验
        LLM->>MM: search_memory(query)
        MM->>MM: RetrievalEngine 多路召回 + 过期过滤
        MM-->>LLM: 语义记忆 + 情节摘要
    end

    alt 第二级: 摘要不够，需要操作细节
        LLM->>MM: search_conversation_traces(keyword)
        MM->>MM: SQLite turns → react_traces → JSONL fallback
        MM-->>LLM: 原始对话记录 (含工具参数/返回值)
    end

    LLM-->>Agent: 响应

    Note over Agent: 后台异步提取
    Agent->>MM: extract_from_turn_v2 (频率控制: ≤5次/session)
```

## 四、记忆类型

### 语义记忆 (SemanticMemory)

```mermaid
classDiagram
    class SemanticMemory {
        +str id
        +MemoryType type
        +MemoryPriority priority
        +str content
        +str subject
        +str predicate
        +float importance_score
        +float confidence
        +datetime expires_at
        +datetime created_at
        +datetime updated_at
    }

    class Episode {
        +str id
        +str session_id
        +str summary
        +str goal
        +str outcome
        +list entities
        +list tools_used
        +datetime started_at
        +datetime ended_at
    }

    class Scratchpad {
        +str current_focus
        +list active_projects
        +datetime updated_at
        +to_markdown()
    }
```

| 类型 | 说明 | 示例 | 默认保留 |
|------|------|------|----------|
| FACT | 事实信息 | "用户的代码目录在 D:\code" | 由 LLM 判断 |
| PREFERENCE | 用户偏好 | "用户喜欢用 Python" | permanent |
| SKILL | 成功模式 | "用 pytest 测试更可靠" | permanent |
| ERROR | 错误教训 | "直接删除文件会导致数据丢失" | 7d |
| RULE | 规则约束 | "禁止虚报执行结果" | 24h (任务规则) / permanent (行为规则) |

### 保留时长机制 (expires_at)

| 优先级 | 保留时长 | 实现方式 |
|--------|----------|----------|
| TRANSIENT | 1 天 | `expires_at = now + 1d`，`end_session` 时清理 |
| SHORT_TERM | 3 天 | `expires_at = now + 3d`，`compute_decay` 低分直接删除 |
| LONG_TERM | 30 天 | `expires_at = now + 30d` |
| PERMANENT | 永不删除 | `expires_at = None` |

v2 提取时 LLM 输出 `duration` 字段 (`permanent|7d|24h|session`)，映射为 `expires_at`。

### 情节记忆 (Episode)

会话结束时由 LLM 生成的交互摘要，包含目标、结果、使用的工具、涉及的实体。
在 `to_markdown()` 中以 "历史操作记录" 标签展示。
不再自动注入 system prompt，仅在 LLM 调用 `search_memory` 时按实体匹配返回。

### 工作记忆 (Scratchpad)

由 compiler `task_description` 驱动的结构化任务状态：
- **当前任务**: 来自 compiler 输出，每轮同步更新
- **近期完成**: 话题切换时自动归档旧任务（带时间戳，最多 5 条）

不再由 LLM 自由生成，`SCRATCHPAD_PROMPT` 已废弃。

## 五、组件关系

```mermaid
flowchart LR
    subgraph core [核心组件]
        Agent[Agent]
        Brain[Brain/LLM]
        Identity[Identity]
    end

    subgraph memory [记忆组件]
        MM[MemoryManager]
        Store[UnifiedStore]
        DB[MemoryStorage SQLite]
        Search[SearchBackend FTS5]
        EX[MemoryExtractor v2]
        RE[RetrievalEngine]
        LC[LifecycleManager]
    end

    subgraph session [会话组件]
        SM[SessionManager]
        S[Session]
    end

    subgraph prompt [提示词组件]
        PA[PromptAssembler]
        PB[PromptBuilder]
    end

    Agent --> MM
    Agent --> Identity
    Agent --> Brain
    Agent --> PA
    PA --> PB
    PB -->|读取 Scratchpad + Core Memory| Store
    MM --> Store
    Store --> DB
    Store --> Search
    MM --> EX
    MM --> RE
    MM --> LC
    Agent --> SM
    SM --> S
    LC -->|每日归纳| Store
```

## 六、提取与注入策略

### 提取策略 (写入端)

| 触发时机 | 方式 | 频率控制 |
|----------|------|----------|
| 用户消息 | `extract_from_turn_v2` (LLM) | ≤5 次/session，≥30 字符 |
| 会话结束 | `generate_episode` (LLM) | 每次会话 1 次 |
| 上下文压缩 | 入队 `extraction_queue` | 去重 (UNIQUE) |
| 凌晨归纳 | `process_unextracted_turns` | 批量处理未提取轮次 |

`extract_quick_facts` (关键词匹配) 已废弃，所有提取由 v2 LLM 完成。

### 注入策略 (读取端) — 渐进式披露

System Prompt 只注入极简信息：
1. **记忆系统自描述** — 告知 LLM 两级搜索机制和使用时机
2. **Scratchpad** — 当前任务 + 近期完成
3. **Core Memory** — MEMORY.md 用户基本信息 + 永久规则

动态记忆 **不再自动注入**，由 LLM 按需两级搜索：

| 级别 | 工具 | 数据源 | 适合场景 |
|------|------|--------|----------|
| 第一级 | `search_memory` | RetrievalEngine 多路召回 → SQLite FTS5 fallback | 偏好/规则/经验/操作摘要 |
| 第二级 | `search_conversation_traces` | SQLite turns → react_traces → JSONL | 操作细节（工具参数、返回值原文） |

检索流程均包含过期记忆过滤 (`expires_at < now()` 自动排除)。

## 七、每日归纳流程

```mermaid
flowchart TB
    START[凌晨触发] --> EXTRACT[批量处理 extraction_queue]
    EXTRACT --> DEDUP[语义去重合并]
    DEDUP --> DECAY[compute_decay 衰减]
    DECAY --> DELETE{importance < 0.1?}
    DELETE -->|是| DEL[直接删除]
    DELETE -->|否 但 < 0.3| DEMOTE[降级为 TRANSIENT]
    DELETE -->|否| KEEP[保留]
    DEL --> EXPIRE[cleanup_expired 清理过期]
    DEMOTE --> EXPIRE
    KEEP --> EXPIRE
    EXPIRE --> ATT[清理过期附件]
    ATT --> REFRESH_MEM[刷新 MEMORY.md]
    REFRESH_MEM --> REFRESH_USER[刷新 USER.md]
    REFRESH_USER --> DONE[完成]
```

## 八、数据目录结构

```
data/memory/
├── memory.db                  # SQLite v2 主存储
│   ├── memories              # 语义记忆 (+ FTS5 索引)
│   ├── episodes              # 情节记忆
│   ├── scratchpad            # 工作记忆
│   ├── conversation_turns    # 对话原文
│   ├── extraction_queue      # 提取重试队列
│   └── attachments           # 文件/媒体记忆
├── memories.json              # v1 兼容（逐步淘汰）
├── conversation_history/      # 原始对话历史 JSONL
│   └── {session_id}.jsonl
└── daily_summaries/           # 每日归纳摘要

identity/
├── MEMORY.md                  # 核心记忆精华
└── USER.md                    # 用户档案（自动生成）
```

## 九、关键文件

| 文件 | 说明 |
|------|------|
| `memory/manager.py` | 记忆管理器，核心协调（提取+保存+检索入口） |
| `memory/storage.py` | SQLite 存储层 (DDL + CRUD) |
| `memory/unified_store.py` | 统一存储封装 (storage + search) |
| `memory/extractor.py` | v2 LLM 提取器 (含 Episode/Scratchpad) |
| `memory/retrieval.py` | 多路召回 + 重排序引擎 |
| `memory/lifecycle.py` | 每日归纳、衰减、清理 |
| `memory/types.py` | SemanticMemory / Episode / Scratchpad 类型定义 |
| `prompt/builder.py` | System Prompt 组装（渐进式披露） |
| `sessions/session.py` | Session 管理（截断摘要优化） |
| `core/agent.py` | Agent 主流程（话题检测 + Scratchpad 同步） |
