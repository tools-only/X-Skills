# 会话管理

会话管理对话状态、追踪上下文使用，并提取长期记忆。

## API 参考

### client.session()

创建新会话或加载已有会话。

**签名**

```python
def session(self, session_id: Optional[str] = None) -> Session
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| session_id | str | 否 | None | 会话 ID。如果为 None 则创建新会话（自动生成 ID） |

**返回值**

| 类型 | 说明 |
|------|------|
| Session | Session 对象 |

**示例：创建新会话**

```python
import openviking as ov

client = ov.OpenViking(path="./data", user="alice")
client.initialize()

# 创建新会话（自动生成 ID）
session = client.session()
print(f"会话 URI: {session.uri}")

client.close()
```

**示例：加载已有会话**

```python
import openviking as ov

client = ov.OpenViking(path="./data", user="alice")
client.initialize()

# 加载已有会话
session = client.session(session_id="abc123")
session.load()
print(f"已加载 {len(session.messages)} 条消息")

client.close()
```

---

### Session.add_message()

向会话添加消息。

**签名**

```python
def add_message(
    self,
    role: str,
    parts: List[Part],
) -> Message
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| role | str | 是 | - | 消息角色："user" 或 "assistant" |
| parts | List[Part] | 是 | - | 消息部分列表（TextPart、ContextPart、ToolPart） |

**返回值**

| 类型 | 说明 |
|------|------|
| Message | 创建的消息对象 |

**Part 类型**

```python
from openviking.message import TextPart, ContextPart, ToolPart

# 文本内容
TextPart(text="你好，有什么可以帮助你的？")

# 上下文引用
ContextPart(
    uri="viking://resources/docs/auth/",
    context_type="resource",  # "resource"、"memory" 或 "skill"
    abstract="认证指南..."
)

# 工具调用
ToolPart(
    tool_id="call_123",
    tool_name="search_web",
    skill_uri="viking://skills/search-web/",
    tool_input={"query": "OAuth 最佳实践"},
    tool_output="",
    tool_status="pending"  # "pending"、"running"、"completed"、"error"
)
```

**示例：文本消息**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# 添加用户消息
session.add_message("user", [
    TextPart(text="如何进行用户认证？")
])

# 添加助手响应
session.add_message("assistant", [
    TextPart(text="你可以使用 OAuth 2.0 进行认证...")
])

client.close()
```

**示例：带上下文引用**

```python
import openviking as ov
from openviking.message import TextPart, ContextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

session.add_message("assistant", [
    TextPart(text="根据文档..."),
    ContextPart(
        uri="viking://resources/docs/auth/",
        context_type="resource",
        abstract="涵盖 OAuth 2.0 的认证指南..."
    )
])

client.close()
```

**示例：带工具调用**

```python
import openviking as ov
from openviking.message import TextPart, ToolPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# 添加带工具调用的消息
msg = session.add_message("assistant", [
    TextPart(text="让我搜索一下..."),
    ToolPart(
        tool_id="call_123",
        tool_name="search_web",
        skill_uri="viking://skills/search-web/",
        tool_input={"query": "OAuth 最佳实践"},
        tool_status="pending"
    )
])

# 稍后更新工具结果
session.update_tool_part(
    message_id=msg.id,
    tool_id="call_123",
    output="找到 5 篇相关文章...",
    status="completed"
)

client.close()
```

---

### Session.used()

追踪对话中实际使用的上下文和技能。

**签名**

```python
def used(
    self,
    contexts: Optional[List[str]] = None,
    skill: Optional[Dict[str, Any]] = None,
) -> None
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| contexts | List[str] | 否 | None | 使用的上下文 URI 列表 |
| skill | Dict | 否 | None | 技能使用信息，包含 uri、input、output、success |

**Skill 字典结构**

```python
{
    "uri": "viking://skills/search-web/",
    "input": "搜索查询",
    "output": "搜索结果...",
    "success": True  # 默认 True
}
```

**示例**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# 搜索相关上下文
results = client.find("认证")

# 在响应中使用上下文
session.add_message("assistant", [
    TextPart(text="根据文档...")
])

# 追踪实际有帮助的上下文
session.used(contexts=[
    "viking://resources/认证文档/"
])

# 追踪技能使用
session.used(skill={
    "uri": "viking://skills/code-search/",
    "input": "搜索认证示例",
    "output": "找到 3 个示例文件",
    "success": True
})

session.commit()

client.close()
```

---

### Session.update_tool_part()

更新工具调用的输出和状态。

**签名**

```python
def update_tool_part(
    self,
    message_id: str,
    tool_id: str,
    output: str,
    status: str = "completed",
) -> None
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| message_id | str | 是 | - | 包含工具调用的消息 ID |
| tool_id | str | 是 | - | 要更新的工具调用 ID |
| output | str | 是 | - | 工具执行输出 |
| status | str | 否 | "completed" | 工具状态："completed" 或 "error" |

**示例**

```python
import openviking as ov
from openviking.message import ToolPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# 添加工具调用
msg = session.add_message("assistant", [
    ToolPart(
        tool_id="call_456",
        tool_name="execute_code",
        skill_uri="viking://skills/code-runner/",
        tool_input={"code": "print('hello')"},
        tool_status="pending"
    )
])

# 执行工具并更新结果
session.update_tool_part(
    message_id=msg.id,
    tool_id="call_456",
    output="hello",
    status="completed"
)

client.close()
```

---

### Session.commit()

提交会话，归档消息并提取长期记忆。

**签名**

```python
def commit(self) -> Dict[str, Any]
```

**返回值**

| 类型 | 说明 |
|------|------|
| Dict | 包含状态和统计信息的提交结果 |

**返回结构**

```python
{
    "session_id": "abc123",
    "status": "committed",
    "memories_extracted": 3,
    "active_count_updated": 5,
    "archived": True,
    "stats": {
        "total_turns": 10,
        "contexts_used": 4,
        "skills_used": 2,
        "memories_extracted": 3
    }
}
```

**提交时发生什么**

1. **归档**：当前消息归档到 `history/archive_N/`
2. **记忆提取**：使用 LLM 提取长期记忆
3. **去重**：新记忆与现有记忆去重
4. **关联**：在记忆和使用的上下文之间创建链接
5. **统计**：更新使用统计

**记忆分类**

| 分类 | 位置 | 说明 |
|------|------|------|
| profile | `user/memories/.overview.md` | 用户档案信息 |
| preferences | `user/memories/preferences/` | 按主题的用户偏好 |
| entities | `user/memories/entities/` | 重要实体（人物、项目） |
| events | `user/memories/events/` | 重要事件 |
| cases | `agent/memories/cases/` | 问题-解决方案案例 |
| patterns | `agent/memories/patterns/` | 交互模式 |

**示例**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# 添加对话
session.add_message("user", [
    TextPart(text="我喜欢深色模式和 vim 快捷键")
])
session.add_message("assistant", [
    TextPart(text="我已记录你对深色模式和 vim 快捷键的偏好。")
])

# 提交会话
result = session.commit()
print(f"状态: {result['status']}")
print(f"提取的记忆数: {result['memories_extracted']}")
print(f"统计: {result['stats']}")

client.close()
```

---

### Session.load()

从存储加载会话数据。

**签名**

```python
def load(self) -> None
```

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 加载已有会话
session = client.session(session_id="existing-session-id")
session.load()

print(f"已加载 {len(session.messages)} 条消息")
for msg in session.messages:
    print(f"  [{msg.role}]: {msg.parts[0].text[:50]}...")

client.close()
```

---

### Session.get_context_for_search()

获取用于搜索查询扩展的会话上下文。

**签名**

```python
def get_context_for_search(
    self,
    query: str,
    max_archives: int = 3,
    max_messages: int = 20
) -> Dict[str, Any]
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| query | str | 是 | - | 用于匹配相关归档的查询 |
| max_archives | int | 否 | 3 | 最多获取的归档数量 |
| max_messages | int | 否 | 20 | 最多获取的最近消息数量 |

**返回值**

| 类型 | 说明 |
|------|------|
| Dict | 包含摘要和最近消息的上下文 |

**返回结构**

```python
{
    "summaries": ["归档1的概览...", "归档2的概览...", ...],
    "recent_messages": [Message, Message, ...]
}
```

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session(session_id="existing-session")
session.load()

context = session.get_context_for_search(
    query="认证",
    max_archives=3,
    max_messages=10
)

print(f"摘要数: {len(context['summaries'])}")
print(f"最近消息数: {len(context['recent_messages'])}")

client.close()
```

---

## Session 属性

| 属性 | 类型 | 说明 |
|------|------|------|
| uri | str | 会话 Viking URI（`viking://session/{session_id}/`） |
| messages | List[Message] | 会话中的当前消息 |
| stats | SessionStats | 会话统计 |
| summary | str | 压缩摘要 |
| usage_records | List[Usage] | 上下文和技能使用记录 |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

session = client.session()

# 访问属性
print(f"URI: {session.uri}")
print(f"消息数: {len(session.messages)}")
print(f"统计: {session.stats}")

client.close()
```

---

## 会话存储结构

```
viking://session/{session_id}/
├── .abstract.md              # L0: 会话概述
├── .overview.md              # L1: 关键决策
├── messages.jsonl            # 当前消息
├── tools/                    # 工具执行
│   └── {tool_id}/
│       └── tool.json
├── .meta.json                # 元数据
├── .relations.json           # 相关上下文
└── history/                  # 归档历史
    ├── archive_001/
    │   ├── messages.jsonl
    │   ├── .abstract.md
    │   └── .overview.md
    └── archive_002/
```

---

## 整体示例

```python
import openviking as ov
from openviking.message import TextPart, ContextPart, ToolPart

# 初始化客户端
client = ov.OpenViking(path="./my_data")
client.initialize()

# 创建新会话
session = client.session()

# 添加用户消息
session.add_message("user", [
    TextPart(text="如何配置 embedding？")
])

# 带会话上下文搜索
results = client.search("embedding 配置", session=session)

# 添加助手响应，带上下文引用
session.add_message("assistant", [
    TextPart(text="根据文档，你可以这样配置 embedding..."),
    ContextPart(
        uri=results.resources[0].uri,
        context_type="resource",
        abstract=results.resources[0].abstract
    )
])

# 追踪实际使用的上下文
session.used(contexts=[results.resources[0].uri])

# 提交会话（归档消息、提取记忆）
result = session.commit()
print(f"提取的记忆数: {result['memories_extracted']}")

client.close()
```

## 最佳实践

### 定期提交

```python
# 在重要交互后提交
if len(session.messages) > 10:
    session.commit()
```

### 追踪实际使用的内容

```python
# 只标记实际有帮助的上下文
if context_was_useful:
    session.used(contexts=[ctx.uri])
```

### 使用会话上下文搜索

```python
# 带对话上下文的搜索结果更好
results = client.search(query, session=session)
```

### 继续前先加载

```python
# 恢复已有会话时始终先加载
session = client.session(session_id="existing-id")
session.load()
```

---

## 相关文档

- [上下文类型](../concepts/context-types.md) - 记忆类型
- [检索](./05-retrieval.md) - 带会话搜索
- [客户端](./01-client.md) - 创建会话
