# 检索

OpenViking 提供两种搜索方法：`find` 用于简单语义搜索，`search` 用于带会话上下文的复杂检索。

## find 与 search 对比

| 方面 | find | search |
|------|------|--------|
| 意图分析 | 否 | 是 |
| 会话上下文 | 否 | 是 |
| 查询扩展 | 否 | 是 |
| 默认topk | 10 | 10 |
| 使用场景 | 简单查询 | 对话式搜索 |

## API 参考

### find()

基本的向量相似度搜索。

**签名**

```python
def find(
    self,
    query: str,
    target_uri: str = "",
    limit: int = 10,
    score_threshold: Optional[float] = None,
    filter: Optional[Dict] = None,
) -> FindResult
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| query | str | 是 | - | 搜索查询字符串 |
| target_uri | str | 否 | "" | 限制搜索到特定 URI 前缀 |
| limit | int | 否 | 10 | 最大结果数 |
| score_threshold | float | 否 | None | 最小相关性分数阈值 |
| filter | Dict | 否 | None | 元数据过滤器 |

**返回值**

| 类型 | 说明 |
|------|------|
| FindResult | 包含上下文的搜索结果 |

**FindResult 结构**

```python
class FindResult:
    memories: List[MatchedContext]   # 记忆上下文
    resources: List[MatchedContext]  # 资源上下文
    skills: List[MatchedContext]     # 技能上下文
    query_plan: Optional[QueryPlan]  # 查询计划（仅 search）
    query_results: Optional[List[QueryResult]]  # 详细结果
    total: int                       # 总数（自动计算）
```

**MatchedContext 结构**

```python
class MatchedContext:
    uri: str                         # Viking URI
    context_type: ContextType        # "resource"、"memory" 或 "skill"
    is_leaf: bool                    # 是否为叶子节点
    abstract: str                    # L0 内容
    category: str                    # 分类
    score: float                     # 相关性分数 (0-1)
    match_reason: str                # 匹配原因
    relations: List[RelatedContext]  # 相关上下文
```

**示例：基本搜索**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.find("如何进行用户认证")

for ctx in results.resources:
    print(f"URI: {ctx.uri}")
    print(f"分数: {ctx.score:.3f}")
    print(f"类型: {ctx.context_type}")
    print(f"摘要: {ctx.abstract[:100]}...")
    print("---")

client.close()
```

**示例：指定目标 URI 搜索**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 仅在资源中搜索
results = client.find(
    "认证",
    target_uri="viking://resources/"
)

# 仅在用户记忆中搜索
results = client.find(
    "偏好",
    target_uri="viking://user/memories/"
)

# 仅在技能中搜索
results = client.find(
    "网络搜索",
    target_uri="viking://skills/"
)

# 在特定项目中搜索
results = client.find(
    "API 端点",
    target_uri="viking://resources/my-project/"
)

client.close()
```

---

### search()

带会话上下文和意图分析的搜索。

**签名**

```python
def search(
    self,
    query: str,
    target_uri: str = "",
    session: Optional[Session] = None,
    limit: int = 3,
    score_threshold: Optional[float] = None,
    filter: Optional[Dict] = None,
) -> FindResult
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| query | str | 是 | - | 搜索查询字符串 |
| target_uri | str | 否 | "" | 限制搜索到特定 URI 前缀 |
| session | Session | 否 | None | 用于上下文感知搜索的会话 |
| limit | int | 否 | 3 | 最大结果数 |
| score_threshold | float | 否 | None | 最小相关性分数阈值 |
| filter | Dict | 否 | None | 元数据过滤器 |

**返回值**

| 类型 | 说明 |
|------|------|
| FindResult | 包含查询计划和上下文的搜索结果 |

**示例：会话感知搜索**

```python
import openviking as ov
from openviking.message import TextPart

client = ov.OpenViking(path="./data")
client.initialize()

# 创建带对话上下文的会话
session = client.session()
session.add_message("user", [
    TextPart(text="我正在用 OAuth 构建登录页面")
])
session.add_message("assistant", [
    TextPart(text="我可以帮你实现 OAuth。")
])

# 搜索理解对话上下文
results = client.search(
    "最佳实践",
    session=session
)

for ctx in results.resources:
    print(f"找到: {ctx.uri}")
    print(f"摘要: {ctx.abstract[:200]}...")

client.close()
```

**示例：不带会话的搜索**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# search 也可以不带会话使用
# 它仍然会对查询进行意图分析
results = client.search(
    "如何实现 OAuth 2.0 授权码流程",
)

for ctx in results.resources:
    print(f"找到: {ctx.uri} (分数: {ctx.score:.3f})")

client.close()
```

---

## 检索流程

```
查询 → 意图分析 → 向量搜索 (L0) → Rerank (L1) → 结果
```

1. **意图分析**（仅 search）：理解查询意图，扩展查询
2. **向量搜索**：使用 Embedding 查找候选
3. **Rerank**：使用内容重新打分提高准确性
4. **结果**：返回 top-k 上下文

## 处理结果

### 渐进式读取内容

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.find("认证")

for ctx in results.resources:
    # 从 L0（摘要）开始 - 已在 ctx.abstract 中
    print(f"摘要: {ctx.abstract}")

    if not ctx.is_leaf:
        # 获取 L1（概览）
        overview = client.overview(ctx.uri)
        print(f"概览: {overview[:500]}...")
    else:
        # 加载 L2（内容）
        content = client.read(ctx.uri)
        print(f"文件内容: {content}")

client.close()
```

### 获取相关资源

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.find("OAuth 实现")

for ctx in results.resources:
    print(f"找到: {ctx.uri}")

    # 获取相关资源
    relations = client.relations(ctx.uri)
    for rel in relations:
        print(f"  相关: {rel['uri']} - {rel['reason']}")

client.close()
```

## 最佳实践

### 使用具体查询

```python
# 好 - 具体的查询
results = client.find("OAuth 2.0 授权码流程实现")

# 效果较差 - 太宽泛
results = client.find("认证")
```

### 限定搜索范围

```python
# 在相关范围内搜索以获得更好的结果
results = client.find(
    "错误处理",
    target_uri="viking://resources/my-project/"
)
```

### 对话中使用会话上下文

```python
# 对话式搜索使用会话
from openviking.message import TextPart

session = client.session()
session.add_message("user", [
    TextPart(text="我正在构建登录页面")
])

# 搜索理解上下文
results = client.search("最佳实践", session=session)
```

### 相关文档

- [资源管理](resources.md) - 资源管理
- [会话管理](sessions.md) - 会话上下文
- [上下文层级](../concepts/context-layers.md) - L0/L1/L2
