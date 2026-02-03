# 技能管理

技能是 Agent 可以调用的能力。本指南介绍如何添加和管理技能。

## API 参考

### add_skill()

添加技能到知识库。

**签名**

```python
def add_skill(
    self,
    data: Any,
    wait: bool = False,
    timeout: float = None,
) -> Dict[str, Any]
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| data | Any | 是 | - | 技能数据（字典、字符串或路径） |
| wait | bool | 否 | False | 是否等待向量化完成 |
| timeout | float | 否 | None | 等待超时时间（秒） |

**支持的数据格式**

1. **字典（技能格式）**：
```python
{
    "name": "skill-name",
    "description": "技能描述",
    "content": "完整的 markdown 内容",
    "allowed_tools": ["Tool1", "Tool2"],  # 可选
    "tags": ["tag1", "tag2"]  # 可选
}
```

2. **字典（MCP 工具格式）** - 自动检测并转换：
```python
{
    "name": "tool_name",
    "description": "工具描述",
    "inputSchema": {
        "type": "object",
        "properties": {...},
        "required": [...]
    }
}
```

3. **字符串（SKILL.md 内容）**：
```python
"""---
name: skill-name
description: 技能描述
---

# 技能内容
"""
```

4. **路径（文件或目录）**：
   - 单个文件：`SKILL.md` 文件路径
   - 目录：包含 `SKILL.md` 的目录路径（包含辅助文件）

**返回值**

| 类型 | 说明 |
|------|------|
| Dict | 包含状态和技能 URI 的结果 |

**返回结构**

```python
{
    "status": "success",
    "uri": "viking://agent/skills/skill-name/",
    "name": "skill-name",
    "auxiliary_files": 0
}
```

**示例：从字典添加技能**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

skill = {
    "name": "search-web",
    "description": "搜索网络获取当前信息",
    "content": """
# search-web

搜索网络获取当前信息。

## 参数
- **query** (string, 必需): 搜索查询
- **limit** (integer, 可选): 最大结果数，默认 10

## 使用场景
当用户需要当前信息时使用。
"""
}

result = client.add_skill(skill)
print(f"已添加: {result['uri']}")

client.close()
```

**示例：从 MCP 工具添加**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# MCP 工具格式会自动检测并转换
mcp_tool = {
    "name": "calculator",
    "description": "执行数学计算",
    "inputSchema": {
        "type": "object",
        "properties": {
            "expression": {
                "type": "string",
                "description": "要计算的数学表达式"
            }
        },
        "required": ["expression"]
    }
}

result = client.add_skill(mcp_tool)
print(f"已添加: {result['uri']}")

client.close()
```

**示例：从 SKILL.md 文件添加**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 从文件路径添加
result = client.add_skill("./skills/search-web/SKILL.md")
print(f"已添加: {result['uri']}")

# 从目录添加（包含辅助文件）
result = client.add_skill("./skills/code-runner/")
print(f"已添加: {result['uri']}")
print(f"辅助文件数: {result['auxiliary_files']}")

client.close()
```

---

## SKILL.md 格式

技能可以使用带 YAML frontmatter 的 SKILL.md 文件定义。

**结构**

```markdown
---
name: skill-name
description: 技能的简要描述
allowed-tools:
  - Tool1
  - Tool2
tags:
  - tag1
  - tag2
---

# 技能名称

完整的 Markdown 格式技能文档。

## 参数
- **param1** (类型, 必需): 描述
- **param2** (类型, 可选): 描述

## 使用场景
何时以及如何使用此技能。

## 示例
技能调用的具体示例。
```

**必需字段**

| 字段 | 类型 | 说明 |
|------|------|------|
| name | str | 技能名称（推荐 kebab-case） |
| description | str | 简要描述 |

**可选字段**

| 字段 | 类型 | 说明 |
|------|------|------|
| allowed-tools | List[str] | 此技能可使用的工具 |
| tags | List[str] | 分类标签 |

---

## 管理技能

### 列出技能

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 列出所有技能
skills = client.ls("viking://agent/skills/")
for skill in skills:
    print(f"{skill['name']}")

# 简单列表（仅名称）
names = client.ls("viking://agent/skills/", simple=True)
print(names)

client.close()
```

### 读取技能内容

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

uri = "viking://agent/skills/search-web/"

# L0: 简要描述
abstract = client.abstract(uri)
print(f"摘要: {abstract}")

# L1: 参数和使用概览
overview = client.overview(uri)
print(f"概览: {overview}")

# L2: 完整技能文档
content = client.read(uri)
print(f"内容: {content}")

client.close()
```

### 搜索技能

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 语义搜索技能
results = client.find(
    "搜索互联网",
    target_uri="viking://agent/skills/",
    limit=5
)

for ctx in results.skills:
    print(f"技能: {ctx.uri}")
    print(f"分数: {ctx.score:.3f}")
    print(f"描述: {ctx.abstract}")
    print("---")

client.close()
```

### 删除技能

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 删除技能
client.rm("viking://agent/skills/old-skill/", recursive=True)

client.close()
```

---

## MCP 转换

OpenViking 自动检测并将 MCP 工具定义转换为技能格式。

**检测**

如果字典包含 `inputSchema` 字段，则视为 MCP 格式：

```python
if "inputSchema" in data:
    # 转换为技能格式
    skill = mcp_to_skill(data)
```

**转换过程**

1. 名称转换为 kebab-case
2. 描述保持不变
3. 从 `inputSchema.properties` 提取参数
4. 从 `inputSchema.required` 标记必需字段
5. 生成 Markdown 内容

**转换示例**

输入（MCP 格式）：
```python
{
    "name": "search_web",
    "description": "搜索网络",
    "inputSchema": {
        "type": "object",
        "properties": {
            "query": {
                "type": "string",
                "description": "搜索查询"
            },
            "limit": {
                "type": "integer",
                "description": "最大结果数"
            }
        },
        "required": ["query"]
    }
}
```

输出（技能格式）：
```python
{
    "name": "search-web",
    "description": "搜索网络",
    "content": """---
name: search-web
description: 搜索网络
---

# search-web

搜索网络

## Parameters

- **query** (string) (required): 搜索查询
- **limit** (integer) (optional): 最大结果数

## Usage

This tool wraps the MCP tool `search-web`. Call this when the user needs functionality matching the description above.
"""
}
```

---

## 技能存储结构

技能存储在 `viking://agent/skills/`：

```
viking://agent/skills/
├── search-web/
│   ├── .abstract.md      # L0: 简要描述
│   ├── .overview.md      # L1: 参数和使用方法
│   ├── SKILL.md          # L2: 完整文档
│   └── [辅助文件]        # 任何附加文件
├── calculator/
│   ├── .abstract.md
│   ├── .overview.md
│   └── SKILL.md
└── ...
```

---

## 最佳实践

### 清晰的描述

```python
# 好 - 具体且可操作
skill = {
    "name": "search-web",
    "description": "使用 Google 搜索网络获取当前信息",
    ...
}

# 不够好 - 太模糊
skill = {
    "name": "search",
    "description": "搜索",
    ...
}
```

### 完整的内容

在技能内容中包含：
- 带类型的清晰参数描述
- 何时使用该技能
- 具体示例
- 边界情况和限制

```python
skill = {
    "name": "search-web",
    "description": "搜索网络获取当前信息",
    "content": """
# search-web

使用 Google 搜索网络获取当前信息。

## 参数
- **query** (string, 必需): 搜索查询。具体的查询效果更好。
- **limit** (integer, 可选): 最大结果数。默认: 10，最大: 100。

## 使用场景
在以下情况使用此技能：
- 用户询问当前事件
- 信息不在知识库中
- 用户明确要求搜索网络

不要在以下情况使用：
- 信息已在资源中可用
- 查询是关于历史事实

## 示例
- "今天天气怎么样？" → search-web(query="今天天气")
- "AI 最新新闻" → search-web(query="AI 新闻 2024", limit=5)

## 限制
- 每小时限制 100 次请求
- 结果可能不包含付费内容
"""
}
```

### 一致的命名

技能名称使用 kebab-case：
- `search-web`（好）
- `searchWeb`（避免）
- `search_web`（避免）

---

## 相关文档

- [上下文类型](../concepts/context-types.md) - 技能概念
- [检索](./05-retrieval.md) - 查找技能
- [会话管理](./04-sessions.md) - 追踪技能使用
