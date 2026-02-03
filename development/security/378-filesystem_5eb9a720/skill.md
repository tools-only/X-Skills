# 文件系统

OpenViking 提供类 Unix 的文件系统操作来管理上下文。

## API 参考

### abstract()

读取 L0 摘要（~100 tokens 摘要）。

**签名**

```python
def abstract(self, uri: str) -> str
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | Viking URI（必须是目录） |

**返回值**

| 类型 | 说明 |
|------|------|
| str | L0 摘要内容（.abstract.md） |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

abstract = client.abstract("viking://resources/docs/")
print(f"摘要: {abstract}")
# 输出: "项目 API 文档，涵盖认证、端点..."

client.close()
```

---

### overview()

读取 L1 概览，针对目录生效。

**签名**

```python
def overview(self, uri: str) -> str
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | Viking URI（必须是目录） |

**返回值**

| 类型 | 说明 |
|------|------|
| str | L1 概览内容（.overview.md） |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

overview = client.overview("viking://resources/docs/")
print(f"概览:\n{overview}")

client.close()
```

---

### read()

读取 L2 完整内容。

**签名**

```python
def read(self, uri: str) -> str
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | Viking URI |

**返回值**

| 类型 | 说明 |
|------|------|
| str | 完整文件内容 |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

content = client.read("viking://resources/docs/api.md")
print(f"内容:\n{content}")

client.close()
```

---

### ls()

列出目录内容。

**签名**

```python
def ls(self, uri: str, **kwargs) -> List[Any]
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | Viking URI |
| simple | bool | 否 | False | 仅返回相对路径 |
| recursive | bool | 否 | False | 递归列出所有子目录 |

**返回值**

| 类型 | 说明 |
|------|------|
| List[Dict] | 条目列表（simple=False 时） |
| List[str] | 路径列表（simple=True 时） |

**条目结构**

```python
{
    "name": "docs",           # 文件/目录名
    "size": 4096,             # 字节大小
    "mode": 16877,            # 文件模式
    "modTime": "2024-01-01T00:00:00Z",  # ISO 时间戳
    "isDir": True,            # 是否为目录
    "uri": "viking://resources/docs/",  # Viking URI
    "meta": {}                # 可选元数据
}
```

**示例：基本列出**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

entries = client.ls("viking://resources/")
for entry in entries:
    type_str = "目录" if entry['isDir'] else "文件"
    print(f"{entry['name']} - {type_str}")

client.close()
```

---

### tree()

获取目录树结构。

**签名**

```python
def tree(self, uri: str) -> List[Dict]
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | Viking URI |

**返回值**

| 类型 | 说明 |
|------|------|
| List[Dict] | 包含 rel_path 的扁平条目列表 |

**条目结构**

```python
[
    {
        "name": "docs",
        "size": 4096,
        "mode": 16877,
        "modTime": "2024-01-01T00:00:00Z",
        "isDir": True,
        "rel_path": "docs/",      # 相对于基础 URI 的路径
        "uri": "viking://resources/docs/"
    },
    ...
]
```

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

entries = client.tree("viking://resources/")
for entry in entries:
    type_str = "目录" if entry['isDir'] else "文件"
    print(f"{entry['rel_path']} - {type_str}")

client.close()
```

---

### rm()

删除文件或目录。

**签名**

```python
def rm(self, uri: str, recursive: bool = False) -> None
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | 要删除的 Viking URI |
| recursive | bool | 否 | False | 递归删除目录 |

**返回值**

| 类型 | 说明 |
|------|------|
| None | - |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 删除单个文件
client.rm("viking://resources/docs/old.md")

# 递归删除目录
client.rm("viking://resources/old-project/", recursive=True)

client.close()
```

---

### mv()

移动文件或目录。

**签名**

```python
def mv(self, from_uri: str, to_uri: str) -> None
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| from_uri | str | 是 | - | 源 Viking URI |
| to_uri | str | 是 | - | 目标 Viking URI |

**返回值**

| 类型 | 说明 |
|------|------|
| None | - |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

client.mv(
    "viking://resources/old-name/",
    "viking://resources/new-name/"
)

client.close()
```

---

### grep()

按模式搜索内容。

**签名**

```python
def grep(
    self,
    uri: str,
    pattern: str,
    case_insensitive: bool = False
) -> Dict
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | 要搜索的 Viking URI |
| pattern | str | 是 | - | 搜索模式（正则表达式） |
| case_insensitive | bool | 否 | False | 忽略大小写 |

**返回值**

| 类型 | 说明 |
|------|------|
| Dict | 包含匹配项的搜索结果 |

**返回结构**

```python
{
    "matches": [
        {
            "uri": "viking://resources/docs/auth.md",
            "line": 15,
            "content": "用户认证由..."
        }
    ],
    "count": 1
}
```

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

results = client.grep(
    "viking://resources/",
    "认证",
    case_insensitive=True
)

print(f"找到 {results['count']} 个匹配")
for match in results['matches']:
    print(f"  {match['uri']}:{match['line']}")
    print(f"    {match['content']}")

client.close()
```

---

### glob()

按模式匹配文件。

**签名**

```python
def glob(self, pattern: str, uri: str = "viking://") -> Dict
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| pattern | str | 是 | - | Glob 模式（如 `**/*.md`） |
| uri | str | 否 | "viking://" | 起始 URI |

**返回值**

| 类型 | 说明 |
|------|------|
| Dict | 匹配的 URI |

**返回结构**

```python
{
    "matches": [
        "viking://resources/docs/api.md",
        "viking://resources/docs/guide.md"
    ],
    "count": 2
}
```

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 查找所有 markdown 文件
results = client.glob("**/*.md", "viking://resources/")
print(f"找到 {results['count']} 个 markdown 文件:")
for uri in results['matches']:
    print(f"  {uri}")

# 查找所有 Python 文件
results = client.glob("**/*.py", "viking://resources/")
print(f"找到 {results['count']} 个 Python 文件")

client.close()
```

---

### link()

创建资源之间的关联。

**签名**

```python
def link(
    self,
    from_uri: str,
    uris: Any,
    reason: str = ""
) -> None
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| from_uri | str | 是 | - | 源 URI |
| uris | str 或 List[str] | 是 | - | 目标 URI |
| reason | str | 否 | "" | 链接原因 |

**返回值**

| 类型 | 说明 |
|------|------|
| None | - |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 单个链接
client.link(
    "viking://resources/docs/auth/",
    "viking://resources/docs/security/",
    reason="认证的安全最佳实践"
)

# 多个链接
client.link(
    "viking://resources/docs/api/",
    [
        "viking://resources/docs/auth/",
        "viking://resources/docs/errors/"
    ],
    reason="相关文档"
)

client.close()
```

---

### relations()

获取资源的关联。

**签名**

```python
def relations(self, uri: str) -> List[Dict[str, Any]]
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | Viking URI |

**返回值**

| 类型 | 说明 |
|------|------|
| List[Dict] | 相关资源列表 |

**返回结构**

```python
[
    {"uri": "viking://resources/docs/security/", "reason": "安全最佳实践"},
    {"uri": "viking://resources/docs/errors/", "reason": "错误处理"}
]
```

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

relations = client.relations("viking://resources/docs/auth/")
for rel in relations:
    print(f"相关: {rel['uri']}")
    print(f"  原因: {rel['reason']}")

client.close()
```

---

### unlink()

删除关联。

**签名**

```python
def unlink(self, from_uri: str, uri: str) -> None
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| from_uri | str | 是 | - | 源 URI |
| uri | str | 是 | - | 要取消链接的目标 URI |

**返回值**

| 类型 | 说明 |
|------|------|
| None | - |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

client.unlink(
    "viking://resources/docs/auth/",
    "viking://resources/docs/security/"
)

client.close()
```

---

## 相关文档

- [Viking URI](../concepts/viking-uri.md) - URI 规范
- [上下文层级](../concepts/context-layers.md) - L0/L1/L2
- [资源管理](resources.md) - 资源管理
