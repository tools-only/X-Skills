# 资源管理

资源是 Agent 可以引用的外部知识。本指南介绍如何添加、管理和检索资源。

## 支持的格式

| 格式 | 扩展名 | 处理方式 |
|------|--------|----------|
| PDF | `.pdf` | 文本和图片提取 |
| Markdown | `.md` | 原生支持 |
| HTML | `.html`, `.htm` | 清洗后文本提取 |
| 纯文本 | `.txt` | 直接导入 |
| JSON/YAML | `.json`, `.yaml`, `.yml` | 结构化解析 |
| 代码 | `.py`, `.js`, `.ts`, `.go`, `.java` 等 | 语法感知解析 |
| 图片 | `.png`, `.jpg`, `.jpeg`, `.gif`, `.webp` | VLM 描述 |
| 视频 | `.mp4`, `.mov`, `.avi` | 帧提取 + VLM |
| 音频 | `.mp3`, `.wav`, `.m4a` | 转录 |
| 文档 | `.docx` | 文本提取 |

## 处理流程

```
输入 → Parser → TreeBuilder → AGFS → SemanticQueue → 向量索引
```

1. **Parser**：根据文件类型提取内容
2. **TreeBuilder**：创建目录结构
3. **AGFS**：存储文件到虚拟文件系统
4. **SemanticQueue**：异步生成 L0/L1
5. **向量索引**：建立语义搜索索引

## API 参考

### add_resource()

添加资源。

**签名**

```python
def add_resource(
    self,
    path: str,
    target: Optional[str] = None,
    reason: str = "",
    instruction: str = "",
    wait: bool = False,
    timeout: float = None,
) -> Dict[str, Any]
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| path | str | 是 | - | 本地文件路径、目录路径或 URL |
| target | str | 否 | None | 目标 Viking URI（必须在 `resources` 作用域） |
| reason | str | 否 | "" | 添加此资源的原因（提升搜索相关性） |
| instruction | str | 否 | "" | 特殊处理指令 |
| wait | bool | 否 | False | 是否等待语义处理完成 |

**返回值**

| 类型 | 说明 |
|------|------|
| Dict[str, Any] | 包含状态和资源信息的结果 |

**返回结构**

```python
{
    "status": "success",           # "success" 或 "error"
    "root_uri": "viking://resources/docs/",  # 根资源 URI
    "source_path": "./docs/",      # 原始源路径
    "errors": [],                  # 错误列表（如有）
    "queue_status": {...}          # 队列状态（仅在 wait=True 时）
}
```

**示例：添加单个文件**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

result = client.add_resource(
    "./documents/guide.md",
    reason="用户指南文档"
)
print(f"已添加到: {result['root_uri']}")

client.wait_processed()
client.close()
```

**示例：从 URL 添加**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

result = client.add_resource(
    "https://example.com/api-docs.md",
    target="viking://resources/external/",
    reason="外部 API 文档"
)

# 等待处理
client.wait_processed()
client.close()
```

**示例：等待处理**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 方式 1：内联等待
result = client.add_resource(
    "./documents/guide.md",
    wait=True
)
print(f"队列状态: {result['queue_status']}")

# 方式 2：单独等待（用于批量处理）
client.add_resource("./file1.md")
client.add_resource("./file2.md")
client.add_resource("./file3.md")

status = client.wait_processed()
print(f"全部处理完成: {status}")

client.close()
```

---

### export_ovpack()

将资源树导出为 `.ovpack` 文件。

**签名**

```python
def export_ovpack(self, uri: str, to: str) -> str
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| uri | str | 是 | - | 要导出的 Viking URI |
| to | str | 是 | - | 目标文件路径 |

**返回值**

| 类型 | 说明 |
|------|------|
| str | 导出文件的路径 |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 导出项目
path = client.export_ovpack(
    "viking://resources/my-project/",
    "./exports/my-project.ovpack"
)
print(f"已导出到: {path}")

client.close()
```

---

### import_ovpack()

导入 `.ovpack` 文件。

**签名**

```python
def import_ovpack(
    self,
    file_path: str,
    parent: str,
    force: bool = False,
    vectorize: bool = True
) -> str
```

**参数**

| 参数 | 类型 | 必填 | 默认值 | 说明 |
|------|------|------|--------|------|
| file_path | str | 是 | - | 本地 `.ovpack` 文件路径 |
| parent | str | 是 | - | 目标父 URI |
| force | bool | 否 | False | 是否覆盖已存在的资源 |
| vectorize | bool | 否 | True | 导入后是否触发向量化 |

**返回值**

| 类型 | 说明 |
|------|------|
| str | 导入资源的根 URI |

**示例**

```python
import openviking as ov

client = ov.OpenViking(path="./data")
client.initialize()

# 导入包
uri = client.import_ovpack(
    "./exports/my-project.ovpack",
    "viking://resources/imported/",
    force=True,
    vectorize=True
)
print(f"已导入到: {uri}")

client.wait_processed()
client.close()
```

---

## 管理资源

### 列出资源

```python
# 列出所有资源
entries = client.ls("viking://resources/")

# 列出详细信息
for entry in entries:
    type_str = "目录" if entry['isDir'] else "文件"
    print(f"{entry['name']} - {type_str}")

# 简单路径列表
paths = client.ls("viking://resources/", simple=True)
# 返回: ["project-a/", "project-b/", "shared/"]

# 递归列出
all_entries = client.ls("viking://resources/", recursive=True)
```

### 读取资源内容

```python
# L0: 摘要
abstract = client.abstract("viking://resources/docs/")

# L1: 概览
overview = client.overview("viking://resources/docs/")

# L2: 完整内容
content = client.read("viking://resources/docs/api.md")
```

### 移动资源

```python
client.mv(
    "viking://resources/old-project/",
    "viking://resources/new-project/"
)
```

### 删除资源

```python
# 删除单个文件
client.rm("viking://resources/docs/old.md")

# 递归删除目录
client.rm("viking://resources/old-project/", recursive=True)
```

### 创建链接

```python
# 链接相关资源
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
```

### 获取关联

```python
relations = client.relations("viking://resources/docs/auth/")
for rel in relations:
    print(f"{rel['uri']}: {rel['reason']}")
```

### 删除链接

```python
client.unlink(
    "viking://resources/docs/auth/",
    "viking://resources/docs/security/"
)
```

## 最佳实践

### 按项目组织

```
viking://resources/
├── project-a/
│   ├── docs/
│   ├── specs/
│   └── references/
├── project-b/
│   └── ...
└── shared/
    └── common-docs/
```

## 相关文档

- [检索](retrieval.md) - 搜索资源
- [文件系统](filesystem.md) - 文件操作
- [上下文类型](../concepts/context-types.md) - 资源概念
