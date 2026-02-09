# file_read 工具使用说明

## 功能概述

`file_read` 工具现在支持两种模式：
1. **单文件模式**：读取单个文件的内容
2. **多文件模式**：一次性读取多个文件的内容

## 单文件模式

### 基本用法

```json
{
  "path": "data/config.json"
}
```

### 指定行范围

```json
{
  "path": "src/main.py",
  "start_line": 10,
  "end_line": 50
}
```

### 自定义选项

```json
{
  "path": "docs/readme.txt",
  "encoding": "utf-8",
  "show_line_numbers": false
}
```

### 返回格式（带行号）

```json
[
  {
    "line": 1,
    "content": "# Configuration File"
  },
  {
    "line": 2,
    "content": "version: 1.0"
  }
]
```

## 多文件模式

### 基本用法

```json
{
  "path": [
    "src/main.py",
    "src/utils.py",
    "data/config.json"
  ]
}
```

### 返回格式

```json
{
  "total_files": 3,
  "success_count": 3,
  "error_count": 0,
  "files": {
    "src/main.py": {
      "status": "success",
      "content": [
        {"line": 1, "content": "import os"},
        {"line": 2, "content": "import sys"}
      ],
      "total_lines": 100
    },
    "src/utils.py": {
      "status": "success",
      "content": [
        {"line": 1, "content": "def helper():"},
        {"line": 2, "content": "    pass"}
      ],
      "total_lines": 50
    },
    "data/config.json": {
      "status": "success",
      "content": [
        {"line": 1, "content": "{"},
        {"line": 2, "content": "  \"version\": \"1.0\""},
        {"line": 3, "content": "}"}
      ],
      "total_lines": 3
    }
  }
}
```

### 部分文件失败的情况

如果某些文件不存在或无法读取：

```json
{
  "total_files": 3,
  "success_count": 2,
  "error_count": 1,
  "files": {
    "src/main.py": {
      "status": "success",
      "content": [...],
      "total_lines": 100
    },
    "src/missing.py": {
      "status": "error",
      "error": "File not found: src/missing.py"
    },
    "src/utils.py": {
      "status": "success",
      "content": [...],
      "total_lines": 50
    }
  },
  "errors": [
    "File not found: src/missing.py"
  ]
}
```

## 参数说明

| 参数 | 类型 | 必需 | 默认值 | 说明 |
|------|------|------|--------|------|
| `path` 或 `file_path` | string \| array | ✅ | - | 单个文件路径或文件路径数组 |
| `start_line` | integer | ❌ | - | 起始行号（从1开始），多文件模式下应用于所有文件 |
| `end_line` | integer | ❌ | - | 结束行号（包含），多文件模式下应用于所有文件 |
| `encoding` | string | ❌ | auto-detect | 文件编码（如 utf-8、gbk） |
| `show_line_numbers` | boolean | ❌ | true | 是否显示行号（JSON格式） |

## 使用场景

### 场景 1：对比多个配置文件

```json
{
  "path": [
    "config/dev.yaml",
    "config/prod.yaml",
    "config/test.yaml"
  ]
}
```

### 场景 2：读取项目的所有主要文件

```json
{
  "path": [
    "README.md",
    "package.json",
    "src/index.ts",
    "src/config.ts"
  ],
  "show_line_numbers": false
}
```

### 场景 3：批量读取代码文件的特定部分

```json
{
  "path": [
    "src/module1.py",
    "src/module2.py",
    "src/module3.py"
  ],
  "start_line": 1,
  "end_line": 20
}
```

## 注意事项

1. **不要读取二进制文件**：如 PDF、Word、图片等，会返回错误
2. **编码自动检测**：如果不指定 encoding，系统会自动检测文件编码
3. **多文件模式**：即使某些文件读取失败，其他文件仍会正常返回
4. **性能考虑**：一次读取大量大文件可能会比较慢，建议合理控制文件数量
5. **参数兼容性**：支持 `path` 和 `file_path` 两种参数名（不同配置文件可能使用不同的参数名）

## LLM 调用示例

作为 LLM Agent，你可以这样调用：

**读取单个文件：**
```json
{
  "tool_name": "file_read",
  "parameters": {
    "path": "data/results.json"
  }
}
```

**同时读取多个相关文件：**
```json
{
  "tool_name": "file_read",
  "parameters": {
    "path": [
      "experiments/exp1/results.json",
      "experiments/exp2/results.json",
      "experiments/exp3/results.json"
    ]
  }
}
```

这样可以一次性获取所有需要的文件内容，提高效率！

