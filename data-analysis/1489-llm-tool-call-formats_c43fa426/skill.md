# LLM 工具调用格式指南

本文档记录了各 LLM 提供者的工具调用格式差异及系统处理方式。

## 概述

不同 LLM 提供者对工具调用（Function Calling / Tool Use）的实现存在差异。虽然大多数遵循 OpenAI 或 Anthropic 的标准格式，但部分模型会返回非标准格式，需要特殊解析处理。

## 当前配置的模型

| 模型 | Provider | API 类型 | 优先级 | 特殊处理 |
|------|----------|---------|--------|---------|
| MiniMax-M2.1 | minimax | anthropic | 0 (最高) | 文本格式解析 |
| qwen3-max | dashscope | openai | 1 | thinking 标签 |
| kimi-k2.5 | moonshot | openai | 2 | 特殊工具格式 |
| claude-opus | yunwu | anthropic | 3 | 标准格式 |

## 各模型详细格式

### 1. Claude (Anthropic 原生)

**API 类型**: `anthropic`

**工具调用格式**: 标准 Anthropic `tool_use` block

```json
{
  "type": "tool_use",
  "id": "toolu_01A09q90qw90lq917835lgs",
  "name": "get_weather",
  "input": {"location": "San Francisco"}
}
```

**无需特殊处理**，系统原生支持。

---

### 2. MiniMax M2.1

**API 类型**: `anthropic` (Anthropic 兼容)  
**基础 URL**: `https://api.minimaxi.com/anthropic`  
**官方文档**: https://platform.minimaxi.com/docs/guides/text-m2-function-call

#### API 选项

MiniMax 提供两种 API：

| API | URL | 特点 |
|-----|-----|------|
| Anthropic 兼容 | `https://api.minimaxi.com/anthropic` | 返回标准 tool_use block |
| OpenAI 兼容 | `https://api.minimaxi.com/v1` | 需设置 `reasoning_split` |

#### Thinking 格式

- **Anthropic API**: thinking 作为独立 block 返回
- **OpenAI API (reasoning_split=True)**: thinking 在 `reasoning_details` 字段
- **OpenAI API (reasoning_split=False)**: thinking 以 `<think>...</think>` 标签包裹在 content 中

#### 特殊文本格式（降级情况）

当 Anthropic 兼容 API 降级时，可能返回文本格式：

```
<minimax:tool_call>
  <invoke name="get_weather">
    <parameter name="location">San Francisco</parameter>
  </invoke>
</minimax:tool_call>
```

**系统处理**：`llm/converters/tools.py` 中的 `parse_text_tool_calls()` 自动解析。

---

### 3. Kimi K2.5 (月之暗面)

**API 类型**: `openai`  
**基础 URL**: `https://api.moonshot.cn/v1`  
**官方文档**: https://platform.moonshot.ai/docs/guide/kimi-k2-5-quickstart

#### 标准格式

Kimi 通过 OpenAI 兼容 API 返回标准 `tool_calls`：

```json
{
  "tool_calls": [{
    "id": "call_xxx",
    "type": "function",
    "function": {
      "name": "get_weather",
      "arguments": "{\"location\": \"Beijing\"}"
    }
  }]
}
```

#### 特殊文本格式（自托管/vLLM 场景）

Kimi K2 原生模型（非 API）可能返回特殊格式：

```
<<|tool_calls_section_begin|>>
<<|tool_call_begin|>>functions.get_weather:0<<|tool_call_argument_begin|>>{"city": "Beijing"}<<|tool_call_end|>>
<<|tool_calls_section_end|>>
```

**格式说明**：
- `functions.{func_name}:{idx}` - 函数 ID，idx 从 0 递增
- 需要手动解析函数名和参数

**系统处理**：`llm/converters/tools.py` 中的 `_parse_kimi_tool_calls()` 自动解析。

#### Thinking 配置

```json
{
  "extra_params": {
    "thinking": {"type": "enabled"}
  }
}
```

---

### 4. DashScope Qwen3 (阿里云通义)

**API 类型**: `openai`  
**基础 URL**: `https://dashscope.aliyuncs.com/compatible-mode/v1`  
**官方文档**: https://help.aliyun.com/zh/model-studio/function-calling

#### 工具调用格式

标准 OpenAI 兼容格式：

```json
{
  "tool_calls": [{
    "id": "call_xxx",
    "type": "function",
    "function": {
      "name": "get_weather",
      "arguments": "{\"location\": \"杭州\"}"
    }
  }]
}
```

#### Thinking 配置

```json
{
  "extra_params": {
    "enable_thinking": true
  }
}
```

Qwen3 的 thinking 可能以 `<think>...</think>` 标签形式出现在 content 中。

---

## 系统处理机制

### 1. 文本格式工具调用解析

位置：`src/openakita/llm/converters/tools.py`

支持的格式：
- `<function_calls>` 通用格式
- `<minimax:tool_call>` MiniMax 格式  
- `<<|tool_calls_section_begin|>>` Kimi K2 格式

```python
def has_text_tool_calls(text: str) -> bool:
    """检测文本中是否包含工具调用"""
    
def parse_text_tool_calls(text: str) -> tuple[str, list[ToolUseBlock]]:
    """解析文本工具调用，返回清理后的文本和工具调用列表"""
```

### 2. Thinking 标签清理

位置：`src/openakita/core/agent.py`

```python
def strip_thinking_tags(text: str) -> str:
    """移除响应中的内部标签"""
```

清理的标签：
- `<thinking>...</thinking>` - Claude extended thinking
- `<think>...</think>` - MiniMax/Qwen thinking
- `<minimax:tool_call>...</minimax:tool_call>` - MiniMax 工具调用
- `<<|tool_calls_section_begin|>>...<<|tool_calls_section_end|>>` - Kimi K2 工具调用

### 3. Provider 层处理

位置：`src/openakita/llm/providers/anthropic.py`

Anthropic provider 在 `_parse_response()` 中：
1. 首先检查标准 `tool_use` block
2. 如果没有，检查文本中是否有工具调用格式
3. 解析文本工具调用并转换为 `ToolUseBlock`

---

## 配置示例

### llm_endpoints.json

```json
{
  "endpoints": [
    {
      "name": "minimax",
      "provider": "minimax",
      "api_type": "anthropic",
      "base_url": "https://api.minimaxi.com/anthropic",
      "api_key_env": "MINIMAX_API_KEY",
      "model": "MiniMax-M2.1",
      "priority": 0
    },
    {
      "name": "kimi",
      "provider": "moonshot",
      "api_type": "openai",
      "base_url": "https://api.moonshot.cn/v1",
      "api_key_env": "KIMI_API_KEY",
      "model": "kimi-k2.5",
      "extra_params": {"thinking": {"type": "enabled"}}
    },
    {
      "name": "dashscope",
      "provider": "dashscope",
      "api_type": "openai",
      "base_url": "https://dashscope.aliyuncs.com/compatible-mode/v1",
      "api_key_env": "DASHSCOPE_API_KEY",
      "model": "qwen3-max-2026-01-23",
      "extra_params": {"enable_thinking": true}
    }
  ]
}
```

---

## 故障排查

### 问题 1: 工具调用以原始 XML 形式显示给用户

**原因**: 文本格式的工具调用未被正确解析

**解决**:
1. 检查 `has_text_tool_calls()` 是否检测到格式
2. 检查 `parse_text_tool_calls()` 解析逻辑
3. 检查 provider 的 `_parse_response()` 是否调用了解析

### 问题 2: Thinking 内容泄露给用户

**原因**: `strip_thinking_tags()` 未清理该格式

**解决**:
1. 添加新的正则表达式匹配模式
2. 确保在发送给用户前调用清理函数

### 问题 3: 模型切换后上下文丢失

**原因**: 不同模型的消息格式不兼容

**解决**:
1. 使用 `task_monitor` 的重试机制
2. 切换时重置上下文到原始用户请求

---

## 参考链接

- [MiniMax 工具使用文档](https://platform.minimaxi.com/docs/guides/text-m2-function-call)
- [Kimi K2 Tool Call Guidance](https://huggingface.co/moonshotai/Kimi-K2-Thinking/blob/main/docs/tool_call_guidance.md)
- [阿里云 DashScope 函数调用](https://help.aliyun.com/zh/model-studio/function-calling)
- [Anthropic Tool Use 文档](https://docs.anthropic.com/en/docs/build-with-claude/tool-use)
- [OpenAI Function Calling](https://platform.openai.com/docs/guides/function-calling)
