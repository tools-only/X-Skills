# CLAUDE.md v5.0

> **核心原则**: 计划 → 确认 → 执行到底 → 验收
> **智能加载**: 只加载必需的文档，保持 context 清洁

---

## 🎯 工作模式（4步）

```
1️⃣ 收到任务 → TodoList 规划
2️⃣ 展示计划 → 用户确认
3️⃣ 执行到底（不问问题）→ 自行决策
4️⃣ 总结验收 → 交付成果
```

**详细说明**: `rules/core/work-mode.md`

---

## 🚨 唯一允许提问的4种情况

1. ❗ **缺少关键凭证** - 数据库密码、API key
2. ❗ **多个对立方案** - 无法从代码库判断
3. ❗ **需求本质矛盾** - 用户要求冲突
4. ❗ **不可逆高风险** - 删除生产数据、强制推送

**其他情况自行决策**: 文件命名、代码风格、依赖版本、UI细节等

**详细说明**: `rules/core/blocking-rules.md`

---

## 🤖 智能 Context 加载

### 系统自动识别任务类型，按需加载文档

**你不需要关心加载什么**，系统会根据你的需求自动选择。

| 任务关键词 | 自动加载的文档 | 预估大小 |
|-----------|---------------|---------|
| 浏览器、自动化、爬虫 | 浏览器自动化指南 | ~15KB |
| 视频、Remotion、动画 | 视频创作指南 | ~25KB |
| 数据、分析、SQL | 数据分析指南 | ~20KB |
| 设计、UI、界面 | 设计指南 | ~30KB |
| 营销、文案、SEO | 营销指南 | ~35KB |
| 开发、代码、功能 | 编码规则 | ~15KB |
| 错误、bug、调试 | 错误目录 | ~12KB |
| 安全、漏洞、审计 | 安全规则 | ~15KB |

**系统说明**: `CONTEXT_MANAGER.md`

---

## 🚀 快速开始（我想要...）

### 不知道用什么工具？

👉 **查看快速决策树**: `index/task-router.md`

- 30秒找到你需要的工具
- 按任务类型分类
- 包含所有能力的快速链接

### 常用能力快速跳转

| 能力 | 快速链接 |
|------|---------|
| 🌐 浏览器自动化 | `capabilities/browser-automation/decision-tree.md` |
| 🎬 视频创作 | `rules/remotion-auto-production.md` |
| 📊 数据分析 | `capabilities/data-analysis-guide.md` |
| 🎨 UI 设计 | `design/DESIGN_MASTER_PERSONA.md` |
| 📝 营销内容 | `vibe-marketing/VIBE_MARKETING_GUIDE.md` |
| 🐛 错误调试 | `errors/top-5-errors.md` |
| 🤖 GPT 专家 | `rules/delegator/triggers.md` |

---

## ⚠️ Top 5 高频错误（快速参考）

| 错误 | 核心问题 | 快速检查 |
|------|---------|---------|
| **E001** | 异步未并行 | 多个异步操作是否用 `Promise.all()`？ |
| **E002** | 轮询无超时 | 轮询是否设置 `maxAttempts`？ |
| **E003** | 错误未抛出 | `catch` 块是否 `throw error`？ |
| **E004** | SQL未用CTE | JOIN 后过滤 → 用 CTE 预过滤大表 |
| **E007** | 资源泄漏 | 所有退出路径都清理资源？ |

**完整错误目录**: `errors/ERROR_CATALOG.md` (E001-E015)

---

## 🧠 核心方法论（长任务）

### 三文件模式

```
task_plan.md     - 任务规划和进度追踪
notes.md         - 研究笔记和发现记录
[deliverable].md - 最终产出物
```

**关键**: 每个重要决策点前重新读取 `task_plan.md`

### 阶段门控

```
Phase 1: 需求理解 → [用户确认 "ready"]
Phase 2: 设计方案 → [确认]
Phase 3: 实现代码
```

**原则**: 永远不进入下一阶段，直到用户明确确认

---

## 🔧 能力库（按需加载）

### MCP Servers（外部数据）

| 任务 | MCP | 文档 |
|-----|-----|------|
| SQL查询 | bytebase | `capabilities/mcp-servers.md` |
| 浏览器 | playwright | `capabilities/browser-automation/` |
| 监控日志 | honeycomb | `capabilities/mcp-servers.md` |

### Skills（自动化任务）

| 类别 | 示例 | 文档 |
|------|------|------|
| Git工作流 | /commit, /create-pr | `capabilities/skills-guide.md` |
| 测试生成 | /write-tests | 同上 |
| UI设计 | ui-ux-pro-max | 同上 |
| 营销 | 24个专业Skills | `capabilities/MARKETING_SKILLS_GUIDE.md` |

**完整清单**: `capabilities/skills-guide.md` (81个Skills)

---

## 🤖 GPT 专家系统

### 可用专家

| 专家 | 专长 | 触发时机 |
|------|------|---------|
| Architect | 系统设计、技术决策 | 架构决策 / 2+失败尝试 |
| Plan Reviewer | 计划验证 | "review this plan" |
| Code Reviewer | 代码质量 | 完成功能后 |
| Security Analyst | 安全审计 | 安全问题 |

**详细文档**: `rules/delegator/`

---

## 📚 完整文档导航

### 索引层（快速查找）
- `index/task-router.md` - 任务路由决策树（30秒找到工具）
- `index/capabilities-index.md` - 能力索引
- `index/tools-index.md` - 工具索引
- `index/error-patterns-index.md` - 错误模式索引

### 规则库（自动加载）
- `rules/core/` - 核心规则（总是加载）
- `rules/domain/` - 领域规则（按需加载）
  - `coding.md`, `testing.md`, `security.md`, `git.md`

### 能力库（按需加载）
- `capabilities/browser-automation/` - 浏览器自动化
- `capabilities/video-creation/` - 视频创作
- `capabilities/data-analysis/` - 数据分析
- `design/` - UI 设计
- `vibe-marketing/` - 营销内容

### 错误库（按需加载）
- `errors/top-5-errors.md` - 高频错误快速参考
- `errors/ERROR_CATALOG.md` - 完整错误目录 (E001-E015)

### 知识库（参考）
- `KNOWLEDGE_MAP.md` - 知识图谱（12个Mermaid图）
- `QUICK_REFERENCE.md` - 一页速查表
- `INDEX.md` - 所有文档的完整索引

---

## 🔧 开发环境

- **OS**: Windows 10.0.26200 | **Shell**: Git Bash
- **路径格式**: Windows (Git Bash 中用正斜杠)
- **当前项目**: 数据分析和自动化（DAA）
- **技术栈**: TypeScript + PostgreSQL (Vercel) + MySQL (my_shell_prod)

---

## 📊 Context 使用优化

### Before (v4.2)
```
CLAUDE.md: 20KB
+ 多个规则文件: 30KB
+ 部分能力文档: 70KB
= 总计: 120KB (60% context)
```

### After (v5.0 - 目标)
```
CLAUDE.md: 5KB
+ 核心规则: 10KB
+ 按需加载 1-2 个文档: 20-30KB
= 总计: 35-45KB (18-23% context)
```

**节省**: ~75KB context (60% → 77%+ free context)

---

## 🎯 Context Engineering 核心

### 智能加载原则

1. **Layer 0** (总是加载): CLAUDE.md + 核心规则 (15KB)
2. **Layer 1** (任务识别): Task Router (3KB)
3. **Layer 2** (按需加载): 相关能力文档 (15-30KB)
4. **Layer 3** (精确匹配): 具体案例/模板 (按需)

### 你需要做什么？

**什么都不需要！** 只需描述你的任务，系统会自动：
- ✅ 识别任务类型
- ✅ 加载相关文档
- ✅ 保持 context 清洁
- ✅ 优化性能

---

## 📝 快速参考

### 我想做 X，应该...

1. **先看** `index/task-router.md` - 30秒找到工具
2. **再看** 具体能力文档 - 深入了解
3. **开始工作** - 系统已自动加载必需内容

### 遇到错误？

1. **先查** `errors/top-5-errors.md` - 5个高频错误
2. **再查** `errors/ERROR_CATALOG.md` - 完整错误库
3. **还不行** - 使用 Debugging Agent

### 需要专家意见？

1. **识别场景** - 见 `rules/delegator/triggers.md`
2. **选择专家** - Architect / Code Reviewer / Security Analyst
3. **自动触发** - 系统会识别并路由

---

## 🚀 开始工作

准备好接收任务了！

**记住**:
- ✅ 快速规划 → 展示 → 确认 → 执行
- ✅ 自行决策 95% 的问题
- ✅ 只在 4 种致命阻塞时提问
- ✅ 系统会自动加载必需文档

---

**版本**: v5.0 (Context Engineering)
**更新**: 2026-02-05
**大小**: ~5KB
**改进**: 从 v4.2 (20KB) 缩减 75%，实现智能按需加载

**升级自**: v4.2 (2026-01-28)
