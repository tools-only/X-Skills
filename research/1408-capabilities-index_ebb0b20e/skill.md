# 能力索引 - Capabilities Index

> **用途**: 快速查找所有可用的能力和文档

---

## 📚 核心能力（按领域分类）

### 🌐 浏览器自动化

| 能力 | 工具 | 文档 | 场景 |
|------|------|------|------|
| **实时操作** | Playwright MCP ⭐ | `capabilities/browser-automation/decision-tree.md` | 对话中操作网页 |
| **批量/脚本** | agent-browser CLI 🚀 | 同上 | >50次操作、定时任务 |
| **网络拦截** | Playwright MCP | 同上 | 录制、修改请求 |

**Skills**: agent-browser, webapp-testing

---

### 🎬 视频创作

| 能力 | 工具 | 文档 | 场景 |
|------|------|------|------|
| **快速模板** | Remotion Templates | `capabilities/REMOTION_TEMPLATES_LIBRARY.md` | 标题、列表、音乐可视化 |
| **自动生成** | Remotion | `rules/remotion-auto-production.md` | 完整视频项目 |
| **创意背景** | Processing | `capabilities/PROCESSING_SKILL.md` | 粒子特效、几何动画 |

**模板数量**: 15个现成模板
**自动化**: 风格匹配、配色、场景设计

---

### 📊 数据分析

| 能力 | 工具 | 文档 | 场景 |
|------|------|------|------|
| **SQL 查询** | bytebase MCP | `capabilities/sql-patterns.md` | MySQL 数据库 |
| **图表生成** | chart MCP | `capabilities/data-viz-guide.md` | 可视化 |
| **监控日志** | honeycomb/guance MCP | `capabilities/mcp-servers.md` | 日志分析 |

**Skills**: 8个数据分析 Skills（Bot毛利率、成本趋势等）

---

### 🎨 UI/UX 设计

| 能力 | 工具 | 文档 | 场景 |
|------|------|------|------|
| **智能设计** | ui-ux-pro-max Skill | `design/DESIGN_MASTER_PERSONA.md` | 完整 UI 设计 |
| **设计风格** | 30种风格库 | `design/UI_DESIGN_STYLES_REFERENCE.md` | 风格选择 |
| **组件库** | theme-factory Skill | `design/design-systems-guide.md` | 设计系统 |

**风格数量**: 30种设计风格
**Skills**: ui-ux-pro-max, frontend-design-offical, canvas-design

---

### 📝 营销内容

| 能力 | 工具 | 文档 | 场景 |
|------|------|------|------|
| **Vibe Marketing** | 研究方法论 | `vibe-marketing/VIBE_MARKETING_GUIDE.md` | 2周研究→1小时 |
| **专业 Skills** | 24个营销 Skills | `capabilities/MARKETING_SKILLS_GUIDE.md` | CRO/SEO/文案 |
| **社交内容** | social-content Skill | 同上 | 社媒内容 |

**Skills 分类**:
- 文案: copywriting, email-sequence
- SEO: seo-audit, schema-markup
- CRO: page-cro, signup-flow-cro, popup-cro
- 策略: pricing-strategy, launch-strategy

---

### 💻 代码开发

| 能力 | 工具 | 文档 | 场景 |
|------|------|------|------|
| **Git 工作流** | /commit, /create-pr | `rules/domain/git.md` | 提交、PR、审查 |
| **测试生成** | /write-tests | `rules/domain/testing.md` | 单元/集成测试 |
| **代码审查** | /code-review | `rules/domain/coding.md` | 质量检查 |

**Plugins（自动激活）**:
- backend-development - 后端架构
- frontend-mobile-development - 前端/移动
- code-documentation - 文档生成

---

### 🔒 安全审计

| 能力 | 工具 | 文档 | 场景 |
|------|------|------|------|
| **安全扫描** | security-scanning Plugin | `rules/domain/security.md` | OWASP Top 10 |
| **GPT 专家** | Security Analyst | `rules/delegator/` | 深度审计 |

**检查项**:
- SQL 注入、XSS、CSRF
- 认证/授权
- 密钥管理
- 速率限制

---

### 🤖 GPT 专家系统

| 专家 | 专长 | 文档 | 触发时机 |
|------|------|------|---------|
| **Architect** | 系统设计 | `rules/delegator/model-selection.md` | 架构决策 |
| **Plan Reviewer** | 计划验证 | 同上 | "review this plan" |
| **Scope Analyst** | 需求分析 | 同上 | 模糊需求 |
| **Code Reviewer** | 代码审查 | 同上 | 完成功能后 |
| **Security Analyst** | 安全审计 | 同上 | 安全问题 |

**详细文档**: `rules/delegator/`

---

## 🔧 特殊能力

### PPT 制作

| 能力 | 输出格式 | 文档 |
|------|---------|------|
| PPT 自动化 | .pptx + HTML + 图片 | `capabilities/PPT_WORKFLOW.md` |

### 创意编程

| 能力 | 特性 | 文档 |
|------|------|------|
| Processing | 6种视觉模式，16种配色 | `capabilities/PROCESSING_SKILL.md` |

### 视频编辑

| 能力 | 特性 | 文档 |
|------|------|------|
| Remotion | 自动风格匹配 | `rules/remotion-auto-production.md` |

---

## 📱 平台集成

### MCP Servers（外部数据访问）

| MCP | 功能 | 文档 |
|-----|------|------|
| **bytebase** | MySQL 查询 | `capabilities/mcp-servers.md` |
| **playwright** | 浏览器自动化 | `capabilities/browser-automation/` |
| **honeycomb** | 监控日志 | `capabilities/mcp-servers.md` |
| **guance** | 观测云 | `capabilities/mcp-servers.md` |
| **chart** | 图表生成 | `capabilities/mcp-servers.md` |
| **stitch** | Google Stitch 项目 | MCP 集成 |
| **firebase** | Firebase 项目管理 | MCP 集成 |

**完整列表**: `capabilities/mcp-servers.md`

---

## 🎯 Skills 完整清单

### Git & 代码质量（8个）
- /commit, /create-pr, /code-review
- /write-tests, /refactor, /debug
- /add-feature, /fix-bug

### 设计（5个）
- ui-ux-pro-max, frontend-design-offical
- canvas-design, theme-factory
- image-editor, image-enhancer

### 营销（24个）
完整列表见 `capabilities/MARKETING_SKILLS_GUIDE.md`

### 数据分析（8个）
位于 `bo-skill-research/shane-skill/data-analysis-agent/skills/`

### 创意（3个）
- nano-banana-pro（图片生成）
- visual-prompt-engineer（Prompt 工程）
- processing-creative（Processing）

### 研究与写作（8个）
- deep-research, literature-review
- citation-validator, synthesizer
- content-research-writer

### 其他（25+个）
- agent-browser（浏览器）
- mcp-builder（MCP 开发）
- skill-creator（Skill 开发）
- notebooklm（知识管理）
- 等等...

**完整清单**: `capabilities/skills-guide.md`（81个Skills）

---

## 🔍 按场景查找

### 我想要自动化...

| 场景 | 推荐能力 | 文档 |
|------|---------|------|
| 网页操作 | Playwright MCP | `browser-automation/` |
| 数据抓取 | agent-browser CLI | 同上 |
| UI 测试 | webapp-testing Skill | `skills-guide.md` |

### 我想要创作...

| 场景 | 推荐能力 | 文档 |
|------|---------|------|
| 产品视频 | Remotion + 模板 | `remotion-auto-production.md` |
| 音乐可视化 | sound-wave 模板 | `REMOTION_TEMPLATES_LIBRARY.md` |
| 创意背景 | Processing | `PROCESSING_SKILL.md` |

### 我想要分析...

| 场景 | 推荐能力 | 文档 |
|------|---------|------|
| 数据库查询 | bytebase MCP | `sql-patterns.md` |
| 日志分析 | honeycomb MCP | `mcp-servers.md` |
| 生成图表 | chart MCP | `data-viz-guide.md` |

### 我想要设计...

| 场景 | 推荐能力 | 文档 |
|------|---------|------|
| UI 界面 | ui-ux-pro-max | `DESIGN_MASTER_PERSONA.md` |
| 选择风格 | 30种风格库 | `UI_DESIGN_STYLES_REFERENCE.md` |
| 组件库 | theme-factory | `design-systems-guide.md` |

---

## 📚 能力文档位置

所有能力文档位于 `~/.claude/capabilities/`：

```
capabilities/
├── browser-automation/
│   └── decision-tree.md
├── PROCESSING_SKILL.md
├── REMOTION_TEMPLATES_LIBRARY.md
├── MARKETING_SKILLS_GUIDE.md
├── skills-guide.md
├── mcp-servers.md
├── data-analysis-guide.md
├── sql-patterns.md
└── ...
```

---

## 🎓 学习路径

### 新手（15分钟）
1. 浏览本索引（5分钟）
2. 查看 `task-router.md`（5分钟）
3. 尝试一个简单任务（5分钟）

### 进阶（1小时）
1. 深入一个领域文档（30分钟）
2. 尝试复杂任务（30分钟）

### 专家（持续）
- 探索所有能力文档
- 组合多个能力解决复杂问题
- 贡献新的 Skills 和能力

---

**版本**: v1.0
**更新**: 2026-02-05
**能力总数**: 100+ 能力，81个 Skills，12个 MCP Servers
