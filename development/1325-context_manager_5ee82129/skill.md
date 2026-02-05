# Context Manager

> **职责**: 智能决定加载哪些文档，保持 context 清洁和高效

---

## 🎯 核心原则

### Progressive Disclosure
- **只加载必需的内容**
- **按任务类型分层加载**
- **避免一次性加载所有文档**

### Context 预算管理
- **总预算**: 200K tokens
- **目标使用**: <40% (80K tokens)
- **留给工作**: >60% (120K tokens)

---

## 🤖 自动加载机制

### 三层加载策略

```
Layer 0: 总是加载 (核心)
  ├─ CLAUDE.md (5KB)
  ├─ rules/core/work-mode.md (3KB)
  └─ rules/core/blocking-rules.md (2KB)
  = 10KB

Layer 1: 按任务类型加载 (领域)
  根据任务选择 1-2 个领域规则
  = 10-15KB

Layer 2: 按需加载 (能力)
  根据具体需求加载 1-3 个能力文档
  = 15-30KB

总计: 35-55KB (17-28% context 使用)
```

---

## 📋 任务类型映射表

### 浏览器自动化

```yaml
task_type: browser-automation
keywords: ["浏览器", "自动化", "爬虫", "测试", "playwright"]
load:
  core: ["work-mode.md", "blocking-rules.md"]
  domain: ["rules/domain/coding.md"]
  capabilities: ["capabilities/browser-automation/decision-tree.md"]
  estimated_size: 15KB
```

### 视频创作

```yaml
task_type: video-creation
keywords: ["视频", "Remotion", "动画", "特效"]
load:
  core: ["work-mode.md", "blocking-rules.md"]
  domain: ["rules/domain/coding.md"]
  capabilities:
    - "rules/remotion-auto-production.md"
    - "capabilities/PROCESSING_SKILL.md"
  estimated_size: 25KB
```

### 数据分析

```yaml
task_type: data-analysis
keywords: ["数据", "分析", "SQL", "报表", "图表"]
load:
  core: ["work-mode.md", "blocking-rules.md"]
  domain: ["rules/domain/coding.md"]
  capabilities:
    - "capabilities/data-analysis-guide.md"
    - "capabilities/sql-patterns.md"
  estimated_size: 20KB
```

### UI 设计

```yaml
task_type: ui-design
keywords: ["设计", "UI", "界面", "风格", "布局"]
load:
  core: ["work-mode.md"]
  domain: []
  capabilities:
    - "design/DESIGN_MASTER_PERSONA.md"
    - "design/UI_DESIGN_STYLES_REFERENCE.md"
  estimated_size: 30KB
```

### 营销内容

```yaml
task_type: marketing
keywords: ["营销", "文案", "SEO", "广告", "策略"]
load:
  core: ["work-mode.md"]
  domain: []
  capabilities:
    - "vibe-marketing/VIBE_MARKETING_GUIDE.md"
    - "capabilities/MARKETING_SKILLS_GUIDE.md"
  estimated_size: 35KB
```

### 代码开发

```yaml
task_type: coding
keywords: ["开发", "代码", "功能", "实现", "bug"]
load:
  core: ["work-mode.md", "blocking-rules.md"]
  domain:
    - "rules/domain/coding.md"
    - "rules/domain/testing.md"
    - "rules/domain/git.md"
  capabilities: []
  estimated_size: 15KB
```

### 错误调试

```yaml
task_type: debugging
keywords: ["错误", "bug", "调试", "失败", "异常"]
load:
  core: ["work-mode.md", "blocking-rules.md"]
  domain: ["rules/domain/coding.md"]
  capabilities:
    - "errors/top-5-errors.md"  # 只加载高频错误
  on_demand:
    - "errors/ERROR_CATALOG.md"  # 需要时再加载完整目录
  estimated_size: 12KB
```

### 安全审计

```yaml
task_type: security
keywords: ["安全", "漏洞", "审计", "权限", "加密"]
load:
  core: ["work-mode.md", "blocking-rules.md"]
  domain: ["rules/domain/security.md"]
  capabilities: ["capabilities/security-best-practices.md"]
  estimated_size: 15KB
```

---

## 🔍 任务识别算法

### 关键词匹配

```typescript
function identify_task_type(user_request: string): string {
  const keywords_map = {
    "browser-automation": ["浏览器", "自动化", "爬虫", "playwright", "测试网页"],
    "video-creation": ["视频", "Remotion", "动画", "特效", "剪辑"],
    "data-analysis": ["数据", "分析", "SQL", "查询", "报表", "图表"],
    "ui-design": ["设计", "UI", "界面", "风格", "布局", "页面"],
    "marketing": ["营销", "文案", "SEO", "广告", "推广", "策略"],
    "coding": ["开发", "代码", "功能", "实现", "编写"],
    "debugging": ["错误", "bug", "调试", "失败", "不工作", "异常"],
    "security": ["安全", "漏洞", "审计", "权限", "加密", "认证"]
  };

  for (const [task_type, keywords] of Object.entries(keywords_map)) {
    for (const keyword of keywords) {
      if (user_request.includes(keyword)) {
        return task_type;
      }
    }
  }

  return "general";  // 默认通用任务
}
```

### 多任务场景

如果识别出多个任务类型，加载它们的并集：

```typescript
function load_for_multiple_tasks(task_types: string[]): LoadPlan {
  const core = ["work-mode.md", "blocking-rules.md"];  // 总是加载

  let domain = new Set<string>();
  let capabilities = new Set<string>();

  for (const task_type of task_types) {
    const plan = get_load_plan(task_type);
    plan.domain.forEach(d => domain.add(d));
    plan.capabilities.forEach(c => capabilities.add(c));
  }

  return {
    core,
    domain: Array.from(domain),
    capabilities: Array.from(capabilities),
    estimated_size: estimate_total_size([...core, ...domain, ...capabilities])
  };
}
```

---

## 📊 加载决策流程

```
用户请求
  ↓
提取关键词
  ↓
识别任务类型
  ↓
┌─────────────────┐
│ 单一任务类型?    │
└─────────────────┘
  ├─ 是 → 加载对应 context
  └─ 否 → 加载多任务 context (取并集)
  ↓
检查 context 大小
  ├─ <50KB → 直接加载
  ├─ 50-80KB → 警告但继续
  └─ >80KB → 要求用户明确优先级
  ↓
加载到内存
  ↓
开始任务
```

---

## 🎯 使用示例

### 示例 1: 简单任务

**用户**: "帮我写一个自动化测试脚本，用 Playwright 打开网页并截图"

**识别**: `browser-automation`

**加载**:
```
core: work-mode.md + blocking-rules.md (10KB)
domain: coding.md (5KB)
capabilities: browser-automation/decision-tree.md (8KB)
总计: 23KB
```

### 示例 2: 复合任务

**用户**: "做一个产品介绍视频，需要数据可视化图表"

**识别**: `video-creation` + `data-analysis`

**加载**:
```
core: work-mode.md + blocking-rules.md (10KB)
domain: coding.md (5KB)
capabilities:
  - remotion-auto-production.md (12KB)
  - data-analysis-guide.md (10KB)
总计: 37KB
```

### 示例 3: 模糊任务

**用户**: "优化一下网站性能"

**识别**: `general` (无明确关键词)

**操作**: 询问用户具体方向
- 前端性能？→ `coding` + `capabilities/frontend-optimization.md`
- 数据库查询？→ `coding` + `capabilities/sql-patterns.md`
- 服务器配置？→ `coding` + `capabilities/devops-guide.md`

---

## ⚙️ 配置

### Context 预算阈值

```yaml
thresholds:
  safe: 50KB      # 安全区，直接加载
  warning: 80KB   # 警告区，提示用户
  limit: 100KB    # 限制区，要求明确优先级
```

### 自动优化策略

```yaml
optimization:
  # 如果超过警告阈值，自动采取措施
  - 移除示例代码块（保留文字说明）
  - 压缩重复内容
  - 只加载"快速参考"部分，不加载详细案例
```

---

## 📝 开发者接口

### 手动指定加载内容

如果自动识别不准确，可以手动指定：

```
@load capabilities/browser-automation/decision-tree.md
@load rules/domain/testing.md
```

### 动态加载

任务执行过程中，如果发现需要额外文档：

```
[发现需要 Remotion 模板库]
→ 动态加载 capabilities/remotion-templates-library.md
→ 继续任务
```

---

## 🔄 持续优化

### 监控指标

- **平均 context 使用**: 目标 <40%
- **加载准确率**: 目标 >90% (加载的文档确实用到了)
- **任务完成率**: 目标 >95% (没有因为缺少文档而失败)

### 反馈循环

```
任务完成后
  ↓
记录实际使用的文档
  ↓
与预加载的文档对比
  ↓
识别未使用的文档 (过度加载)
识别缺失的文档 (加载不足)
  ↓
更新映射表
```

---

**版本**: v1.0
**更新**: 2026-02-05
**状态**: Active
