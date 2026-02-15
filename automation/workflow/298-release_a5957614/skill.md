# OpenAkita 发版流程

本文档记录 OpenAkita 的标准发版操作流程。**所有贡献者（包括 AI 辅助工具）发版时必须遵循此流程。**

---

## 分支策略

```
main (开发分支)
 │── 新功能开发、大版本迭代
 │── tag 从此分支打出的发布标记为 pre-release（开发版）
 │
 ├── v1.21.x (定稿/稳定分支)
 │    │── 从 main 分出，只接受 bugfix
 │    │── tag v1.21.9, v1.21.10... → 正式稳定版
 │    └── 用户默认下载和自动更新走这个渠道
 │
 ├── v1.22.x (下一个稳定分支)
 │    └── ...
 └── ...
```

### 规则

| 分支 | 用途 | tag 发布后 | latest.json |
|---|---|---|---|
| `main` | 开发、新功能 | 标记为 **pre-release** | 更新 `latest-dev.json` |
| `v{x}.{y}.x` | 定稿稳定、只修 bug | 标记为 **正式版** | 更新 `latest.json`（用户可见） |

### 创建稳定分支

当 `main` 上的功能趋于稳定，准备发一个正式版时：

```bash
# 从 main 创建稳定分支
git checkout main
git checkout -b v1.22.x
git push origin v1.22.x

# 在稳定分支上打 tag
python scripts/version.py set 1.22.0
git add . && git commit -m "chore: bump version to v1.22.0"
git tag v1.22.0
git push origin v1.22.x v1.22.0
```

之后 `main` 继续开发 v1.23.0+，`v1.22.x` 只接收 cherry-pick 过来的 bugfix。

### 在稳定分支上修 bug

```bash
git checkout v1.22.x
# 修复代码...
python scripts/version.py set 1.22.1
git add . && git commit -m "fix: 修复 xxx 问题"
git tag v1.22.1
git push origin v1.22.x v1.22.1
```

---

## 发版步骤（适用于所有分支）

### Step 1: 前置准备

1. 确认代码在目标分支上（`main` 或 `v{x}.{y}.x`）
2. 确认 CI 绿灯（CI 同时保护 `main` 和 `v*.x` 分支）
3. 运行版本号同步：
   ```bash
   python scripts/version.py set <new_version>
   ```
4. 提交版本号变更并推送

### Step 2: 预验证（可选但推荐）

- 在 GitHub Actions 页面手动触发 **"Release Dry Run"** workflow
- 等待三个平台（Windows / Linux / macOS）全部构建成功

### Step 3: 打 tag 触发正式发布

```bash
git tag v<version>
git push origin v<version>
```

### Step 4: 等待 CI 完成

- Release 会先以 **Draft** 状态创建（用户不可见）
- CI 自动检测 tag 来源：
  - 来自 `v*.x` 稳定分支 → 发布为正式版，更新 `latest.json`
  - 来自 `main` → 发布为 pre-release，更新 `latest-dev.json`
- 全平台 Desktop 构建成功后，Release 自动发布

---

## 如果 CI 构建失败

**绝不要为了修 CI 而递增版本号**（如 v1.21.9 → v1.21.10）。

### 方法 1: 强制更新 tag（推荐）

```bash
# 修复代码（在对应分支上）
git add . && git commit -m "fix: 修复构建问题"
git push origin <当前分支>

# 强制更新 tag
git tag -f v<version>
git push origin v<version> --force
```

### 方法 2: 手动重跑 workflow

GitHub Actions 页面 → 找到失败的 Release workflow → **"Re-run all jobs"**

### 方法 3: workflow_dispatch 手动触发

手动触发 Release workflow，指定 `tag` 和 `ref` 参数。

---

## 注意事项

- **Draft Release 对用户不可见**，放心重跑
- **PyPI 发布是幂等的**（`twine upload --skip-existing`）
- 发布成功后，CI 自动更新 `latest.json`（稳定版）或 `latest-dev.json`（开发版）
- 应用内自动更新默认检查 `latest.json`（稳定渠道），用户只会收到稳定版更新
- 官网下载页默认展示 `latest.json` 中的稳定版

---

## Workflow 文件说明

| Workflow | 文件 | 触发 | 用途 |
|---|---|---|---|
| CI | `.github/workflows/ci.yml` | push/PR to `main`、`v*.x` | lint + test + 完整 Tauri 构建检查 |
| Release Dry Run | `.github/workflows/release-dryrun.yml` | 手动触发 | 全平台构建预验证，不发布 |
| Release | `.github/workflows/release.yml` | tag push `v*` | Draft → Build → 渠道检测 → Publish |
