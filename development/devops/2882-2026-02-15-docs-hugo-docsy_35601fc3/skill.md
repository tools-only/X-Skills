# Docs Hugo Docsy Migration Implementation Plan
> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Convert `docs/` into a Hugo site using Docsy, with concise docs focused on introduction, installation (docker-compose + docker-swarm), and usage, while reserving a multilingual structure with English only.

**Architecture:** Hugo site rooted at `docs/` using Hugo Modules to import Docsy. Content lives under `content/en/` with a clear section structure. Static assets move to `static/images/`. Legacy Markdown in `docs/` is removed or relocated after migration to avoid duplicate/ambiguous sources.

**Tech Stack:** Hugo (extended), Docsy (Hugo module), Markdown, YAML config.

---

### Task 1: Initialize Hugo + Docsy site scaffolding

**Files:**
- Create: `docs/hugo.yaml`
- Create: `docs/go.mod`
- Create: `docs/content/en/_index.md`
- Create: `docs/content/en/docs/_index.md`

**Step 1: Run a failing build to confirm no config exists**

Run: `cd docs && hugo`
Expected: FAIL with error similar to "unable to find config file".

**Step 2: Add Hugo configuration with Docsy module**

Create `docs/hugo.yaml`:
```yaml
baseURL: "https://docs.example.com/"
title: "IntentKit"
languageCode: "en-us"
defaultContentLanguage: "en"
defaultContentLanguageInSubdir: true

module:
  imports:
    - path: github.com/google/docsy

languages:
  en:
    languageName: "English"
    contentDir: "content/en"
    title: "IntentKit"

menu:
  main:
    - name: "Docs"
      url: "/docs/"
      weight: 10
    - name: "GitHub"
      url: "https://github.com/crestalnetwork/intentkit"
      weight: 20

params:
  ui:
    sidebar_menu_compact: true
    breadcrumb_disable: false
  links:
    developer:
      - name: "GitHub"
        url: "https://github.com/crestalnetwork/intentkit"
        icon: "fab fa-github"
        desc: "Source code"
```

Create `docs/go.mod`:
```go
module github.com/crestalnetwork/intentkit/docs

go 1.22
```

Create `docs/content/en/_index.md`:
```markdown
---
title: "IntentKit"
linkTitle: "IntentKit"
---

IntentKit is an agent system and runtime for building, managing, and integrating autonomous agents.
```

Create `docs/content/en/docs/_index.md`:
```markdown
---
title: "Documentation"
linkTitle: "Docs"
weight: 1
---
```

**Step 3: Run Hugo module tidy and rebuild**

Run: `cd docs && hugo mod tidy`
Expected: Go module downloads Docsy and creates `go.sum`.

Run: `cd docs && hugo`
Expected: PASS, site builds into `docs/public/`.

**Step 4: Commit**

```bash
git add docs/hugo.yaml docs/go.mod docs/go.sum docs/content/en/_index.md docs/content/en/docs/_index.md
git commit -m "feat: initialize hugo docsy site"
```

---

### Task 2: Create introduction, installation, and usage docs (English only)

**Files:**
- Create: `docs/content/en/docs/introduction.md`
- Create: `docs/content/en/docs/installation/_index.md`
- Create: `docs/content/en/docs/installation/docker-compose.md`
- Create: `docs/content/en/docs/installation/docker-swarm.md`
- Create: `docs/content/en/docs/usage/_index.md`
- Create: `docs/content/en/docs/usage/getting-started.md`

**Step 1: Run a failing content check**

Run: `cd docs && hugo --printPathWarnings`
Expected: FAIL or warnings about missing sections (no installation/usage pages).

**Step 2: Add the introduction page**

Create `docs/content/en/docs/introduction.md`:
```markdown
---
title: "Introduction"
weight: 1
---

IntentKit provides a modular agent system with an API server, autonomous runner, scheduler, and integrations.
Use it to create and run agents, manage skills, and integrate them into existing products.
```

**Step 3: Add installation section with docker-compose**

Create `docs/content/en/docs/installation/_index.md`:
```markdown
---
title: "Installation"
weight: 2
---
```

Create `docs/content/en/docs/installation/docker-compose.md`:
```markdown
---
title: "Docker Compose"
weight: 1
---

Prerequisites:
- Docker and Docker Compose

Clone the repo and start the stack:

```bash
cp example.env .env
docker compose up -d
```

The API becomes available at:
- http://localhost:8000
```

**Step 4: Add installation section with docker-swarm**

Create `docs/content/en/docs/installation/docker-swarm.md`:
```markdown
---
title: "Docker Swarm"
weight: 2
---

Prerequisites:
- Docker with Swarm mode enabled

Initialize swarm and deploy:

```bash
docker swarm init
docker stack deploy -c docs/deployment/docker-compose.yml intentkit
```

Check services:

```bash
docker service ls
```
```

**Step 5: Add usage section**

Create `docs/content/en/docs/usage/_index.md`:
```markdown
---
title: "Usage"
weight: 3
---
```

Create `docs/content/en/docs/usage/getting-started.md`:
```markdown
---
title: "Getting Started"
weight: 1
---

1. Start the API server with Docker Compose or Swarm.
2. Use the Agent API at `http://localhost:8000/redoc` to explore endpoints.
3. Create an agent configuration and interact with it via the API.
```

**Step 6: Rebuild and verify**

Run: `cd docs && hugo --minify`
Expected: PASS, `docs/public/` updated without errors.

**Step 7: Commit**

```bash
git add docs/content/en/docs
git commit -m "feat: add intro install usage docs"
```

---

### Task 3: Migrate assets and remove legacy Markdown

**Files:**
- Move: `docs/images/intentkit_banner.png` -> `docs/static/images/intentkit_banner.png`
- Move: `docs/agent.webp` -> `docs/static/images/agent.webp`
- Move: `docs/arch.jpg` -> `docs/static/images/arch.jpg`
- Remove: legacy Markdown under `docs/` that is superseded
- Modify: `README.md`

**Step 1: Run a failing reference check**

Run: `cd docs && rg "images/" content`
Expected: FAIL (no references updated yet).

**Step 2: Move images to Hugo static directory**

```bash
mkdir -p docs/static/images
mv docs/images/intentkit_banner.png docs/static/images/intentkit_banner.png
mv docs/agent.webp docs/static/images/agent.webp
mv docs/arch.jpg docs/static/images/arch.jpg
```

**Step 3: Remove legacy Markdown to avoid duplicate sources**

Remove:
- `docs/README.md`
- `docs/agent.md`
- `docs/agent_api.md`
- `docs/architecture.md`
- `docs/configuration.md`
- `docs/discord.md`
- `docs/llm.md`
- `docs/openai_compatible.md`
- `docs/how_to/readme.md`
- `docs/how_to/clean_memory.md`
- `docs/skills/*.md`
- `docs/contributing/*.md`

**Step 4: Update root README to point to the new docs**

Modify `README.md` to reference the Hugo docs path (e.g., `docs/` or the hosted site URL).

**Step 5: Rebuild and verify**

Run: `cd docs && hugo --minify`
Expected: PASS with no broken paths.

**Step 6: Commit**

```bash
git add docs README.md
git commit -m "refactor: migrate docs to hugo docsy"
```
