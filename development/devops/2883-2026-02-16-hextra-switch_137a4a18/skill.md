# Switch Docs Theme to Hextra Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace Docsy with Hextra for the Hugo docs site, clean Docsy-specific artifacts, and apply common Hextra configuration suited to IntentKit.

**Architecture:** Keep Hugo site rooted at `docs/`, preserve `content/en/` structure, and swap the theme module to Hextra. Update site configuration to Hextra’s menu, navbar, and edit-link patterns, and remove Docsy/PostCSS dependencies and Docker build steps that were only needed for Docsy.

**Tech Stack:** Hugo (extended), Hextra (Hugo module), YAML config, Markdown.

---

### Task 1: Capture current Docsy artifacts and verify build baseline

**Files:**
- Inspect: `docs/hugo.yaml`
- Inspect: `docs/go.mod`
- Inspect: `docs/go.sum`
- Inspect: `docs/Dockerfile`
- Inspect: `docs/package.json`
- Inspect: `docs/postcss.config.js`

**Step 1: Confirm current theme and artifacts**

Run: `sed -n '1,200p' docs/hugo.yaml`
Expected: Docsy module import and Docsy params present.

**Step 2: Run a baseline build**

Run: `cd docs && hugo --minify`
Expected: PASS, current Docsy build succeeds.

**Step 3: Commit**

```bash
git add -A
git commit -m "chore: snapshot docsy docs state"
```

---

### Task 2: Switch Hugo module to Hextra and update configuration

**Files:**
- Modify: `docs/hugo.yaml`
- Modify: `docs/go.mod`
- Modify: `docs/go.sum`

**Step 1: Update Hugo config for Hextra**

Replace Docsy module import with Hextra and update params to common Hextra patterns (menu search item, GitHub icon, navbar settings, edit link). Use Hextra configuration guidance for menu items and navbar logo/title options. [Navigation config](https://imfing.github.io/hextra/docs/guide/configuration/), [Logo and title config](https://imfing.github.io/hextra/docs/guide/configuration/).

Target config (example):
```yaml
baseURL: "https://docs.example.com/"
title: "IntentKit"
languageCode: "en-us"
defaultContentLanguage: "en"
defaultContentLanguageInSubdir: true

module:
  imports:
    - path: github.com/imfing/hextra

languages:
  en:
    languageName: "English"
    contentDir: "content/en"
    title: "IntentKit"

menu:
  main:
    - name: "Docs"
      pageRef: /docs
      weight: 1
    - name: "Search"
      weight: 2
      params:
        type: search
    - name: "Theme"
      weight: 3
      params:
        type: theme-toggle
        label: true
    - name: "GitHub"
      weight: 4
      url: "https://github.com/crestalnetwork/intentkit"
      params:
        icon: github

params:
  navbar:
    displayTitle: true
    displayLogo: false
  editURL:
    enable: true
    base: "https://github.com/crestalnetwork/intentkit/edit/main/docs/content"
```

**Step 2: Update go.mod / go.sum**

Run: `cd docs && hugo mod tidy`
Expected: go.sum updated with Hextra module, Docsy removed.

**Step 3: Commit**

```bash
git add docs/hugo.yaml docs/go.mod docs/go.sum
git commit -m "feat: switch docs theme to hextra"
```

---

### Task 3: Remove Docsy/PostCSS-specific artifacts and update Docker build

**Files:**
- Delete: `docs/package.json`
- Delete: `docs/package-lock.json`
- Delete: `docs/postcss.config.js`
- Modify: `docs/Dockerfile`

**Step 1: Remove PostCSS files**

Delete the Docsy-only PostCSS artifacts and any leftover lockfile.

**Step 2: Simplify Dockerfile build**

Remove node/npm install steps and rely on Hugo build only. Example:
```dockerfile
FROM gohugoio/hugo:0.155.3-ext-alpine AS builder
RUN apk add --no-cache git
WORKDIR /src
COPY . .
RUN hugo --minify

FROM caddy:2-alpine
WORKDIR /srv
COPY --from=builder /src/public /srv
EXPOSE 80
CMD ["caddy", "file-server", "--root", "/srv", "--listen", ":80"]
```

**Step 3: Verify build**

Run: `cd docs && hugo --minify`
Expected: PASS with Hextra output.

**Step 4: Commit**

```bash
git add docs/Dockerfile docs/package.json docs/package-lock.json docs/postcss.config.js
git commit -m "chore: remove docsy build artifacts"
```

---

### Task 4: Optional Hextra look-and-feel polish (project-specific)

**Files:**
- Create: `docs/assets/css/custom.css`

**Step 1: Add a small brand color and width tweak**

Add a minimal custom CSS file to adjust primary color and keep layout width aligned with product docs. Use Hextra custom CSS conventions. [Custom CSS](https://imfing.github.io/hextra/docs/advanced/customization/).

```css
:root {
  --primary-hue: 220deg;
  --primary-saturation: 80%;
  --primary-lightness: 45%;
}
```

**Step 2: Verify build**

Run: `cd docs && hugo --minify`
Expected: PASS.

**Step 3: Commit**

```bash
git add docs/assets/css/custom.css
git commit -m "style: add hextra brand color"
```
