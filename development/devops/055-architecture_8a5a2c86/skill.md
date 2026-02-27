# Architecture Reference

Complete structural blueprint for TanStack Start on Cloudflare Workers. Claude uses this to generate all project files from scratch.

## Project File Tree

```
PROJECT_NAME/
├── src/
│   ├── routes/
│   │   ├── __root.tsx              # Root layout (HTML shell, theme, CSS import)
│   │   ├── index.tsx               # Landing / auth redirect
│   │   ├── login.tsx               # Login page
│   │   ├── register.tsx            # Register page
│   │   ├── _authed.tsx             # Auth guard layout route
│   │   ├── _authed/
│   │   │   ├── dashboard.tsx       # Dashboard with stat cards
│   │   │   ├── items.tsx           # Items list table
│   │   │   ├── items.$id.tsx       # Edit item
│   │   │   └── items.new.tsx       # Create item
│   │   └── api/
│   │       └── auth/
│   │           └── $.ts            # better-auth API catch-all
│   ├── components/
│   │   ├── ui/                     # shadcn/ui components (auto-installed)
│   │   ├── app-sidebar.tsx         # Navigation sidebar
│   │   ├── theme-toggle.tsx        # Light/dark/system toggle
│   │   ├── user-nav.tsx            # User dropdown menu
│   │   └── stat-card.tsx           # Dashboard stat card
│   ├── db/
│   │   ├── schema.ts               # Drizzle schema (all tables)
│   │   └── index.ts                # Drizzle client factory
│   ├── lib/
│   │   ├── auth.server.ts          # better-auth server config
│   │   ├── auth.client.ts          # better-auth React hooks
│   │   └── utils.ts                # cn() helper for shadcn/ui
│   ├── server/
│   │   └── functions.ts            # Server functions (CRUD, auth checks)
│   ├── styles/
│   │   └── app.css                 # Tailwind v4 + shadcn/ui CSS variables
│   ├── router.tsx                  # TanStack Router configuration
│   ├── client.tsx                  # Client entry (hydrateRoot)
│   ├── ssr.tsx                     # SSR entry
│   └── routeTree.gen.ts            # Auto-generated route tree (do not edit)
├── drizzle/                        # Generated migrations
├── public/                         # Static assets (favicon, etc.)
├── vite.config.ts                  # Vite + Cloudflare + TanStack Start config
├── wrangler.jsonc                  # Cloudflare Workers config
├── drizzle.config.ts               # Drizzle Kit config
├── tsconfig.json                   # TypeScript config
├── package.json                    # Dependencies + scripts
├── .dev.vars                       # Local env vars (NOT committed)
└── .gitignore
```

## Dependencies

### Runtime

```json
{
  "react": "^19.0.0",
  "react-dom": "^19.0.0",
  "@tanstack/react-router": "^1.120.0",
  "@tanstack/react-start": "^1.120.0",
  "drizzle-orm": "^0.38.0",
  "better-auth": "^1.2.0",
  "zod": "^3.24.0",
  "class-variance-authority": "^0.7.0",
  "clsx": "^2.1.0",
  "tailwind-merge": "^3.0.0",
  "lucide-react": "^0.480.0"
}
```

### Dev

```json
{
  "@cloudflare/vite-plugin": "^1.0.0",
  "@tailwindcss/vite": "^4.0.0",
  "@vitejs/plugin-react": "^4.4.0",
  "tailwindcss": "^4.0.0",
  "typescript": "^5.7.0",
  "drizzle-kit": "^0.30.0",
  "wrangler": "^4.0.0",
  "tw-animate-css": "^1.2.0"
}
```

### package.json scripts

```json
{
  "dev": "vite",
  "build": "vite build",
  "preview": "vite preview",
  "deploy": "wrangler deploy",
  "db:generate": "drizzle-kit generate",
  "db:migrate:local": "wrangler d1 migrations apply PROJECT_NAME-db --local",
  "db:migrate:remote": "wrangler d1 migrations apply PROJECT_NAME-db --remote"
}
```

## Vite Configuration

Plugin order matters. Cloudflare must be first.

```typescript
// vite.config.ts
import { defineConfig } from "vite";
import { cloudflare } from "@cloudflare/vite-plugin";
import { tanstackStart } from "@tanstack/react-start/plugin/vite";
import tailwindcss from "@tailwindcss/vite";
import viteReact from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [
    cloudflare({ viteEnvironment: { name: "ssr" } }),
    tailwindcss(),
    tanstackStart(),
    viteReact(),
  ],
});
```

## Wrangler Configuration

```jsonc
// wrangler.jsonc
{
  "$schema": "node_modules/wrangler/config-schema.json",
  "name": "PROJECT_NAME",
  "compatibility_date": "2025-04-01",
  "compatibility_flags": ["nodejs_compat"],
  "main": "@tanstack/react-start/server-entry",
  "account_id": "ACCOUNT_ID",
  "d1_databases": [
    {
      "binding": "DB",
      "database_name": "PROJECT_NAME-db",
      "database_id": "DATABASE_ID",
      "migrations_dir": "drizzle"
    }
  ]
}
```

Key points:
- `main` MUST be `"@tanstack/react-start/server-entry"` — this is the Nitro server entry
- `nodejs_compat` is required (NOT `node_compat`)
- `migrations_dir` tells wrangler where Drizzle outputs migrations
- Add `account_id` to avoid interactive prompts during deploy

## TypeScript Configuration

```json
{
  "compilerOptions": {
    "target": "ES2022",
    "module": "ESNext",
    "moduleResolution": "bundler",
    "jsx": "react-jsx",
    "strict": true,
    "skipLibCheck": true,
    "esModuleInterop": true,
    "resolveJsonModule": true,
    "isolatedModules": true,
    "noEmit": true,
    "paths": {
      "@/*": ["./src/*"]
    },
    "types": ["@cloudflare/workers-types/2023-07-01"]
  },
  "include": ["src/**/*", "vite.config.ts"]
}
```

## Environment Variables

### Auth (required)

| Variable | Purpose | Example |
|----------|---------|---------|
| `BETTER_AUTH_SECRET` | Session encryption key | `openssl rand -hex 32` |
| `BETTER_AUTH_URL` | Auth base URL | `http://localhost:3000` (dev) |
| `TRUSTED_ORIGINS` | Allowed origins (comma-separated) | `http://localhost:3000,https://app.workers.dev` |

### Google OAuth (optional)

| Variable | Purpose |
|----------|---------|
| `GOOGLE_CLIENT_ID` | From Google Cloud Console |
| `GOOGLE_CLIENT_SECRET` | From Google Cloud Console |

### .dev.vars template

```
BETTER_AUTH_SECRET=<generated-hex-32>
BETTER_AUTH_URL=http://localhost:3000
TRUSTED_ORIGINS=http://localhost:3000
# GOOGLE_CLIENT_ID=
# GOOGLE_CLIENT_SECRET=
```
