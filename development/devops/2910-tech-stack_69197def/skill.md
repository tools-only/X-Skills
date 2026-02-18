# Tech Stack Reference

Complete dependency list for vite-flare-starter.

## Core Stack

| Technology | Version | Purpose |
|-----------|---------|---------|
| React | 19 | UI framework |
| Vite | 6.x | Build tool + dev server |
| Hono | 4.x | Backend API framework (on Workers) |
| Cloudflare Workers | — | Runtime + deployment |
| D1 | — | Edge SQLite database |
| Drizzle ORM | 0.38+ | Type-safe database queries |
| better-auth | latest | Authentication (Google OAuth + email/password) |
| Tailwind CSS | v4 | Styling |
| shadcn/ui | latest | UI component library |
| TanStack Query | v5 | Server state management |

## Runtime Dependencies

| Package | Purpose |
|---------|---------|
| `react`, `react-dom` | UI rendering |
| `hono` | API routes and middleware |
| `better-auth` | Auth (server + client) |
| `drizzle-orm` | Database ORM |
| `@tanstack/react-query` | Data fetching + caching |
| `react-hook-form` | Form handling |
| `@hookform/resolvers` | Form validation |
| `zod` | Schema validation |
| `react-router-dom` | Client-side routing |
| `@radix-ui/*` | Headless UI primitives (shadcn/ui) |
| `tailwind-merge` | Tailwind class merging |
| `class-variance-authority` | Variant styling |
| `clsx` | Conditional classnames |
| `lucide-react` | Icon library |
| `sonner` | Toast notifications |
| `date-fns` | Date formatting |

## Dev Dependencies

| Package | Purpose |
|---------|---------|
| `@cloudflare/vite-plugin` | Workers + Vite integration |
| `wrangler` | Cloudflare CLI |
| `drizzle-kit` | Database migrations |
| `typescript` | Type checking |
| `@tailwindcss/vite` | Tailwind v4 Vite plugin |
| `vitest` | Testing |
| `eslint` | Linting |

## Cloudflare Bindings

| Binding | Type | Purpose |
|---------|------|---------|
| `DB` | D1 Database | Primary application database |
| `AVATARS` | R2 Bucket | User avatar storage |
| `FILES` | R2 Bucket | General file uploads |
| `AI` | Workers AI | AI model inference |

## Project Structure

```
src/
├── client/                 # React frontend
│   ├── components/         # UI components
│   ├── hooks/              # Custom hooks + TanStack Query
│   ├── pages/              # Route pages
│   ├── lib/                # Utilities (auth client, etc.)
│   └── main.tsx            # App entry point
├── server/                 # Hono backend
│   ├── index.ts            # Worker entry point
│   ├── routes/             # API routes
│   ├── middleware/          # Auth, CORS, etc.
│   └── db/                 # Drizzle schema + queries
└── shared/                 # Shared types between client/server
```

## Key Commands

| Command | Purpose |
|---------|---------|
| `pnpm dev` | Start local dev server |
| `pnpm build` | Production build |
| `pnpm deploy` | Deploy to Cloudflare |
| `pnpm db:migrate:local` | Apply migrations locally |
| `pnpm db:migrate:remote` | Apply migrations to production |
| `pnpm db:generate` | Generate migration from schema changes |
| `pnpm test` | Run tests |
| `pnpm lint` | Run linter |
