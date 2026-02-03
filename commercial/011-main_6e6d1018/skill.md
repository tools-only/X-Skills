# Database & ORM Skill

**Drizzle ORM + Vercel Postgres + Vercel KV (Redis)** für ManufacturingInsideAnalyzer

## Wann aktiviert?

- Keywords: database, Drizzle, ORM, schema, migration, Vercel Postgres, Redis, KV
- Du arbeitest mit DB-Dateien: `schema.ts`, `migrations/`, `drizzle.config.ts`

## ⚠️ Wichtig: Vercel Postgres ≠ Vercel KV

- **Vercel Postgres**: Relational Database (SQL) → Drizzle ORM
- **Vercel KV**: Redis Key-Value Store → `@vercel/kv` SDK (kein ORM!)

## 1. Drizzle ORM Setup (Vercel Postgres)

### Installation

```bash
npm install drizzle-orm @vercel/postgres
npm install -D drizzle-kit
```

### Konfiguration

**drizzle.config.ts:**
```ts
import { defineConfig } from 'drizzle-kit';

export default defineConfig({
  schema: './src/db/schema.ts',
  out: './drizzle',
  dialect: 'postgresql',
  dbCredentials: {
    url: process.env.POSTGRES_URL!
  }
});
```

**.env:**
```
POSTGRES_URL=postgres://...
```

### Schema Definition

```ts
// src/db/schema.ts
import { pgTable, serial, text, timestamp, integer, boolean } from 'drizzle-orm/pg-core';

export const users = pgTable('users', {
  id: serial('id').primaryKey(),
  email: text('email').notNull().unique(),
  name: text('name').notNull(),
  createdAt: timestamp('created_at').defaultNow()
});

export const analyses = pgTable('analyses', {
  id: serial('id').primaryKey(),
  userId: integer('user_id').references(() => users.id),
  fileName: text('file_name').notNull(),
  result: text('result'), // JSON als String
  createdAt: timestamp('created_at').defaultNow(),
  completed: boolean('completed').default(false)
});

// Types exportieren
export type User = typeof users.$inferSelect;
export type NewUser = typeof users.$inferInsert;
```

### Database Client

```ts
// src/db/index.ts
import { drizzle } from 'drizzle-orm/vercel-postgres';
import { sql } from '@vercel/postgres';
import * as schema from './schema';

export const db = drizzle(sql, { schema });
```

## 2. Queries mit Drizzle

### Select (Read)

```ts
import { db } from './db';
import { users, analyses } from './db/schema';
import { eq, and, desc } from 'drizzle-orm';

// Alle Users
const allUsers = await db.select().from(users);

// User by ID
const user = await db.select()
  .from(users)
  .where(eq(users.id, 123));

// Mit Relations (Join)
const userWithAnalyses = await db.select()
  .from(users)
  .leftJoin(analyses, eq(users.id, analyses.userId))
  .where(eq(users.id, 123));

// Komplexe Query
const recentAnalyses = await db.select()
  .from(analyses)
  .where(and(
    eq(analyses.completed, true),
    eq(analyses.userId, userId)
  ))
  .orderBy(desc(analyses.createdAt))
  .limit(10);
```

### Insert (Create)

```ts
// Single Insert
const newUser = await db.insert(users).values({
  email: 'test@example.com',
  name: 'Test User'
}).returning();

// Bulk Insert
await db.insert(analyses).values([
  { userId: 1, fileName: 'file1.csv' },
  { userId: 1, fileName: 'file2.csv' }
]);
```

### Update

```ts
await db.update(analyses)
  .set({ completed: true, result: JSON.stringify(data) })
  .where(eq(analyses.id, analysisId));
```

### Delete

```ts
await db.delete(analyses)
  .where(eq(analyses.id, analysisId));
```

## 3. Migrations

### Create Migration

```bash
# Nach Schema-Änderungen
npx drizzle-kit generate
```

Erzeugt SQL-File in `drizzle/`:
```sql
CREATE TABLE IF NOT EXISTS "users" (
  "id" serial PRIMARY KEY NOT NULL,
  "email" text NOT NULL,
  "name" text NOT NULL
);
```

### Apply Migration

```bash
# Lokal
npx drizzle-kit push

# Production (Vercel)
# In build command: npm run build && npx drizzle-kit push
```

### Seeding

```ts
// scripts/seed.ts
import { db } from '../src/db';
import { users } from '../src/db/schema';

async function seed() {
  await db.insert(users).values([
    { email: 'admin@example.com', name: 'Admin' },
    { email: 'user@example.com', name: 'User' }
  ]);
  console.log('Seeded!');
}

seed();
```

```bash
npx tsx scripts/seed.ts
```

## 4. Vercel KV (Redis)

**Für:** Caching, Sessions, Rate Limiting

### Setup

```bash
npm install @vercel/kv
```

### Usage

```ts
import { kv } from '@vercel/kv';

// Set
await kv.set('user:123', { name: 'John', email: 'john@example.com' });

// Get
const user = await kv.get<User>('user:123');

// Expire (TTL)
await kv.set('session:abc', sessionData, { ex: 3600 }); // 1 hour

// Increment (Rate Limiting)
const count = await kv.incr('rate:user@example.com');

// Delete
await kv.del('user:123');

// Hash (für komplexe Objekte)
await kv.hset('user:123', { name: 'John', age: 30 });
const age = await kv.hget('user:123', 'age');
```

## 5. ManufacturingInsideAnalyzer Patterns

### Rate Limiting (Vercel KV)

```ts
async function checkAnalysisLimit(email: string): Promise<boolean> {
  const key = `analysis-limit:${email}`;
  const current = await kv.get<number>(key) || 0;

  if (current >= 5) return false; // Max 5 per 24h

  await kv.set(key, current + 1, { ex: 86400 }); // 24h TTL
  return true;
}
```

### Session Management

```ts
import { kv } from '@vercel/kv';

export async function createSession(userId: number): Promise<string> {
  const sessionId = crypto.randomUUID();

  await kv.set(
    `session:${sessionId}`,
    { userId, createdAt: Date.now() },
    { ex: 3600 } // 1 hour
  );

  return sessionId;
}

export async function getSession(sessionId: string) {
  return await kv.get(`session:${sessionId}`);
}
```

### Caching Analysis Results

```ts
// Cache für teure AI-Analysen
async function analyzeWithCache(fileHash: string, data: any) {
  const cached = await kv.get(`analysis:${fileHash}`);
  if (cached) return cached;

  const result = await performExpensiveAnalysis(data);

  await kv.set(`analysis:${fileHash}`, result, { ex: 86400 });
  return result;
}
```

## 6. Frankfurt Region (DSGVO!)

**Vercel Postgres:**
```json
// vercel.json
{
  "regions": ["fra1"]
}
```

**Vercel KV:**
Automatisch in Frankfurt wenn Vercel Project in `fra1`.

## 7. Troubleshooting

**Connection Error?**
→ Prüfe `POSTGRES_URL` in `.env` und Vercel Dashboard

**Migration fehlgeschlagen?**
→ `npx drizzle-kit drop` (⚠️ Löscht Daten!) dann `push`

**KV not working?**
→ Erstelle KV Store im Vercel Dashboard, verlinke mit Project

## Ressourcen

- [Drizzle ORM Docs](https://orm.drizzle.team/)
- [Vercel Postgres](https://vercel.com/docs/storage/vercel-postgres)
- [Vercel KV](https://vercel.com/docs/storage/vercel-kv)
