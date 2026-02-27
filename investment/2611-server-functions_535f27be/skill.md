# Server Functions Reference

Patterns for `createServerFn` on Cloudflare Workers with D1, Drizzle, and better-auth.

## Core Pattern

```typescript
import { createServerFn } from "@tanstack/react-start";
import { getDb } from "../db";

export const getItems = createServerFn({ method: "GET" }).handler(async () => {
  const db = getDb(); // Per-request Drizzle client
  return db.select().from(items).all();
});
```

**CRITICAL**: Create the Drizzle client inside the handler, NOT at module level. The `env` import from `cloudflare:workers` is only available during request handling.

## Env Access on Cloudflare

```typescript
// ✅ CORRECT — inside handler or called from handler
import { env } from "cloudflare:workers";

export const getItems = createServerFn({ method: "GET" }).handler(async () => {
  const db = drizzle(env.DB, { schema }); // env available here
  return db.select().from(items).all();
});

// ❌ WRONG — module level
import { env } from "cloudflare:workers";
const db = drizzle(env.DB, { schema }); // env is undefined here!
```

## Input Validation

Use `.inputValidator()` with Zod or inline validator:

```typescript
import { z } from "zod";

export const createItem = createServerFn({ method: "POST" })
  .inputValidator(
    z.object({
      name: z.string().min(1),
      description: z.string().optional(),
    })
  )
  .handler(async ({ data }) => {
    const db = getDb();
    const id = crypto.randomUUID();
    await db.insert(items).values({ id, ...data, createdAt: Date.now() });
    return { id };
  });
```

## Protected Server Functions

Check auth session inside the handler. Throw redirect for unauthenticated users:

```typescript
import { redirect } from "@tanstack/react-router";
import { getAuth } from "../lib/auth.server";

async function requireSession(request?: Request) {
  const auth = getAuth();
  const session = await auth.api.getSession({ headers: request?.headers ?? new Headers() });
  if (!session) {
    throw redirect({ to: "/login" });
  }
  return session;
}

// Reusable auth check as a server function
export const getSessionFn = createServerFn({ method: "GET" }).handler(
  async ({ request }) => {
    const auth = getAuth();
    const session = await auth.api.getSession({
      headers: request.headers,
    });
    return session;
  }
);

export const getItems = createServerFn({ method: "GET" }).handler(
  async ({ request }) => {
    const session = await requireSession(request);
    const db = getDb();
    return db
      .select()
      .from(items)
      .where(eq(items.userId, session.user.id))
      .all();
  }
);
```

## Route Loader Pattern

Use server functions in route `loader` for data fetching:

```typescript
import { createFileRoute } from "@tanstack/react-router";
import { getItems } from "../../server/functions";

export const Route = createFileRoute("/_authed/items")({
  loader: () => getItems(),
  component: ItemsPage,
});

function ItemsPage() {
  const items = Route.useLoaderData();

  return (
    <div>
      {items.map((item) => (
        <div key={item.id}>{item.name}</div>
      ))}
    </div>
  );
}
```

## Route beforeLoad (Auth Guard)

Use `beforeLoad` on layout routes to protect groups of routes:

```typescript
import { createFileRoute, redirect } from "@tanstack/react-router";
import { getSessionFn } from "../server/functions";

export const Route = createFileRoute("/_authed")({
  beforeLoad: async () => {
    const session = await getSessionFn();
    if (!session) {
      throw redirect({ to: "/login" });
    }
    return { session };
  },
});
```

Child routes access the session via `useRouteContext()`:

```typescript
function DashboardPage() {
  const { session } = Route.useRouteContext();
  return <h1>Welcome, {session.user.name}!</h1>;
}
```

## Mutation + Router Invalidation

After mutations, invalidate the router to refetch loader data:

```typescript
function CreateItemForm() {
  const router = useRouter();

  const handleSubmit = async (data: NewItem) => {
    await createItem({ data });
    router.invalidate(); // Refetches all active loaders
    router.navigate({ to: "/items" });
  };

  return <form onSubmit={...}>...</form>;
}
```

## Type Safety

Use Drizzle's type inference for server function input/output:

```typescript
import { items } from "../db/schema";
import type { InferSelectModel, InferInsertModel } from "drizzle-orm";

type Item = InferSelectModel<typeof items>;
type NewItem = InferInsertModel<typeof items>;

export const getItems = createServerFn({ method: "GET" }).handler(
  async (): Promise<Item[]> => {
    const db = getDb();
    return db.select().from(items).all();
  }
);
```

## Error Handling

Server functions can throw errors that propagate to the client:

```typescript
export const getItem = createServerFn({ method: "GET" })
  .inputValidator(z.object({ id: z.string() }))
  .handler(async ({ data }) => {
    const db = getDb();
    const item = await db
      .select()
      .from(items)
      .where(eq(items.id, data.id))
      .get();

    if (!item) {
      throw new Error("Item not found");
    }

    return item;
  });
```

For auth failures, always use `throw redirect()` — not error responses.
