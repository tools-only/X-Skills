---
name: nextjs-tanstack-stack
description: Use when building Next.js applications with TanStack ecosystem (Table, Query, Form, Virtual), Zustand state management, or implementing data-intensive dashboards, virtualized tables, forms with validation, or performance optimization. Triggers on "Next.js App Router", "TanStack", "data table", "virtualization", "memoization", "server components", "client state", "Zustand".
---

# Next.js + TanStack Stack Implementation

Production patterns for Next.js App Router applications using TanStack ecosystem and Zustand, aligned with boring-over-clever philosophy.

## Stack Coverage

| Library | Version | Purpose |
|---------|---------|---------|
| Next.js | 14+ | App Router, RSC |
| TanStack Query | v5 | Server state |
| TanStack Table | v8 | Data grids |
| TanStack Form | v1 | Form handling |
| TanStack Virtual | v3 | List virtualization |
| Zustand | v4+ | Client state |

## When to Use This Skill

**Use for:**
- Data-intensive dashboards and tables
- Forms with complex validation and field arrays
- Large list rendering with virtualization
- Server/client component architecture
- Integrating server and client state

**Delegate to other skills:**
- Visual design decisions → `frontend-design`
- UX planning and wireframes → `ux-designer`
- Browser automation testing → `webapp-testing`

---

## Core Architecture Patterns

### 1. Server vs Client Component Boundary

```
Server Component (default):
├─ Data fetching at edge
├─ No interactivity needed
├─ SEO-critical content
└─ Secret/env access

Client Component ('use client'):
├─ Event handlers (onClick, onChange)
├─ Browser APIs (localStorage, window)
├─ useState, useEffect, useRef
├─ TanStack hooks (Query/Form/Table)
└─ Zustand stores
```

**Composition Pattern:**
```tsx
// app/markets/page.tsx (Server Component)
import { MarketTable } from './market-table';
import { getMarkets } from '@/lib/api';

export default async function MarketsPage() {
  const initialData = await getMarkets(); // Server fetch
  return <MarketTable initialData={initialData} />;
}

// app/markets/market-table.tsx (Client Component)
'use client';
import { useQuery } from '@tanstack/react-query';
import { getMarkets } from '@/lib/api';

interface Props {
  initialData: Awaited<ReturnType<typeof getMarkets>>;
}

export function MarketTable({ initialData }: Props) {
  const { data } = useQuery({
    queryKey: ['markets'],
    queryFn: getMarkets,
    initialData,
  });
  return <Table data={data} />;
}
```

### 2. Provider Architecture

```tsx
// lib/providers.tsx
'use client';

import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { useState, type ReactNode } from 'react';

export function Providers({ children }: { children: ReactNode }) {
  // useState prevents new client on every render
  const [queryClient] = useState(() => new QueryClient({
    defaultOptions: {
      queries: {
        staleTime: 5 * 1000,
        retry: (failureCount, error) => {
          if ((error as { status?: number }).status === 404) return false;
          return failureCount < 3;
        },
      },
      mutations: {
        onError: (error) => {
          console.error('[Mutation Error]', { error, timestamp: new Date().toISOString() });
        },
      },
    },
  }));

  return (
    <QueryClientProvider client={queryClient}>
      {children}
    </QueryClientProvider>
  );
}
```

### 3. TanStack Query Patterns

**Query Key Factory:**
```tsx
// lib/query-keys.ts
export const queryKeys = {
  markets: {
    all: ['markets'] as const,
    detail: (id: string) => ['markets', id] as const,
    filtered: (filters: MarketFilters) => ['markets', 'filtered', filters] as const,
  },
  trades: {
    byMarket: (marketId: string) => ['trades', marketId] as const,
  },
} as const;
```

**Query Hook with Select:**
```tsx
// hooks/use-markets.ts
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';
import { queryKeys } from '@/lib/query-keys';

export function useMarkets(filters?: MarketFilters) {
  return useQuery({
    queryKey: filters ? queryKeys.markets.filtered(filters) : queryKeys.markets.all,
    queryFn: () => fetchMarkets(filters),
    select: (data) => data.filter((m) => m.active), // Derive in query, not render
  });
}

export function useUpdateMarket() {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: updateMarket,
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.markets.all });
    },
    onError: (error, variables) => {
      console.error('[UpdateMarket Failed]', { error, marketId: variables.id });
    },
  });
}
```

### 4. TanStack Table + Virtualization

```tsx
'use client';

import { useVirtualizer } from '@tanstack/react-virtual';
import {
  useReactTable,
  getCoreRowModel,
  getSortedRowModel,
  type ColumnDef,
  type SortingState,
} from '@tanstack/react-table';
import { useRef, useState, useMemo } from 'react';

interface VirtualizedTableProps<T> {
  data: T[];
  columns: ColumnDef<T>[];
  estimateRowHeight?: number;
}

export function VirtualizedTable<T>({
  data,
  columns,
  estimateRowHeight = 48,
}: VirtualizedTableProps<T>) {
  const parentRef = useRef<HTMLDivElement>(null);
  const [sorting, setSorting] = useState<SortingState>([]);

  const table = useReactTable({
    data,
    columns,
    state: { sorting },
    onSortingChange: setSorting,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
  });

  const { rows } = table.getRowModel();

  const virtualizer = useVirtualizer({
    count: rows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => estimateRowHeight,
    overscan: 10,
  });

  const virtualRows = virtualizer.getVirtualItems();

  return (
    <div ref={parentRef} className="h-[600px] overflow-auto">
      <table className="w-full">
        <thead className="sticky top-0 bg-white z-10">
          {table.getHeaderGroups().map((headerGroup) => (
            <tr key={headerGroup.id}>
              {headerGroup.headers.map((header) => (
                <th
                  key={header.id}
                  onClick={header.column.getToggleSortingHandler()}
                  className="cursor-pointer select-none"
                >
                  {flexRender(header.column.columnDef.header, header.getContext())}
                </th>
              ))}
            </tr>
          ))}
        </thead>
        <tbody>
          <tr style={{ height: `${virtualizer.getTotalSize()}px` }}>
            <td colSpan={columns.length} className="relative">
              {virtualRows.map((virtualRow) => {
                const row = rows[virtualRow.index];
                return (
                  <div
                    key={row.id}
                    className="absolute w-full flex"
                    style={{
                      top: virtualRow.start,
                      height: `${virtualRow.size}px`,
                    }}
                  >
                    {row.getVisibleCells().map((cell) => (
                      <div key={cell.id} className="flex-1">
                        {flexRender(cell.column.columnDef.cell, cell.getContext())}
                      </div>
                    ))}
                  </div>
                );
              })}
            </td>
          </tr>
        </tbody>
      </table>
    </div>
  );
}
```

### 5. Zustand Store Patterns

**Slice Pattern with Actions Namespace:**
```tsx
// stores/market-filters.ts
import { create } from 'zustand';
import { devtools, persist } from 'zustand/middleware';

interface MarketFiltersState {
  search: string;
  category: string | null;
  sortBy: 'volume' | 'price' | 'name';
  actions: {
    setSearch: (search: string) => void;
    setCategory: (category: string | null) => void;
    setSortBy: (sortBy: 'volume' | 'price' | 'name') => void;
    reset: () => void;
  };
}

const initialState = {
  search: '',
  category: null,
  sortBy: 'volume' as const,
};

export const useMarketFilters = create<MarketFiltersState>()(
  devtools(
    persist(
      (set) => ({
        ...initialState,
        actions: {
          setSearch: (search) => set({ search }),
          setCategory: (category) => set({ category }),
          setSortBy: (sortBy) => set({ sortBy }),
          reset: () => set(initialState),
        },
      }),
      { name: 'market-filters' }
    ),
    { name: 'MarketFilters' }
  )
);

// Atomic selectors prevent unnecessary re-renders
export const selectSearch = (state: MarketFiltersState) => state.search;
export const selectCategory = (state: MarketFiltersState) => state.category;
export const selectSortBy = (state: MarketFiltersState) => state.sortBy;
export const selectActions = (state: MarketFiltersState) => state.actions;
```

**Usage with Atomic Selectors:**
```tsx
function SearchInput() {
  const search = useMarketFilters(selectSearch);
  const { setSearch } = useMarketFilters(selectActions);

  return <input value={search} onChange={(e) => setSearch(e.target.value)} />;
}
```

### 6. Error Boundary Pattern

```tsx
// components/error-boundary.tsx
'use client';

import { Component, type ReactNode, type ErrorInfo } from 'react';

interface Props {
  children: ReactNode;
  fallback: (error: Error, reset: () => void) => ReactNode;
  onError?: (error: Error, info: ErrorInfo) => void;
}

interface State {
  error: Error | null;
}

export class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  componentDidCatch(error: Error, info: ErrorInfo) {
    console.error('[ErrorBoundary]', {
      error: error.message,
      stack: error.stack,
      componentStack: info.componentStack,
    });
    this.props.onError?.(error, info);
  }

  reset = () => {
    this.setState({ error: null });
  };

  render() {
    if (this.state.error) {
      return this.props.fallback(this.state.error, this.reset);
    }
    return this.props.children;
  }
}
```

**With TanStack Query:**
```tsx
import { QueryErrorResetBoundary } from '@tanstack/react-query';

function MarketDashboard() {
  return (
    <QueryErrorResetBoundary>
      {({ reset }) => (
        <ErrorBoundary
          fallback={(error, localReset) => (
            <div>
              <p>Error: {error.message}</p>
              <button onClick={() => { reset(); localReset(); }}>
                Retry
              </button>
            </div>
          )}
        >
          <MarketTable />
        </ErrorBoundary>
      )}
    </QueryErrorResetBoundary>
  );
}
```

### 7. TanStack Form Pattern

```tsx
'use client';

import { useForm } from '@tanstack/react-form';
import { zodValidator } from '@tanstack/zod-form-adapter';
import { z } from 'zod';
import { useCreateMarket } from '@/hooks/use-markets';

const marketSchema = z.object({
  name: z.string().min(1, 'Required').max(100),
  description: z.string().min(10, 'At least 10 characters'),
  outcomes: z.array(z.object({
    label: z.string().min(1),
    probability: z.number().min(0).max(100),
  })).min(2, 'At least 2 outcomes'),
});

type MarketFormData = z.infer<typeof marketSchema>;

export function CreateMarketForm() {
  const createMutation = useCreateMarket();

  const form = useForm({
    defaultValues: {
      name: '',
      description: '',
      outcomes: [{ label: '', probability: 50 }, { label: '', probability: 50 }],
    } satisfies MarketFormData,
    validatorAdapter: zodValidator(),
    validators: {
      onChange: marketSchema,
    },
    onSubmit: async ({ value }) => {
      await createMutation.mutateAsync(value);
    },
  });

  return (
    <form
      onSubmit={(e) => {
        e.preventDefault();
        form.handleSubmit();
      }}
    >
      <form.Field name="name">
        {(field) => (
          <div>
            <input
              value={field.state.value}
              onChange={(e) => field.handleChange(e.target.value)}
              onBlur={field.handleBlur}
            />
            {field.state.meta.errors.length > 0 && (
              <span className="text-red-500">{field.state.meta.errors[0]}</span>
            )}
          </div>
        )}
      </form.Field>

      <button
        type="submit"
        disabled={form.state.isSubmitting || !form.state.canSubmit}
      >
        {form.state.isSubmitting ? 'Creating...' : 'Create Market'}
      </button>
    </form>
  );
}
```

---

## Performance Optimization Checklist

1. **Memoization**: `useMemo` for derived data, `useCallback` for handlers passed to children
2. **Selector Pattern**: Extract atomic selectors from Zustand stores
3. **Query Deduplication**: Consistent query keys prevent duplicate requests
4. **Virtualization**: TanStack Virtual for lists > 100 items
5. **Code Splitting**: `dynamic()` for route-level, `lazy()` for component-level
6. **Suspense Boundaries**: Wrap data-fetching at feature boundaries

---

## Best Practices

1. Keep Server Components for static content and initial data fetching
2. Colocate Query hooks with components that use them
3. Use factory functions for query keys (`queryKeys.markets.detail(id)`)
4. Separate Zustand actions into `actions` namespace
5. Apply error boundaries at feature boundaries, not globally
6. Prefer `select` in useQuery over deriving in render
7. Use TypeScript strict mode with explicit return types
8. Test hooks with React Testing Library + MSW

---

## Common Pitfalls

- **'use client' everywhere**: Breaks RSC benefits; only add when needed
- **New QueryClient on render**: Use `useState(() => new QueryClient())`
- **Full store subscription**: Use atomic selectors, not `useStore((s) => s)`
- **Mixed state boundaries**: Keep server/client state clearly separated
- **Silent mutation errors**: Always log with context in `onError`
- **Over-memoizing**: Don't memoize primitives or static arrays

---

## Resources

### Reference Files
- [`references/tanstack-query-patterns.md`](references/tanstack-query-patterns.md) — Query factories, caching, optimistic updates
- [`references/tanstack-table-patterns.md`](references/tanstack-table-patterns.md) — Column defs, sorting, filtering, virtualization
- [`references/tanstack-form-patterns.md`](references/tanstack-form-patterns.md) — Validation, field arrays, async submission
- [`references/zustand-patterns.md`](references/zustand-patterns.md) — Store slicing, middleware, persistence
- [`references/nextjs-app-router.md`](references/nextjs-app-router.md) — Layouts, streaming, parallel routes
- [`references/performance-patterns.md`](references/performance-patterns.md) — Profiling, bundle analysis, optimization

### Example Files
- [`examples/virtualized-data-table.tsx`](examples/virtualized-data-table.tsx) — Table + Virtual for 10k+ rows
- [`examples/query-with-error-boundary.tsx`](examples/query-with-error-boundary.tsx) — Query + ErrorBoundary composition
- [`examples/form-with-validation.tsx`](examples/form-with-validation.tsx) — TanStack Form + Zod
- [`examples/zustand-query-sync.tsx`](examples/zustand-query-sync.tsx) — Filter state synced with Query
