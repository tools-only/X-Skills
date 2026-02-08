# TanStack Query Patterns

Production patterns for server state management with TanStack Query v5.

## Query Key Factories

### Typed Factory Pattern

```tsx
// lib/query-keys.ts
export const queryKeys = {
  markets: {
    all: ['markets'] as const,
    lists: () => [...queryKeys.markets.all, 'list'] as const,
    list: (filters: MarketFilters) => [...queryKeys.markets.lists(), filters] as const,
    details: () => [...queryKeys.markets.all, 'detail'] as const,
    detail: (id: string) => [...queryKeys.markets.details(), id] as const,
  },
  trades: {
    all: ['trades'] as const,
    byMarket: (marketId: string) => [...queryKeys.trades.all, 'market', marketId] as const,
    byUser: (userId: string) => [...queryKeys.trades.all, 'user', userId] as const,
  },
  user: {
    all: ['user'] as const,
    current: () => [...queryKeys.user.all, 'current'] as const,
    settings: () => [...queryKeys.user.all, 'settings'] as const,
  },
} as const;
```

### Usage

```tsx
// Invalidate all market queries
queryClient.invalidateQueries({ queryKey: queryKeys.markets.all });

// Invalidate only market lists (not details)
queryClient.invalidateQueries({ queryKey: queryKeys.markets.lists() });

// Invalidate specific market
queryClient.invalidateQueries({ queryKey: queryKeys.markets.detail('abc123') });
```

---

## Query Hook Patterns

### Basic Query with Error Handling

```tsx
import { useQuery } from '@tanstack/react-query';
import { queryKeys } from '@/lib/query-keys';

export function useMarket(id: string) {
  return useQuery({
    queryKey: queryKeys.markets.detail(id),
    queryFn: async () => {
      const response = await fetch(`/api/markets/${id}`);
      if (!response.ok) {
        throw new Error(`Failed to fetch market: ${response.status}`);
      }
      return response.json() as Promise<Market>;
    },
    staleTime: 30 * 1000, // 30 seconds
    retry: (failureCount, error) => {
      // Don't retry on 404
      if ((error as { status?: number }).status === 404) return false;
      return failureCount < 3;
    },
  });
}
```

### Query with Select (Derived Data)

```tsx
export function useActiveMarkets() {
  return useQuery({
    queryKey: queryKeys.markets.all,
    queryFn: fetchAllMarkets,
    select: (data) => data.filter((m) => m.status === 'active'),
  });
}

// More efficient: transform once, not every render
export function useMarketNames() {
  return useQuery({
    queryKey: queryKeys.markets.all,
    queryFn: fetchAllMarkets,
    select: (data) => data.map((m) => ({ id: m.id, name: m.name })),
  });
}
```

### Dependent Queries

```tsx
export function useUserTrades(userId: string | undefined) {
  const userQuery = useUser(userId);

  return useQuery({
    queryKey: queryKeys.trades.byUser(userId!),
    queryFn: () => fetchUserTrades(userId!),
    enabled: !!userId && !!userQuery.data, // Only fetch when user is loaded
  });
}
```

### Parallel Queries

```tsx
import { useQueries } from '@tanstack/react-query';

export function useMarketDetails(ids: string[]) {
  return useQueries({
    queries: ids.map((id) => ({
      queryKey: queryKeys.markets.detail(id),
      queryFn: () => fetchMarket(id),
      staleTime: 60 * 1000,
    })),
  });
}
```

---

## Infinite Queries (Pagination)

### Cursor-Based

```tsx
import { useInfiniteQuery } from '@tanstack/react-query';

interface PageResult {
  items: Market[];
  nextCursor: string | null;
}

export function useInfiniteMarkets() {
  return useInfiniteQuery({
    queryKey: queryKeys.markets.lists(),
    queryFn: async ({ pageParam }) => {
      const response = await fetch(`/api/markets?cursor=${pageParam ?? ''}`);
      return response.json() as Promise<PageResult>;
    },
    initialPageParam: null as string | null,
    getNextPageParam: (lastPage) => lastPage.nextCursor,
    getPreviousPageParam: (firstPage) => null, // If bidirectional
  });
}

// Usage
function MarketList() {
  const { data, fetchNextPage, hasNextPage, isFetchingNextPage } = useInfiniteMarkets();

  const allMarkets = data?.pages.flatMap((page) => page.items) ?? [];

  return (
    <>
      {allMarkets.map((market) => <MarketCard key={market.id} market={market} />)}
      {hasNextPage && (
        <button onClick={() => fetchNextPage()} disabled={isFetchingNextPage}>
          {isFetchingNextPage ? 'Loading...' : 'Load More'}
        </button>
      )}
    </>
  );
}
```

### Offset-Based

```tsx
export function useInfiniteMarkets(pageSize = 20) {
  return useInfiniteQuery({
    queryKey: ['markets', 'infinite', pageSize],
    queryFn: async ({ pageParam }) => {
      const response = await fetch(`/api/markets?offset=${pageParam}&limit=${pageSize}`);
      return response.json() as Promise<{ items: Market[]; total: number }>;
    },
    initialPageParam: 0,
    getNextPageParam: (lastPage, allPages) => {
      const loadedCount = allPages.flatMap((p) => p.items).length;
      return loadedCount < lastPage.total ? loadedCount : undefined;
    },
  });
}
```

---

## Mutations

### Basic Mutation

```tsx
import { useMutation, useQueryClient } from '@tanstack/react-query';

export function useCreateMarket() {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: async (data: CreateMarketInput) => {
      const response = await fetch('/api/markets', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(data),
      });
      if (!response.ok) throw new Error('Failed to create market');
      return response.json() as Promise<Market>;
    },
    onSuccess: (newMarket) => {
      // Invalidate list queries
      queryClient.invalidateQueries({ queryKey: queryKeys.markets.lists() });

      // Optionally set the new item in cache
      queryClient.setQueryData(queryKeys.markets.detail(newMarket.id), newMarket);
    },
    onError: (error, variables) => {
      console.error('[CreateMarket Failed]', { error, input: variables });
    },
  });
}
```

### Optimistic Updates

```tsx
export function useToggleFavorite() {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: async ({ marketId, isFavorite }: { marketId: string; isFavorite: boolean }) => {
      const response = await fetch(`/api/markets/${marketId}/favorite`, {
        method: isFavorite ? 'DELETE' : 'POST',
      });
      if (!response.ok) throw new Error('Failed to toggle favorite');
      return response.json();
    },
    onMutate: async ({ marketId, isFavorite }) => {
      // Cancel outgoing refetches
      await queryClient.cancelQueries({ queryKey: queryKeys.markets.detail(marketId) });

      // Snapshot previous value
      const previousMarket = queryClient.getQueryData<Market>(
        queryKeys.markets.detail(marketId)
      );

      // Optimistically update
      if (previousMarket) {
        queryClient.setQueryData(queryKeys.markets.detail(marketId), {
          ...previousMarket,
          isFavorite: !isFavorite,
        });
      }

      return { previousMarket };
    },
    onError: (error, variables, context) => {
      // Rollback on error
      if (context?.previousMarket) {
        queryClient.setQueryData(
          queryKeys.markets.detail(variables.marketId),
          context.previousMarket
        );
      }
    },
    onSettled: (_, __, { marketId }) => {
      // Refetch to ensure sync
      queryClient.invalidateQueries({ queryKey: queryKeys.markets.detail(marketId) });
    },
  });
}
```

### Mutation with Retry

```tsx
export function useUpdateMarket() {
  return useMutation({
    mutationFn: updateMarket,
    retry: 3,
    retryDelay: (attemptIndex) => Math.min(1000 * 2 ** attemptIndex, 30000),
  });
}
```

---

## Caching Strategies

### Placeholder Data

```tsx
// Show stale data immediately while fetching fresh data
export function useMarket(id: string) {
  const queryClient = useQueryClient();

  return useQuery({
    queryKey: queryKeys.markets.detail(id),
    queryFn: () => fetchMarket(id),
    placeholderData: () => {
      // Try to find this market in the list cache
      const markets = queryClient.getQueryData<Market[]>(queryKeys.markets.all);
      return markets?.find((m) => m.id === id);
    },
  });
}
```

### Initial Data from Server Component

```tsx
// app/markets/[id]/page.tsx (Server Component)
export default async function MarketPage({ params }: { params: { id: string } }) {
  const market = await getMarket(params.id);
  return <MarketDetail initialData={market} />;
}

// components/market-detail.tsx (Client Component)
'use client';

export function MarketDetail({ initialData }: { initialData: Market }) {
  const { data } = useQuery({
    queryKey: queryKeys.markets.detail(initialData.id),
    queryFn: () => fetchMarket(initialData.id),
    initialData,
    staleTime: 60 * 1000, // Consider initial data fresh for 1 minute
  });

  return <div>{data.name}</div>;
}
```

### Prefetching

```tsx
// Prefetch on hover
function MarketCard({ market }: { market: Market }) {
  const queryClient = useQueryClient();

  const prefetchDetails = () => {
    queryClient.prefetchQuery({
      queryKey: queryKeys.markets.detail(market.id),
      queryFn: () => fetchMarket(market.id),
      staleTime: 60 * 1000,
    });
  };

  return (
    <Link
      href={`/markets/${market.id}`}
      onMouseEnter={prefetchDetails}
      onFocus={prefetchDetails}
    >
      {market.name}
    </Link>
  );
}
```

---

## Error Handling

### Global Error Handler

```tsx
const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      retry: (failureCount, error) => {
        // Don't retry on 4xx errors
        const status = (error as { status?: number }).status;
        if (status && status >= 400 && status < 500) return false;
        return failureCount < 3;
      },
    },
    mutations: {
      onError: (error, variables, context) => {
        // Global error logging
        console.error('[Mutation Error]', {
          error,
          variables,
          timestamp: new Date().toISOString(),
        });

        // Toast notification
        toast.error('Operation failed. Please try again.');
      },
    },
  },
});
```

### Query Error Boundary

```tsx
import { QueryErrorResetBoundary } from '@tanstack/react-query';
import { ErrorBoundary } from '@/components/error-boundary';

function MarketSection() {
  return (
    <QueryErrorResetBoundary>
      {({ reset }) => (
        <ErrorBoundary
          fallback={(error, localReset) => (
            <div className="p-4 bg-red-50 rounded">
              <p className="text-red-800">{error.message}</p>
              <button
                onClick={() => {
                  reset();
                  localReset();
                }}
                className="mt-2 px-4 py-2 bg-red-600 text-white rounded"
              >
                Retry
              </button>
            </div>
          )}
        >
          <MarketList />
        </ErrorBoundary>
      )}
    </QueryErrorResetBoundary>
  );
}
```

---

## Suspense Integration

### With Suspense Boundary

```tsx
import { useSuspenseQuery } from '@tanstack/react-query';
import { Suspense } from 'react';

function MarketDetails({ id }: { id: string }) {
  const { data } = useSuspenseQuery({
    queryKey: queryKeys.markets.detail(id),
    queryFn: () => fetchMarket(id),
  });

  return <div>{data.name}</div>;
}

// Usage
function MarketPage({ id }: { id: string }) {
  return (
    <Suspense fallback={<MarketSkeleton />}>
      <MarketDetails id={id} />
    </Suspense>
  );
}
```

---

## DevTools

```tsx
import { ReactQueryDevtools } from '@tanstack/react-query-devtools';

export function Providers({ children }: { children: React.ReactNode }) {
  return (
    <QueryClientProvider client={queryClient}>
      {children}
      <ReactQueryDevtools initialIsOpen={false} />
    </QueryClientProvider>
  );
}
```

---

## Testing Patterns

### Mock Query Client

```tsx
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render } from '@testing-library/react';

function createTestQueryClient() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false },
      mutations: { retry: false },
    },
  });
}

function renderWithClient(ui: React.ReactElement) {
  const testQueryClient = createTestQueryClient();
  return render(
    <QueryClientProvider client={testQueryClient}>
      {ui}
    </QueryClientProvider>
  );
}
```

### MSW Integration

```tsx
import { setupServer } from 'msw/node';
import { http, HttpResponse } from 'msw';

const server = setupServer(
  http.get('/api/markets/:id', ({ params }) => {
    return HttpResponse.json({ id: params.id, name: 'Test Market' });
  })
);

beforeAll(() => server.listen());
afterEach(() => server.resetHandlers());
afterAll(() => server.close());
```
