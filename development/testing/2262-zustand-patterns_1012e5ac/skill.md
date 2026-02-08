# Zustand Patterns

Client state management patterns with Zustand v4+, aligned with locality and composability principles.

## Store Design Principles

1. **Actions Namespace**: Group mutations under `actions` to separate state from behavior
2. **Atomic Selectors**: Create granular selectors to prevent unnecessary re-renders
3. **Slice Pattern**: Split large stores into feature-based slices
4. **Immutable Updates**: Never mutate state directly; Zustand uses immer-like syntax

---

## Basic Store Pattern

### With Actions Namespace

```tsx
// stores/counter.ts
import { create } from 'zustand';

interface CounterState {
  count: number;
  actions: {
    increment: () => void;
    decrement: () => void;
    reset: () => void;
    setCount: (count: number) => void;
  };
}

export const useCounter = create<CounterState>((set) => ({
  count: 0,
  actions: {
    increment: () => set((state) => ({ count: state.count + 1 })),
    decrement: () => set((state) => ({ count: state.count - 1 })),
    reset: () => set({ count: 0 }),
    setCount: (count) => set({ count }),
  },
}));

// Atomic selectors
export const selectCount = (state: CounterState) => state.count;
export const selectActions = (state: CounterState) => state.actions;
```

### Usage with Selectors

```tsx
function Counter() {
  // Only re-renders when count changes
  const count = useCounter(selectCount);
  // Actions are stable, never cause re-renders
  const { increment, decrement } = useCounter(selectActions);

  return (
    <div>
      <span>{count}</span>
      <button onClick={increment}>+</button>
      <button onClick={decrement}>-</button>
    </div>
  );
}
```

---

## Feature Store Pattern

### Market Filters Store

```tsx
// stores/market-filters.ts
import { create } from 'zustand';
import { devtools, persist } from 'zustand/middleware';

type SortField = 'volume' | 'price' | 'name' | 'change24h';
type SortDirection = 'asc' | 'desc';

interface MarketFiltersState {
  // State
  search: string;
  category: string | null;
  sortField: SortField;
  sortDirection: SortDirection;
  showFavoritesOnly: boolean;

  // Computed (derived in selectors, not stored)
  // Actions
  actions: {
    setSearch: (search: string) => void;
    setCategory: (category: string | null) => void;
    setSort: (field: SortField, direction?: SortDirection) => void;
    toggleFavoritesOnly: () => void;
    reset: () => void;
  };
}

const initialState = {
  search: '',
  category: null,
  sortField: 'volume' as SortField,
  sortDirection: 'desc' as SortDirection,
  showFavoritesOnly: false,
};

export const useMarketFilters = create<MarketFiltersState>()(
  devtools(
    persist(
      (set) => ({
        ...initialState,
        actions: {
          setSearch: (search) => set({ search }, false, 'setSearch'),
          setCategory: (category) => set({ category }, false, 'setCategory'),
          setSort: (field, direction) =>
            set(
              (state) => ({
                sortField: field,
                sortDirection: direction ?? (state.sortField === field
                  ? state.sortDirection === 'asc' ? 'desc' : 'asc'
                  : 'desc'),
              }),
              false,
              'setSort'
            ),
          toggleFavoritesOnly: () =>
            set(
              (state) => ({ showFavoritesOnly: !state.showFavoritesOnly }),
              false,
              'toggleFavoritesOnly'
            ),
          reset: () => set(initialState, false, 'reset'),
        },
      }),
      {
        name: 'market-filters',
        partialize: (state) => ({
          // Only persist these fields
          sortField: state.sortField,
          sortDirection: state.sortDirection,
          showFavoritesOnly: state.showFavoritesOnly,
        }),
      }
    ),
    { name: 'MarketFilters' }
  )
);

// Atomic selectors
export const selectSearch = (state: MarketFiltersState) => state.search;
export const selectCategory = (state: MarketFiltersState) => state.category;
export const selectSort = (state: MarketFiltersState) => ({
  field: state.sortField,
  direction: state.sortDirection,
});
export const selectShowFavoritesOnly = (state: MarketFiltersState) =>
  state.showFavoritesOnly;
export const selectActions = (state: MarketFiltersState) => state.actions;
```

---

## Slice Pattern (Large Stores)

### Combined Store with Slices

```tsx
// stores/slices/ui-slice.ts
import { type StateCreator } from 'zustand';

export interface UISlice {
  sidebarOpen: boolean;
  theme: 'light' | 'dark';
  uiActions: {
    toggleSidebar: () => void;
    setTheme: (theme: 'light' | 'dark') => void;
  };
}

export const createUISlice: StateCreator<UISlice> = (set) => ({
  sidebarOpen: true,
  theme: 'light',
  uiActions: {
    toggleSidebar: () => set((state) => ({ sidebarOpen: !state.sidebarOpen })),
    setTheme: (theme) => set({ theme }),
  },
});

// stores/slices/user-slice.ts
export interface UserSlice {
  user: User | null;
  isAuthenticated: boolean;
  userActions: {
    setUser: (user: User | null) => void;
    logout: () => void;
  };
}

export const createUserSlice: StateCreator<UserSlice> = (set) => ({
  user: null,
  isAuthenticated: false,
  userActions: {
    setUser: (user) => set({ user, isAuthenticated: !!user }),
    logout: () => set({ user: null, isAuthenticated: false }),
  },
});

// stores/app-store.ts
import { create } from 'zustand';
import { createUISlice, type UISlice } from './slices/ui-slice';
import { createUserSlice, type UserSlice } from './slices/user-slice';

type AppStore = UISlice & UserSlice;

export const useAppStore = create<AppStore>()((...args) => ({
  ...createUISlice(...args),
  ...createUserSlice(...args),
}));

// Selectors
export const selectSidebarOpen = (state: AppStore) => state.sidebarOpen;
export const selectTheme = (state: AppStore) => state.theme;
export const selectUser = (state: AppStore) => state.user;
export const selectIsAuthenticated = (state: AppStore) => state.isAuthenticated;
```

---

## Middleware Patterns

### Persist Middleware

```tsx
import { persist, createJSONStorage } from 'zustand/middleware';

export const useSettings = create<SettingsState>()(
  persist(
    (set) => ({
      // ...state and actions
    }),
    {
      name: 'app-settings',
      storage: createJSONStorage(() => localStorage),
      partialize: (state) => ({
        // Only persist specific fields
        theme: state.theme,
        language: state.language,
      }),
      version: 1,
      migrate: (persistedState, version) => {
        if (version === 0) {
          // Migration from v0 to v1
          return { ...persistedState, newField: 'default' };
        }
        return persistedState as SettingsState;
      },
    }
  )
);
```

### Devtools Middleware

```tsx
import { devtools } from 'zustand/middleware';

export const useStore = create<StoreState>()(
  devtools(
    (set) => ({
      count: 0,
      actions: {
        // Third argument is the action name for devtools
        increment: () => set((s) => ({ count: s.count + 1 }), false, 'increment'),
      },
    }),
    {
      name: 'MyStore',
      enabled: process.env.NODE_ENV === 'development',
    }
  )
);
```

### Immer Middleware

```tsx
import { immer } from 'zustand/middleware/immer';

interface State {
  items: Item[];
  actions: {
    addItem: (item: Item) => void;
    updateItem: (id: string, updates: Partial<Item>) => void;
    removeItem: (id: string) => void;
  };
}

export const useItems = create<State>()(
  immer((set) => ({
    items: [],
    actions: {
      addItem: (item) =>
        set((state) => {
          state.items.push(item); // Direct mutation with immer
        }),
      updateItem: (id, updates) =>
        set((state) => {
          const item = state.items.find((i) => i.id === id);
          if (item) Object.assign(item, updates);
        }),
      removeItem: (id) =>
        set((state) => {
          const index = state.items.findIndex((i) => i.id === id);
          if (index !== -1) state.items.splice(index, 1);
        }),
    },
  }))
);
```

---

## Integration with TanStack Query

### Syncing Filters with Query

```tsx
// hooks/use-filtered-markets.ts
import { useQuery } from '@tanstack/react-query';
import { useMarketFilters, selectSearch, selectCategory, selectSort } from '@/stores/market-filters';

export function useFilteredMarkets() {
  const search = useMarketFilters(selectSearch);
  const category = useMarketFilters(selectCategory);
  const sort = useMarketFilters(selectSort);

  // Query key includes all filter values
  return useQuery({
    queryKey: ['markets', 'filtered', { search, category, sort }],
    queryFn: () => fetchMarkets({ search, category, ...sort }),
    // Debounce search by keeping stale data while typing
    placeholderData: (previousData) => previousData,
  });
}
```

### Debounced Search

```tsx
import { useDeferredValue } from 'react';

function MarketSearch() {
  const search = useMarketFilters(selectSearch);
  const { setSearch } = useMarketFilters(selectActions);

  // Defer the value for query
  const deferredSearch = useDeferredValue(search);

  const { data, isLoading } = useQuery({
    queryKey: ['markets', 'search', deferredSearch],
    queryFn: () => searchMarkets(deferredSearch),
    enabled: deferredSearch.length >= 2,
  });

  return (
    <input
      value={search}
      onChange={(e) => setSearch(e.target.value)}
      placeholder="Search markets..."
    />
  );
}
```

---

## URL Sync Pattern

```tsx
// hooks/use-url-sync.ts
import { useEffect } from 'react';
import { useRouter, useSearchParams } from 'next/navigation';
import { useMarketFilters, selectCategory, selectActions } from '@/stores/market-filters';

export function useUrlSync() {
  const router = useRouter();
  const searchParams = useSearchParams();
  const category = useMarketFilters(selectCategory);
  const { setCategory } = useMarketFilters(selectActions);

  // Sync URL → Store on mount
  useEffect(() => {
    const urlCategory = searchParams.get('category');
    if (urlCategory && urlCategory !== category) {
      setCategory(urlCategory);
    }
  }, []);

  // Sync Store → URL on change
  useEffect(() => {
    const params = new URLSearchParams(searchParams);
    if (category) {
      params.set('category', category);
    } else {
      params.delete('category');
    }
    router.replace(`?${params.toString()}`, { scroll: false });
  }, [category]);
}
```

---

## Testing Patterns

### Reset Store Between Tests

```tsx
// stores/test-utils.ts
import { useMarketFilters } from './market-filters';

export function resetAllStores() {
  useMarketFilters.setState(useMarketFilters.getInitialState());
}

// In tests
beforeEach(() => {
  resetAllStores();
});
```

### Mock Store

```tsx
import { create } from 'zustand';

function createMockStore(initialState: Partial<MarketFiltersState>) {
  return create<MarketFiltersState>()((set) => ({
    search: '',
    category: null,
    ...initialState,
    actions: {
      setSearch: (search) => set({ search }),
      setCategory: (category) => set({ category }),
      reset: () => set({ search: '', category: null }),
    },
  }));
}
```

---

## Anti-Patterns to Avoid

### DON'T: Subscribe to Entire Store

```tsx
// BAD: Re-renders on ANY state change
function Component() {
  const state = useStore(); // Subscribes to everything
  return <div>{state.count}</div>;
}

// GOOD: Only re-renders when count changes
function Component() {
  const count = useStore(selectCount);
  return <div>{count}</div>;
}
```

### DON'T: Create Selectors Inline

```tsx
// BAD: New selector on every render
function Component() {
  const count = useStore((state) => state.count);
}

// GOOD: Stable selector reference
const selectCount = (state: State) => state.count;
function Component() {
  const count = useStore(selectCount);
}
```

### DON'T: Store Derived State

```tsx
// BAD: Storing computed values
interface State {
  items: Item[];
  totalCount: number; // Derived from items.length
}

// GOOD: Derive in selector
interface State {
  items: Item[];
}
const selectTotalCount = (state: State) => state.items.length;
```
