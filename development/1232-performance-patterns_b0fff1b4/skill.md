# Performance Patterns

React and Next.js performance optimization techniques.

## Memoization Guidelines

### When to Use useMemo

```tsx
// ✅ GOOD: Expensive computation
const sortedData = useMemo(() => {
  return data.sort((a, b) => b.value - a.value);
}, [data]);

// ✅ GOOD: Object passed to memoized child
const filters = useMemo(() => ({
  search: searchTerm,
  category: selectedCategory,
}), [searchTerm, selectedCategory]);

// ❌ BAD: Simple value
const doubled = useMemo(() => count * 2, [count]); // Just use: count * 2

// ❌ BAD: Stable reference
const emptyArray = useMemo(() => [], []); // Just use: const emptyArray = []
```

### When to Use useCallback

```tsx
// ✅ GOOD: Handler passed to memoized child
const handleClick = useCallback((id: string) => {
  setSelected(id);
}, []);

// ✅ GOOD: Handler in dependency array
const fetchData = useCallback(async () => {
  const result = await api.getData();
  setData(result);
}, []);

useEffect(() => {
  fetchData();
}, [fetchData]);

// ❌ BAD: Handler not passed to child or in deps
function Component() {
  // This doesn't need useCallback
  const handleClick = () => console.log('clicked');
  return <button onClick={handleClick}>Click</button>;
}
```

### React.memo

```tsx
// ✅ GOOD: Expensive render, stable props
const ExpensiveList = memo(function ExpensiveList({ items }: Props) {
  return (
    <ul>
      {items.map((item) => (
        <ComplexItem key={item.id} item={item} />
      ))}
    </ul>
  );
});

// With custom comparison
const MarketCard = memo(
  function MarketCard({ market }: { market: Market }) {
    return <div>{market.name}</div>;
  },
  (prevProps, nextProps) => {
    // Only re-render if these fields change
    return (
      prevProps.market.id === nextProps.market.id &&
      prevProps.market.price === nextProps.market.price
    );
  }
);

// ❌ BAD: Cheap render
const Label = memo(({ text }: { text: string }) => <span>{text}</span>);
// Just use: function Label({ text }) { return <span>{text}</span>; }
```

---

## Code Splitting

### Route-Level Splitting (Automatic)

```tsx
// Next.js App Router automatically code-splits each route
// app/markets/page.tsx → separate chunk
// app/trades/page.tsx → separate chunk
```

### Component-Level Lazy Loading

```tsx
import dynamic from 'next/dynamic';

// Lazy load heavy component
const Chart = dynamic(() => import('@/components/chart'), {
  loading: () => <ChartSkeleton />,
  ssr: false, // Client-only (no SSR)
});

// Lazy load with named export
const AdvancedTable = dynamic(
  () => import('@/components/tables').then((mod) => mod.AdvancedTable),
  { loading: () => <TableSkeleton /> }
);

// Usage
function Dashboard() {
  const [showChart, setShowChart] = useState(false);

  return (
    <div>
      <button onClick={() => setShowChart(true)}>Show Chart</button>
      {showChart && <Chart data={data} />}
    </div>
  );
}
```

### Conditional Imports

```tsx
// Load heavy library only when needed
async function exportToPDF() {
  const { jsPDF } = await import('jspdf');
  const doc = new jsPDF();
  // ...
}
```

---

## Virtualization

### TanStack Virtual for Lists

```tsx
import { useVirtualizer } from '@tanstack/react-virtual';

function VirtualList({ items }: { items: Item[] }) {
  const parentRef = useRef<HTMLDivElement>(null);

  const virtualizer = useVirtualizer({
    count: items.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 50, // Estimated row height
    overscan: 5, // Extra items to render
  });

  return (
    <div ref={parentRef} className="h-[400px] overflow-auto">
      <div
        style={{
          height: `${virtualizer.getTotalSize()}px`,
          position: 'relative',
        }}
      >
        {virtualizer.getVirtualItems().map((virtualItem) => {
          const item = items[virtualItem.index];
          return (
            <div
              key={item.id}
              style={{
                position: 'absolute',
                top: 0,
                left: 0,
                width: '100%',
                height: `${virtualItem.size}px`,
                transform: `translateY(${virtualItem.start}px)`,
              }}
            >
              {item.name}
            </div>
          );
        })}
      </div>
    </div>
  );
}
```

### Window Virtualization

```tsx
import { useWindowVirtualizer } from '@tanstack/react-virtual';

function WindowVirtualList({ items }: { items: Item[] }) {
  const virtualizer = useWindowVirtualizer({
    count: items.length,
    estimateSize: () => 100,
    overscan: 5,
  });

  return (
    <div style={{ height: `${virtualizer.getTotalSize()}px`, position: 'relative' }}>
      {virtualizer.getVirtualItems().map((virtualItem) => (
        <div
          key={virtualItem.key}
          style={{
            position: 'absolute',
            top: 0,
            left: 0,
            width: '100%',
            transform: `translateY(${virtualItem.start}px)`,
          }}
        >
          <ItemCard item={items[virtualItem.index]} />
        </div>
      ))}
    </div>
  );
}
```

---

## Image Optimization

```tsx
import Image from 'next/image';

// Optimized image with lazy loading
<Image
  src="/hero.jpg"
  alt="Hero image"
  width={1200}
  height={600}
  priority={false} // Lazy load (default)
  placeholder="blur"
  blurDataURL="data:image/..." // Low-res placeholder
/>

// Above the fold - load immediately
<Image
  src="/logo.png"
  alt="Logo"
  width={120}
  height={40}
  priority // Load immediately
/>

// Fill container
<div className="relative h-48 w-full">
  <Image
    src="/cover.jpg"
    alt="Cover"
    fill
    sizes="(max-width: 768px) 100vw, 50vw"
    className="object-cover"
  />
</div>
```

---

## Bundle Analysis

### Setup

```bash
# Install
npm install @next/bundle-analyzer

# next.config.js
const withBundleAnalyzer = require('@next/bundle-analyzer')({
  enabled: process.env.ANALYZE === 'true',
});

module.exports = withBundleAnalyzer({
  // config
});

# Run analysis
ANALYZE=true npm run build
```

### Common Optimizations

```tsx
// ❌ BAD: Import entire library
import _ from 'lodash';
const result = _.debounce(fn, 300);

// ✅ GOOD: Import specific function
import debounce from 'lodash/debounce';
const result = debounce(fn, 300);

// ❌ BAD: Import entire icon library
import { FaHome, FaUser } from 'react-icons/fa';

// ✅ GOOD: Import from specific path
import FaHome from 'react-icons/fa/FaHome';
import FaUser from 'react-icons/fa/FaUser';
```

---

## React DevTools Profiler

### Identifying Re-renders

```tsx
// Wrap components to track renders
function useRenderCount(componentName: string) {
  const renderCount = useRef(0);
  renderCount.current += 1;

  useEffect(() => {
    console.log(`${componentName} rendered ${renderCount.current} times`);
  });
}

// In component
function MarketCard({ market }: Props) {
  useRenderCount('MarketCard');
  // ...
}
```

### Profiler Component

```tsx
import { Profiler, type ProfilerOnRenderCallback } from 'react';

const onRenderCallback: ProfilerOnRenderCallback = (
  id,
  phase,
  actualDuration,
  baseDuration,
  startTime,
  commitTime
) => {
  console.log({
    id,
    phase,
    actualDuration: `${actualDuration.toFixed(2)}ms`,
    baseDuration: `${baseDuration.toFixed(2)}ms`,
  });
};

function App() {
  return (
    <Profiler id="MarketList" onRender={onRenderCallback}>
      <MarketList />
    </Profiler>
  );
}
```

---

## Deferred Loading

### Intersection Observer

```tsx
import { useInView } from 'react-intersection-observer';

function LazySection({ children }: { children: React.ReactNode }) {
  const { ref, inView } = useInView({
    triggerOnce: true,
    threshold: 0.1,
  });

  return (
    <div ref={ref}>
      {inView ? children : <Skeleton />}
    </div>
  );
}

// Usage
function Dashboard() {
  return (
    <div>
      <HeroSection /> {/* Loads immediately */}

      <LazySection>
        <HeavyChartSection /> {/* Loads when scrolled into view */}
      </LazySection>

      <LazySection>
        <DataTableSection />
      </LazySection>
    </div>
  );
}
```

### useDeferredValue

```tsx
import { useDeferredValue, useMemo } from 'react';

function SearchResults({ query }: { query: string }) {
  // Deferred value allows UI to stay responsive
  const deferredQuery = useDeferredValue(query);
  const isStale = query !== deferredQuery;

  const results = useMemo(() => {
    return filterItems(items, deferredQuery);
  }, [deferredQuery]);

  return (
    <div style={{ opacity: isStale ? 0.7 : 1 }}>
      {results.map((item) => (
        <ResultItem key={item.id} item={item} />
      ))}
    </div>
  );
}
```

---

## Web Vitals Monitoring

```tsx
// app/layout.tsx
import { SpeedInsights } from '@vercel/speed-insights/next';
import { Analytics } from '@vercel/analytics/react';

export default function RootLayout({ children }) {
  return (
    <html>
      <body>
        {children}
        <SpeedInsights />
        <Analytics />
      </body>
    </html>
  );
}
```

### Custom Reporting

```tsx
// lib/vitals.ts
import { onCLS, onFCP, onFID, onLCP, onTTFB } from 'web-vitals';

type MetricHandler = (metric: { name: string; value: number }) => void;

export function reportWebVitals(onPerfEntry: MetricHandler) {
  onCLS(onPerfEntry);
  onFCP(onPerfEntry);
  onFID(onPerfEntry);
  onLCP(onPerfEntry);
  onTTFB(onPerfEntry);
}

// Usage
reportWebVitals((metric) => {
  console.log(metric);
  // Send to analytics
  analytics.track('Web Vital', {
    name: metric.name,
    value: Math.round(metric.value),
  });
});
```

---

## Common Anti-Patterns

### Avoid Inline Objects in Props

```tsx
// ❌ BAD: New object on every render
<Table filters={{ search: query, sort: 'name' }} />

// ✅ GOOD: Memoized object
const filters = useMemo(() => ({ search: query, sort: 'name' }), [query]);
<Table filters={filters} />
```

### Avoid Inline Functions as Props

```tsx
// ❌ BAD: New function on every render
<Button onClick={() => handleClick(item.id)} />

// ✅ GOOD: Stable reference
const handleItemClick = useCallback((id: string) => {
  handleClick(id);
}, [handleClick]);

<Button onClick={() => handleItemClick(item.id)} />

// ✅ BETTER: Pass data via data attribute
<Button data-id={item.id} onClick={handleButtonClick} />

function handleButtonClick(e: React.MouseEvent<HTMLButtonElement>) {
  const id = e.currentTarget.dataset.id;
  handleClick(id);
}
```

### Avoid State in Parent When Child Owns It

```tsx
// ❌ BAD: Parent manages child state
function Parent() {
  const [isOpen, setIsOpen] = useState(false);
  return <Modal isOpen={isOpen} onClose={() => setIsOpen(false)} />;
}

// ✅ GOOD: Child manages own state
function Parent() {
  return <Modal trigger={<Button>Open</Button>} />;
}

function Modal({ trigger, children }) {
  const [isOpen, setIsOpen] = useState(false);
  return (
    <>
      <span onClick={() => setIsOpen(true)}>{trigger}</span>
      {isOpen && <ModalContent onClose={() => setIsOpen(false)}>{children}</ModalContent>}
    </>
  );
}
```
