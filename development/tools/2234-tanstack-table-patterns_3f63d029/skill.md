# TanStack Table Patterns

Comprehensive patterns for building performant data tables with TanStack Table v8.

## Column Definition Patterns

### Basic Typed Columns

```tsx
import { createColumnHelper, type ColumnDef } from '@tanstack/react-table';

interface Market {
  id: string;
  name: string;
  volume: number;
  price: number;
  change24h: number;
  updatedAt: Date;
}

const columnHelper = createColumnHelper<Market>();

export const columns: ColumnDef<Market>[] = [
  columnHelper.accessor('name', {
    header: 'Market',
    cell: (info) => info.getValue(),
  }),
  columnHelper.accessor('volume', {
    header: 'Volume',
    cell: (info) => formatCurrency(info.getValue()),
    sortingFn: 'basic',
  }),
  columnHelper.accessor('price', {
    header: 'Price',
    cell: (info) => `$${info.getValue().toFixed(2)}`,
  }),
  columnHelper.accessor('change24h', {
    header: '24h Change',
    cell: (info) => {
      const value = info.getValue();
      const color = value >= 0 ? 'text-green-500' : 'text-red-500';
      return <span className={color}>{value.toFixed(2)}%</span>;
    },
  }),
  columnHelper.accessor('updatedAt', {
    header: 'Updated',
    cell: (info) => formatRelativeTime(info.getValue()),
    sortingFn: 'datetime',
  }),
];
```

### Computed Columns

```tsx
columnHelper.display({
  id: 'actions',
  header: () => <span className="sr-only">Actions</span>,
  cell: ({ row }) => (
    <DropdownMenu>
      <DropdownMenuItem onClick={() => handleEdit(row.original)}>
        Edit
      </DropdownMenuItem>
      <DropdownMenuItem onClick={() => handleDelete(row.original.id)}>
        Delete
      </DropdownMenuItem>
    </DropdownMenu>
  ),
}),
```

### Column Grouping

```tsx
export const groupedColumns: ColumnDef<Market>[] = [
  columnHelper.group({
    id: 'identification',
    header: 'Identification',
    columns: [
      columnHelper.accessor('id', { header: 'ID' }),
      columnHelper.accessor('name', { header: 'Name' }),
    ],
  }),
  columnHelper.group({
    id: 'metrics',
    header: 'Metrics',
    columns: [
      columnHelper.accessor('volume', { header: 'Volume' }),
      columnHelper.accessor('price', { header: 'Price' }),
    ],
  }),
];
```

---

## Sorting Patterns

### Controlled Sorting State

```tsx
import { useState } from 'react';
import { useReactTable, getSortedRowModel, type SortingState } from '@tanstack/react-table';

function SortableTable({ data, columns }) {
  const [sorting, setSorting] = useState<SortingState>([
    { id: 'volume', desc: true }, // Default sort
  ]);

  const table = useReactTable({
    data,
    columns,
    state: { sorting },
    onSortingChange: setSorting,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
  });

  return (
    <table>
      <thead>
        {table.getHeaderGroups().map((headerGroup) => (
          <tr key={headerGroup.id}>
            {headerGroup.headers.map((header) => (
              <th
                key={header.id}
                onClick={header.column.getToggleSortingHandler()}
                className="cursor-pointer select-none"
              >
                {flexRender(header.column.columnDef.header, header.getContext())}
                {{
                  asc: ' ↑',
                  desc: ' ↓',
                }[header.column.getIsSorted() as string] ?? null}
              </th>
            ))}
          </tr>
        ))}
      </thead>
      {/* body */}
    </table>
  );
}
```

### Multi-Column Sort

```tsx
const table = useReactTable({
  // ...
  enableMultiSort: true,
  maxMultiSortColCount: 3,
});
```

### Custom Sorting Function

```tsx
columnHelper.accessor('status', {
  header: 'Status',
  sortingFn: (rowA, rowB, columnId) => {
    const order = ['active', 'pending', 'closed'];
    const a = order.indexOf(rowA.getValue(columnId));
    const b = order.indexOf(rowB.getValue(columnId));
    return a - b;
  },
}),
```

---

## Filtering Patterns

### Global Filter

```tsx
import { useState } from 'react';
import { getFilteredRowModel, type FilterFn } from '@tanstack/react-table';

const fuzzyFilter: FilterFn<Market> = (row, columnId, filterValue) => {
  const value = row.getValue(columnId);
  return String(value).toLowerCase().includes(String(filterValue).toLowerCase());
};

function FilterableTable({ data, columns }) {
  const [globalFilter, setGlobalFilter] = useState('');

  const table = useReactTable({
    data,
    columns,
    state: { globalFilter },
    onGlobalFilterChange: setGlobalFilter,
    globalFilterFn: fuzzyFilter,
    getCoreRowModel: getCoreRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  });

  return (
    <>
      <input
        value={globalFilter}
        onChange={(e) => setGlobalFilter(e.target.value)}
        placeholder="Search all columns..."
      />
      <table>{/* ... */}</table>
    </>
  );
}
```

### Column Filters

```tsx
import { useState } from 'react';
import { type ColumnFiltersState } from '@tanstack/react-table';

function ColumnFilterTable({ data, columns }) {
  const [columnFilters, setColumnFilters] = useState<ColumnFiltersState>([]);

  const table = useReactTable({
    data,
    columns,
    state: { columnFilters },
    onColumnFiltersChange: setColumnFilters,
    getCoreRowModel: getCoreRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  });

  // Per-column filter input
  const nameColumn = table.getColumn('name');
  return (
    <input
      value={(nameColumn?.getFilterValue() as string) ?? ''}
      onChange={(e) => nameColumn?.setFilterValue(e.target.value)}
      placeholder="Filter by name..."
    />
  );
}
```

### Faceted Filters (Unique Values)

```tsx
import { getFacetedRowModel, getFacetedUniqueValues } from '@tanstack/react-table';

const table = useReactTable({
  // ...
  getFacetedRowModel: getFacetedRowModel(),
  getFacetedUniqueValues: getFacetedUniqueValues(),
});

// Get unique values for a column
const categoryColumn = table.getColumn('category');
const uniqueCategories = categoryColumn?.getFacetedUniqueValues();
// Map<categoryValue, count>
```

---

## Row Selection

### Checkbox Selection

```tsx
import { useState } from 'react';
import { type RowSelectionState } from '@tanstack/react-table';

function SelectableTable({ data, columns }) {
  const [rowSelection, setRowSelection] = useState<RowSelectionState>({});

  const selectColumn: ColumnDef<Market> = {
    id: 'select',
    header: ({ table }) => (
      <input
        type="checkbox"
        checked={table.getIsAllRowsSelected()}
        indeterminate={table.getIsSomeRowsSelected()}
        onChange={table.getToggleAllRowsSelectedHandler()}
      />
    ),
    cell: ({ row }) => (
      <input
        type="checkbox"
        checked={row.getIsSelected()}
        disabled={!row.getCanSelect()}
        onChange={row.getToggleSelectedHandler()}
      />
    ),
  };

  const table = useReactTable({
    data,
    columns: [selectColumn, ...columns],
    state: { rowSelection },
    onRowSelectionChange: setRowSelection,
    getCoreRowModel: getCoreRowModel(),
    enableRowSelection: true,
  });

  const selectedRows = table.getSelectedRowModel().rows;
  // selectedRows.map(row => row.original)
}
```

### Conditional Selection

```tsx
const table = useReactTable({
  // ...
  enableRowSelection: (row) => row.original.status !== 'locked',
});
```

---

## Virtualization Integration

### Basic Virtual Rows

```tsx
import { useVirtualizer } from '@tanstack/react-virtual';
import { useRef } from 'react';

function VirtualTable({ data, columns }) {
  const parentRef = useRef<HTMLDivElement>(null);

  const table = useReactTable({
    data,
    columns,
    getCoreRowModel: getCoreRowModel(),
  });

  const { rows } = table.getRowModel();

  const virtualizer = useVirtualizer({
    count: rows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 48, // Row height estimate
    overscan: 10, // Render 10 extra rows above/below viewport
  });

  return (
    <div ref={parentRef} className="h-[600px] overflow-auto">
      <div style={{ height: `${virtualizer.getTotalSize()}px`, position: 'relative' }}>
        {virtualizer.getVirtualItems().map((virtualRow) => {
          const row = rows[virtualRow.index];
          return (
            <div
              key={row.id}
              style={{
                position: 'absolute',
                top: 0,
                left: 0,
                width: '100%',
                height: `${virtualRow.size}px`,
                transform: `translateY(${virtualRow.start}px)`,
              }}
            >
              {row.getVisibleCells().map((cell) => (
                <div key={cell.id}>
                  {flexRender(cell.column.columnDef.cell, cell.getContext())}
                </div>
              ))}
            </div>
          );
        })}
      </div>
    </div>
  );
}
```

### Dynamic Row Heights

```tsx
const virtualizer = useVirtualizer({
  count: rows.length,
  getScrollElement: () => parentRef.current,
  estimateSize: () => 48,
  measureElement: (element) => element.getBoundingClientRect().height,
});
```

### Virtual Columns (Horizontal)

```tsx
const columnVirtualizer = useVirtualizer({
  horizontal: true,
  count: columns.length,
  getScrollElement: () => parentRef.current,
  estimateSize: (index) => columns[index].size ?? 150,
  overscan: 5,
});
```

---

## Pagination

### Client-Side Pagination

```tsx
import { useState } from 'react';
import { getPaginationRowModel, type PaginationState } from '@tanstack/react-table';

function PaginatedTable({ data, columns }) {
  const [pagination, setPagination] = useState<PaginationState>({
    pageIndex: 0,
    pageSize: 20,
  });

  const table = useReactTable({
    data,
    columns,
    state: { pagination },
    onPaginationChange: setPagination,
    getCoreRowModel: getCoreRowModel(),
    getPaginationRowModel: getPaginationRowModel(),
  });

  return (
    <>
      <table>{/* ... */}</table>
      <div className="flex items-center gap-2">
        <button
          onClick={() => table.previousPage()}
          disabled={!table.getCanPreviousPage()}
        >
          Previous
        </button>
        <span>
          Page {table.getState().pagination.pageIndex + 1} of {table.getPageCount()}
        </span>
        <button
          onClick={() => table.nextPage()}
          disabled={!table.getCanNextPage()}
        >
          Next
        </button>
        <select
          value={pagination.pageSize}
          onChange={(e) => table.setPageSize(Number(e.target.value))}
        >
          {[10, 20, 50, 100].map((size) => (
            <option key={size} value={size}>
              Show {size}
            </option>
          ))}
        </select>
      </div>
    </>
  );
}
```

### Server-Side Pagination

```tsx
function ServerPaginatedTable() {
  const [pagination, setPagination] = useState({ pageIndex: 0, pageSize: 20 });

  const { data, isLoading } = useQuery({
    queryKey: ['markets', pagination],
    queryFn: () => fetchMarkets({
      page: pagination.pageIndex,
      limit: pagination.pageSize,
    }),
  });

  const table = useReactTable({
    data: data?.items ?? [],
    columns,
    pageCount: data?.totalPages ?? -1,
    state: { pagination },
    onPaginationChange: setPagination,
    getCoreRowModel: getCoreRowModel(),
    manualPagination: true, // Disable client-side pagination
  });
}
```

---

## Column Visibility

```tsx
import { useState } from 'react';
import { type VisibilityState } from '@tanstack/react-table';

function TableWithColumnToggle({ data, columns }) {
  const [columnVisibility, setColumnVisibility] = useState<VisibilityState>({
    id: false, // Hide ID column by default
  });

  const table = useReactTable({
    data,
    columns,
    state: { columnVisibility },
    onColumnVisibilityChange: setColumnVisibility,
    getCoreRowModel: getCoreRowModel(),
  });

  return (
    <>
      <div className="flex gap-2">
        {table.getAllLeafColumns().map((column) => (
          <label key={column.id} className="flex items-center gap-1">
            <input
              type="checkbox"
              checked={column.getIsVisible()}
              onChange={column.getToggleVisibilityHandler()}
            />
            {column.id}
          </label>
        ))}
      </div>
      <table>{/* ... */}</table>
    </>
  );
}
```

---

## Column Resizing

```tsx
const table = useReactTable({
  data,
  columns,
  columnResizeMode: 'onChange', // or 'onEnd'
  getCoreRowModel: getCoreRowModel(),
});

// In header cell
<th
  style={{ width: header.getSize() }}
>
  {flexRender(header.column.columnDef.header, header.getContext())}
  <div
    onMouseDown={header.getResizeHandler()}
    onTouchStart={header.getResizeHandler()}
    className={`resizer ${header.column.getIsResizing() ? 'isResizing' : ''}`}
  />
</th>
```

---

## Performance Tips

1. **Memoize columns**: Define outside component or wrap in `useMemo`
2. **Virtualize large datasets**: Use TanStack Virtual for 100+ rows
3. **Avoid inline functions in cells**: Use `useCallback` or define outside
4. **Use `manualSorting`/`manualFiltering`** for server-side operations
5. **Batch state updates**: Combine sorting + filtering changes
