# Radix UI Components Skill

**Für DresdenAIInsights** - 40+ Radix Primitives, Accessibility-First UI Components

## Wann aktiviert?

- Du arbeitest mit Radix UI (`@radix-ui/*`)
- Keywords: Dialog, Dropdown, Popover, Tooltip, Accordion
- Accessibility-Fragen (ARIA, Keyboard Navigation)

## Core Principles

### 1. Composition Pattern

Radix UI nutzt **Composition statt Props**:

```tsx
// ✅ Radix Way - Composable
import * as Dialog from '@radix-ui/react-dialog';

<Dialog.Root>
  <Dialog.Trigger>Open</Dialog.Trigger>
  <Dialog.Portal>
    <Dialog.Overlay />
    <Dialog.Content>
      <Dialog.Title>Titel</Dialog.Title>
      <Dialog.Description>Beschreibung</Dialog.Description>
      <Dialog.Close>Close</Dialog.Close>
    </Dialog.Content>
  </Dialog.Portal>
</Dialog.Root>

// ❌ Nicht Radix-Like
<Dialog title="Titel" description="..." onClose={...} />
```

### 2. Unstyled by Default

Radix liefert nur Logik & Accessibility - **du** stylest:

```tsx
<Dialog.Content className="fixed inset-0 bg-white rounded-lg shadow-xl">
  {/* Tailwind, CSS-in-JS, oder plain CSS */}
</Dialog.Content>
```

## Häufigste Components

### Dialog (Modal)

```tsx
import * as Dialog from '@radix-ui/react-dialog';

export function ContactModal() {
  return (
    <Dialog.Root>
      <Dialog.Trigger asChild>
        <button className="btn-primary">Kontakt</button>
      </Dialog.Trigger>

      <Dialog.Portal>
        <Dialog.Overlay className="fixed inset-0 bg-black/50" />

        <Dialog.Content className="fixed top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 bg-white p-6 rounded-lg">
          <Dialog.Title className="text-xl font-bold">
            Kontaktformular
          </Dialog.Title>

          <Dialog.Description className="text-gray-600 mt-2">
            Sende uns eine Nachricht
          </Dialog.Description>

          {/* Form content */}

          <Dialog.Close asChild>
            <button className="absolute top-2 right-2">×</button>
          </Dialog.Close>
        </Dialog.Content>
      </Dialog.Portal>
    </Dialog.Root>
  );
}
```

**Accessibility:** ✅ Auto-Fokus, Esc-Taste, Trap-Focus, ARIA-Labels

### Dropdown Menu

```tsx
import * as DropdownMenu from '@radix-ui/react-dropdown-menu';

<DropdownMenu.Root>
  <DropdownMenu.Trigger>
    Menu
  </DropdownMenu.Trigger>

  <DropdownMenu.Portal>
    <DropdownMenu.Content className="bg-white shadow-lg rounded">
      <DropdownMenu.Item className="px-4 py-2 hover:bg-gray-100">
        Profil
      </DropdownMenu.Item>
      <DropdownMenu.Item className="px-4 py-2">
        Einstellungen
      </DropdownMenu.Item>

      <DropdownMenu.Separator className="h-px bg-gray-200" />

      <DropdownMenu.Item className="px-4 py-2 text-red-600">
        Logout
      </DropdownMenu.Item>
    </DropdownMenu.Content>
  </DropdownMenu.Portal>
</DropdownMenu.Root>
```

**Keyboard:** ↑/↓ Navigation, Enter zum Auswählen, Esc zum Schließen

### Accordion

```tsx
import * as Accordion from '@radix-ui/react-accordion';

<Accordion.Root type="single" collapsible>
  <Accordion.Item value="item-1">
    <Accordion.Header>
      <Accordion.Trigger className="flex justify-between w-full">
        Was ist Manufacturing Inside Analyzer?
        <span>▼</span>
      </Accordion.Trigger>
    </Accordion.Header>

    <Accordion.Content className="overflow-hidden data-[state=open]:animate-slideDown">
      Ein KI-gestütztes Tool zur Analyse von Produktionsdaten...
    </Accordion.Content>
  </Accordion.Item>

  <Accordion.Item value="item-2">
    {/* ... */}
  </Accordion.Item>
</Accordion.Root>
```

### Tooltip

```tsx
import * as Tooltip from '@radix-ui/react-tooltip';

<Tooltip.Provider>
  <Tooltip.Root>
    <Tooltip.Trigger asChild>
      <button>Hover me</button>
    </Tooltip.Trigger>

    <Tooltip.Portal>
      <Tooltip.Content className="bg-black text-white px-2 py-1 rounded text-sm">
        Helpful information
        <Tooltip.Arrow className="fill-black" />
      </Tooltip.Content>
    </Tooltip.Portal>
  </Tooltip.Root>
</Tooltip.Provider>
```

## DresdenAIInsights Spezifika

### Theme Integration

```tsx
// components/ui/Dialog.tsx
import * as DialogPrimitive from '@radix-ui/react-dialog';
import { cn } from '@/lib/utils';

export const DialogContent = React.forwardRef<
  React.ElementRef<typeof DialogPrimitive.Content>,
  React.ComponentPropsWithoutRef<typeof DialogPrimitive.Content>
>(({ className, children, ...props }, ref) => (
  <DialogPrimitive.Portal>
    <DialogPrimitive.Overlay className="fixed inset-0 bg-black/80 backdrop-blur-sm" />
    <DialogPrimitive.Content
      ref={ref}
      className={cn(
        "fixed left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2",
        "w-full max-w-lg bg-gradient-to-br from-gray-900 to-gray-800",
        "border border-cyan-500/20 rounded-lg shadow-2xl",
        className
      )}
      {...props}
    >
      {children}
    </DialogPrimitive.Content>
  </DialogPrimitive.Portal>
));
```

### 40+ Components Wrapper

```
src/components/ui/
├── accordion.tsx
├── alert-dialog.tsx
├── dialog.tsx
├── dropdown-menu.tsx
├── popover.tsx
├── tooltip.tsx
└── ... (37 weitere)
```

**Pattern:** Wrapped Radix Primitive mit Theme-Styles

## Accessibility Best Practices

### ARIA Labels

```tsx
<Dialog.Content aria-describedby="dialog-desc">
  <Dialog.Title>Titel</Dialog.Title>
  <Dialog.Description id="dialog-desc">
    Beschreibung
  </Dialog.Description>
</Dialog.Content>
```

### asChild Pattern

```tsx
// ✅ Verhindert zusätzliche DOM-Nodes
<Dialog.Trigger asChild>
  <button>Open</button>
</Dialog.Trigger>

// ❌ Erzeugt <button><button>...</button></button>
<Dialog.Trigger>
  <button>Open</button>
</Dialog.Trigger>
```

### Fokus-Management

```tsx
// Auto-Fokus auf erstes Input
<Dialog.Content>
  <input autoFocus />
</Dialog.Content>

// Oder spezifisch
import { useRef, useEffect } from 'react';

const inputRef = useRef<HTMLInputElement>(null);
useEffect(() => inputRef.current?.focus(), []);

<input ref={inputRef} />
```

## Performance-Tipps

### Lazy Loading

```tsx
import { lazy, Suspense } from 'react';

const HeavyDialog = lazy(() => import('./HeavyDialog'));

<Suspense fallback={<Spinner />}>
  <HeavyDialog />
</Suspense>
```

### Portal Optimization

```tsx
// Ein Portal pro App, nicht pro Component
<Tooltip.Provider delayDuration={200}>
  <App />
</Tooltip.Provider>
```

## Troubleshooting

**Dialog schließt nicht?**
→ Prüfe `<Dialog.Root>` umschließt alles

**Styles nicht sichtbar?**
→ Radix ist unstyled - du musst `className` hinzufügen

**Keyboard Navigation funktioniert nicht?**
→ Prüfe `aria-` Attributes und `role`

## Ressourcen

- [Radix UI Docs](https://www.radix-ui.com/primitives)
- [Accessibility Guide](https://www.radix-ui.com/primitives/docs/overview/accessibility)
