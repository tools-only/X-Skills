# React + Vite Development Skill

Dieser Skill hilft bei der Entwicklung von React-Anwendungen mit Vite Build-System. Er ist optimiert für deine Projekte: **digitalTwin**, **DresdenAIInsights**, und **ManufacturingInsideAnalyzer**.

## Wann wird dieser Skill aktiviert?

- Du arbeitest mit `.tsx`-Dateien in `components/` oder `src/`
- Du erwähnst Keywords wie: component, React, Vite, useState, useEffect
- Du möchtest UI-Komponenten erstellen oder debuggen

## Kern-Prinzipien

### 1. Komponenten-Struktur

**Single Responsibility**: Jede Komponente hat EINEN klaren Zweck.

**Beispiel - GUT:**
```tsx
// components/UserProfile.tsx
export function UserProfile({ userId }: { userId: string }) {
  const { user, loading } = useUser(userId);

  if (loading) return <Spinner />;
  return <div>{user.name}</div>;
}
```

**Beispiel - SCHLECHT:**
```tsx
// Macht zu viel: Fetching + Rendering + Validation
export function UserProfile() {
  const [user, setUser] = useState();
  useEffect(() => { fetch(...) }, []);
  const validate = () => { ... };
  return <div>...</div>;
}
```

### 2. Props über Drilling

**Problem**: Props durch 5+ Komponenten durchreichen
**Lösung**: Context API oder Custom Hooks

**Beispiel:**
```tsx
// hooks/useTheme.ts
export function useTheme() {
  return useContext(ThemeContext);
}

// components/Button.tsx
export function Button() {
  const { theme } = useTheme(); // Kein Props-Drilling!
  return <button className={theme}>Click</button>;
}
```

### 3. Performance-Optimierung

**Nur wenn nötig memoizen:**
- `useMemo`: Teure Berechnungen (Filter, Transformationen)
- `useCallback`: Callbacks in dep-Arrays
- `React.memo`: Komponenten mit vielen Re-Renders

**Beispiel:**
```tsx
function DataTable({ items }: { items: Item[] }) {
  // GUT: Teurer Filter wird gecached
  const filtered = useMemo(
    () => items.filter(i => i.score > 0.5),
    [items]
  );

  return <Table data={filtered} />;
}
```

### 4. TypeScript-First

**Alle Props typisieren:**
```tsx
interface ButtonProps {
  variant: 'primary' | 'secondary';
  onClick: () => void;
  children: React.ReactNode;
  disabled?: boolean;
}

export function Button({ variant, onClick, children, disabled }: ButtonProps) {
  return <button onClick={onClick}>{children}</button>;
}
```

### 5. Error Boundaries

**Jede kritische Komponente braucht Fehlerbehandlung:**
```tsx
<ErrorBoundary fallback={<ErrorMessage />}>
  <CriticalComponent />
</ErrorBoundary>
```

## Projektspezifische Guidelines

### digitalTwin (React 19 + Gemini AI)

- **AI-Responses immer validieren**: Gemini-Output kann fehlerhaft sein
- **Streaming UI**: Zeige Loading-States bei AI-Calls
- **Camera-Permissions**: Immer Fallback für denied permissions
- **Mock-Daten**: Dev-Mode sollte ohne Camera funktionieren

**Beispiel:**
```tsx
function CameraView() {
  const [stream, setStream] = useState<MediaStream | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    navigator.mediaDevices.getUserMedia({ video: true })
      .then(setStream)
      .catch(() => setError('Camera access denied'));
  }, []);

  if (error) return <MockCameraView />;
  return <video srcObject={stream} />;
}
```

### DresdenAIInsights (React 18.3 + Radix UI)

- **40+ Radix UI Components**: Nutze vorhandene Primitives
- **Accessibility**: ARIA-Labels für alle interaktiven Elemente
- **Particle Effects**: Performance-Monitor bei Animationen
- **IONOS-Deployment**: Statische Builds nur (kein SSR)

**Beispiel Radix:**
```tsx
import * as Dialog from '@radix-ui/react-dialog';

function ContactModal() {
  return (
    <Dialog.Root>
      <Dialog.Trigger>Kontakt</Dialog.Trigger>
      <Dialog.Portal>
        <Dialog.Overlay />
        <Dialog.Content>
          <Dialog.Title>Kontaktformular</Dialog.Title>
          {/* Content */}
        </Dialog.Content>
      </Dialog.Portal>
    </Dialog.Root>
  );
}
```

### ManufacturingInsideAnalyzer (React 19 + Vercel)

- **DSGVO-Compliance**: Siehe DSGVO-Skill
- **Smart Sampling**: UI muss große Datasets (1567×591) handhaben
- **Vercel Edge Functions**: API-Calls via `/api/` Routes
- **10-Sekunden-Timeout**: Zeige Progress bei langen Analysen

**Beispiel:**
```tsx
function AnalysisProgress({ fileSize }: { fileSize: number }) {
  const [progress, setProgress] = useState(0);
  const timeoutMs = 10000;

  useEffect(() => {
    const interval = setInterval(() => {
      setProgress(p => Math.min(p + 10, 95));
    }, timeoutMs / 10);
    return () => clearInterval(interval);
  }, []);

  return <ProgressBar value={progress} max={100} />;
}
```

## Vite-spezifische Best Practices

### 1. Environment Variables

**Nur `VITE_`-Prefix wird exposed:**
```ts
// .env
VITE_API_URL=https://api.example.com
SECRET_KEY=xxx  // ❌ Wird NICHT in Bundle inkludiert

// code
const apiUrl = import.meta.env.VITE_API_URL; // ✅
```

### 2. Code Splitting

**Route-basiertes Lazy Loading:**
```tsx
import { lazy, Suspense } from 'react';

const Dashboard = lazy(() => import('./pages/Dashboard'));

function App() {
  return (
    <Suspense fallback={<Spinner />}>
      <Routes>
        <Route path="/dashboard" element={<Dashboard />} />
      </Routes>
    </Suspense>
  );
}
```

### 3. Asset Handling

**Bilder optimieren:**
```tsx
// ❌ Nicht optimal
<img src="/logo.png" />

// ✅ Vite optimiert automatisch
import logoUrl from './assets/logo.png';
<img src={logoUrl} alt="Logo" />
```

## Testing

**Jede Komponente braucht mindestens:**
1. Smoke Test (rendert ohne Fehler)
2. User Interaction Test (Clicks, Inputs)
3. Edge Case Test (Leere Daten, Fehler)

**Beispiel mit Vitest:**
```tsx
import { render, screen } from '@testing-library/react';
import { Button } from './Button';

test('renders button', () => {
  render(<Button>Click me</Button>);
  expect(screen.getByText('Click me')).toBeInTheDocument();
});

test('calls onClick', async () => {
  const onClick = vi.fn();
  render(<Button onClick={onClick}>Click</Button>);
  await userEvent.click(screen.getByText('Click'));
  expect(onClick).toHaveBeenCalledOnce();
});
```

## Deployment-Checklist

Vor jedem Deployment:
- [ ] `npm run build` erfolgreich
- [ ] TypeScript Errors: `npx tsc --noEmit`
- [ ] Keine Console.logs im Production-Code
- [ ] Environment Variables korrekt gesetzt
- [ ] Asset-Pfade relativ (kein `C:\...`)
- [ ] Bundle Size < 500KB (gzipped)

## Ressourcen

Für tiefere Themen siehe:
- [styling.md](resources/styling.md) - TailwindCSS + CSS-in-JS
- [state-management.md](resources/state-management.md) - Context vs. Zustand
- [performance.md](resources/performance.md) - Profiling + Optimization
- [accessibility.md](resources/accessibility.md) - ARIA + Keyboard Navigation

## Layout Anti-Patterns (KRITISCH - Lesson Learned 2025-11-18)

### Problem: Überschüssige Nested Divs

**Root Cause**: Unbewusste Verschachtelung führt zu Layout-Problemen und schwer debugbaren Style-Konflikten.

**Goldene Regel**: **Maximal 2 nested divs** pro logische UI-Section.

**SCHLECHT:**
```tsx
// ❌ 5 Ebenen verschachtelt - unnötig komplex
function LoginGate({ children }: { children: React.ReactNode }) {
  return (
    <div className="outer-wrapper">
      <div className="container">
        <div className="inner-container">
          <div className="content-wrapper">
            <div className="content">
              {children}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
```

**GUT:**
```tsx
// ✅ Maximal 2 Ebenen - einfach und debuggbar
function LoginGate({ children }: { children: React.ReactNode }) {
  return (
    <div className="min-h-screen bg-slate-900 flex items-center justify-center p-4">
      <div className="w-full max-w-md">
        {children}
      </div>
    </div>
  );
}
```

**Warum?**
- **Debugging**: Weniger Ebenen → einfacher zu inspizieren
- **Performance**: Weniger DOM-Nodes → schneller Rendering
- **Maintainability**: Klare Struktur → leichter zu verstehen
- **CSS Specificity**: Weniger Verschachtelung → weniger Konflikte

**Inspection Technique** (Chrome DevTools):
```bash
# Open DevTools → Elements Tab
# Click component → Check nesting depth
# If >2 levels of plain divs (no semantic tags) → REFACTOR!
```

---

### Problem: Conflicting Styles durch Layout-Ebenen

**Symptom**: Styles werden überschrieben, obwohl Tailwind-Classes korrekt sind.

**Root Cause**: Übergeordnete Wrapper haben konkurrierende Layout-Properties.

**Beispiel - Conflict:**
```tsx
// ❌ Parent hat flex, Child hat grid → Konflikt!
function App() {
  return (
    <div className="flex flex-col min-h-screen"> {/* flex container */}
      <main className="grid grid-cols-2 gap-4">  {/* grid versucht sich durchzusetzen */}
        {/* Layout bricht! */}
      </main>
    </div>
  );
}

// ✅ Konsistentes Layout-System
function App() {
  return (
    <div className="flex flex-col min-h-screen">
      <main className="flex-1 p-4">  {/* flex child */}
        <div className="grid grid-cols-2 gap-4"> {/* grid INNERHALB von flex */}
          {/* Layout funktioniert! */}
        </div>
      </main>
    </div>
  );
}
```

**Debugging Checklist**:
- [ ] Prüfe Parent-Layout (flex vs grid vs block)
- [ ] Nur EINE Layout-Strategie pro Hierarchie-Ebene
- [ ] Verwende semantische Tags (`<main>`, `<section>`, `<article>`) statt generischer `<div>`
- [ ] Tailwind DevTools: Inspect applied classes

---

### Problem: CSP (Content Security Policy) Violations

**Symptom**: Production-Build funktioniert, aber Styles fehlen oder werden blockiert.

**Root Cause**: External CDNs (z.B. Tailwind CDN) werden von CSP blockiert.

**SCHLECHT:**
```html
<!-- ❌ NIEMALS in Production! -->
<head>
  <script src="https://cdn.tailwindcss.com"></script>
</head>
```

**GUT:**
```bash
# ✅ Lokale Installation
npm install -D tailwindcss postcss autoprefixer
npx tailwindcss init -p

# tailwind.config.js
export default {
  content: ['./index.html', './src/**/*.{js,ts,jsx,tsx}'],
  theme: { extend: {} },
  plugins: [],
}

# postcss.config.js
export default {
  plugins: {
    tailwindcss: {},
    autoprefixer: {},
  },
}
```

**CSP Header Verification**:
```bash
# Check production CSP
curl -I https://your-production-url.com | grep -i content-security-policy

# Expected (without CDN):
# Content-Security-Policy: default-src 'self'; script-src 'self' 'unsafe-inline'
```

**Pre-Commit Hook**:
```bash
# .husky/pre-commit
#!/bin/sh

# Check for CDN links in HTML/TSX files
if grep -r "cdn\." src/ index.html 2>/dev/null; then
  echo "❌ ERROR: External CDN detected!"
  echo "Remove CDN links and use local packages instead."
  exit 1
fi
```

**Regel**:
- **Nie** externe CDNs in Production (Tailwind, jQuery, Bootstrap)
- **Immer** lokale npm packages verwenden
- **CSP-Header** testen nach jedem Deployment
- **Pre-commit hooks** für CDN-Detection

---

## Theme Switching Patterns (Dark Mode)

### Implementation mit Context + localStorage

**Pattern**: Theme-State in Context, Persistence in localStorage

```tsx
// contexts/ThemeContext.tsx
import { createContext, useContext, useState, useEffect, ReactNode } from 'react';

type Theme = 'light' | 'dark';

interface ThemeContextType {
  theme: Theme;
  toggleTheme: () => void;
}

const ThemeContext = createContext<ThemeContextType | undefined>(undefined);

export function ThemeProvider({ children }: { children: ReactNode }) {
  const [theme, setTheme] = useState<Theme>(() => {
    // Read from localStorage on mount
    const stored = localStorage.getItem('theme');
    return (stored === 'dark' || stored === 'light') ? stored : 'light';
  });

  useEffect(() => {
    // Persist to localStorage on change
    localStorage.setItem('theme', theme);

    // Update document class for Tailwind dark: selector
    document.documentElement.classList.toggle('dark', theme === 'dark');
  }, [theme]);

  const toggleTheme = () => {
    setTheme(prev => prev === 'light' ? 'dark' : 'light');
  };

  return (
    <ThemeContext.Provider value={{ theme, toggleTheme }}>
      {children}
    </ThemeContext.Provider>
  );
}

export function useTheme() {
  const context = useContext(ThemeContext);
  if (!context) {
    throw new Error('useTheme must be used within ThemeProvider');
  }
  return context;
}
```

**Usage in Components**:
```tsx
// components/ThemeToggle.tsx
import { useTheme } from '../contexts/ThemeContext';

export function ThemeToggle() {
  const { theme, toggleTheme } = useTheme();

  return (
    <button
      onClick={toggleTheme}
      className="p-2 rounded-lg bg-slate-200 dark:bg-slate-700"
      aria-label={`Switch to ${theme === 'light' ? 'dark' : 'light'} mode`}
    >
      {theme === 'light' ? '🌙' : '☀️'}
    </button>
  );
}
```

**Tailwind Config**:
```js
// tailwind.config.js
export default {
  content: ['./index.html', './src/**/*.{js,ts,jsx,tsx}'],
  darkMode: 'class', // ✅ Use 'dark' class on <html>
  theme: {
    extend: {
      colors: {
        // Custom dark mode colors
        dark: {
          bg: '#0f172a',      // slate-900
          surface: '#1e293b', // slate-800
          border: '#334155',  // slate-700
        }
      }
    }
  }
}
```

**CSS Variables Alternative** (für komplexe Themes):
```css
/* index.css */
:root {
  --color-bg: #ffffff;
  --color-text: #000000;
  --color-primary: #3b82f6;
}

:root.dark {
  --color-bg: #0f172a;
  --color-text: #f1f5f9;
  --color-primary: #60a5fa;
}

/* Usage in components */
.button {
  background-color: var(--color-primary);
  color: var(--color-text);
}
```

**System Preference Detection**:
```tsx
// Detect OS dark mode preference
const [theme, setTheme] = useState<Theme>(() => {
  const stored = localStorage.getItem('theme');
  if (stored) return stored as Theme;

  // Fallback to system preference
  return window.matchMedia('(prefers-color-scheme: dark)').matches
    ? 'dark'
    : 'light';
});

// Listen to system preference changes
useEffect(() => {
  const mediaQuery = window.matchMedia('(prefers-color-scheme: dark)');

  const handleChange = (e: MediaQueryListEvent) => {
    if (!localStorage.getItem('theme')) {
      setTheme(e.matches ? 'dark' : 'light');
    }
  };

  mediaQuery.addEventListener('change', handleChange);
  return () => mediaQuery.removeEventListener('change', handleChange);
}, []);
```

**Regel**:
- **localStorage** für User-Präferenz (überschreibt System)
- **System Preference** als Fallback
- **document.documentElement.classList** für Tailwind `dark:` selector
- **CSS Variables** für komplexe Theme-Systeme mit vielen Farben

---

## Design System Components (Accessibility-First)

### Dropdown mit Accessibility

**Problem**: White text on white background (unreadable)

**Root Cause**: Browser default styles für `<option>` ignorieren `color` property.

**SCHLECHT:**
```tsx
// ❌ Options erben NICHT parent color/background!
<select className="bg-slate-800 text-white">
  <option value="1">Option 1</option>
  <option value="2">Option 2</option>
</select>
// Resultat: Weißer Text auf weißem Hintergrund im Dropdown!
```

**GUT:**
```tsx
// ✅ Explizite Styles auf <option> Elementen
<select className="bg-slate-800 text-white border border-slate-600 rounded px-3 py-2">
  <option value="1" className="bg-slate-800 text-white">Option 1</option>
  <option value="2" className="bg-slate-800 text-white">Option 2</option>
  <option value="3" className="bg-slate-800 text-white">Option 3</option>
</select>
```

**Best Practice - Reusable Component**:
```tsx
// components/Select.tsx
interface SelectOption {
  value: string;
  label: string;
}

interface SelectProps {
  options: SelectOption[];
  value: string;
  onChange: (value: string) => void;
  label: string;
  disabled?: boolean;
}

export function Select({ options, value, onChange, label, disabled }: SelectProps) {
  return (
    <div className="flex flex-col gap-2">
      <label
        htmlFor={label}
        className="text-sm font-medium text-slate-200"
      >
        {label}
      </label>
      <select
        id={label}
        value={value}
        onChange={(e) => onChange(e.target.value)}
        disabled={disabled}
        className="bg-slate-800 text-white border border-slate-600 rounded px-3 py-2
                   focus:outline-none focus:ring-2 focus:ring-blue-500
                   disabled:opacity-50 disabled:cursor-not-allowed"
        aria-label={label}
      >
        {options.map(({ value, label }) => (
          <option
            key={value}
            value={value}
            className="bg-slate-800 text-white" // ✅ CRITICAL!
          >
            {label}
          </option>
        ))}
      </select>
    </div>
  );
}
```

**Usage**:
```tsx
<Select
  label="Analysis Type"
  value={selectedType}
  onChange={setSelectedType}
  options={[
    { value: 'quality', label: 'Quality Analysis' },
    { value: 'maintenance', label: 'Predictive Maintenance' },
    { value: 'performance', label: 'Performance Optimization' },
  ]}
/>
```

**Accessibility Checklist**:
- [ ] `<label>` mit `for` / `htmlFor` verknüpft
- [ ] `aria-label` für Screen Readers
- [ ] `disabled` state mit opacity + cursor
- [ ] `focus:ring` für Keyboard-Navigation sichtbar
- [ ] Explizite `className` auf JEDEM `<option>`

---

### File Upload Component

**Pattern**: Drag & Drop + File Input mit Preview

```tsx
// components/FileUpload.tsx
import { useState, useRef } from 'react';

interface FileUploadProps {
  accept: string;
  maxSizeMB: number;
  onFileSelected: (file: File) => void;
}

export function FileUpload({ accept, maxSizeMB, onFileSelected }: FileUploadProps) {
  const [dragActive, setDragActive] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  const handleFile = (file: File) => {
    // Validate file size
    const maxSizeBytes = maxSizeMB * 1024 * 1024;
    if (file.size > maxSizeBytes) {
      setError(`File too large. Max size: ${maxSizeMB}MB`);
      return;
    }

    // Validate file type
    const fileExtension = file.name.split('.').pop()?.toLowerCase();
    const acceptedExtensions = accept.split(',').map(ext => ext.trim().replace('.', ''));

    if (!fileExtension || !acceptedExtensions.includes(fileExtension)) {
      setError(`Invalid file type. Accepted: ${accept}`);
      return;
    }

    setError(null);
    onFileSelected(file);
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setDragActive(false);

    const file = e.dataTransfer.files[0];
    if (file) handleFile(file);
  };

  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) handleFile(file);
  };

  return (
    <div className="w-full">
      <div
        className={`
          border-2 border-dashed rounded-lg p-8 text-center cursor-pointer
          transition-colors duration-200
          ${dragActive ? 'border-blue-500 bg-blue-50' : 'border-slate-600 bg-slate-800'}
          hover:border-blue-400
        `}
        onDragEnter={() => setDragActive(true)}
        onDragLeave={() => setDragActive(false)}
        onDragOver={(e) => e.preventDefault()}
        onDrop={handleDrop}
        onClick={() => inputRef.current?.click()}
      >
        <input
          ref={inputRef}
          type="file"
          accept={accept}
          onChange={handleChange}
          className="hidden"
          aria-label="File upload"
        />

        <div className="flex flex-col items-center gap-2">
          <svg className="w-12 h-12 text-slate-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M15 13l-3-3m0 0l-3 3m3-3v12" />
          </svg>
          <p className="text-slate-200">
            Drag & drop file here, or <span className="text-blue-400 underline">browse</span>
          </p>
          <p className="text-sm text-slate-400">
            Accepted: {accept} (Max {maxSizeMB}MB)
          </p>
        </div>
      </div>

      {error && (
        <p className="mt-2 text-sm text-red-400" role="alert">
          {error}
        </p>
      )}
    </div>
  );
}
```

**Usage**:
```tsx
<FileUpload
  accept=".csv,.xlsx,.json"
  maxSizeMB={5}
  onFileSelected={(file) => {
    console.log('Selected file:', file.name);
    // Process file...
  }}
/>
```

---

### Progress Indicator (Long-Running Operations)

**Pattern**: Fake progress with realistic timing

```tsx
// components/AnalysisProgress.tsx
import { useState, useEffect } from 'react';

interface AnalysisProgressProps {
  timeoutMs: number;  // Expected completion time (e.g., 10000)
  onComplete?: () => void;
}

export function AnalysisProgress({ timeoutMs, onComplete }: AnalysisProgressProps) {
  const [progress, setProgress] = useState(0);

  useEffect(() => {
    // Incremental progress updates (10 steps)
    const interval = setInterval(() => {
      setProgress(prev => {
        const next = prev + 10;
        if (next >= 100) {
          clearInterval(interval);
          onComplete?.();
          return 100;
        }
        return next;
      });
    }, timeoutMs / 10);

    return () => clearInterval(interval);
  }, [timeoutMs, onComplete]);

  return (
    <div className="w-full">
      <div className="flex justify-between mb-2">
        <span className="text-sm text-slate-400">Analyzing data...</span>
        <span className="text-sm text-slate-200">{progress}%</span>
      </div>

      <div className="w-full h-2 bg-slate-700 rounded-full overflow-hidden">
        <div
          className="h-full bg-blue-500 transition-all duration-300 ease-out"
          style={{ width: `${progress}%` }}
          role="progressbar"
          aria-valuenow={progress}
          aria-valuemin={0}
          aria-valuemax={100}
        />
      </div>

      {progress >= 95 && (
        <p className="mt-2 text-sm text-yellow-400">
          Finalizing analysis...
        </p>
      )}
    </div>
  );
}
```

**Usage with API Call**:
```tsx
function AnalysisView() {
  const [analyzing, setAnalyzing] = useState(false);
  const [result, setResult] = useState<string | null>(null);

  const handleAnalyze = async (data: any) => {
    setAnalyzing(true);

    try {
      const response = await fetch('/api/analyze', {
        method: 'POST',
        body: JSON.stringify(data),
      });

      const result = await response.json();
      setResult(result.analysis);
    } catch (error) {
      console.error('Analysis failed:', error);
    } finally {
      setAnalyzing(false);
    }
  };

  if (analyzing) {
    return <AnalysisProgress timeoutMs={10000} />;
  }

  if (result) {
    return <div>{result}</div>;
  }

  return <button onClick={() => handleAnalyze(data)}>Analyze</button>;
}
```

**Regel**:
- **Fake progress** ist OK für UX (besser als Spinner)
- **95% cap** bis API tatsächlich fertig (verhindert 100% ohne Ergebnis)
- **aria-* attributes** für Accessibility
- **Optimistic UI**: Zeige Fortschritt, auch wenn Backend-Timing variiert

---

## Support

Bei Fragen zu spezifischen Patterns:
1. Prüfe project-spezifische `AGENTS.md`
2. Siehe vorhandene Komponenten als Beispiele
3. Frage nach Code-Review via Agent

---

## Lesson Learned Origins

**Layout Anti-Patterns**: Incident 2025-11-18 (Nested divs caused style conflicts)
**CSP Violations**: Incident 2025-11-18 (Tailwind CDN blocked in production)
**Dropdown Accessibility**: Bug 2025-11-19 (White text on white background)
**Theme Switching**: Feature 2025-11-12 (OAuth login with dark mode)
