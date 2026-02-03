# Accessibility & WCAG 2.2 Skill

**WCAG 2.2 Compliance** für barrierefreie Web-Apps (fokussiert auf DresdenAIInsights)

## Wann aktiviert?

- Keywords: accessibility, a11y, WCAG, ARIA, screen reader, keyboard navigation
- Du arbeitest an UI-Components
- Accessibility-Fragen

## WCAG 2.2 - Die 4 Prinzipien

### 1. Perceivable (Wahrnehmbar)
Informationen müssen für alle Sinne zugänglich sein.

### 2. Operable (Bedienbar)
UI muss bedienbar sein (Maus, Keyboard, Touch).

### 3. Understandable (Verständlich)
Inhalte und Bedienung müssen verständlich sein.

### 4. Robust
Kompatibel mit aktuellen und zukünftigen Technologien.

---

## Level AA Compliance (Mindest-Standard)

### 1. Text Alternatives (1.1.1)

**Regel:** Alle Nicht-Text-Inhalte brauchen Alternativen.

```tsx
// ✅ GUT
<img src="logo.png" alt="Dresden AI Insights Logo" />
<button aria-label="Menü öffnen">☰</button>

// ❌ SCHLECHT
<img src="logo.png" /> // Missing alt
<button>☰</button> // Icon ohne Label
```

### 2. Keyboard Access (2.1.1)

**Regel:** Alle Funktionen per Keyboard bedienbar.

```tsx
// ✅ Keyboard-navigierbar
<button onClick={handleClick}>Action</button>

// ❌ Nur Maus
<div onClick={handleClick}>Action</div> // Nicht tab-bar!

// ✅ Mit Keyboard-Support
<div
  role="button"
  tabIndex={0}
  onClick={handleClick}
  onKeyDown={(e) => e.key === 'Enter' && handleClick()}
>
  Action
</div>
```

### 3. Focus Visible (2.4.7)

**Regel:** Fokus muss sichtbar sein.

```css
/* ✅ Fokus-Indikator */
button:focus-visible {
  outline: 2px solid #0066cc;
  outline-offset: 2px;
}

/* ❌ NIEMALS Outline entfernen ohne Alternative! */
button:focus {
  outline: none; /* Gefahr! */
}
```

### 4. Color Contrast (1.4.3)

**Mindest-Kontrast:**
- Normal-Text: **4.5:1**
- Großer Text (18pt+): **3:1**

```tsx
// ✅ Ausreichend Kontrast
<p style={{ color: '#333', background: '#fff' }}>Text</p> // 12.6:1

// ❌ Zu wenig Kontrast
<p style={{ color: '#999', background: '#fff' }}>Text</p> // 2.8:1
```

**Tool:** [WebAIM Contrast Checker](https://webaim.org/resources/contrastchecker/)

### 5. Semantic HTML (1.3.1)

```tsx
// ✅ Semantisch
<nav>
  <ul>
    <li><a href="/">Home</a></li>
    <li><a href="/about">About</a></li>
  </ul>
</nav>

<main>
  <article>
    <h1>Titel</h1>
    <p>Inhalt...</p>
  </article>
</main>

// ❌ Nicht-semantisch
<div className="nav">
  <div><span onClick={...}>Home</span></div>
</div>
```

---

## ARIA (Accessible Rich Internet Applications)

### Wann ARIA nutzen?

**Regel #1:** ARIA nur wenn natives HTML nicht reicht!

```tsx
// ❌ Unnötiges ARIA
<button role="button">Click</button> // Button ist schon Button!

// ✅ Natives HTML bevorzugen
<button>Click</button>
```

### ARIA Roles

```tsx
// Custom Components
<div role="dialog" aria-labelledby="dialog-title">
  <h2 id="dialog-title">Bestätigung</h2>
  <p>Möchten Sie fortfahren?</p>
</div>

// Navigation
<div role="navigation" aria-label="Hauptnavigation">
  {/* Nav items */}
</div>

// Status Updates
<div role="status" aria-live="polite">
  Datei wird hochgeladen...
</div>
```

### ARIA States & Properties

```tsx
// Expanded State
<button aria-expanded={isOpen} aria-controls="menu-id">
  Menu
</button>
<div id="menu-id" hidden={!isOpen}>
  {/* Menu items */}
</div>

// Disabled
<button aria-disabled="true">Submit</button>

// Required
<input aria-required="true" />

// Invalid
<input aria-invalid={hasError} aria-describedby="error-msg" />
{hasError && <span id="error-msg">E-Mail ungültig</span>}
```

---

## Keyboard Navigation

### Tab Order

```tsx
// Standard: Natürliche Reihenfolge (top → bottom, left → right)
<input tabIndex={0} /> // Normale Tab-Reihenfolge
<button tabIndex={-1}>Skip from tab</button> // Fokussierbar via JS, nicht via Tab
<div tabIndex={0}>Fokussierbar</div> // Macht div tab-bar

// ❌ NIEMALS positive tabIndex verwenden!
<input tabIndex={1} /> // Verwirrt Nutzer!
```

### Keyboard Event Handling

```tsx
function AccessibleButton({ onClick, children }) {
  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      onClick();
    }
  };

  return (
    <div
      role="button"
      tabIndex={0}
      onClick={onClick}
      onKeyDown={handleKeyDown}
    >
      {children}
    </div>
  );
}
```

### Focus Management

```tsx
import { useRef, useEffect } from 'react';

function Dialog({ isOpen }) {
  const firstFocusRef = useRef<HTMLButtonElement>(null);

  useEffect(() => {
    if (isOpen) {
      firstFocusRef.current?.focus(); // Auto-Fokus beim Öffnen
    }
  }, [isOpen]);

  return (
    <div role="dialog">
      <button ref={firstFocusRef}>Erste Aktion</button>
      <button>Zweite Aktion</button>
    </div>
  );
}
```

---

## Testing

### Automated Testing

```bash
npm install -D @axe-core/react
```

```tsx
// index.tsx (Development nur!)
if (process.env.NODE_ENV !== 'production') {
  import('@axe-core/react').then((axe) => {
    axe.default(React, ReactDOM, 1000);
  });
}
```

### Manual Testing

**Keyboard-Only Test:**
1. Verstecke die Maus
2. Nutze nur Tab, Enter, Pfeiltasten
3. Ist alles erreichbar und bedienbar?

**Screen Reader Test:**
- **Windows:** NVDA (gratis)
- **macOS:** VoiceOver (integriert)
- **Linux:** Orca

**Empfohlene Browser/Screen Reader Combos:**
- NVDA + Firefox
- JAWS + Chrome
- VoiceOver + Safari

---

## DresdenAIInsights Spezifika

### Particle Animations

```tsx
// ✅ Reduzierte Animation für vestibular disorders
function ParticleBackground() {
  const prefersReducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  if (prefersReducedMotion) {
    return <StaticBackground />;
  }

  return <AnimatedParticles />;
}
```

### Color Scheme

```tsx
// Dark Theme mit ausreichend Kontrast
const theme = {
  bg: '#0a0a0a',       // Dunkel
  text: '#ffffff',     // Weiß (21:1 Kontrast!)
  accent: '#00d9ff',   // Cyan (8.6:1 auf Dunkel)
  muted: '#9ca3af'     // Grau (4.6:1)
};
```

### Interactive Elements

```tsx
// Radix UI hat Accessibility built-in
<Dialog.Root>
  <Dialog.Trigger>Open</Dialog.Trigger>
  <Dialog.Content>
    {/* Automatisch: Focus-Trap, Esc-Close, ARIA-Labels */}
  </Dialog.Content>
</Dialog.Root>
```

---

## WCAG 2.2 Neue Success Criteria

### Focus Appearance (2.4.11)

Fokus-Indikator muss:
- Mindestens 2 CSS pixels dick sein
- Kontrast ≥ 3:1 zu unfokussiertem State

```css
button:focus-visible {
  outline: 2px solid #0066cc; /* ✅ 2px */
  outline-offset: 2px;
}
```

### Dragging Movements (2.5.7)

Alternative zu Drag & Drop anbieten:

```tsx
// ✅ Drag UND Buttons
<DraggableItem>
  <button onClick={moveUp}>▲</button>
  <button onClick={moveDown}>▼</button>
</DraggableItem>
```

### Target Size (2.5.8)

Tap-Targets mindestens **24×24 CSS pixels**:

```css
button {
  min-width: 44px; /* iOS/Android Empfehlung */
  min-height: 44px;
  padding: 12px 24px;
}
```

---

## Compliance Checklist

- [ ] Alle Bilder haben `alt`-Text
- [ ] Fokus-Indikatoren sichtbar (2px, 3:1 Kontrast)
- [ ] Farbkontrast ≥ 4.5:1 (Normal), ≥ 3:1 (Groß)
- [ ] Alle Funktionen per Keyboard erreichbar
- [ ] ARIA-Labels wo nötig
- [ ] Semantisches HTML verwendet
- [ ] Heading-Hierarchie korrekt (h1 → h2 → h3)
- [ ] Forms haben Labels
- [ ] Error Messages sind klar
- [ ] Animationen respektieren `prefers-reduced-motion`
- [ ] Screen Reader Test durchgeführt
- [ ] axe DevTools keine Errors

---

## Ressourcen

- [WCAG 2.2 Guideline](https://www.w3.org/WAI/WCAG22/quickref/)
- [axe DevTools](https://www.deque.com/axe/devtools/)
- [WebAIM](https://webaim.org/)
- [A11y Project](https://www.a11yproject.com/)
