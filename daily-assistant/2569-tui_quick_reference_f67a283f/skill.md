# TUI Component Quick Reference

**Building a new TUI component? Start here.**

---

## 🚀 Quick Start Template

```typescript
import { Container, Text, truncateToWidth } from "@mariozechner/pi-tui";
import type { Theme } from "@mariozechner/pi-coding-agent";
import { 
  getVisibleWidth, 
  calculateBorderFill, 
  validateLineWidth,
  renderTopBorder,
  renderBottomBorder
} from "./tui-validation.js";

export class MyComponent {
  private theme: Theme;
  private cachedWidth?: number;
  private cachedLines?: string[];

  constructor(theme: Theme) {
    this.theme = theme;
  }

  render(width: number): string[] {
    // Cache check
    if (this.cachedLines && this.cachedWidth === width) {
      return this.cachedLines;
    }

    const lines: string[] = [];

    // Top border
    const title = this.theme.fg("accent", this.theme.bold("My Component"));
    const topBorder = renderTopBorder(title, width, (s) => this.theme.fg("border", s));
    validateLineWidth(topBorder, width, "MyComponent", "top border");
    lines.push(topBorder);

    // Empty line
    const emptyLine = this.theme.fg("border", `│${" ".repeat(width - 2)}│`);
    validateLineWidth(emptyLine, width, "MyComponent", "empty line");
    lines.push(emptyLine);

    // Content
    const content = "Your content here";
    const contentLine = this.renderContentLine(content, width);
    lines.push(contentLine);

    // Empty line
    lines.push(emptyLine);

    // Bottom border
    const bottomBorder = renderBottomBorder(width, (s) => this.theme.fg("border", s));
    validateLineWidth(bottomBorder, width, "MyComponent", "bottom border");
    lines.push(bottomBorder);

    // Cache results
    this.cachedWidth = width;
    this.cachedLines = lines;
    return lines;
  }

  private renderContentLine(content: string, width: number): string {
    const innerWidth = width - 4; // 2 border + 2 spaces
    const truncated = truncateToWidth(content, innerWidth);
    const contentWidth = getVisibleWidth(truncated);
    const padding = Math.max(0, innerWidth - contentWidth);
    
    const line = this.theme.fg("border", "│") +
      "  " +
      truncated +
      " ".repeat(padding) +
      this.theme.fg("border", "│");
    
    validateLineWidth(line, width, "MyComponent", "content");
    return line;
  }

  invalidate(): void {
    this.cachedWidth = undefined;
    this.cachedLines = undefined;
  }
}
```

---

## 📏 The Golden Rules

### 1. ALWAYS use `truncateToWidth()` on final lines

```typescript
// ❌ BAD
lines.push(content);

// ✅ GOOD
lines.push(truncateToWidth(content, width));
```

### 2. NEVER hardcode border calculations

```typescript
// ❌ BAD
"┌─ Title " + "─".repeat(width - 15) + "┐"

// ✅ GOOD
const titleWidth = getVisibleWidth(title);
const fill = calculateBorderFill(width, titleWidth);
"┌─ " + title + " " + "─".repeat(fill) + "┐"
```

### 3. ALWAYS validate critical lines

```typescript
// ❌ BAD
lines.push(border);

// ✅ GOOD
validateLineWidth(border, width, "MyComponent", "top border");
lines.push(border);
```

### 4. Use `getVisibleWidth()` not `.length`

```typescript
// ❌ BAD
const width = title.length;  // Includes ANSI codes!

// ✅ GOOD
const width = getVisibleWidth(title);  // Strips ANSI
```

---

## 🎨 Common Patterns

### Bordered Panel

```typescript
// Top border
const title = this.theme.fg("accent", "Panel Title");
const topBorder = renderTopBorder(title, width, (s) => this.theme.fg("border", s));
validateLineWidth(topBorder, width, "Component", "top border");
lines.push(topBorder);

// Content lines
for (const item of items) {
  const line = renderContentLine(item, width, (s) => this.theme.fg("border", s));
  validateLineWidth(line, width, "Component", "content");
  lines.push(line);
}

// Bottom border
const bottomBorder = renderBottomBorder(width, (s) => this.theme.fg("border", s));
validateLineWidth(bottomBorder, width, "Component", "bottom border");
lines.push(bottomBorder);
```

### Simple List (No Borders)

```typescript
for (const item of items) {
  const styled = this.theme.fg("dim", `• ${item}`);
  const line = truncateToWidth(styled, width);
  validateLineWidth(line, width, "Component", "list item");
  lines.push(line);
}
```

### Progress Bar

```typescript
import { renderProgressBar } from "./progress-bar.js";

const barWidth = Math.min(20, width - 10); // Leave room for label
const bar = renderProgressBar(progress, barWidth, theme);
const label = `Progress: ${progress}%`;
const line = truncateToWidth(`${label} [${bar}]`, width);
validateLineWidth(line, width, "Component", "progress");
lines.push(line);
```

---

## ⚠️ Common Mistakes

### Mistake #1: Hardcoded Border Math

```typescript
// ❌ WRONG
"┌─ " + title + " " + "─".repeat(width - 22) + "┐"
```

**Why it fails:** Assumes title is always exactly 17 visible characters.

**Fix:**
```typescript
// ✅ CORRECT
const titleWidth = getVisibleWidth(title);
const fill = calculateBorderFill(width, titleWidth);
"┌─ " + title + " " + "─".repeat(fill) + "┐"
```

---

### Mistake #2: Using `.length` on Styled Strings

```typescript
// ❌ WRONG
const title = this.theme.fg("accent", "Title");
const width = title.length;  // Includes \x1b[31m...\x1b[0m codes!
```

**Why it fails:** ANSI escape codes add invisible characters.

**Fix:**
```typescript
// ✅ CORRECT
const title = this.theme.fg("accent", "Title");
const width = getVisibleWidth(title);  // Strips ANSI codes
```

---

### Mistake #3: No Validation

```typescript
// ❌ WRONG
lines.push(border);
```

**Why it fails:** If border is too wide, Pi crashes with no debug info.

**Fix:**
```typescript
// ✅ CORRECT
validateLineWidth(border, width, "MyComponent", "top border");
lines.push(border);
```

Then check `~/.pi/agent/tui-width-errors.log` if crashes occur.

---

### Mistake #4: Forgetting to Cache

```typescript
// ❌ BAD
render(width: number): string[] {
  // Recalculates everything on every render
  const lines = [];
  // ... expensive rendering ...
  return lines;
}
```

**Why it fails:** Performance - renders hundreds of times per second.

**Fix:**
```typescript
// ✅ GOOD
render(width: number): string[] {
  if (this.cachedLines && this.cachedWidth === width) {
    return this.cachedLines;
  }
  
  const lines = [];
  // ... rendering ...
  
  this.cachedWidth = width;
  this.cachedLines = lines;
  return lines;
}

invalidate(): void {
  this.cachedWidth = undefined;
  this.cachedLines = undefined;
}
```

---

## 🧪 Testing Checklist

Before committing a new TUI component:

- [ ] Import `tui-validation.ts` utilities
- [ ] Use `calculateBorderFill()` for all borders
- [ ] Use `getVisibleWidth()` not `.length`
- [ ] Call `validateLineWidth()` on all critical lines
- [ ] Use `truncateToWidth()` on final output
- [ ] Implement caching (`cachedLines`, `cachedWidth`)
- [ ] Test in narrow terminal (80 cols)
- [ ] Test in wide terminal (200+ cols)
- [ ] Test with themed/colored content
- [ ] Check `~/.pi/agent/tui-width-errors.log` for violations

---

## 📖 API Reference

### `tui-validation.ts`

```typescript
// Width measurement
getVisibleWidth(str: string): number
  // Strips ANSI codes, returns actual visible character count

// Border calculations
calculateBorderFill(width: number, titleWidth: number): number
  // Returns fill needed: width - titleWidth - 5

// Rendering helpers
renderTopBorder(title: string, width: number, theme: (s: string) => string): string
  // Renders: ┌─ title ──────────┐

renderBottomBorder(width: number, theme: (s: string) => string): string
  // Renders: └─────────────────┘

renderEmptyLine(width: number, theme: (s: string) => string): string
  // Renders: │                 │

renderContentLine(content: string, width: number, theme: (s: string) => string): string
  // Renders: │  content        │

// Validation
validateLineWidth(line: string, maxWidth: number, component: string, context: string): void
  // Logs violations to ~/.pi/agent/tui-width-errors.log
  // Does NOT throw - just logs for debugging
```

---

## 🔗 Related Documentation

- **Pi TUI API:** `/opt/homebrew/lib/node_modules/@mariozechner/pi-coding-agent/docs/tui.md`
- **Dex TUI Reference:** `06-Resources/Dex_System/Pi_TUI_Reference.md`
- **Bug fix details:** `System/PRDs/tui-width-bug-fix.md`
- **Examples:** Look at `progress-indicator.ts`, `task-board.ts`, `daily-plan-wizard.ts`

---

## 💡 Pro Tips

1. **Start with a working example** - Copy `progress-indicator.ts` structure
2. **Validate early** - Add validation from day 1, not after bugs appear
3. **Test narrow first** - 80-column terminals expose issues fastest
4. **Use theme callbacks** - `(s) => theme.fg("border", s)` makes colors consistent
5. **Cache aggressively** - Rendering is expensive, cache whenever width unchanged

---

**Questions?** Check `System/PRDs/tui-width-fix-summary.md` for debugging tips.
