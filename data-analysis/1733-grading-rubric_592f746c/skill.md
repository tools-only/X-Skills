# Design Grading Rubric

Detailed scoring criteria for evaluating web UI design quality from a world-class designer perspective.

## Scoring Philosophy

**Goal**: Identify designs that are distinctive, intentional, and memorable.

**Key Principle**: Bold maximalism and refined minimalism both work. The differentiator is **intentionality**, not intensity. Generic, safe, predictable design always fails.

---

## Dimension 1: Typography (20%)

### Score 9-10: Exceptional
- Distinctive, memorable font pairing (display + body)
- Perfect typographic hierarchy (H1 > H2 > H3 > body clearly distinct)
- Intentional font weights creating rhythm (300, 400, 600, 800)
- Line-height optimized per context (1.2 headings, 1.6 body)
- Letter-spacing refined for large text
- Type scale following mathematical ratio (1.25 or 1.333)

### Score 7-8: Good
- Non-generic fonts with clear character
- Clear hierarchy with 3+ distinct levels
- Good weight variation
- Readable line-height
- Minor refinement opportunities

### Score 5-6: Acceptable
- Somewhat generic but not default fonts
- Basic hierarchy present but weak
- Limited weight usage (2 weights)
- Acceptable but not optimized line-height

### Score 3-4: Poor
- Default/system fonts (Inter, Roboto, Arial, sans-serif)
- Weak or inconsistent hierarchy
- Single font weight throughout
- Line-height issues affecting readability

### Score 1-2: Failing
- Browser defaults visible
- No hierarchy
- Unreadable text
- Accessibility failures

### Typography Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| `font-family: Inter` | -3 | Direct CSS match |
| `font-family: Roboto` | -3 | Direct CSS match |
| `font-family: Arial` | -3 | Direct CSS match |
| `font-family: system-ui` | -2 | Lazy default |
| Single `font-weight` value | -2 | No variation |
| `line-height: 1` or `1.2` for body | -1 | Too tight |
| No `font-weight` on headings | -1 | Missing emphasis |

---

## Dimension 2: Color & Contrast (15%)

### Score 9-10: Exceptional
- Unique, memorable color palette
- Perfect WCAG AAA contrast (7:1+)
- Intentional color roles (primary, secondary, accent)
- Color creates mood/emotion
- Masterful use of neutrals
- Dark/light modes both excellent

### Score 7-8: Good
- Non-generic color choices
- WCAG AA compliant (4.5:1)
- Clear color system
- Good accent usage
- Minor contrast issues in non-critical areas

### Score 5-6: Acceptable
- Somewhat generic but functional
- Mostly AA compliant
- Basic color differentiation
- Missing accent or too many accents

### Score 3-4: Poor
- Generic palettes (blue/purple gradient cliches)
- Some contrast failures
- Competing colors without hierarchy
- No clear system

### Score 1-2: Failing
- Major contrast failures
- Unreadable combinations
- No intentional palette

### Color Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| Purple-to-blue gradient on white | -3 | AI aesthetic cliche |
| Contrast < 4.5:1 on text | -3 | WCAG AA failure |
| Contrast < 3:1 on large text | -2 | WCAG failure |
| 5+ competing accent colors | -2 | Visual chaos |
| No accent color at all | -1 | Flat, boring |
| Pure black (#000) on white (#fff) | -1 | Harsh, unrefined |

---

## Dimension 3: Layout & Spacing (20%)

### Score 9-10: Exceptional
- Unexpected, memorable composition
- Intentional asymmetry or perfect symmetry (never accidental)
- Grid-breaking elements that create interest
- Perfect rhythm in spacing
- Generous negative space used meaningfully
- Visual flow guides the eye

### Score 7-8: Good
- Interesting layout choices
- Good spacing rhythm
- Some unexpected elements
- Clear visual hierarchy through space
- Minor alignment issues

### Score 5-6: Acceptable
- Standard but clean layouts
- Consistent spacing
- Predictable but functional
- No visual interest

### Score 3-4: Poor
- Cookie-cutter grid layouts
- Inconsistent spacing
- Cramped or too sparse
- No rhythm or intention

### Score 1-2: Failing
- Broken layouts
- Overlapping content
- Unusable spacing

### Layout Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| 12-column grid with no breaks | -2 | Predictable |
| Equal spacing everywhere (all 16px) | -2 | No rhythm |
| Cards in perfect 3x3 grid | -2 | Cookie-cutter |
| No margin on container | -2 | Cramped |
| Centered everything | -1 | Safe/boring |
| No max-width on text | -1 | Hard to read |

---

## Dimension 4: Motion & Interaction (15%)

### Score 9-10: Exceptional
- Meaningful, purposeful animations
- Perfect timing and easing
- Micro-interactions that delight
- Page transitions that guide
- Motion reinforces hierarchy
- No gratuitous animation

### Score 7-8: Good
- CSS transitions on interactive elements
- Good hover states
- Some entrance animations
- Appropriate duration (150-300ms)

### Score 5-6: Acceptable
- Basic hover states
- Simple transitions
- Nothing broken but nothing special

### Score 3-4: Poor
- Missing hover states
- No transitions (instant state changes)
- Jarring interactions

### Score 1-2: Failing
- Broken interactions
- Inaccessible motion
- No feedback on actions

### Motion Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| No `transition` on buttons | -2 | Jarring clicks |
| No hover state changes | -2 | Dead feeling |
| `transition: all` (lazy) | -1 | Performance issue |
| Duration > 500ms | -1 | Feels slow |
| Duration < 100ms | -1 | Too abrupt |
| No focus states | -2 | Accessibility |

---

## Dimension 5: Visual Polish (15%)

### Score 9-10: Exceptional
- Rich visual depth and texture
- Masterful use of shadows
- Subtle gradients and overlays
- Refined borders and dividers
- Attention to every detail
- Cohesive visual language

### Score 7-8: Good
- Good use of shadows for depth
- Some texture or gradient
- Clean borders
- Polished feel overall

### Score 5-6: Acceptable
- Basic shadows present
- Functional but flat
- Standard borders
- Missing refinement

### Score 3-4: Poor
- Flat solid backgrounds
- No shadows or depth
- Harsh borders or none
- Unpolished details

### Score 1-2: Failing
- Visually broken
- Clashing styles
- No cohesion

### Polish Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| `background: #fff` (pure flat) | -2 | No atmosphere |
| No `box-shadow` on cards | -2 | Flat, no depth |
| `border: 1px solid #ccc` | -1 | Generic, lazy |
| `border-radius: 0` everywhere | -1 | Harsh |
| No background variation | -1 | Monotonous |

---

## Dimension 6: Accessibility (15%)

### Score 9-10: Exceptional
- WCAG AAA compliant
- Perfect keyboard navigation
- Excellent screen reader support
- Skip links and landmarks
- Reduced motion support
- Focus visible on all elements

### Score 7-8: Good
- WCAG AA compliant
- Keyboard accessible
- Proper focus indicators
- Good semantic HTML
- ARIA where needed

### Score 5-6: Acceptable
- Mostly AA compliant
- Basic keyboard support
- Some focus indicators
- Minor issues

### Score 3-4: Poor
- Contrast failures
- Missing focus states
- Keyboard traps
- Poor semantics

### Score 1-2: Failing
- Major accessibility barriers
- Unusable without mouse
- Critical WCAG failures

### Accessibility Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| No visible focus indicator | -3 | :focus-visible missing |
| `outline: none` without replacement | -3 | Focus removed |
| Non-semantic buttons (`div onclick`) | -2 | Inaccessible |
| Images without alt text | -2 | Screen reader fail |
| Low contrast text | -3 | WCAG failure |
| No skip link | -1 | Keyboard navigation |

---

## Overall Score Calculation

```
Overall = (Typography × 0.20) +
          (Color × 0.15) +
          (Layout × 0.20) +
          (Motion × 0.15) +
          (Polish × 0.15) +
          (Accessibility × 0.15)
```

## Pass/Fail Thresholds

| Dimension | Pass Threshold | Critical Threshold |
|-----------|----------------|-------------------|
| Typography | >= 6 | < 4 blocks overall pass |
| Color | >= 6 | < 4 blocks overall pass |
| Layout | >= 6 | < 4 blocks overall pass |
| Motion | >= 5 | < 3 blocks overall pass |
| Polish | >= 5 | < 3 blocks overall pass |
| Accessibility | >= 7 | < 5 blocks overall pass |

**Overall Pass**: Score >= 8.0 AND no dimension below critical threshold

---

## Quick Reference Card

### What Makes a 10/10?
- Typography: Distinctive fonts, perfect hierarchy
- Color: Unique palette, AAA contrast
- Layout: Memorable composition, intentional space
- Motion: Purposeful, delightful micro-interactions
- Polish: Rich depth, refined details
- Accessibility: Exceeds WCAG AA

### What Makes a 5/10?
- Typography: Generic but readable
- Color: Functional, AA compliant
- Layout: Standard grid, predictable
- Motion: Basic hover states
- Polish: Flat but clean
- Accessibility: Mostly compliant

### What Makes a 2/10?
- Typography: Browser defaults, no hierarchy
- Color: Contrast failures, no system
- Layout: Broken or unusable
- Motion: No feedback, jarring
- Polish: Visually broken
- Accessibility: Major barriers
