# UX Grading Rubric

Detailed scoring criteria for evaluating web application user experience. Focus on usability, efficiency, and user satisfaction—not visual aesthetics.

## Scoring Philosophy

**Goal**: Identify interfaces where users can accomplish tasks with minimal friction, confusion, or frustration.

**Key Principle**: Good UX is invisible. Users should focus on their goals, not on figuring out the interface. If users have to think about how to use something, that's a UX failure.

---

## Dimension 1: Usability (20%)

*Can users accomplish their goals with minimal effort?*

### Score 9-10: Exceptional
- Primary actions immediately obvious
- Tasks completable in minimum steps
- Clear, action-oriented labels (verbs)
- Consistent patterns throughout
- Zero guesswork required
- Works for first-time and returning users

### Score 7-8: Good
- Primary actions visible and clear
- Tasks completable with minor friction
- Labels mostly descriptive
- Mostly consistent patterns
- Occasional confusion for new users

### Score 5-6: Acceptable
- Primary actions findable but not obvious
- Tasks completable with effort
- Some vague labels ("Submit", "Continue")
- Inconsistencies between sections
- Users need trial and error

### Score 3-4: Poor
- Primary actions hidden or unclear
- Tasks require many unnecessary steps
- Labels ambiguous or confusing
- Major inconsistencies
- Users frequently get lost

### Score 1-2: Failing
- Primary actions impossible to find
- Tasks cannot be completed
- Labels meaningless
- No consistent patterns
- Users abandon in frustration

### Usability Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| "Submit" or "Click Here" buttons | -3 | Vague action text |
| Hidden primary CTA | -3 | Main action not visible |
| > 3 clicks to primary action | -2 | Deep navigation |
| Inconsistent button styles | -2 | Buttons look different across pages |
| No clear starting point | -2 | New users don't know what to do |
| Too many choices (> 7 options) | -1 | Choice overload |

---

## Dimension 2: Information Architecture (15%)

*Is content organized logically and easy to find?*

### Score 9-10: Exceptional
- Information hierarchy crystal clear
- Navigation matches user mental models
- Labels use user language, not internal jargon
- Critical info above the fold
- Scanning-friendly layout
- Search works perfectly

### Score 7-8: Good
- Clear hierarchy with minor issues
- Navigation mostly intuitive
- Labels mostly clear
- Important info visible
- Reasonable scanning
- Search functional

### Score 5-6: Acceptable
- Hierarchy present but weak
- Navigation requires learning
- Some jargon or unclear labels
- Some scrolling for key info
- Dense layout
- Basic search

### Score 3-4: Poor
- Flat or confusing hierarchy
- Navigation unintuitive
- Internal jargon prevalent
- Important info buried
- Wall-of-text layout
- Search broken or missing

### Score 1-2: Failing
- No discernible organization
- Users can't navigate
- Labels meaningless
- Critical info hidden
- Unusable information density

### Information Architecture Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| > 3 levels of navigation | -3 | Too deep |
| Internal jargon labels | -2 | Non-user language |
| No section headings | -2 | Can't scan |
| Important info below fold | -2 | Must scroll to key content |
| No breadcrumbs in deep nav | -1 | Can't tell where you are |
| Inconsistent categorization | -1 | Same item in multiple categories |

---

## Dimension 3: User Flows (20%)

*Are paths to goals efficient with clear progress?*

### Score 9-10: Exceptional
- Shortest possible path to goals
- Clear progress indicators (steps, %)
- Obvious next action at every point
- Easy to go back or cancel
- Success state unmistakable
- Works for all user paths

### Score 7-8: Good
- Efficient paths with minimal extra steps
- Progress shown in multi-step flows
- Next action usually obvious
- Back/cancel available
- Clear completion

### Score 5-6: Acceptable
- Paths have some unnecessary steps
- Basic progress indication
- Next action requires some thought
- Back possible but not obvious
- Completion could be clearer

### Score 3-4: Poor
- Many unnecessary steps
- No progress indication
- Unclear what to do next
- Difficult to go back
- Uncertain if task completed

### Score 1-2: Failing
- Impossible to complete flows
- Users stuck mid-flow
- Dead ends
- Can't go back
- No completion feedback

### User Flow Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| Dead end (no next action) | -3 | Page with no obvious next step |
| No progress in multi-step | -3 | Missing step indicator |
| Can't go back | -2 | No back button or working browser back |
| Unnecessary steps | -2 | Steps that don't add value |
| Unclear completion | -2 | No success confirmation |
| No cancel option | -1 | Trapped in flow |

---

## Dimension 4: Affordances & Signifiers (15%)

*Do interactive elements look interactive?*

### Score 9-10: Exceptional
- Every interactive element obviously clickable
- Non-interactive elements clearly not clickable
- Icons paired with labels
- Distinct hover/focus/active states
- Disabled state obvious
- Touch targets appropriately sized (44px+)

### Score 7-8: Good
- Most interactive elements clearly styled
- Good distinction from non-interactive
- Key icons labeled
- Hover states present
- Disabled distinguishable

### Score 5-6: Acceptable
- Buttons identifiable but not compelling
- Some confusion possible
- Icons sometimes unclear
- Basic hover states
- Disabled somewhat visible

### Score 3-4: Poor
- Buttons look like text
- Links not underlined
- Icons without labels
- Missing hover states
- Can't tell disabled from enabled

### Score 1-2: Failing
- Can't identify what's clickable
- No visual distinction for interactive
- Icons meaningless
- No state changes
- Users click everywhere trying to find actions

### Affordance Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| Text-only buttons (no border/bg) | -3 | Buttons look like regular text |
| Links without underline/color | -2 | Text links not styled |
| Icons without text labels | -2 | Ambiguous meaning |
| No hover state change | -2 | No feedback on hover |
| Disabled same as enabled | -2 | Can't tell state |
| Touch target < 44px | -1 | Hard to tap on mobile |

---

## Dimension 5: Feedback & Status (15%)

*Does the system respond to user actions?*

### Score 9-10: Exceptional
- Instant feedback on every action
- Clear loading states for async operations
- Success confirmations for important actions
- Undo available for reversible actions
- Real-time validation
- Progress for long operations

### Score 7-8: Good
- Feedback on most actions
- Loading states present
- Key actions confirmed
- Some undo options
- Validation at submission

### Score 5-6: Acceptable
- Feedback on primary actions
- Basic loading indication
- Minimal confirmations
- Limited undo
- Form validation at submit only

### Score 3-4: Poor
- Feedback often missing
- No loading states
- No confirmations
- No undo
- Errors only on submit

### Score 1-2: Failing
- No feedback on actions
- App seems frozen during loading
- Users don't know if actions worked
- Mistakes permanent
- Users refresh page to check status

### Feedback Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| No loading spinner/state | -3 | Button doesn't change during async |
| No success confirmation | -2 | Action completes silently |
| Silent failures | -3 | Errors not shown |
| No undo for destructive actions | -2 | Delete is permanent |
| No progress on long operations | -1 | Just a spinner for 30+ seconds |
| Stale data shown | -1 | Page doesn't refresh |

---

## Dimension 6: Error Handling (15%)

*Are errors prevented, and when they occur, are they helpful?*

### Score 9-10: Exceptional
- Errors prevented through constraints and defaults
- Inline validation before submission
- Error messages specific and actionable
- User input preserved on errors
- Clear recovery path suggested
- Neutral, helpful tone

### Score 7-8: Good
- Key errors prevented
- Validation on important fields
- Error messages explain problem
- Input mostly preserved
- Recovery usually possible

### Score 5-6: Acceptable
- Some error prevention
- Basic validation
- Error messages present but vague
- Input sometimes lost
- Recovery unclear

### Score 3-4: Poor
- No error prevention
- Validation only on submit
- Generic "Error occurred" messages
- Input lost on error
- No recovery guidance

### Score 1-2: Failing
- Errors frequent and unprevented
- No validation
- No error messages (silent failures)
- All input lost
- Users stuck with no guidance

### Error Handling Anti-Patterns

| Pattern | Score Impact | Detection |
|---------|--------------|-----------|
| "An error occurred" (generic) | -3 | No specific message |
| No inline validation | -2 | All errors on submit |
| Error clears form input | -3 | User loses work |
| No suggested fix | -2 | Just shows error, no recovery |
| Blames user ("Invalid input") | -1 | Accusatory tone |
| Red text without icon | -1 | Accessibility issue |

---

## Overall Score Calculation

```
Overall = (Usability × 0.20) +
          (Info Architecture × 0.15) +
          (User Flows × 0.20) +
          (Affordances × 0.15) +
          (Feedback × 0.15) +
          (Error Handling × 0.15)
```

## Pass/Fail Thresholds

| Dimension | Pass Threshold | Critical Threshold |
|-----------|----------------|-------------------|
| Usability | >= 6 | < 4 blocks overall pass |
| Info Architecture | >= 5 | < 3 blocks overall pass |
| User Flows | >= 6 | < 4 blocks overall pass |
| Affordances | >= 5 | < 3 blocks overall pass |
| Feedback | >= 5 | < 3 blocks overall pass |
| Error Handling | >= 5 | < 3 blocks overall pass |

**Overall Pass**: Score >= 8.0 AND no dimension below critical threshold

---

## Quick Reference Card

### What Makes a 10/10 UX?
- **Usability**: Tasks effortless, zero confusion
- **Info Architecture**: Everything exactly where expected
- **User Flows**: Shortest path, clear progress, obvious next steps
- **Affordances**: Every interactive element obviously interactive
- **Feedback**: Instant response on every action
- **Error Handling**: Errors prevented, or if they happen, helpful and recoverable

### What Makes a 5/10 UX?
- **Usability**: Tasks completable with friction
- **Info Architecture**: Findable with effort
- **User Flows**: Extra steps, unclear progress
- **Affordances**: Mostly identifiable interactive elements
- **Feedback**: Basic loading and confirmations
- **Error Handling**: Generic messages, some input preserved

### What Makes a 2/10 UX?
- **Usability**: Can't figure out how to do things
- **Info Architecture**: Lost, can't find anything
- **User Flows**: Stuck, dead ends, can't go back
- **Affordances**: No idea what's clickable
- **Feedback**: Feels broken, no response
- **Error Handling**: Silent failures, lost work, no guidance

---

## Common UX Fixes by Dimension

### Quick Usability Wins
- Change "Submit" to specific action ("Create Account", "Send Message")
- Make primary CTA visually prominent (size, color, position)
- Add helpful placeholder text

### Quick Info Architecture Wins
- Add section headings
- Move important info above the fold
- Use breadcrumbs for deep navigation

### Quick User Flow Wins
- Add step indicator (Step 2 of 3)
- Add "Back" and "Cancel" buttons
- Add success confirmation page/toast

### Quick Affordance Wins
- Add background/border to buttons
- Underline links or use button style
- Add text labels to icons

### Quick Feedback Wins
- Add `disabled` state + spinner to buttons during async
- Add toast/snackbar for success
- Add skeleton loaders

### Quick Error Handling Wins
- Add specific error messages
- Add inline validation on blur
- Preserve form input on error
