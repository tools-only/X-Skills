# Burndown Detection Patterns

Consolidated patterns from `/deslop` (25 patterns) and `/qa` (6 dimensions) for the `/burndown` skill.

## Pattern Categories

| Category | Source | Agent | Count |
|----------|--------|-------|-------|
| Code Slop | /deslop | Code Slop Hunter | 25 patterns |
| Architecture | /qa | Architecture Auditor | 15+ patterns |
| Scalability | /qa | Debt Classifier | 10+ patterns |
| Maintainability | /qa | Debt Classifier | 10+ patterns |

---

## Category 1: Code Slop Patterns (from /deslop)

### Python Patterns (1-9)

**Pattern 1: Unnecessary Defensive Checks**
```python
# SLOP: Type already guarantees this
def process(data: dict[str, Any]) -> None:
    if data is None:
        return
    if not isinstance(data, dict):
        raise TypeError("Expected dict")

# CLEAN: Trust the type signature
def process(data: dict[str, Any]) -> None:
    # ... actual logic
```
- Severity: Medium
- Fix: Remove redundant checks that duplicate type hints

**Pattern 2: Verbose Variable Declarations**
```python
# SLOP
result = some_function()
return result

# CLEAN
return some_function()
```
- Severity: Low
- Fix: Inline single-use variables

**Pattern 3: Over-Commented Code**
```python
# SLOP
# Initialize the user service
user_service = UserService()
# Get the user by ID
user = user_service.get_by_id(user_id)

# CLEAN
user_service = UserService()
user = user_service.get_by_id(user_id)
```
- Severity: Low
- Fix: Remove comments that restate obvious operations

**Pattern 4: Unnecessary Exception Handling**
```python
# SLOP
try:
    result = calculate_total(items)
except Exception as e:
    logger.error(f"Error: {e}")
    raise

# CLEAN
result = calculate_total(items)  # Let exceptions propagate
```
- Severity: Medium
- Fix: Remove try/except that just logs and re-raises

**Pattern 5: Verbose String Formatting**
```python
# SLOP
message = "User {} has {} items".format(user.name, len(items))

# CLEAN
message = f"User {user.name} has {len(items)} items"
```
- Severity: Low
- Fix: Use f-strings consistently

**Pattern 6: Redundant Truthiness Patterns**
```python
# SLOP
if len(items) > 0:
if value == True:
return True if condition else False

# CLEAN
if items:
if value:
return condition
```
- Severity: Low
- Fix: Use Pythonic truthiness

**Pattern 7: Overly Verbose Dict/List Operations**
```python
# SLOP
if key in data.keys():
for i in range(len(items)):

# CLEAN
if key in data:
for i, item in enumerate(items):
```
- Severity: Low
- Fix: Use Pythonic idioms

**Pattern 8: Unnecessary Type Conversions**
```python
# SLOP
str(f"Hello {name}")
list([x for x in items])

# CLEAN
f"Hello {name}"
[x for x in items]
```
- Severity: Low
- Fix: Remove redundant conversions

**Pattern 9: Import Slop**
```python
# SLOP
from typing import Any, Dict, List, Optional, Union  # Only using Optional
from module import function1
from module import function2  # Should be combined

# CLEAN
from typing import Optional
from module import function1, function2
```
- Severity: Low
- Fix: Remove unused imports, consolidate

---

### React/TypeScript Patterns (10-20)

**Pattern 10: Unnecessary Null Guards**
```typescript
// SLOP: Type guarantees user.profile exists
const name = user?.profile?.name ?? "";

// CLEAN
const name = user.profile.name;
```
- Severity: Medium
- Fix: Trust TypeScript types

**Pattern 11: Unnecessary Memoization**
```typescript
// SLOP
const fullName = useMemo(() => `${first} ${last}`, [first, last]);
const handleClick = useCallback(() => setOpen(true), []);

// CLEAN
const fullName = `${first} ${last}`;
const handleClick = () => setOpen(true);
```
- Severity: Medium
- Fix: Remove useMemo/useCallback for cheap operations

**Pattern 12: Verbose State Declarations**
```typescript
// SLOP
const [isLoading, setIsLoading] = useState<boolean>(false);

// CLEAN
const [isLoading, setIsLoading] = useState(false);
```
- Severity: Low
- Fix: Let TypeScript infer types

**Pattern 13: Over-Defensive Props Handling**
```typescript
// SLOP
const Component = ({ items }: Props) => {
    if (!items) return null;
    if (!Array.isArray(items)) return null;
    // ...
};

// CLEAN
const Component = ({ items }: Props) => {
    // Trust your types
};
```
- Severity: Medium
- Fix: Remove runtime checks that duplicate compile-time guarantees

**Pattern 14: Unnecessary Fragment Wrappers**
```typescript
// SLOP
return (
    <>
        <Component />
    </>
);

// CLEAN
return <Component />;
```
- Severity: Low
- Fix: Remove fragments with single child

**Pattern 15: Verbose Event Handler Patterns**
```typescript
// SLOP
onClick={(e) => handleClick(e)}

// CLEAN
onClick={handleClick}
```
- Severity: Low
- Fix: Pass handler directly

**Pattern 16: Type Assertion Abuse**
```typescript
// SLOP
const data = response as any;
const user = data as unknown as User;
// @ts-ignore

// AVOID
```
- Severity: High
- Fix: Fix types properly instead of asserting

**Pattern 17: Inconsistent Async Patterns**
```typescript
// SLOP
async function loadData() {
    return await fetchData();  // Unnecessary await
}

// CLEAN
async function loadData() {
    return fetchData();
}
```
- Severity: Low
- Fix: Remove unnecessary await before return

**Pattern 18: Verbose Conditional Rendering**
```typescript
// SLOP
{condition ? <Component /> : null}

// CLEAN
{condition && <Component />}
```
- Severity: Low
- Fix: Use short-circuit rendering

**Pattern 19: Unnecessary Spread Operations**
```typescript
// SLOP
<Component {...{ prop1, prop2 }} />

// CLEAN
<Component prop1={prop1} prop2={prop2} />
```
- Severity: Low
- Fix: Use explicit props

**Pattern 20: Effect Cleanup Slop**
```typescript
// SLOP
useEffect(() => {
    doSomething();
    return () => {};  // Empty cleanup
}, []);

// CLEAN
useEffect(() => {
    doSomething();
}, []);
```
- Severity: Low
- Fix: Remove empty cleanup functions

---

### Cross-Language Patterns (21-25)

**Pattern 21: Comment Slop**
- `// TODO: implement` with no specifics
- Comments that restate function names
- Commented-out code blocks
- Section dividers: `// ========== HELPERS ==========`
- Severity: Low
- Fix: Remove or make actionable

**Pattern 22: Naming Inconsistencies**
```python
# SLOP: Mixed conventions
def getUserById():  # camelCase
def get_user_by_email():  # snake_case
```
- Severity: Medium
- Fix: Enforce consistent naming per language

**Pattern 23: Dead Code**
- Unreachable code after return
- Conditions always true/false
- Unused parameters
- Unused imports
- Functions never called
- Severity: Medium
- Fix: Delete

**Pattern 24: Logging/Console Slop**
```typescript
// SLOP
console.log("Component mounted");
console.log("Props:", props);
```
- Severity: Low
- Fix: Remove debug logging

**Pattern 25: Error Message Slop**
```python
# SLOP
raise ValueError("An error occurred while processing the data")

# CLEAN
raise ValueError(f"Invalid user_id: {user_id}")
```
- Severity: Low
- Fix: Make error messages specific

---

## Category 2: Architecture Patterns (from /qa)

### Size & Complexity

**A1: Oversized Files**
| Threshold | Severity |
|-----------|----------|
| > 500 lines | Critical |
| > 300 lines | High |
| > 200 lines (React) | Medium |

Fix: Split into focused modules

**A2: Long Functions**
- > 50 lines: High
- Cyclomatic complexity > 10: High
- Nesting depth > 4: Medium

Fix: Extract helper functions

**A3: Too Many Exports**
- > 10 public exports: Medium

Fix: Split module or create barrel file

---

### Separation of Concerns

**A4: Frontend Layer Violations**
- Components containing business logic + data fetching + presentation
- Components directly calling APIs
- UI state mixed with server state
- Severity: High

Fix: Separate into hooks/services/components

**A5: Backend Layer Violations**
- Route handlers with business logic
- Services directly accessing DB
- Domain objects knowing about HTTP
- Severity: High

Fix: Introduce service/repository layers

**A6: Cross-Domain Imports**
- Feature importing from unrelated feature
- Shared folders with domain-specific code
- Severity: Medium

Fix: Extract to shared abstractions

---

### Coupling & Abstraction

**A7: Circular Dependencies**
- Any import cycle
- Severity: Critical

Fix: Extract interface, dependency inversion

**A8: Deep Imports**
- `../../../other/internal/private`
- Severity: Medium

Fix: Use public API / barrel exports

**A9: Leaky Abstractions**
- DB schemas exposed across boundaries
- API response shapes in domain code
- Severity: High

Fix: Map to domain types at boundaries

**A10: Over-Abstraction**
- Interfaces with only one implementation
- Unnecessary indirection
- Severity: Medium

Fix: Inline or simplify

**A11: Hardcoded Values**
- URLs, credentials, feature flags in code
- Magic numbers without constants
- Severity: High

Fix: Move to configuration

**A12: Duplication**
- 2+ instances of copy-pasted logic
- Near-duplicate functions
- Severity: Medium

Fix: Extract shared function

---

## Category 3: Scalability Patterns (from /qa)

**S1: N+1 Query Pattern**
```python
# SLOP
for user in users:
    orders = get_orders(user.id)  # Query in loop!
```
- Severity: Critical
- Fix: Use prefetch_related / eager loading

**S2: Missing Pagination**
- List endpoints without limits
- Queries returning unbounded results
- Severity: High
- Fix: Add pagination parameters

**S3: Sequential Awaits**
```typescript
// SLOP
const users = await getUsers();
const orders = await getOrders();

// CLEAN
const [users, orders] = await Promise.all([getUsers(), getOrders()]);
```
- Severity: Medium
- Fix: Parallelize with Promise.all / asyncio.gather

**S4: Missing Timeouts**
- External calls without timeout
- Severity: High
- Fix: Add timeout configuration

**S5: Unbounded Data Structures**
- Caches without TTL
- Lists that grow without limit
- Severity: High
- Fix: Add limits and eviction

**S6: Stale Cache Reads**
- Writes without cache invalidation
- Severity: Medium
- Fix: Invalidate on write

**S7: Unbounded Rendering**
- Long lists without virtualization
- Severity: Medium (frontend)
- Fix: Use react-virtual or similar

---

## Category 4: Maintainability Patterns (from /qa)

**M1: Consistency Violations**
- Mixed patterns for same problem
- Inconsistent naming conventions
- Inconsistent file organization
- Severity: Medium
- Fix: Standardize on one pattern

**M2: Dead Code (Exhaustive)**
- Unused exports
- Unused functions/classes
- Unreachable branches
- Unused dependencies
- Dead feature folders
- Severity: Medium
- Fix: Delete

**M3: Type Safety Gaps**
- Python: `Any`, missing types, `# type: ignore`
- TypeScript: `any`, `as`, `@ts-ignore`
- Severity: High
- Fix: Add proper types

**M4: Error Handling Gaps**
- Empty catch blocks
- Generic exception catching
- Missing error boundaries
- Swallowed errors
- Severity: High
- Fix: Handle or propagate properly

**M5: Testability Blockers**
- Hard dependencies on singletons
- Side effects in constructors
- Global mutable state
- Severity: Medium
- Fix: Dependency injection

---

## Severity Summary

| Severity | Action Required | Examples |
|----------|-----------------|----------|
| **Critical** | MUST fix to pass checkpoint | Circular deps, N+1, 500+ line files |
| **High** | Should fix | Type assertions, layer violations, hardcoding |
| **Medium** | Nice to fix | Memoization abuse, consistency, duplication |
| **Low** | Optional | Verbose patterns, comments, formatting |
