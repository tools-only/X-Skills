# Contextvar Remediation Skill

**Use this skill when**: Converting explicit user_id parameters to ambient context using Python's contextvars.

**Critical Success Factors**:
1. Trust but verify - Always validate conclusions against the decision guide
2. Understand the architecture - Know how services are instantiated (singleton vs per-request)
3. Be systematic - Update services first, then callsites
4. Clean code - Omit optional parameters rather than passing None explicitly

---

## Investigation Phase

### Step 1: Understand Current Usage Patterns

Before making changes, **thoroughly investigate** how the code currently works:

```bash
# Find all explicit user_id usage
grep -r "user_id: Optional" target_directory/
grep -r "user_id: str" target_directory/
grep -r "user_id=" target_directory/ | wc -l
```

**Key Questions to Answer**:
- How are services instantiated? (Singleton factory? Per-request? Per-user?)
- Where is `set_current_user_id()` already being called?
- Which code runs during active requests vs background tasks?
- Does a decision guide exist? (e.g., CONTEXTVAR_DECISION_GUIDE.md)

### Step 2: Read and Internalize the Decision Guide

**DO NOT SKIP THIS STEP**. Read the project's contextvar decision guide completely.

**Critical patterns to identify**:
- ❌ **Anti-pattern**: `user_id: Optional[str] = None` on every method (method-level dual-mode)
- ✅ **Correct**: Constructor-only dual-mode OR pure ambient context
- ❌ **Anti-pattern**: Setting context then still passing user_id explicitly
- ✅ **Correct**: Set context once, let it propagate

**The guide will tell you**:
- Which pattern to use (constructor vs pure ambient)
- How to handle background tasks vs active flows
- Special cases (database connections, infrastructure boundaries)

### Step 3: Categorize Code by Context Type

Classify each file/method into:

1. **Active Request Flow**: HTTP endpoints, WebSocket handlers, tool execution
   - **Strategy**: Use pure ambient context
   - **Entry point**: Auth layer sets `set_current_user_id()`

2. **Background Tasks**: Scheduled jobs, batch processing, admin scripts
   - **Strategy**: Explicit iteration with context setting per user
   - **Pattern**: `for user in users: set_current_user_id(user_id); service.method()`

3. **Internal Services**: Shared infrastructure called by both above
   - **Strategy**: Depends on instantiation pattern (see Step 4)

4. **Infrastructure Boundary**: Database connections, cache clients
   - **Strategy**: Keep explicit (may have connection pooling concerns)

### Step 4: Identify Service Instantiation Pattern

**Critical decision point**: How are services created?

#### Pattern A: Singleton Factory (created at startup)
```python
# main.py
factory = ServiceFactory()  # Created once, no user context yet
```

**Implication**: Constructor-only dual-mode **DOESN'T WORK** because:
- Services are created before any user context exists
- Constructor can't call `get_current_user_id()` (would fail)
- Methods would still need to fetch context anyway

**Solution**: Pure ambient context in all methods.

#### Pattern B: Per-Request Instantiation
```python
# request handler
service = SomeService(user_id=get_current_user_id())
```

**Implication**: Constructor-only dual-mode **CAN WORK** because:
- Services are created during request with context available
- Constructor sets `self.user_id` once
- Methods use `self.user_id` (no parameters needed)

**Solution**: Constructor-only pattern acceptable.

#### Pattern C: Mixed Usage
Some methods called from both active and background contexts.

**Solution**: Background tasks must call `set_current_user_id()` before service calls.

---

## Implementation Phase

### Step 5: Update Service Layer (Services First, Callsites Later)

**Order matters**: Update service method signatures BEFORE updating callsites to see what breaks.

#### For Pure Ambient Context Pattern:

```python
# Before
class MemoryService:
    def find_memories(self, query: str, user_id: Optional[str] = None):
        memories = self.db.get_memories(query, user_id=user_id)
        return memories

# After
class MemoryService:
    def find_memories(self, query: str):
        memories = self.db.get_memories(query)  # db layer handles context
        return memories
```

**Key changes**:
1. Remove `user_id: Optional[str] = None` from method signature
2. Remove `user_id` from docstring Args section
3. Remove `user_id=...` from internal service calls
4. **Important**: Just omit the parameter, don't pass `user_id=None`

#### For Constructor-Only Pattern:

```python
# Before
class MemoryService:
    def find_memories(self, query: str, user_id: Optional[str] = None):
        memories = self.db.get_memories(query, user_id=user_id)
        return memories

# After
class MemoryService:
    def __init__(self, ..., user_id: Optional[str] = None):
        self.user_id = user_id or get_current_user_id()

    def find_memories(self, query: str):
        memories = self.db.get_memories(query, user_id=self.user_id)
        return memories
```

#### Methods That Need Explicit Context Reads

Some methods need user_id for logging or business logic:

```python
def run_full_process(self) -> Dict[str, int]:
    """Run process for current user."""
    from utils.user_context import get_current_user_id
    user_id = get_current_user_id()  # For logging only

    stats = self.do_work()
    logger.info(f"Completed for user {user_id}: {stats}")
    return stats
```

**Only call `get_current_user_id()` when actually needed**, not speculatively.

### Step 6: Update Infrastructure/Database Layer

**The infrastructure boundary typically keeps Optional[str]**:

```python
class DatabaseLayer:
    def _resolve_user_id(self, user_id: Optional[str] = None) -> str:
        """Resolve user_id from explicit param or ambient context."""
        if user_id is not None:
            return user_id
        from utils.user_context import get_current_user_id
        return get_current_user_id()

    def get_memories(self, query: str, user_id: Optional[str] = None):
        resolved_user_id = self._resolve_user_id(user_id)
        # Use resolved_user_id for actual query
```

**Why keep Optional here?**
- Backward compatibility during migration
- Explicit control for background tasks
- Database connection pooling may require explicit user_id

**Services call with omitted parameter**:
```python
# Service layer
memories = self.db.get_memories(query)  # Not user_id=None, just omit it
```

### Step 7: Update Callsites Systematically

**Use replace_all for mechanical changes**:

```python
# Find all callsites
grep -r "self\.service\.method.*user_id=" .

# Update pattern
# Before: service.find_memories(query, user_id=user_id)
# After:  service.find_memories(query)
```

**For background tasks, ensure context is set**:

```python
# Before
for user in users:
    user_id = str(user["id"])
    results = service.process(user_id=user_id)

# After
for user in users:
    user_id = str(user["id"])
    set_current_user_id(user_id)  # Set ambient context
    results = service.process()    # Uses ambient context
```

**Critical**: Background tasks MUST call `set_current_user_id()` at the start of each user iteration.

### Step 8: Clean Up Explicit None Passing

After initial updates, search for explicit `user_id=None`:

```python
# Verbose (works but cluttered)
memories = self.db.get_memories(query, user_id=None)

# Clean (preferred)
memories = self.db.get_memories(query)
```

**Why clean it up?**
- If the parameter has a default of `None`, omitting it is clearer
- Shows intent: "use the default behavior" vs "explicitly passing None"
- Reduces visual noise

---

## Validation Phase

### Step 9: Verify Background Tasks Still Work

**Test that scheduled jobs handle user iteration correctly**:

```python
# Check each background task file
grep -A 10 "set_current_user_id" scheduled_tasks.py

# Verify pattern:
# 1. Loop over users
# 2. Call set_current_user_id(user_id)
# 3. Call service methods (no user_id arg)
```

**If set_current_user_id is missing**, background tasks will fail with:
```
RuntimeError: No user context set. Ensure authentication is properly initialized.
```

This is **correct fail-fast behavior** - it surfaces the bug immediately.

### Step 10: Verify Active Flows Inherit Context

**Check that request handlers set context at entry**:

```python
# API/WebSocket handlers should have:
@route("/endpoint")
def handler(request):
    user_id = authenticate(request)
    set_current_user_id(user_id)  # Set once at entry

    # All downstream service calls inherit context automatically
    result = service.do_work()
    return result
```

**The beauty of ambient context**:
- Set once at the entry point
- Propagates through entire call stack automatically
- No manual threading needed

### Step 11: Count Your Wins

**Measure the improvement**:

```bash
# Before remediation
grep -r "user_id: Optional" | wc -l
grep -r "user_id=" | wc -l

# After remediation
grep -r "user_id: Optional" | wc -l   # Should be much lower
grep -r "user_id=" | wc -l             # Should be much lower
```

**Report the impact**:
- X method signatures simplified
- Y explicit parameters removed
- Z callsites updated
- Module now uses consistent ambient context pattern

---

## Common Pitfalls and Solutions

### Pitfall 1: Trusting Initial Analysis Without Verification

**Symptom**: Agent or assistant concludes "NO CHANGES NEEDED" or claims current pattern is optimal.

**Solution**:
1. Read the decision guide yourself
2. Check if current pattern matches guide recommendations
3. Look for anti-patterns explicitly called out in guide
4. **Trust but verify** - even if analysis seems thorough

**Example from real work**:
```
Agent: "lt_memory uses Optional[str] throughout - this is the gold standard"
Reality: Decision guide explicitly says "Method-level dual-mode is almost always wrong"
```

### Pitfall 2: Not Understanding Service Instantiation

**Symptom**: Applying constructor-only pattern to singleton services created at startup.

**Solution**:
1. Check where services are created (main.py, factory, etc.)
2. Determine if user context exists at creation time
3. If services are singletons: Use pure ambient in methods
4. If per-request: Constructor-only pattern OK

**Example**:
```python
# ❌ Wrong: Singleton service with constructor dual-mode
class Service:
    def __init__(self, user_id: Optional[str] = None):
        self.user_id = user_id or get_current_user_id()  # FAILS - no context at startup!

# ✅ Right: Methods fetch from ambient when needed
class Service:
    def find_data(self, query: str):
        # Context set by caller (request handler or background loop)
        return self.db.query(query)  # db layer reads ambient context
```

### Pitfall 3: Removing user_id from Background Task Methods

**Symptom**: Background tasks that iterate over users fail because methods can't identify which user.

**Solution**:
- Background tasks **call `set_current_user_id()`** before calling service methods
- Service methods use ambient context (no explicit parameter needed)
- Each iteration sets context for that user

**Pattern**:
```python
# Scheduled job
def run_for_all_users():
    for user in get_users():
        set_current_user_id(str(user["id"]))  # ← Critical
        service.process()  # ← No user_id parameter
```

### Pitfall 4: Passing user_id=None Explicitly

**Symptom**: After remediation, code has many `user_id=None` arguments.

**Solution**: Omit the parameter entirely if default is `None`.

```python
# ❌ Verbose
result = db.query(sql, user_id=None)

# ✅ Clean
result = db.query(sql)
```

### Pitfall 5: Forgetting Import Statements

**Symptom**: `NameError: name 'get_current_user_id' is not defined`

**Solution**: Add import to files that now need to read context:

```python
from utils.user_context import get_current_user_id
```

**Files that typically need the import**:
- Services that log user_id
- Methods that call `get_current_user_id()` for business logic
- Background task files that set context

---

## Red Flags to Watch For

### 🚩 Red Flag: "Optional[str] = None on Every Method"

If every method in a service accepts `user_id: Optional[str] = None`, this is **method-level dual-mode** - an anti-pattern.

**Why it's wrong**:
- Pollutes every method signature unnecessarily
- Doesn't enforce context setting (silent failures)
- Creates ambiguity: use parameter or ambient context?

**Solution**: Remove from all methods, use pure ambient context.

### 🚩 Red Flag: Setting Context Then Passing Explicitly

```python
set_current_user_id(user_id)
result = service.method(user_id=user_id)  # ← Redundant!
```

**Why it's wrong**: You've already set ambient context, just trust it.

**Solution**: Remove the explicit parameter.

### 🚩 Red Flag: Convenience Wrapper Functions

```python
def find_memories_for_current_user(query: str):
    user_id = get_current_user_id()
    return service.find_memories(query, user_id)
```

**Why it's wrong**: If the underlying service needs wrappers, the service API is wrong.

**Solution**: Fix the service to use ambient context, delete wrappers.

### 🚩 Red Flag: Request Body Contains user_id

```python
class RequestModel(BaseModel):
    user_id: str  # ← Security risk!
    data: Dict
```

**Why it's wrong**: Users shouldn't specify their own user_id - it comes from authentication.

**Solution**: Remove from request model, use authenticated session's user_id.

---

## Decision Tree

```
Is this code called during active user requests?
├─ YES: Use pure ambient context
│  └─ Remove user_id from all method signatures
│
└─ NO: Is this a background task that processes multiple users?
   ├─ YES: Keep iteration explicit, set context per user
   │  └─ for user in users: set_current_user_id(user_id); service.method()
   │
   └─ NO: Is this infrastructure boundary (DB, cache)?
      ├─ YES: Keep Optional[str], implement _resolve_user_id() helper
      └─ NO: Re-evaluate - might be misclassified
```

---

## Expected Outcomes

After successful remediation:

✅ **Service layer**: No `user_id` parameters on public methods
✅ **Background tasks**: Explicit `set_current_user_id()` before service calls
✅ **Active flows**: Context set once at entry point, propagates automatically
✅ **Infrastructure**: Keeps Optional[str] with fallback to ambient context
✅ **Tests**: Updated fixtures to set context before calling services
✅ **Callsites**: All `user_id=...` arguments removed from service calls
✅ **Fail-fast**: Missing context raises immediately rather than silently failing

**Code Impact**:
- Fewer parameters (cleaner signatures)
- Less parameter threading (simpler call chains)
- Consistent pattern (one way to handle user context)
- Better security (can't pass wrong user_id)

---

## Checklist for Skill Execution

Before considering remediation complete:

- [ ] Read project's contextvar decision guide completely
- [ ] Identified all explicit user_id usage patterns (grep counts)
- [ ] Determined service instantiation pattern (singleton vs per-request)
- [ ] Chose correct remediation pattern (pure ambient vs constructor-only)
- [ ] Updated service layer method signatures (removed user_id params)
- [ ] Updated infrastructure layer (added _resolve_user_id if needed)
- [ ] Updated all callsites (removed user_id arguments)
- [ ] Cleaned up explicit user_id=None passing
- [ ] Verified background tasks set context before service calls
- [ ] Verified active flows inherit context from auth layer
- [ ] Added imports where needed (get_current_user_id, set_current_user_id)
- [ ] Tested that missing context fails fast with clear error
- [ ] Counted the impact (X signatures, Y parameters, Z callsites improved)
- [ ] Verified no red flags remain (check list above)

---

## Summary

**The Core Principle**:
- Active request flows: Set context once at entry, trust propagation
- Background tasks: Set context per user in iteration loop
- Services: Pure ambient context (no user_id parameters)
- Infrastructure: Optional with fallback to ambient

**The Critical Insight**:
If you need singleton services but apply constructor-only pattern, you've created an unsolvable problem. Constructor runs before context exists. This is why you must understand instantiation timing.

**The Success Formula**:
1. Read the guide (don't skip)
2. Understand the architecture (singleton vs per-request)
3. Choose the right pattern (pure ambient for singletons)
4. Update systematically (services first, callsites second)
5. Verify behavior (background tasks set context, active flows inherit)
6. Clean up (omit parameters, don't pass None)

**Trust the Process**: Contextvars work. The runtime enforces isolation. If you set context before calling services, ambient context **will** be available. Don't hedge with Optional parameters "just in case" - that defeats the purpose.
