# Debugging Rubric - Failure Pattern Taxonomy

Reference for identifying and resolving common failure patterns in web applications.

## Severity Levels

| Level | Description | Action |
|-------|-------------|--------|
| **P0 - Critical** | Service completely down | Fix immediately, all other work stops |
| **P1 - High** | Major feature broken | Fix within current iteration |
| **P2 - Medium** | Minor feature broken | Fix if time permits |
| **P3 - Low** | Cosmetic/UX issue | Log for later |

## Failure Categories

### 1. Authentication Failures

#### Symptoms
- 401 Unauthorized responses
- Sudden logout
- "Session expired" errors
- Redirect loops to login page

#### Common Causes
| Pattern | Root Cause | Fix |
|---------|------------|-----|
| 401 on all requests | Token expired/missing | Check token refresh logic |
| 401 after redirect | Auth header lost in redirect | Use 308 instead of 307, or preserve headers |
| Logout cascade | 401 triggers `clearToken()` | Add retry before clearing |
| CORS + auth | Preflight fails, no auth header | Fix CORS config |

#### Log Patterns
```
# Azure/Backend
"Authentication failed" | "Invalid token" | "Token expired"
"JWT decode error" | "Signature verification failed"

# Frontend Console
"401" | "Unauthorized" | "clearToken" | "logout"
```

#### Diagnostic Steps
1. Check if token exists in request headers
2. Verify token hasn't expired (decode JWT, check `exp`)
3. Check for redirect responses (307/308) losing headers
4. Verify CORS allows Authorization header

---

### 2. Network/Connectivity Failures

#### Symptoms
- Connection refused
- Timeout errors
- CORS errors
- WebSocket disconnection

#### Common Causes
| Pattern | Root Cause | Fix |
|---------|------------|-----|
| Connection refused | Service not running | Check deployment status |
| Timeout | Service overloaded or hung | Check resource limits, scaling |
| CORS error | Missing/wrong CORS config | Add origin to allowed list |
| WS disconnect | Firewall/proxy blocking | Check WebSocket URL config |

#### Log Patterns
```
# Azure/Backend
"Connection refused" | "ETIMEDOUT" | "ECONNRESET"
"Max retries exceeded" | "Connection pool exhausted"

# Frontend Console
"Failed to fetch" | "NetworkError" | "CORS" | "blocked"
"WebSocket connection failed" | "wss://" | "ws://"
```

#### Diagnostic Steps
1. Check service health endpoints directly
2. Verify DNS resolution
3. Check firewall/security group rules
4. Verify WebSocket URL isn't falling back to localhost

---

### 3. Configuration/Environment Failures

#### Symptoms
- Feature works locally but not in staging/prod
- Requests going to wrong URL
- Missing environment variable errors
- Fallback to localhost URLs

#### Common Causes
| Pattern | Root Cause | Fix |
|---------|------------|-----|
| localhost in prod | Empty env var → fallback | Set env var properly |
| Wrong API URL | Env var mismatch | Verify NEXT_PUBLIC_* vars |
| Feature flag off | Config not deployed | Check deployment config |
| Secret missing | Not in production env | Add to Azure/deployment |

#### Log Patterns
```
# Azure/Backend
"KeyError" | "Missing environment variable" | "Config not found"
"undefined" (in URL construction)

# Frontend Console
"localhost:3000" | "localhost:8000" (in production!)
"undefined" | "null" (in API calls)
```

#### Diagnostic Steps
1. Check Network tab for localhost URLs
2. Grep code for fallback patterns: `|| 'http://localhost'`
3. Verify env vars are set: `printenv | grep NEXT_PUBLIC`
4. Check for undefined in URL construction

---

### 4. Data/State Failures

#### Symptoms
- Stale data displayed
- "Cannot read property of undefined"
- Type errors in data processing
- Empty lists when data expected

#### Common Causes
| Pattern | Root Cause | Fix |
|---------|------------|-----|
| undefined.property | Null check missing | Add optional chaining |
| Stale data | Cache not invalidated | Force cache refresh |
| Type mismatch | API response changed | Update type definitions |
| Empty result | Query filter too strict | Check query parameters |

#### Log Patterns
```
# Backend
"TypeError" | "KeyError" | "AttributeError"
"NoneType has no attribute" | "undefined is not a function"

# Frontend Console
"Cannot read property" | "undefined" | "null"
"TypeError: x is not a function"
```

#### Diagnostic Steps
1. Check API response shape in Network tab
2. Verify TypeScript types match API
3. Add console.log before error point
4. Check for race conditions in data loading

---

### 5. Deployment/Infrastructure Failures

#### Symptoms
- Container won't start
- Health check fails
- Memory/CPU limits exceeded
- Cold start timeouts

#### Common Causes
| Pattern | Root Cause | Fix |
|---------|------------|-----|
| Container crash | OOM or unhandled exception | Check limits, add error handling |
| Health fail | App not ready on startup | Increase startup probe timeout |
| Resource exhausted | Undersized container | Increase limits |
| Cold start | Too slow to start | Optimize startup, add warmup |

#### Log Patterns
```
# Azure
"OOMKilled" | "Container terminated"
"Liveness probe failed" | "Readiness probe failed"
"Resource limit exceeded" | "CPU throttled"
```

#### Diagnostic Steps
1. Check container logs: `az containerapp logs show`
2. Check resource metrics in Azure portal
3. Verify health endpoint is fast (< 10s)
4. Check for blocking startup operations

---

### 6. Database/External Service Failures

#### Symptoms
- Slow responses
- "Connection pool exhausted"
- Intermittent 500 errors
- Data inconsistency

#### Common Causes
| Pattern | Root Cause | Fix |
|---------|------------|-----|
| Pool exhausted | Connections not released | Add connection timeout |
| Slow queries | Missing index or N+1 | Optimize query, add index |
| Inconsistency | Race condition | Add transaction/locking |
| Intermittent 500 | External service flaky | Add retry logic |

#### Log Patterns
```
# Backend
"Connection pool exhausted" | "Too many connections"
"Query timeout" | "Lock wait timeout"
"External service error" | "Retry attempt"
```

#### Diagnostic Steps
1. Check database metrics (connections, query time)
2. Enable query logging to find slow queries
3. Check external service status pages
4. Add request tracing (LogFire spans)

---

## Correlation Checklist

When multiple services have errors, use this checklist to find root cause:

### Time Correlation
- [ ] When did errors start? Check deployment history
- [ ] Do errors correlate with traffic spikes?
- [ ] Are errors periodic (cron jobs, batch processes)?

### Service Correlation
- [ ] Which service failed first?
- [ ] Are downstream services failing due to upstream?
- [ ] Is there a shared dependency (database, cache)?

### Change Correlation
- [ ] What was deployed recently?
- [ ] Were environment variables changed?
- [ ] Were external services updated?

### Pattern Correlation
- [ ] Same error across multiple services = shared dependency
- [ ] Cascade of 401s = auth service issue
- [ ] Timeouts everywhere = network/infrastructure issue

## Quick Reference: Error → Action

| Error Type | First Action |
|------------|--------------|
| 401 Unauthorized | Check token, check redirects |
| 403 Forbidden | Check permissions, CORS |
| 404 Not Found | Check routes, URL construction |
| 500 Internal | Check backend logs, stack trace |
| 502 Bad Gateway | Check upstream service health |
| 503 Unavailable | Check container status, scaling |
| Connection refused | Check service is running |
| CORS error | Check CORS config, credentials |
| undefined/null | Check data flow, add null checks |
| Timeout | Check service load, query performance |
