---
name: langfuse-incident-runbook
description: |
  Troubleshoot and respond to Langfuse-related incidents and outages.
  Use when experiencing Langfuse outages, debugging production issues,
  or responding to LLM observability incidents.
  Trigger with phrases like "langfuse incident", "langfuse outage",
  "langfuse down", "langfuse production issue", "langfuse troubleshoot".
allowed-tools: Read, Write, Edit, Bash(curl:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Incident Runbook

## Overview
Step-by-step procedures for responding to Langfuse-related incidents.

## Prerequisites
- Access to Langfuse dashboard
- Application logs access
- Metrics/monitoring dashboards
- Escalation contacts

## Incident Severity Levels

| Severity | Description | Response Time | Escalation |
|----------|-------------|---------------|------------|
| P1 | Complete outage, no traces | 15 min | Immediate |
| P2 | Degraded, partial data loss | 1 hour | 4 hours |
| P3 | Slow/delayed traces | 4 hours | Next business day |
| P4 | Minor issues, workaround exists | 24 hours | Best effort |

## Quick Diagnostics

### Step 1: Initial Assessment (2 minutes)

```bash
#!/bin/bash
# quick-diagnosis.sh

echo "=== Langfuse Quick Diagnosis ==="
echo "Time: $(date)"
echo ""

# 1. Check Langfuse status
echo "1. Langfuse Status:"
curl -s https://status.langfuse.com/api/v2/status.json | jq '.status.description'

# 2. Check API connectivity
echo ""
echo "2. API Connectivity:"
curl -s -o /dev/null -w "HTTP %{http_code} in %{time_total}s\n" \
  https://cloud.langfuse.com/api/public/health

# 3. Check authentication
echo ""
echo "3. Auth Test:"
AUTH=$(echo -n "$LANGFUSE_PUBLIC_KEY:$LANGFUSE_SECRET_KEY" | base64)
curl -s -o /dev/null -w "HTTP %{http_code}\n" \
  -H "Authorization: Basic $AUTH" \
  "https://cloud.langfuse.com/api/public/traces?limit=1"

# 4. Check application health
echo ""
echo "4. Application Metrics:"
curl -s http://localhost:3000/api/metrics | grep langfuse | head -5
```

### Step 2: Determine Incident Type

| Symptom | Likely Cause | Go To |
|---------|--------------|-------|
| No traces appearing | SDK not flushing | Section A |
| 401/403 errors | Authentication issue | Section B |
| High latency | Network/rate limits | Section C |
| Missing data | Partial failures | Section D |
| Complete outage | Langfuse service issue | Section E |

---

## Section A: Traces Not Appearing

### Symptoms
- Dashboard shows no new traces
- No errors in application logs
- Application functioning normally

### Diagnosis Steps

```typescript
// 1. Verify SDK is enabled
console.log("Langfuse enabled:", process.env.LANGFUSE_ENABLED !== "false");
console.log("Environment:", process.env.NODE_ENV);

// 2. Check for pending events
// Add this to your code temporarily
const langfuse = getLangfuse();
console.log("Pending events:", langfuse.pendingItems?.length || "unknown");

// 3. Force flush and check for errors
try {
  await langfuse.flushAsync();
  console.log("Flush successful");
} catch (error) {
  console.error("Flush failed:", error);
}
```

### Resolution Steps

1. **Check shutdown handlers**
   ```typescript
   // Ensure shutdown is registered
   process.on("beforeExit", async () => {
     await langfuse.shutdownAsync();
   });
   ```

2. **Reduce batch size temporarily**
   ```typescript
   const langfuse = new Langfuse({
     flushAt: 1,        // Immediate flush
     flushInterval: 1000,
   });
   ```

3. **Enable debug logging**
   ```bash
   DEBUG=langfuse* npm start
   ```

---

## Section B: Authentication Errors

### Symptoms
- 401 Unauthorized errors
- 403 Forbidden errors
- "Invalid API key" messages

### Diagnosis Steps

```bash
# 1. Verify environment variables
echo "Public key starts with: ${LANGFUSE_PUBLIC_KEY:0:10}"
echo "Secret key is set: ${LANGFUSE_SECRET_KEY:+yes}"
echo "Host: ${LANGFUSE_HOST:-https://cloud.langfuse.com}"

# 2. Test credentials directly
curl -v -X GET \
  -H "Authorization: Basic $(echo -n "$LANGFUSE_PUBLIC_KEY:$LANGFUSE_SECRET_KEY" | base64)" \
  "${LANGFUSE_HOST:-https://cloud.langfuse.com}/api/public/traces?limit=1"
```

### Resolution Steps

1. **Verify keys match project**
   - Go to Langfuse Dashboard > Settings > API Keys
   - Ensure keys are from the correct project
   - Check keys haven't been revoked

2. **Check for key rotation**
   - If keys were recently rotated, update all environments
   - Verify secret manager has latest values

3. **Verify host URL**
   - Cloud: `https://cloud.langfuse.com`
   - Self-hosted: Your instance URL (no trailing slash)

---

## Section C: High Latency / Timeouts

### Symptoms
- Slow API responses
- Request timeouts
- 429 Rate limit errors

### Diagnosis Steps

```typescript
// Check flush timing
const start = Date.now();
await langfuse.flushAsync();
console.log(`Flush took ${Date.now() - start}ms`);

// Check batch sizes
console.log("Current batch size:", langfuse.pendingItems?.length);
```

```bash
# Network latency test
curl -w "DNS: %{time_namelookup}s\nConnect: %{time_connect}s\nTLS: %{time_appconnect}s\nTotal: %{time_total}s\n" \
  -o /dev/null -s https://cloud.langfuse.com/api/public/health
```

### Resolution Steps

1. **For rate limits**
   ```typescript
   // Increase batching
   const langfuse = new Langfuse({
     flushAt: 50,
     flushInterval: 10000,
   });
   ```

2. **For network issues**
   - Check firewall rules allow outbound HTTPS
   - Verify DNS resolution
   - Consider using a closer region (self-hosted)

3. **Implement circuit breaker**
   ```typescript
   class CircuitBreaker {
     private failures = 0;
     private lastFailure?: Date;
     private readonly threshold = 5;
     private readonly resetMs = 60000;

     async execute<T>(operation: () => Promise<T>): Promise<T | null> {
       if (this.isOpen()) {
         console.warn("Circuit breaker open, skipping Langfuse");
         return null;
       }

       try {
         const result = await operation();
         this.reset();
         return result;
       } catch (error) {
         this.recordFailure();
         throw error;
       }
     }

     private isOpen(): boolean {
       if (this.failures < this.threshold) return false;
       if (!this.lastFailure) return false;
       return Date.now() - this.lastFailure.getTime() < this.resetMs;
     }

     private recordFailure() {
       this.failures++;
       this.lastFailure = new Date();
     }

     private reset() {
       this.failures = 0;
       this.lastFailure = undefined;
     }
   }
   ```

---

## Section D: Missing/Partial Data

### Symptoms
- Some traces appear, others don't
- Missing spans or generations
- Incomplete trace data

### Diagnosis Steps

```typescript
// Check for errors in trace operations
const trace = langfuse.trace({ name: "test" });
console.log("Trace ID:", trace.id);

const span = trace.span({ name: "test-span" });
console.log("Span ID:", span.id);

// Verify end() is called
span.end({ output: { test: true } });
console.log("Span ended");

await langfuse.flushAsync();
console.log("Flushed");
```

### Resolution Steps

1. **Ensure all spans are ended**
   ```typescript
   const span = trace.span({ name: "operation" });
   try {
     return await doWork();
   } finally {
     span.end(); // Always end in finally
   }
   ```

2. **Check for exceptions swallowing**
   ```typescript
   try {
     await langfuse.flushAsync();
   } catch (error) {
     console.error("Langfuse flush error:", error);
     // Don't swallow - log for debugging
   }
   ```

---

## Section E: Langfuse Service Outage

### Symptoms
- status.langfuse.com shows issues
- All API calls failing
- Multiple users affected

### Immediate Actions

1. **Check status page**: https://status.langfuse.com

2. **Enable fallback mode**
   ```typescript
   // Graceful degradation
   const langfuse = new Langfuse({
     enabled: false, // Disable during outage
   });
   ```

3. **Queue events locally**
   ```typescript
   // Store events to file during outage
   const pendingEvents: any[] = [];

   function queueEvent(event: any) {
     pendingEvents.push({
       ...event,
       timestamp: new Date().toISOString(),
     });

     if (pendingEvents.length > 1000) {
       // Write to file
       fs.writeFileSync(
         `langfuse-backup-${Date.now()}.json`,
         JSON.stringify(pendingEvents)
       );
       pendingEvents.length = 0;
     }
   }
   ```

4. **Monitor for recovery**
   ```bash
   # Watch status
   watch -n 30 'curl -s https://status.langfuse.com/api/v2/status.json | jq .status'
   ```

---

## Post-Incident Checklist

- [ ] Verify traces are appearing in dashboard
- [ ] Check no data was lost during incident
- [ ] Review error rates returning to normal
- [ ] Update incident documentation
- [ ] Schedule post-mortem if P1/P2
- [ ] Update runbook with learnings

## Escalation Contacts

| Level | Contact | When |
|-------|---------|------|
| L1 | On-call engineer | All incidents |
| L2 | Platform team lead | P1/P2 unresolved 30min |
| L3 | Langfuse support | Service-side issues |

## Resources
- [Langfuse Status](https://status.langfuse.com)
- [Langfuse Discord](https://langfuse.com/discord)
- [Langfuse GitHub Issues](https://github.com/langfuse/langfuse/issues)

## Next Steps
For data export and retention, see `langfuse-data-handling`.
