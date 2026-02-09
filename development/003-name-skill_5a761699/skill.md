---
name: twinmind-incident-runbook
description: |
  Incident response procedures for TwinMind integration failures.
  Use when experiencing outages, debugging production issues,
  or responding to alerts related to TwinMind.
  Trigger with phrases like "twinmind incident", "twinmind outage",
  "twinmind down", "twinmind emergency", "twinmind runbook".
allowed-tools: Read, Bash(curl:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Incident Runbook

## Overview
Procedures for diagnosing and resolving TwinMind integration incidents.

## Prerequisites
- Access to monitoring dashboards
- Production environment credentials
- On-call rotation contacts
- TwinMind support contact

## Incident Classification

| Severity | Description | Response Time | Examples |
|----------|-------------|---------------|----------|
| P1 - Critical | Complete service outage | 15 minutes | All transcriptions failing |
| P2 - High | Major feature degraded | 1 hour | Summaries not generating |
| P3 - Medium | Minor functionality impacted | 4 hours | Slow transcription |
| P4 - Low | Cosmetic or edge case | 24 hours | Occasional timeout |

## Initial Response Checklist

```markdown
## Incident Response - First 5 Minutes

### 1. Acknowledge Alert
- [ ] Acknowledge in PagerDuty/Opsgenie
- [ ] Join incident Slack channel
- [ ] Note start time: __________

### 2. Quick Health Check
- [ ] Check TwinMind status page: https://status.twinmind.com
- [ ] Check our service health endpoint
- [ ] Check recent deployments
- [ ] Check infrastructure status

### 3. Initial Assessment
- [ ] Identify affected functionality
- [ ] Estimate user impact
- [ ] Classify severity (P1/P2/P3/P4)

### 4. Communication
- [ ] Update status page (if P1/P2)
- [ ] Notify stakeholders
- [ ] Create incident ticket
```

## Diagnostic Commands

### Check TwinMind API Status

```bash
# Quick health check
curl -s -o /dev/null -w "%{http_code}" \
  -H "Authorization: Bearer $TWINMIND_API_KEY" \
  https://api.twinmind.com/v1/health

# Detailed health check
curl -s -H "Authorization: Bearer $TWINMIND_API_KEY" \
  https://api.twinmind.com/v1/health | jq

# Check status page API
curl -s https://status.twinmind.com/api/v2/status.json | jq '.status'
```

### Check Our Service Health

```bash
# Application health
curl -s http://localhost:8080/health | jq

# TwinMind-specific health
curl -s http://localhost:8080/health/twinmind | jq

# Check recent errors in logs
kubectl logs -l app=twinmind-service --tail=100 | grep -i error

# Check metrics endpoint
curl -s http://localhost:8080/metrics | grep twinmind_errors
```

### Check Rate Limits

```bash
# Get current rate limit status
curl -I -H "Authorization: Bearer $TWINMIND_API_KEY" \
  https://api.twinmind.com/v1/health 2>/dev/null | grep -i ratelimit

# Check our rate limit tracking
curl -s http://localhost:8080/metrics | grep rate_limit
```

## Common Incident Scenarios

### Scenario 1: All Transcriptions Failing

**Symptoms:**
- 100% error rate on transcriptions
- `twinmind_transcriptions_total{status="error"}` spiking

**Diagnosis Steps:**

```bash
# 1. Check TwinMind API availability
curl -v -H "Authorization: Bearer $TWINMIND_API_KEY" \
  https://api.twinmind.com/v1/health

# 2. Check error types in logs
kubectl logs -l app=twinmind-service --tail=200 | \
  grep -E "(error|Error|ERROR)" | \
  awk '{print $NF}' | sort | uniq -c | sort -rn

# 3. Test a simple transcription
curl -X POST \
  -H "Authorization: Bearer $TWINMIND_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{"audio_url":"https://example.com/test.mp3"}' \
  https://api.twinmind.com/v1/transcribe

# 4. Check if specific to our API key
# Try with a different key if available
```

**Resolution Steps:**

1. If TwinMind API is down:
   - Update status page
   - Notify users
   - Wait for TwinMind resolution
   - Monitor status.twinmind.com

2. If our API key is invalid:
   - Check key expiration
   - Regenerate key if needed
   - Update secrets manager
   - Restart services

3. If network issue:
   - Check DNS resolution
   - Check firewall rules
   - Verify egress connectivity

---

### Scenario 2: High Latency

**Symptoms:**
- P95 latency > 5 seconds
- Users reporting slow transcriptions
- `twinmind_api_latency_seconds` histogram showing high values

**Diagnosis Steps:**

```bash
# 1. Check latency metrics
curl -s http://localhost:8080/metrics | \
  grep twinmind_api_latency

# 2. Check for queuing
curl -s http://localhost:8080/metrics | \
  grep queue

# 3. Check TwinMind processing times
# Look at recent transcription durations
curl -s http://localhost:8080/metrics | \
  grep twinmind_transcription_duration

# 4. Check resource usage
kubectl top pods -l app=twinmind-service
```

**Resolution Steps:**

1. If TwinMind is slow:
   - Check status.twinmind.com for degradation notice
   - Consider switching to faster model (ear-2)
   - Queue non-urgent requests

2. If our service is overloaded:
   - Scale up replicas
   - Increase rate limiting
   - Queue requests

3. If network latency:
   - Check regional connectivity
   - Consider caching DNS
   - Verify no proxy issues

---

### Scenario 3: Rate Limiting

**Symptoms:**
- 429 errors increasing
- `twinmind_errors_total{error_type="RATE_LIMITED"}` > 0
- Rate limit remaining gauge near 0

**Diagnosis Steps:**

```bash
# 1. Check current rate limit status
curl -I -H "Authorization: Bearer $TWINMIND_API_KEY" \
  https://api.twinmind.com/v1/health 2>/dev/null | grep -i ratelimit

# 2. Check request volume
curl -s http://localhost:8080/metrics | \
  grep twinmind_api_requests_total

# 3. Identify heavy consumers
kubectl logs -l app=twinmind-service --tail=1000 | \
  grep -i "transcribe" | \
  awk '{print $4}' | sort | uniq -c | sort -rn | head -20
```

**Resolution Steps:**

1. Immediate:
   - Enable request queue if not already
   - Reject non-critical requests
   - Implement stricter client-side limiting

2. Short-term:
   - Request rate limit increase from TwinMind
   - Upgrade to higher tier
   - Distribute load across multiple API keys

3. Long-term:
   - Implement smarter caching
   - Batch requests where possible
   - Optimize client request patterns

---

### Scenario 4: Authentication Failures

**Symptoms:**
- 401 errors
- "Invalid API key" messages
- All authenticated requests failing

**Diagnosis Steps:**

```bash
# 1. Verify API key is set
echo $TWINMIND_API_KEY | head -c 20

# 2. Test API key directly
curl -H "Authorization: Bearer $TWINMIND_API_KEY" \
  https://api.twinmind.com/v1/me

# 3. Check if key was rotated
# Check TwinMind dashboard for key history

# 4. Verify secrets manager
aws secretsmanager get-secret-value --secret-id twinmind-api-key | jq
```

**Resolution Steps:**

1. If key expired/revoked:
   - Generate new API key in TwinMind dashboard
   - Update secrets manager
   - Restart affected services

2. If key leaked:
   - Immediately revoke key
   - Generate new key
   - Audit access logs
   - Update secrets
   - Review security practices

## Escalation Path

```
Level 1: On-Call Engineer (0-15 min)
   ↓
Level 2: Team Lead (15-30 min)
   ↓
Level 3: Engineering Manager (30-60 min)
   ↓
Level 4: VP Engineering (60+ min)
```

### TwinMind Support Contact

- **Email:** support@twinmind.com
- **Enterprise Support:** enterprise-support@twinmind.com
- **Status Page:** https://status.twinmind.com
- **Community:** https://community.twinmind.com

## Post-Incident Checklist

```markdown
## Post-Incident Review

### Immediate (Within 24 hours)
- [ ] Confirm issue fully resolved
- [ ] Update status page to resolved
- [ ] Notify stakeholders of resolution
- [ ] Document timeline in incident ticket

### Follow-up (Within 1 week)
- [ ] Schedule post-mortem meeting
- [ ] Write incident report
- [ ] Identify root cause
- [ ] Create action items
- [ ] Update runbook if needed

### Long-term
- [ ] Implement preventive measures
- [ ] Update monitoring/alerting
- [ ] Share learnings with team
- [ ] Update documentation
```

## Incident Report Template

```markdown
# Incident Report: [TITLE]

**Date:** YYYY-MM-DD
**Duration:** HH:MM - HH:MM (X hours Y minutes)
**Severity:** P1/P2/P3/P4
**Impact:** [Number of affected users/requests]

## Summary
[Brief description of what happened]

## Timeline
- HH:MM - Alert triggered
- HH:MM - On-call acknowledged
- HH:MM - [Key actions taken]
- HH:MM - Issue resolved

## Root Cause
[Detailed explanation of why the incident occurred]

## Resolution
[How the issue was fixed]

## Action Items
- [ ] [Preventive measure 1]
- [ ] [Monitoring improvement]
- [ ] [Documentation update]

## Lessons Learned
[Key takeaways for the team]
```

## Output
- Incident classification guide
- Diagnostic commands
- Common scenario playbooks
- Escalation procedures
- Post-incident checklist
- Report template

## Resources
- [TwinMind Status Page](https://status.twinmind.com)
- [TwinMind Support](https://twinmind.com/support)
- [PagerDuty Best Practices](https://response.pagerduty.com/)

## Next Steps
For data handling procedures, see `twinmind-data-handling`.
