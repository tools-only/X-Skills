---
skill_name: api-docs
skill_category: documentation
description: API documentation patterns with endpoint specs, auth flows, and error handling
allowed_tools: [Read, Edit, Glob, Grep, Write]
token_estimate: 600
version: 1.0
last_updated: 2026-01-13
owner: Claude Copilot
status: active

tags: [api, documentation, openapi, rest, endpoints, anti-pattern, best-practice, validation]
related_skills: [tutorial-patterns]

trigger_files: ["**/api/**", "**/routes/**", "**/controllers/**", "**/*.openapi.*", "**/swagger.*"]
trigger_keywords: [api docs, endpoint documentation, api reference, swagger, openapi, rest api]

quality_keywords: [anti-pattern, pattern, validation, best-practice, api-first, versioning]
---

# API Documentation

Patterns for creating clear, complete API documentation that developers can use without guessing.

## Purpose

- Provide consistent structure for API endpoint documentation
- Ensure all critical information (auth, errors, examples) is included
- Enable developers to make successful API calls on first try

---

## Core Patterns

### Pattern 1: Endpoint Specification

**When to use:** Documenting any REST API endpoint.

**Implementation:**
```markdown
## `METHOD /path/:param`

Brief description of what this endpoint does.

### Authentication
Bearer token required in Authorization header.

### Parameters
| Name | Type | In | Required | Description |
|------|------|-----|----------|-------------|
| id | string | path | Yes | Resource identifier |
| limit | integer | query | No | Max results (default: 20) |

### Request Body
```json
{
  "field": "value"
}
```

### Response
```json
{
  "id": "abc123",
  "created_at": "2024-01-01T00:00:00Z"
}
```

### Errors
| Code | Meaning | Resolution |
|------|---------|------------|
| 400 | Invalid request | Check request body format |
| 401 | Unauthorized | Include valid Bearer token |
| 404 | Not found | Verify resource ID exists |

### Example
```bash
curl -X GET "https://api.example.com/resource/abc123" \
  -H "Authorization: Bearer TOKEN"
```
```

**Benefits:**
- Developers find information predictably
- All edge cases documented upfront
- Copy-paste examples reduce errors

### Pattern 2: Error Response Documentation

**When to use:** Documenting error responses for any endpoint.

**Implementation:**
```markdown
### Error Response Format
All errors return this structure:
```json
{
  "error": {
    "code": "RESOURCE_NOT_FOUND",
    "message": "User with ID 'xyz' not found",
    "details": {
      "resource": "user",
      "id": "xyz"
    }
  }
}
```

### Common Error Codes
| Code | HTTP Status | Description | Action |
|------|-------------|-------------|--------|
| `INVALID_REQUEST` | 400 | Malformed request body | Validate JSON structure |
| `UNAUTHORIZED` | 401 | Missing or invalid token | Refresh authentication |
| `FORBIDDEN` | 403 | Insufficient permissions | Check user roles |
| `NOT_FOUND` | 404 | Resource doesn't exist | Verify ID |
| `RATE_LIMITED` | 429 | Too many requests | Wait and retry |
```

**Benefits:**
- Consistent error handling in client code
- Actionable resolution guidance
- Reduces support burden

---

## Anti-Patterns

### Anti-Pattern 1: Missing Authentication Details

| Aspect | Description |
|--------|-------------|
| **WHY** | Developers waste time debugging 401 errors; increases support tickets |
| **DETECTION** | Endpoint docs without "Authentication" section; no example with auth header |
| **FIX** | Always include auth method, header format, and token example |

**Bad Example:**
```markdown
## GET /api/users

Returns list of users.

### Response
```json
[{"id": "1", "name": "John"}]
```
```

**Good Example:**
```markdown
## GET /api/users

Returns list of users.

### Authentication
Requires Bearer token in Authorization header.
```
Authorization: Bearer <your-api-key>
```

### Response
```json
[{"id": "1", "name": "John"}]
```

### Example
```bash
curl -X GET "https://api.example.com/api/users" \
  -H "Authorization: Bearer sk_live_abc123"
```
```

### Anti-Pattern 2: Undocumented Error States

| Aspect | Description |
|--------|-------------|
| **WHY** | Clients can't handle errors gracefully; leads to poor UX and debugging chaos |
| **DETECTION** | Docs only show success responses; no error codes or troubleshooting |
| **FIX** | Document every possible error with code, cause, and resolution |

**Bad Example:**
```markdown
### Response
```json
{"status": "success", "data": {...}}
```
```

**Good Example:**
```markdown
### Success Response (200)
```json
{"status": "success", "data": {...}}
```

### Error Responses
| Status | Code | When | Resolution |
|--------|------|------|------------|
| 400 | `INVALID_EMAIL` | Email format wrong | Use valid email |
| 409 | `USER_EXISTS` | Email already registered | Use login endpoint |
| 422 | `WEAK_PASSWORD` | Password too simple | 8+ chars, mixed case |
```

### Anti-Pattern 3: Example-Free Documentation

| Aspect | Description |
|--------|-------------|
| **WHY** | Forces developers to guess; increases time to first successful call |
| **DETECTION** | No curl/code examples; no sample request/response pairs |
| **FIX** | Include complete, copy-paste-ready examples for every endpoint |

**Bad Example:**
```markdown
## POST /api/users
Creates a new user with the provided data.
```

**Good Example:**
```markdown
## POST /api/users
Creates a new user.

### Request
```bash
curl -X POST "https://api.example.com/api/users" \
  -H "Authorization: Bearer TOKEN" \
  -H "Content-Type: application/json" \
  -d '{"email": "user@example.com", "name": "John Doe"}'
```

### Response
```json
{
  "id": "usr_abc123",
  "email": "user@example.com",
  "name": "John Doe",
  "created_at": "2024-01-15T10:30:00Z"
}
```
```

---

## Validation Checklist

### Pre-Documentation
- [ ] Endpoint actually exists and works
- [ ] Understand all parameters and their validation
- [ ] Tested all error cases

### Documentation
- [ ] Authentication section included
- [ ] All parameters documented with types
- [ ] Request body example provided (if applicable)
- [ ] Success response example included
- [ ] Error responses documented with resolution
- [ ] Complete curl example provided

### Post-Documentation
- [ ] Example actually works when copied
- [ ] All status codes verified against implementation
- [ ] Cross-linked to related endpoints

---

## Related Resources

- Related skills: `skill_get("tutorial-patterns")`
- OpenAPI Specification: https://swagger.io/specification/

---

## Changelog

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-01-13 | Initial version with anti-patterns |
