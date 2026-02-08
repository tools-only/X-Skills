# Audit Logging

MCP Gateway Registry provides comprehensive audit logging for compliance, security monitoring, and operational visibility. All API requests and MCP server access events are logged to MongoDB/DocumentDB with automatic retention management.

![Audit Log Viewer](img/audit-log.png)

## Overview

Audit logging captures two types of events:

1. **Registry API Access** - All REST API requests to the Registry (`/api/*`, `/v0.1/*`)
2. **MCP Server Access** - All MCP protocol requests proxied through the Gateway

Sensitive data such as authentication tokens, session cookies, and passwords are never logged. Credentials are masked to show only the last 6 characters as a hint for debugging.

## Security and Privacy

### Data That Is NOT Logged

The following sensitive data is explicitly excluded from audit logs:

- **Authentication tokens** (Bearer tokens, JWT tokens)
- **Session cookies** (Cookie header values)
- **Passwords** (form fields, query parameters)
- **API keys** (full values)
- **Refresh tokens**
- **Authorization header values**

### Data Masking

When credential hints are logged for debugging purposes, they are automatically masked:

- Full token: `eyJhbGciOiJSUzI1NiIsInR5...` becomes `***zI1Ni`
- Tokens shorter than 6 characters become `***`

Query parameters with sensitive names (token, password, key, secret, api_key, etc.) are automatically masked.

## Event Schemas

### Registry API Access Event

Logged for every REST API request to the Registry.

```json
{
  "timestamp": "2026-02-06T10:30:00.000Z",
  "log_type": "registry_api_access",
  "version": "1.0",
  "request_id": "abc123-def456-...",
  "correlation_id": null,
  "identity": {
    "username": "john.doe@example.com",
    "auth_method": "oauth2",
    "provider": "keycloak",
    "groups": ["mcp-registry-admin", "developers"],
    "scopes": ["registry-admins"],
    "is_admin": true,
    "credential_type": "session_cookie",
    "credential_hint": "***abc123"
  },
  "request": {
    "method": "POST",
    "path": "/api/servers",
    "query_params": {},
    "client_ip": "192.168.1.100",
    "forwarded_for": "10.0.0.1",
    "user_agent": "Mozilla/5.0...",
    "content_length": 1024
  },
  "response": {
    "status_code": 201,
    "duration_ms": 45.32,
    "content_length": 512
  },
  "action": {
    "operation": "create",
    "resource_type": "server",
    "resource_id": "my-mcp-server",
    "description": "Create new MCP server"
  },
  "authorization": {
    "decision": "ALLOW",
    "required_permission": "servers:write",
    "evaluated_scopes": ["registry-admins"]
  }
}
```

### MCP Server Access Event

Logged for every MCP protocol request proxied through the Gateway.

```json
{
  "timestamp": "2026-02-06T10:30:00.000Z",
  "log_type": "mcp_server_access",
  "version": "1.0",
  "request_id": "xyz789-...",
  "correlation_id": null,
  "identity": {
    "username": "ai-agent@example.com",
    "auth_method": "jwt_bearer",
    "provider": "keycloak",
    "groups": [],
    "scopes": ["mcp-server-cloudflare-docs"],
    "is_admin": false,
    "credential_type": "bearer_token",
    "credential_hint": "***def456"
  },
  "mcp_server": {
    "name": "cloudflare-docs",
    "path": "/cloudflare-docs",
    "version": "1.0.0",
    "proxy_target": "http://internal-mcp-server:8080/mcp"
  },
  "mcp_request": {
    "method": "tools/call",
    "tool_name": "search_docs",
    "resource_uri": null,
    "mcp_session_id": "session-123",
    "transport": "streamable-http",
    "jsonrpc_id": "1"
  },
  "mcp_response": {
    "status": "success",
    "duration_ms": 123.45,
    "error_code": null,
    "error_message": null
  }
}
```

## Data Fields Reference

### Identity Fields

| Field | Description |
|-------|-------------|
| `username` | Username or identifier of the requester |
| `auth_method` | Authentication method: `oauth2`, `traditional`, `jwt_bearer`, `anonymous` |
| `provider` | Identity provider: `cognito`, `entra_id`, `keycloak` |
| `groups` | Groups the user belongs to |
| `scopes` | OAuth scopes granted to the user |
| `is_admin` | Whether the user has admin privileges |
| `credential_type` | Type of credential: `session_cookie`, `bearer_token`, `none` |
| `credential_hint` | Masked hint of the credential (last 6 chars only) |

### Action Fields (Registry API only)

| Field | Description |
|-------|-------------|
| `operation` | Operation type: `create`, `read`, `update`, `delete`, `list`, `toggle`, `rate`, `login`, `logout`, `search` |
| `resource_type` | Resource type: `server`, `agent`, `auth`, `federation`, `health`, `search` |
| `resource_id` | Identifier of the resource being acted upon |
| `description` | Human-readable description of the action |

### MCP Request Fields (MCP Access only)

| Field | Description |
|-------|-------------|
| `method` | JSON-RPC method name: `tools/call`, `tools/list`, `resources/read`, `resources/list`, etc. |
| `tool_name` | Name of the tool being called (for `tools/call` method) |
| `resource_uri` | URI of the resource being accessed (for `resources/read` method) |
| `mcp_session_id` | MCP session identifier |
| `transport` | Transport protocol: `streamable-http`, `sse`, `stdio` |
| `jsonrpc_id` | JSON-RPC request ID |

## Data Retention

Audit logs are automatically expired using MongoDB/DocumentDB TTL (Time-To-Live) indexes.

### Default Retention

- **Default retention period**: 7 days
- **TTL index field**: `timestamp`

### Configuring Retention

Set the `AUDIT_LOG_MONGODB_TTL_DAYS` environment variable to customize retention:

```bash
# Keep logs for 30 days
export AUDIT_LOG_MONGODB_TTL_DAYS=30

# Keep logs for 90 days (compliance requirement)
export AUDIT_LOG_MONGODB_TTL_DAYS=90
```

The TTL index is created when running the DocumentDB initialization script:

```bash
./scripts/init-documentdb.sh
```

### Important Notes

- TTL indexes run approximately once per minute in MongoDB/DocumentDB
- Documents may persist slightly longer than the TTL value
- Changing the TTL requires dropping and recreating the index with `--recreate` flag
- For compliance requirements, consider also streaming logs to a long-term archive

## Storage

### MongoDB Collection

Audit events are stored in the `audit_events_{namespace}` collection with the following indexes:

| Index | Purpose |
|-------|---------|
| `request_id` (unique) | Fast lookup by request ID |
| `identity.username` + `timestamp` | Query by user over time range |
| `action.operation` + `timestamp` | Query by operation type over time range |
| `action.resource_type` + `timestamp` | Query by resource type over time range |
| `timestamp` (TTL) | Automatic expiration after configured days |

### Storage Sizing

Typical event sizes:
- Registry API event: ~1-2 KB
- MCP Server Access event: ~1-2 KB

Estimated storage (without compression):
- 1,000 requests/day for 7 days: ~14 MB
- 10,000 requests/day for 30 days: ~600 MB
- 100,000 requests/day for 90 days: ~18 GB

## Viewing Audit Logs

### Admin UI

Administrators can view audit logs in the Registry UI:

1. Navigate to **Settings** > **Audit** > **Audit Logs**
2. Select log stream: **Registry API** or **MCP Access**
3. Apply filters (time range, username, operation, status)
4. Click any row to view full event details
5. Export filtered results as JSONL or CSV

### API Access

Query audit events programmatically:

```bash
# Get recent Registry API events
curl -H "Authorization: Bearer $TOKEN" \
  "https://registry.example.com/api/audit/events?stream=registry_api&limit=50"

# Get MCP access events for a specific user
curl -H "Authorization: Bearer $TOKEN" \
  "https://registry.example.com/api/audit/events?stream=mcp_access&username=john.doe"

# Export events as JSONL
curl -H "Authorization: Bearer $TOKEN" \
  "https://registry.example.com/api/audit/export?stream=registry_api&format=jsonl"
```

### MongoDB/DocumentDB Direct Query

```javascript
// Find all events for a user in the last 24 hours
db.audit_events_default.find({
  "identity.username": "john.doe@example.com",
  "timestamp": { $gte: new Date(Date.now() - 24*60*60*1000) }
}).sort({ timestamp: -1 })

// Count events by operation type
db.audit_events_default.aggregate([
  { $match: { log_type: "registry_api_access" } },
  { $group: { _id: "$action.operation", count: { $sum: 1 } } },
  { $sort: { count: -1 } }
])

// Find failed MCP requests
db.audit_events_default.find({
  "log_type": "mcp_server_access",
  "mcp_response.status": "error"
})
```

## Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `AUDIT_LOG_ENABLED` | `true` | Enable/disable audit logging |
| `AUDIT_LOG_MONGODB_TTL_DAYS` | `7` | Log retention period in days |

### Non-Blocking Design

Audit logging is designed to never impact request processing:

- Logging happens asynchronously after the response is sent
- Failures in audit logging are logged as warnings but don't fail requests
- High-volume scenarios use batched writes (if enabled)

## Compliance Considerations

### SOC 2 / ISO 27001

Audit logs support compliance requirements by capturing:

- **Who**: User identity with auth method and provider
- **What**: Operation performed with resource details
- **When**: Precise UTC timestamp
- **Where**: Client IP and forwarded-for headers
- **Outcome**: Success/failure status with error details

### GDPR

- User identifiers (usernames) are logged for accountability
- No PII beyond usernames is captured
- Logs can be exported and deleted per data subject requests
- TTL-based retention supports data minimization

### Additional Recommendations

For production compliance deployments:

1. Stream audit logs to a SIEM (Splunk, Datadog, etc.) for long-term retention
2. Set up alerts for suspicious patterns (failed auths, privilege escalation)
3. Regularly review admin actions in the audit log
4. Document your retention policy and ensure TTL matches it

## Troubleshooting

### Logs Not Appearing

1. Verify audit logging is enabled: `AUDIT_LOG_ENABLED=true`
2. Check MongoDB connection: ensure the Registry can write to the database
3. Look for warnings in Registry logs: `grep "audit" registry.log`

### TTL Not Working

1. Verify the TTL index exists: `db.audit_events_default.getIndexes()`
2. Note that MongoDB TTL runs approximately every 60 seconds
3. Documents may persist up to 60 seconds beyond their expiration time

### Missing Events

1. Check if the request completed (cancelled requests may not be logged)
2. Verify the log stream filter matches the event type
3. For MCP access, ensure the path is not an API path (starts with `/api/`)
