# Scopes Management

This document describes the scope configuration file format used by the MCP Gateway Registry for fine-grained access control.

## Overview

Scopes define what resources (MCP servers, agents) users can access and what actions they can perform. The registry uses JSON-based scope configuration files that can be loaded during initialization or managed via the CLI.

## Scope Configuration File Format

### Example Files

- `scripts/registry-admins.json` - Bootstrap admin scope loaded during database initialization
- `cli/examples/public-mcp-users.json` - Example scope for users with limited access

### Complete Field Reference

```json
{
  "_id": "scope-name",
  "scope_name": "scope-name",
  "description": "Human-readable description of this scope",
  "group_mappings": ["group-name-1", "group-uuid-2"],
  "server_access": [
    {
      "server": "server-name",
      "methods": ["initialize", "tools/list", "tools/call"],
      "tools": ["tool-name-1", "tool-name-2"]
    },
    {
      "agents": {
        "actions": [
          {"action": "list_agents", "resources": ["/agent-path"]},
          {"action": "get_agent", "resources": ["/agent-path"]}
        ]
      }
    }
  ],
  "ui_permissions": {
    "list_agents": ["all"],
    "get_agent": ["/specific-agent"],
    "publish_agent": [],
    "list_service": ["all"],
    "toggle_service": ["service-name"]
  },
  "create_in_idp": true
}
```

## Field Descriptions

### Top-Level Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `_id` | string | Yes | Unique identifier for the scope document in MongoDB. Should match `scope_name`. |
| `scope_name` | string | No | Human-readable scope name. If omitted, `_id` is used. |
| `description` | string | No | Description explaining the purpose of this scope. |
| `group_mappings` | array | Yes | List of IdP group names or IDs that map to this scope. |
| `server_access` | array | Yes | List of MCP server access rules and agent action permissions. |
| `ui_permissions` | object | No | UI-level permissions for the registry web interface. |
| `create_in_idp` | boolean | No | When true, the CLI will create the group in the IdP (Keycloak/Entra). |

### group_mappings Field

The `group_mappings` array contains IdP group identifiers that should be mapped to this scope. When a user authenticates, their IdP groups are matched against these mappings to determine their effective scopes.

**Important for Entra ID:**
- Entra ID uses Group Object IDs (GUIDs), not group names
- You must include the Group Object ID from Azure Portal > Groups > Overview
- Example: `"5f605d68-06bc-4208-b992-bb378eee12c5"`

**For Keycloak:**
- Use the group name as defined in Keycloak
- Example: `"public-mcp-users"`

**Example with both:**
```json
{
  "group_mappings": [
    "public-mcp-users",
    "5f605d68-06bc-4208-b992-bb378eee12c5"
  ]
}
```

This means users in either the Keycloak group `public-mcp-users` OR the Entra ID group with Object ID `5f605d68-06bc-4208-b992-bb378eee12c5` will receive this scope.

### server_access Field

The `server_access` array defines what MCP servers and methods users can access. Each entry can be either a server access rule or an agent actions block.

#### Server Access Rule

```json
{
  "server": "server-name-or-wildcard",
  "methods": ["method-1", "method-2"],
  "tools": ["tool-name-or-wildcard"]
}
```

| Field | Description |
|-------|-------------|
| `server` | Server name or `"*"` for all servers |
| `methods` | List of allowed MCP methods (see below) |
| `tools` | List of allowed tool names or `["*"]` for all tools |

**Standard MCP Methods:**
- `initialize` - Initialize MCP session
- `notifications/initialized` - Session initialized notification
- `ping` - Health check
- `tools/list` - List available tools
- `tools/call` - Execute a tool
- `resources/list` - List available resources
- `resources/templates/list` - List resource templates
- `GET`, `POST`, `PUT`, `DELETE` - HTTP methods for REST API access

**Example - Full MCP access to specific servers:**
```json
{
  "server": "context7",
  "methods": [
    "initialize",
    "notifications/initialized",
    "ping",
    "tools/list",
    "tools/call",
    "resources/list",
    "resources/templates/list"
  ],
  "tools": ["*"]
}
```

**Example - Wildcard access (admin):**
```json
{
  "server": "*",
  "methods": ["all"],
  "tools": ["all"]
}
```

#### Agent Actions Block

Agent actions define what operations users can perform on A2A agents.

```json
{
  "agents": {
    "actions": [
      {"action": "action-name", "resources": ["/agent-path-1", "/agent-path-2"]}
    ]
  }
}
```

**Available Agent Actions:**

| Action | Description | API Endpoint |
|--------|-------------|--------------|
| `list_agents` | View agents in listings | `GET /api/agents` |
| `get_agent` | View agent details | `GET /api/agents/{path}` |
| `publish_agent` | Register new agents | `POST /api/agents/register` |
| `modify_agent` | Update existing agents | `PUT /api/agents/{path}` |
| `delete_agent` | Remove agents | `DELETE /api/agents/{path}` |

**Resource Patterns:**
- `/agent-name` - Specific agent path (e.g., `/flight-booking`)
- `all` - All agents (wildcard access)

**Example - Limited agent access:**
```json
{
  "agents": {
    "actions": [
      {"action": "list_agents", "resources": ["/flight-booking", "/code-reviewer"]},
      {"action": "get_agent", "resources": ["/flight-booking", "/code-reviewer"]}
    ]
  }
}
```

**Example - Full agent admin access:**
```json
{
  "agents": {
    "actions": [
      {"action": "list_agents", "resources": ["all"]},
      {"action": "get_agent", "resources": ["all"]},
      {"action": "publish_agent", "resources": ["all"]},
      {"action": "modify_agent", "resources": ["all"]},
      {"action": "delete_agent", "resources": ["all"]}
    ]
  }
}
```

### ui_permissions Field

UI permissions control what actions users can perform in the web interface and REST API for service/agent management.

```json
{
  "ui_permissions": {
    "permission_name": ["resource-1", "resource-2"]
  }
}
```

**Available UI Permissions:**

| Permission | Description | Applies To |
|------------|-------------|------------|
| `list_agents` | View agents in UI | Agent paths or `"all"` |
| `get_agent` | View agent details | Agent paths or `"all"` |
| `publish_agent` | Register new agents via UI | Agent paths or `"all"` |
| `modify_agent` | Edit agents via UI | Agent paths or `"all"` |
| `delete_agent` | Delete agents via UI | Agent paths or `"all"` |
| `list_service` | View MCP servers in UI | Server names or `"all"` |
| `register_service` | Register new MCP servers | Server names or `"all"` |
| `health_check_service` | Run health checks | Server names or `"all"` |
| `toggle_service` | Enable/disable servers | Server names or `"all"` |
| `modify_service` | Edit server configurations | Server names or `"all"` |

**Example - Read-only access:**
```json
{
  "ui_permissions": {
    "list_service": ["all"],
    "list_agents": ["/flight-booking"],
    "get_agent": ["/flight-booking"]
  }
}
```

**Example - Full admin access:**
```json
{
  "ui_permissions": {
    "list_agents": ["all"],
    "get_agent": ["all"],
    "publish_agent": ["all"],
    "modify_agent": ["all"],
    "delete_agent": ["all"],
    "list_service": ["all"],
    "register_service": ["all"],
    "health_check_service": ["all"],
    "toggle_service": ["all"],
    "modify_service": ["all"]
  }
}
```

## Complete Examples

### Admin Scope (registry-admins.json)

Full access to all servers, agents, and UI functions:

```json
{
  "_id": "registry-admins",
  "group_mappings": ["registry-admins"],
  "server_access": [
    {
      "server": "*",
      "methods": ["all"],
      "tools": ["all"]
    },
    {
      "agents": {
        "actions": [
          {"action": "list_agents", "resources": ["all"]},
          {"action": "get_agent", "resources": ["all"]},
          {"action": "publish_agent", "resources": ["all"]},
          {"action": "modify_agent", "resources": ["all"]},
          {"action": "delete_agent", "resources": ["all"]}
        ]
      }
    }
  ],
  "ui_permissions": {
    "list_agents": ["all"],
    "get_agent": ["all"],
    "publish_agent": ["all"],
    "modify_agent": ["all"],
    "delete_agent": ["all"],
    "list_service": ["all"],
    "register_service": ["all"],
    "health_check_service": ["all"],
    "toggle_service": ["all"],
    "modify_service": ["all"]
  }
}
```

### Limited User Scope (public-mcp-users.json)

Access to specific MCP servers and one agent:

```json
{
  "scope_name": "public-mcp-users",
  "description": "Users with access to public MCP servers and flight-booking agent",
  "server_access": [
    {
      "server": "context7",
      "methods": [
        "initialize",
        "notifications/initialized",
        "ping",
        "tools/list",
        "tools/call",
        "resources/list",
        "resources/templates/list"
      ],
      "tools": ["*"]
    },
    {
      "server": "api",
      "methods": ["initialize", "GET", "POST", "servers", "agents", "search"],
      "tools": []
    },
    {
      "agents": {
        "actions": [
          {"action": "list_agents", "resources": ["/flight-booking"]},
          {"action": "get_agent", "resources": ["/flight-booking"]}
        ]
      }
    }
  ],
  "group_mappings": [
    "public-mcp-users",
    "5f605d68-06bc-4208-b992-bb378eee12c5"
  ],
  "ui_permissions": {
    "list_service": ["all"],
    "list_agents": ["/flight-booking"],
    "get_agent": ["/flight-booking"]
  },
  "create_in_idp": true
}
```

## Managing Scopes

### Using the CLI

Import a scope from JSON file:
```bash
uv run python api/registry_management.py \
  --token-file .token \
  --registry-url https://registry.example.com \
  import-group cli/examples/public-mcp-users.json
```

List all scopes:
```bash
uv run python api/registry_management.py \
  --token-file .token \
  --registry-url https://registry.example.com \
  list-groups
```

### Bootstrap Admin Scope

The `registry-admins` scope is automatically loaded during database initialization:
- **Local (MongoDB CE)**: `docker compose up mongodb-init`
- **Production (DocumentDB)**: `./terraform/aws-ecs/scripts/run-documentdb-init.sh`

### Server Path Variations

When defining server access, you may need to include path variations to handle different URL patterns:

```json
{
  "server_access": [
    {"server": "context7", "methods": [...], "tools": ["*"]},
    {"server": "/context7", "methods": [...], "tools": ["*"]},
    {"server": "/context7/", "methods": [...], "tools": ["*"]}
  ]
}
```

This ensures access works regardless of whether the server is accessed as:
- `context7`
- `/context7`
- `/context7/`

## Entra ID Integration

When using Microsoft Entra ID (Azure AD) as the identity provider:

1. **Create a group in Azure Portal:**
   - Navigate to Azure Portal > Azure Active Directory > Groups
   - Create a new Security group
   - Note the Group Object ID (GUID)

2. **Add the Object ID to group_mappings:**
   ```json
   {
     "group_mappings": [
       "my-keycloak-group",
       "12345678-1234-1234-1234-123456789012"
     ]
   }
   ```

3. **Assign users to the Azure AD group:**
   - Users in this group will receive the scope permissions when they authenticate

4. **Configure Entra ID app to include groups in tokens:**
   - In the App Registration, configure the `groups` claim
   - Set `groupMembershipClaims` to `"SecurityGroup"` in the manifest

## Troubleshooting

### User Not Getting Expected Permissions

1. Check group membership in IdP (Keycloak/Entra)
2. Verify `group_mappings` includes the correct group name/ID
3. Check registry logs for scope mapping messages
4. Use the debug endpoint: `GET /api/debug/user-context`

### Scope Not Found

1. Ensure the scope was imported: `list-groups` command
2. Check MongoDB collection: `mcp_scopes_default`
3. Re-run database initialization if bootstrap scope missing

### Entra ID Groups Not Working

1. Verify Group Object ID (not display name) is in `group_mappings`
2. Check that `groupMembershipClaims` is configured in app manifest
3. Verify user is assigned to the group in Azure Portal
4. Check that optional claims include `groups` in ID token
