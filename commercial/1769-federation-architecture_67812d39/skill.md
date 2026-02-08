# Federation Architecture Design

This document describes the peer-to-peer federation architecture for MCP Gateway Registry instances.

## Overview

Federation enables multiple MCP Gateway Registry instances to share servers and agents across organizational boundaries. Unlike external registry integration (Anthropic, ASOR), peer-to-peer federation connects registry instances that run the same codebase, enabling bidirectional synchronization with fine-grained control.

## Architecture Principles

### Peer-to-Peer Symmetric Design

The codebase has **no concept of "hub" or "spoke"**. Every registry instance runs identical code with identical capabilities:

- Any registry can be both an exporter and an importer simultaneously
- There is no role flag, no hierarchy, no hardcoded topology
- Terms like "Hub" and "LOB" are purely organizational labels describing deployment choices

```
Symmetric Architecture:

  Registry A  <------>  Registry B
       |                    |
       v                    v
  Registry C  <------>  Registry D

Any registry can pull from any other registry.
Any registry can export to any other registry.
```

### Common Deployment Patterns

**Hub Pulls from LOBs (Centralized Visibility)**

A central IT team maintains a Hub Registry. Each Line of Business (LOB) maintains their own registry. The Hub pulls from LOBs to provide centralized visibility.

```
         LOB-A Registry   LOB-B Registry   LOB-C Registry
                  \           |            /
                   \          |           /
                    Hub Registry (Central IT)
                    "Show me everything across all LOBs"
```

**LOBs Pull from Hub (Inheritance)**

LOBs inherit shared tools and agents from a central Hub.

```
                    Hub Registry (Central IT)
                   /          |           \
                  /           |            \
         LOB-A Registry   LOB-B Registry   LOB-C Registry
         "Inherit shared tools from central"
```

**Mesh Topology**

For organizations with peer relationships between registries.

```
         Registry A  <------>  Registry B
              ^                     ^
              |                     |
              v                     v
         Registry C  <------>  Registry D
```

## Data Model

### Peer Registry Configuration

Each peer is configured with the following attributes:

| Field | Type | Description |
|-------|------|-------------|
| `peer_id` | string | Unique identifier for the peer (e.g., "lob-a", "hub") |
| `name` | string | Human-readable name (e.g., "LOB-A Registry") |
| `endpoint` | string | Base URL of the peer registry |
| `enabled` | boolean | Whether sync is enabled for this peer |
| `sync_mode` | enum | `all`, `whitelist`, or `tag_filter` |
| `whitelist_servers` | string[] | Server paths to include (if sync_mode=whitelist) |
| `whitelist_agents` | string[] | Agent paths to include (if sync_mode=whitelist) |
| `tag_filters` | string[] | Tags to filter by (if sync_mode=tag_filter) |
| `sync_interval_minutes` | int | Interval for scheduled sync (0 = manual only) |
| `federation_token` | string | Static token for authenticating to this peer (encrypted) |
| `expected_client_id` | string | OAuth2 client ID expected from this peer (for peer identification) |

### Sync Modes

**All Mode**
```json
{
  "peer_id": "lob-a",
  "sync_mode": "all"
}
```
Imports all public servers and agents from the peer.

**Whitelist Mode**
```json
{
  "peer_id": "lob-a",
  "sync_mode": "whitelist",
  "whitelist_servers": ["/critical-tool", "/shared-service"],
  "whitelist_agents": ["/data-analyst-agent"]
}
```
Imports only the specified servers and agents.

**Tag Filter Mode**
```json
{
  "peer_id": "lob-a",
  "sync_mode": "tag_filter",
  "tag_filters": ["production", "shared"]
}
```
Imports servers and agents that have any of the specified tags.

### Visibility Control

Servers and agents have a `visibility` field that controls federation export:

| Visibility | Behavior |
|------------|----------|
| `internal` | Not exported via federation (default) |
| `public` | Exported to all authenticated peers |
| `group-restricted` | Exported only to peers in specified groups |

## Authentication

### Static Token Authentication (Recommended)

The primary authentication method for federation uses static tokens. This is IdP-agnostic and works regardless of whether registries use Keycloak, Entra ID, Cognito, or no identity provider at all.

**On the Exporting Registry:**

```bash
# .env on the exporting registry
FEDERATION_STATIC_TOKEN_AUTH_ENABLED=true
FEDERATION_STATIC_TOKEN=<generated-secret-key>
```

**On the Importing Registry:**

```json
{
  "peer_id": "lob-a",
  "endpoint": "https://lob-a-registry.corp.com",
  "federation_token": "<token-from-lob-a>"
}
```

The token is encrypted using Fernet symmetric encryption before storage in MongoDB/DocumentDB.

### Encryption at Rest

Secrets stored in peer configurations (federation tokens, OAuth client secrets) are encrypted using Fernet (AES-128-CBC):

```bash
# Generate encryption key (one-time)
python -c "from cryptography.fernet import Fernet; print(Fernet.generate_key().decode())"

# Set in environment
FEDERATION_ENCRYPTION_KEY=<generated-fernet-key>
```

This approach is:
- Database-agnostic (same code for MongoDB CE and DocumentDB)
- Simple (one env var holds the symmetric key)
- No extra dependencies on database side

### OAuth2 Client Credentials (Alternative)

For organizations requiring short-lived tokens and JWT audit trails, OAuth2 client credentials flow is supported as an alternative:

```json
{
  "peer_id": "lob-external",
  "endpoint": "https://external-registry.example.com",
  "auth_config": {
    "token_endpoint": "https://keycloak.example.com/realms/mcp-gateway/protocol/openid-connect/token",
    "client_id": "federation-hub-m2m",
    "client_secret": "<encrypted>",
    "scope": null
  }
}
```

This works with any OAuth2-compliant provider (Keycloak, Entra ID, Cognito).

## Sync Process

### Sync Triggers

Sync can be triggered in three ways:

1. **Manual via API**: `POST /api/peers/{peer_id}/sync`
2. **Manual via UI**: Click "Sync Now" in the Settings UI
3. **Scheduled**: Background scheduler checks every 60 seconds and triggers sync for peers where `sync_interval_minutes > 0` and the interval has elapsed

### Sync Flow

```
1. Load peer configuration from MongoDB
2. Decrypt authentication credentials
3. Authenticate to peer's federation export API
4. Fetch servers and agents based on sync_mode filters
5. Apply visibility filtering on the exporting side
6. Store synced items with federation metadata:
   - is_federated: true
   - source_peer_id: "lob-a"
   - upstream_path: "/original/path"
   - is_read_only: true
7. Update sync status (last_sync, generation number)
8. Mark items not present in sync as orphaned
```

### Federation Metadata

Synced servers and agents carry metadata indicating their federated origin:

```json
{
  "name": "External Tool",
  "path": "/lob-a/external-tool",
  "sync_metadata": {
    "is_federated": true,
    "source_peer_id": "lob-a",
    "upstream_path": "/external-tool",
    "last_synced_at": "2026-02-05T10:30:00Z",
    "is_read_only": true
  }
}
```

### Path Namespacing

Synced items are namespaced under their peer ID to prevent collisions:

| Original Path (on peer) | Synced Path (on importer) |
|------------------------|---------------------------|
| `/my-tool` | `/lob-a/my-tool` |
| `/data-agent` | `/lob-a/data-agent` |

### Orphan Detection

When a server or agent is removed from the upstream peer, the sync process:

1. Detects the item is no longer present in the export
2. Increments the generation number
3. Items from previous generations are marked as orphaned
4. Admins can review orphaned items before removal

## API Endpoints

### Federation Export API

Endpoints exposed by registries for peers to pull data:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/federation/health` | GET | Health check (unauthenticated) |
| `/api/v1/federation/servers` | GET | Export servers for sync |
| `/api/v1/federation/agents` | GET | Export agents for sync |

### Peer Management API

Endpoints for managing peer configurations:

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/peers` | GET | List all configured peers |
| `/api/peers` | POST | Add a new peer |
| `/api/peers/{peer_id}` | GET | Get peer configuration |
| `/api/peers/{peer_id}` | PUT | Update peer configuration |
| `/api/peers/{peer_id}` | DELETE | Remove a peer |
| `/api/peers/{peer_id}/sync` | POST | Trigger sync for a peer |
| `/api/peers/{peer_id}/status` | GET | Get sync status |
| `/api/peers/{peer_id}/enable` | POST | Enable a peer |
| `/api/peers/{peer_id}/disable` | POST | Disable a peer |
| `/api/peers/sync` | POST | Sync all enabled peers |

## Security Considerations

### Principle of Least Privilege

The `FEDERATION_STATIC_TOKEN` grants access only to federation-scoped endpoints:

| Accessible | Not Accessible |
|------------|----------------|
| `/api/v1/federation/*` | `/api/servers/*` |
| `/api/peers/*` | `/api/agents/*` |
| | `/api/admin/*` |
| | `/v0.1/*` |

### Token Security

1. **Encryption at rest**: Tokens are Fernet-encrypted in MongoDB
2. **Transport security**: All federation traffic uses HTTPS
3. **Token rotation**: Admin API available for rotating tokens without restart
4. **Revocation**: Immediate revocation via admin endpoint

### Read-Only Federated Items

Synced servers and agents are marked as read-only:
- Cannot be modified through the local registry API
- Cannot be deleted (controlled by upstream)
- Must be managed at the source registry

## Comparison with External Registry Federation

| Aspect | Peer-to-Peer Federation | External Registry (Anthropic, ASOR) |
|--------|------------------------|-------------------------------------|
| Protocol | Same codebase, symmetric | Third-party APIs |
| Direction | Bidirectional | Import only |
| Authentication | Static token or OAuth2 | Provider-specific |
| Sync control | Fine-grained (whitelist, tags) | Configuration-based |
| Visibility | Configurable per item | All-or-nothing |
| Path handling | Namespaced by peer_id | Tagged by source |

## Related Documentation

- [Federation Operational Guide](../federation-operational-guide.md) - Setup and operations
- [Federation Guide](../federation.md) - External registry integration (Anthropic, ASOR)
- [Static Token Auth](../static-token-auth.md) - Static token authentication
- [Storage Architecture](storage-architecture-mongodb-documentdb.md) - Database design
