# Static Token Auth for Registry API

## Overview

Static Token Auth allows Registry API endpoints (`/api/*`, `/v0.1/*`) to be accessed using a static API key instead of IdP-based JWT validation. This is useful in trusted network environments where configuring a full identity provider (Keycloak, Entra ID, Cognito) is not practical.

MCP Gateway endpoints are **not affected** and continue to require full IdP authentication.

## When to Use

- CI/CD pipelines that register or query MCP servers and agents
- CLI tooling in trusted network environments
- Development and testing environments without an IdP
- Automated scripts that interact with the Registry API

## Configuration

Two environment variables control this feature:

| Variable | Description | Default |
|----------|-------------|---------|
| `REGISTRY_STATIC_TOKEN_AUTH_ENABLED` | Enable static token auth for Registry API | `false` |
| `REGISTRY_API_TOKEN` | Static API key that clients must send as a Bearer token | (empty) |

Both must be set for the feature to activate. If `REGISTRY_STATIC_TOKEN_AUTH_ENABLED=true` but `REGISTRY_API_TOKEN` is empty, the auth server logs an error and falls back to standard IdP JWT validation.

### Generate a Token

```bash
python3 -c "import secrets; print(secrets.token_urlsafe(32))"
```

### Docker Compose

Add to `.env`:

```bash
REGISTRY_STATIC_TOKEN_AUTH_ENABLED=true
REGISTRY_API_TOKEN=your-generated-token
```

These are passed to the auth server container via `docker-compose.yml`.

### AWS ECS (Terraform)

Set in `terraform.tfvars`:

```hcl
registry_static_token_auth_enabled = true
registry_api_token                 = "your-generated-token"
```

Alternatively, set the token via environment variable to avoid storing it in a file:

```bash
export TF_VAR_registry_api_token="your-generated-token"
```

## Usage

Clients send the static API key as a Bearer token in the Authorization header:

```bash
curl -H "Authorization: Bearer your-generated-token" \
  http://localhost:7860/api/services/list
```

Using the Registry CLI:

```bash
# Save token to a file
echo "your-generated-token" > .network-trusted-token

# Use with registry_management.py
uv run python api/registry_management.py \
  --registry-url http://localhost:7860 \
  --token-file .network-trusted-token \
  list
```

## How It Works

1. The auth server checks if `REGISTRY_STATIC_TOKEN_AUTH_ENABLED` is true and the request path matches `/api/*` or `/v0.1/*`
2. If the request has a session cookie (browser/UI), the bypass is skipped and normal session auth is used
3. If no session cookie is present, the Bearer token is validated against `REGISTRY_API_TOKEN`
4. On success, the request proceeds with a `network-trusted` identity that has full admin permissions on Registry API endpoints

## Security Considerations

- The static API key is a shared secret. Treat it like a password.
- Rotate the token periodically by updating `REGISTRY_API_TOKEN` and restarting the auth server.
- This feature does not affect MCP Gateway endpoints, which always require IdP authentication.
- Use network-level controls (VPC, security groups, firewall rules) in addition to the static token.
- For production deployments with external access, prefer IdP-based authentication.
