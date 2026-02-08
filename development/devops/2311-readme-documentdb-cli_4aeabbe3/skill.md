# DocumentDB CLI Tools

Command-line tools for inspecting and managing DocumentDB collections in the MCP Gateway Registry.

## Overview

The DocumentDB CLI provides commands to:
- List all collections in the database
- Inspect collection schemas and statistics
- Count documents in collections
- Search and query documents
- View sample documents

## Files

- [`manage-documentdb.py`](../../../scripts/manage-documentdb.py) - Python script that performs DocumentDB operations
- [`run-documentdb-cli.sh`](run-documentdb-cli.sh) - Shell wrapper that runs the Python script inside an ECS task with VPC access

## Usage

### Prerequisites

- AWS credentials configured
- DocumentDB endpoint stored in SSM Parameter Store at `/mcp-gateway/documentdb/endpoint`
- DocumentDB credentials stored in Secrets Manager at `mcp-gateway/documentdb/credentials`
- ECS cluster and task definition deployed

### Commands

#### List All Collections

```bash
./terraform/aws-ecs/scripts/run-documentdb-cli.sh list
```

**Output:**
```
Found 6 collections in database 'mcp_registry'
================================================================================

Collection: mcp_agents_default
  Documents: 12
  Size: 0.05 MB

Collection: mcp_embeddings_1536_default
  Documents: 156
  Size: 2.34 MB

Collection: mcp_scopes_default
  Documents: 8
  Size: 0.02 MB

Collection: mcp_servers_default
  Documents: 24
  Size: 0.15 MB
```

#### Inspect Collection Schema and Stats

```bash
./terraform/aws-ecs/scripts/run-documentdb-cli.sh inspect mcp_servers_default
```

**Output:**
```
Collection: mcp_servers_default
================================================================================

Document Count: 24

--- Collection Statistics ---
Size: 0.15 MB
Storage Size: 0.25 MB
Total Index Size: 0.08 MB
Average Object Size: 6234 bytes

--- Indexes ---

Index: _id_
  Keys: {
    "_id": 1
  }

Index: path_1
  Keys: {
    "path": 1
  }
  Unique: True

--- Sample Document Schema ---
{
  "_id": "ObjectId",
  "path": "str",
  "name": "str",
  "enabled": "bool",
  "description": "str",
  "created_at": "datetime",
  "updated_at": "datetime",
  "metadata": "dict"
}
```

#### Count Documents

```bash
./terraform/aws-ecs/scripts/run-documentdb-cli.sh count mcp_servers_default
```

**Output:**
```
Collection: mcp_servers_default
Document Count: 24
```

#### Search Documents (List)

```bash
# Show first 10 documents (default)
./terraform/aws-ecs/scripts/run-documentdb-cli.sh search mcp_servers_default

# Show first 20 documents
./terraform/aws-ecs/scripts/run-documentdb-cli.sh search mcp_servers_default 20
```

**Output:**
```
Collection: mcp_servers_default
Showing 10 documents (limit: 10)
================================================================================

--- Document 1 ---
{
  "_id": "507f1f77bcf86cd799439011",
  "path": "/currenttime",
  "name": "CurrentTime Server",
  "enabled": true,
  "description": "Returns current time in various formats",
  "created_at": "2024-01-15T10:30:00Z",
  ...
}

--- Document 2 ---
...
```

#### Sample Document

```bash
./terraform/aws-ecs/scripts/run-documentdb-cli.sh sample mcp_servers_default
```

**Output:**
```
Collection: mcp_servers_default
Sample Document:
================================================================================
{
  "_id": "507f1f77bcf86cd799439011",
  "path": "/currenttime",
  "name": "CurrentTime Server",
  "enabled": true,
  "description": "Returns current time in various formats",
  ...
}
```

#### Query with Filter

```bash
# Find enabled servers
./terraform/aws-ecs/scripts/run-documentdb-cli.sh query mcp_servers_default '{"enabled": true}'

# Find server by path
./terraform/aws-ecs/scripts/run-documentdb-cli.sh query mcp_servers_default '{"path": "/currenttime"}'

# Query with limit
./terraform/aws-ecs/scripts/run-documentdb-cli.sh query mcp_servers_default '{"enabled": true}' 5
```

**Output:**
```
Collection: mcp_servers_default
Filter: {"enabled": true}
Found 18 documents (limit: 10)
================================================================================

--- Document 1 ---
{
  "_id": "507f1f77bcf86cd799439011",
  "path": "/currenttime",
  "enabled": true,
  ...
}
```

## Environment Variables

- `DOCUMENTDB_HOST` - Override DocumentDB endpoint (optional, read from SSM if not set)
- `AWS_REGION` - AWS region (default: us-east-1)

## How It Works

1. The shell script reads DocumentDB connection details from AWS services:
   - Endpoint from SSM Parameter Store (`/mcp-gateway/documentdb/endpoint`)
   - Credentials from Secrets Manager (`mcp-gateway/documentdb/credentials`)
   - VPC configuration from the registry ECS service

2. It launches an ECS Fargate task using the `mcp-gateway-v2-registry` task definition

3. The task runs inside the VPC with network access to DocumentDB

4. The Python script executes the requested command and outputs results

5. Logs are retrieved from CloudWatch and displayed

## Common Collections

- `mcp_servers_default` - MCP server registrations
- `mcp_agents_default` - Agent registrations
- `mcp_scopes_default` - Authorization scope definitions
- `mcp_embeddings_1536_default` - Vector embeddings for semantic search
- `mcp_groups_default` - Group definitions
- `mcp_security_scans_default` - Security scan results

## Troubleshooting

### No logs found

If you see "No logs found", the task may have failed to start. Check:
1. Task definition exists: `aws ecs describe-task-definition --task-definition mcp-gateway-v2-registry`
2. Network configuration is correct
3. DocumentDB credentials are valid

### Connection timeout

If the task hangs or times out:
1. Verify security groups allow traffic to DocumentDB on port 27017
2. Verify task is running in the same VPC as DocumentDB
3. Check DocumentDB cluster status

### Invalid filter JSON

For query commands, ensure the filter is valid JSON:
```bash
# Correct
./run-documentdb-cli.sh query mcp_servers_default '{"enabled": true}'

# Incorrect (missing quotes around JSON)
./run-documentdb-cli.sh query mcp_servers_default {"enabled": true}
```

## Direct Python Usage (Local)

If you have direct network access to DocumentDB (e.g., VPN, bastion host):

```bash
# Set environment variables
export DOCUMENTDB_HOST=your-cluster.docdb.amazonaws.com
export DOCUMENTDB_USERNAME=admin
export DOCUMENTDB_PASSWORD=yourpassword
export DOCUMENTDB_DATABASE=mcp_registry
export DOCUMENTDB_USE_TLS=true
export DOCUMENTDB_TLS_CA_FILE=/path/to/global-bundle.pem

# Run commands directly
cd scripts
python manage-documentdb.py list
python manage-documentdb.py inspect --collection mcp_servers_default
python manage-documentdb.py search --collection mcp_servers_default --limit 5
```

## See Also

- [OpenSearch CLI](run-aoss-cli.sh) - Similar tool for OpenSearch Serverless indexes
- [DocumentDB Initialization](run-documentdb-init.sh) - Initialize DocumentDB indexes and load scopes
- [View CloudWatch Logs](view-cloudwatch-logs.sh) - View ECS service logs
