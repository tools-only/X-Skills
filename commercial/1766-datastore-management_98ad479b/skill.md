# Datastore Management Guide

**Last Updated:** January 3, 2026
**Applies to:** All storage backends (MongoDB CE, AWS DocumentDB)

---

## Table of Contents

1. [Overview](#overview)
2. [Local Development - MongoDB CE](#local-development---mongodb-ce)
3. [Production - AWS DocumentDB via ECS](#production---aws-documentdb-via-ecs)
4. [Common Operations](#common-operations)
5. [Troubleshooting](#troubleshooting)

---

## Overview

The MCP Gateway Registry uses MongoDB-compatible datastores for storage:

- **Local Development:** MongoDB Community Edition 8.2 in Docker
- **Production:** AWS DocumentDB (MongoDB-compatible managed service)

This guide explains how to access and manage datastores in both environments.

---

## Local Development - MongoDB CE

### Prerequisites

- MongoDB container running: `docker compose ps mongodb` shows "healthy"
- No authentication required (configured for local dev simplicity)

### Accessing the Datastore (mongosh)

#### Method 1: Direct Docker Exec (Recommended)

```bash
# Connect to MongoDB shell
docker exec -it mcp-mongodb mongosh

# You should see:
# Current Mongosh Log ID: ...
# Connecting to: mongodb://127.0.0.1:27017/?directConnection=true
# ...
# rs0 [direct: primary] test>
```

#### Method 2: Connect from Host Machine

If you have `mongosh` installed locally:

```bash
mongosh mongodb://localhost:27017/mcp_registry
```

### Basic Datastore Operations

Once connected to mongosh:

```javascript
// Switch to the registry database
use mcp_registry

// List all collections
show collections
// Expected output:
//   mcp_agents_default
//   mcp_embeddings_1536_default
//   mcp_federation_config_default
//   mcp_scopes_default
//   mcp_security_scans_default
//   mcp_servers_default

// Check replica set status
rs.status()

// View database statistics
db.stats()
```

### Viewing Collection Contents

#### List All Servers

```javascript
// Count total servers
db.mcp_servers_default.countDocuments()

// View all servers (formatted)
db.mcp_servers_default.find().pretty()

// View specific server by path
db.mcp_servers_default.findOne({ path: "/servers/financial-data" })

// List only server names and paths
db.mcp_servers_default.find(
  {},
  { "manifest.serverInfo.name": 1, path: 1, _id: 0 }
)
```

#### List All Agents

```javascript
// Count total agents
db.mcp_agents_default.countDocuments()

// View all agents
db.mcp_agents_default.find().pretty()

// Find agents by tag
db.mcp_agents_default.find({ tags: "finance" }).pretty()
```

#### View Vector Embeddings

```javascript
// Count embeddings
db.mcp_embeddings_1536_default.countDocuments()

// View embedding metadata (without the large vector array)
db.mcp_embeddings_1536_default.find(
  {},
  {
    path: 1,
    entity_type: 1,
    name: 1,
    embedding_metadata: 1,
    indexed_at: 1,
    _id: 0
  }
).pretty()

// Check specific embedding
db.mcp_embeddings_1536_default.findOne({ path: "/servers/financial-data" })
```

#### View Scopes

```javascript
// List all scopes
db.mcp_scopes_default.find().pretty()

// Find server scopes
db.mcp_scopes_default.find({ scope_type: "server_scope" }).pretty()

// Find group mappings
db.mcp_scopes_default.find({ scope_type: "group_mapping" }).pretty()
```

#### View Security Scans

```javascript
// Count security scans
db.mcp_security_scans_default.countDocuments()

// View latest scans
db.mcp_security_scans_default.find().sort({ scan_timestamp: -1 }).limit(5).pretty()

// Find scans for specific server
db.mcp_security_scans_default.find({ server_path: "/servers/financial-data" }).pretty()
```

### Collection Indexes

```javascript
// View indexes on servers collection
db.mcp_servers_default.getIndexes()

// View indexes on embeddings collection
db.mcp_embeddings_1536_default.getIndexes()

// Check index usage stats
db.mcp_servers_default.aggregate([{ $indexStats: {} }])
```

### Query Performance Analysis

```javascript
// Explain query execution plan
db.mcp_servers_default.find({ path: "/servers/financial-data" }).explain("executionStats")

// Find slow operations (if profiling enabled)
db.system.profile.find({ millis: { $gt: 100 } }).sort({ ts: -1 }).limit(5).pretty()
```

### Exiting mongosh

```javascript
// Exit the shell
exit
```

Or press `Ctrl+D`

---

## Production - AWS DocumentDB via ECS

### Prerequisites

- AWS ECS cluster running with DocumentDB
- ECS exec permissions configured
- `manage-documentdb.py` script available in registry container

### Accessing DocumentDB via ECS Exec

#### Step 1: SSH into Registry Container

```bash
# From your local machine
# Use the ecs-ssh helper script
cd terraform/aws-ecs
./scripts/ecs-ssh.sh registry

# Or manually with AWS CLI
aws ecs execute-command \
  --cluster mcp-gateway-ecs-cluster \
  --task <task-id> \
  --container registry \
  --interactive \
  --command "/bin/bash"
```

#### Step 2: Activate Python Virtual Environment

```bash
# Inside the ECS container
source .venv/bin/activate

# Verify Python environment
which python
# Should show: /app/.venv/bin/python
```

#### Step 3: Run DocumentDB Management Script

The `manage-documentdb.py` script provides commands for managing collections and querying data.

##### List All Collections

```bash
python scripts/manage-documentdb.py list
```

##### Inspect a Collection

```bash
# Show collection schema and indexes
python scripts/manage-documentdb.py inspect --collection mcp_servers_default
```

##### Count Documents

```bash
# Count all documents in collection
python scripts/manage-documentdb.py count --collection mcp_servers_default
```

##### Search Documents

```bash
# List documents with optional limit
python scripts/manage-documentdb.py search --collection mcp_servers_default --limit 5

# Search specific collection
python scripts/manage-documentdb.py search --collection mcp_agents_default --limit 10
```

##### View Sample Document

```bash
# Show one sample document from collection
python scripts/manage-documentdb.py sample --collection mcp_servers_default
```

##### Query with Filter

```bash
# Query with MongoDB filter syntax
python scripts/manage-documentdb.py query \
  --collection mcp_servers_default \
  --filter '{"path": "/servers/financial-data"}'

# Query enabled servers
python scripts/manage-documentdb.py query \
  --collection mcp_servers_default \
  --filter '{"enabled": true}'

# Query by tags
python scripts/manage-documentdb.py query \
  --collection mcp_servers_default \
  --filter '{"tags": "finance"}'
```

##### View Embeddings

```bash
# Sample embedding document (shows structure without large vector array)
python scripts/manage-documentdb.py sample --collection mcp_embeddings_1536_default

# Count total embeddings
python scripts/manage-documentdb.py count --collection mcp_embeddings_1536_default
```

**Note:** The script automatically reads connection parameters from environment variables in the ECS container (`DOCUMENTDB_HOST`, `DOCUMENTDB_USERNAME`, `DOCUMENTDB_PASSWORD`, etc.).

---

## Common Operations

### Checking Datastore Health

#### Local (MongoDB CE)

```javascript
// In mongosh
db.serverStatus()
db.stats()
rs.status()
```

#### Production (DocumentDB)

```bash
# Use count command to verify connection and check collections
python scripts/manage-documentdb.py list

# Check specific collection
python scripts/manage-documentdb.py count --collection mcp_servers_default
```

### Searching for Specific Documents

#### Local (MongoDB CE)

```javascript
// Search servers by tag
db.mcp_servers_default.find({ tags: "finance" })

// Search by partial name match
db.mcp_servers_default.find({
  "manifest.serverInfo.name": /financial/i
})

// Complex query with multiple conditions
db.mcp_servers_default.find({
  enabled: true,
  tags: { $in: ["finance", "data"] }
})
```

#### Production (DocumentDB)

```bash
# Search by tags
python scripts/manage-documentdb.py query \
  --collection mcp_servers_default \
  --filter '{"tags": "finance"}'

# Query enabled servers
python scripts/manage-documentdb.py query \
  --collection mcp_servers_default \
  --filter '{"enabled": true}'
```

### Viewing Recent Activity

#### Local (MongoDB CE)

```javascript
// Recent server registrations
db.mcp_servers_default.find().sort({ registered_at: -1 }).limit(5)

// Recent embeddings
db.mcp_embeddings_1536_default.find().sort({ indexed_at: -1 }).limit(5)

// Recent security scans
db.mcp_security_scans_default.find().sort({ scan_timestamp: -1 }).limit(5)
```

#### Production (DocumentDB)

```bash
# View recent servers (sorted by registration)
python scripts/manage-documentdb.py search \
  --collection mcp_servers_default \
  --limit 5

# View recent embeddings
python scripts/manage-documentdb.py search \
  --collection mcp_embeddings_1536_default \
  --limit 5
```

### Backup and Export

#### Local (MongoDB CE)

##### Option 1: Binary Backup (mongodump) - Recommended for Full Backups

```bash
# Export entire database (BSON format - preserves data types)
docker exec mcp-mongodb mongodump \
  --db=mcp_registry \
  --out=/tmp/mongodb-backup

# Copy backup from container to host
docker cp mcp-mongodb:/tmp/mongodb-backup ./mongodb-backup-$(date +%Y%m%d)

# Restore from backup (if needed)
docker cp ./mongodb-backup-20260103 mcp-mongodb:/tmp/restore-backup
docker exec mcp-mongodb mongorestore \
  --db=mcp_registry \
  /tmp/restore-backup/mcp_registry
```

##### Option 2: JSON Export (mongoexport) - Human-Readable, Portable

```bash
# Export specific collection to JSON (one document per line)
docker exec mcp-mongodb mongoexport \
  --db=mcp_registry \
  --collection=mcp_servers_default \
  --out=/tmp/servers.json

# Copy to host
docker cp mcp-mongodb:/tmp/servers.json ./servers-backup-$(date +%Y%m%d).json

# Pretty-print JSON (optional, for readability)
docker exec mcp-mongodb mongoexport \
  --db=mcp_registry \
  --collection=mcp_servers_default \
  --jsonArray \
  --pretty \
  --out=/tmp/servers-pretty.json

# Import from JSON (if needed)
docker cp ./servers-backup-20260103.json mcp-mongodb:/tmp/import-servers.json
docker exec mcp-mongodb mongoimport \
  --db=mcp_registry \
  --collection=mcp_servers_default \
  --file=/tmp/import-servers.json
```

##### Export All Collections

```bash
# Export all collections to JSON
COLLECTIONS="mcp_servers_default mcp_agents_default mcp_scopes_default mcp_embeddings_1536_default mcp_security_scans_default mcp_federation_config_default"

for collection in $COLLECTIONS; do
  echo "Exporting $collection..."
  docker exec mcp-mongodb mongoexport \
    --db=mcp_registry \
    --collection=$collection \
    --out=/tmp/${collection}.json
  docker cp mcp-mongodb:/tmp/${collection}.json ./${collection}-$(date +%Y%m%d).json
done
```

#### Production (DocumentDB)

##### Option 1: AWS Automated Backups (Recommended)

AWS DocumentDB provides automated continuous backups with point-in-time recovery:

```bash
# Create manual snapshot (from local machine with AWS CLI)
aws docdb create-db-cluster-snapshot \
  --db-cluster-snapshot-identifier mcp-registry-manual-$(date +%Y%m%d) \
  --db-cluster-identifier mcp-registry-prod

# List available snapshots
aws docdb describe-db-cluster-snapshots \
  --db-cluster-identifier mcp-registry-prod

# Restore from snapshot (creates new cluster)
aws docdb restore-db-cluster-from-snapshot \
  --db-cluster-identifier mcp-registry-restored \
  --snapshot-identifier mcp-registry-manual-20260103 \
  --engine docdb
```

##### Option 2: Binary Backup with mongodump (from ECS Container)

```bash
# SSH into ECS container
cd terraform/aws-ecs
./scripts/ecs-ssh.sh registry
source .venv/bin/activate

# Export entire database to BSON
mongodump \
  --host=$DOCUMENTDB_HOST \
  --port=27017 \
  --username=$DOCUMENTDB_USERNAME \
  --password=$DOCUMENTDB_PASSWORD \
  --ssl \
  --sslCAFile=/app/global-bundle.pem \
  --db=mcp_registry \
  --out=/tmp/documentdb-backup

# Upload to S3
BACKUP_DATE=$(date +%Y%m%d-%H%M%S)
aws s3 cp /tmp/documentdb-backup \
  s3://mcp-gateway-backups/documentdb-backup-${BACKUP_DATE}/ \
  --recursive

# Cleanup temporary files
rm -rf /tmp/documentdb-backup

echo "Backup uploaded to: s3://mcp-gateway-backups/documentdb-backup-${BACKUP_DATE}/"
```

##### Option 3: JSON Export of Specific Collections (from ECS Container)

```bash
# SSH into ECS container
cd terraform/aws-ecs
./scripts/ecs-ssh.sh registry
source .venv/bin/activate

# Export specific collection to JSON
mongoexport \
  --host=$DOCUMENTDB_HOST \
  --port=27017 \
  --username=$DOCUMENTDB_USERNAME \
  --password=$DOCUMENTDB_PASSWORD \
  --ssl \
  --sslCAFile=/app/global-bundle.pem \
  --db=mcp_registry \
  --collection=mcp_servers_default \
  --out=/tmp/servers-export.json

# Upload to S3
aws s3 cp /tmp/servers-export.json \
  s3://mcp-gateway-backups/exports/servers-$(date +%Y%m%d).json

# Cleanup
rm /tmp/servers-export.json
```

##### Restore from S3 Backup

```bash
# SSH into ECS container
cd terraform/aws-ecs
./scripts/ecs-ssh.sh registry
source .venv/bin/activate

# Download backup from S3
aws s3 cp s3://mcp-gateway-backups/documentdb-backup-20260103-120000/ \
  /tmp/restore-backup/ \
  --recursive

# Restore using mongorestore
mongorestore \
  --host=$DOCUMENTDB_HOST \
  --port=27017 \
  --username=$DOCUMENTDB_USERNAME \
  --password=$DOCUMENTDB_PASSWORD \
  --ssl \
  --sslCAFile=/app/global-bundle.pem \
  --db=mcp_registry \
  /tmp/restore-backup/mcp_registry

# Cleanup
rm -rf /tmp/restore-backup
```

**Important Notes:**

- **mongodump/mongorestore**: Binary format (BSON), preserves all data types including binary data and dates
- **mongoexport/mongoimport**: JSON format, human-readable but may lose type information
- **For production**: Use AWS automated backups for disaster recovery, manual exports for data migration
- **S3 bucket**: Replace `mcp-gateway-backups` with your actual S3 bucket name
- **Embeddings**: Vector embeddings are large; consider excluding from exports if not needed:
  ```bash
  mongodump --excludeCollection=mcp_embeddings_1536_default ...
  ```

---

## Troubleshooting

### Cannot Connect to MongoDB (Local)

**Problem:** `docker exec -it mcp-mongodb mongosh` fails

**Solutions:**

```bash
# Check if container is running
docker compose ps mongodb

# Check container logs
docker compose logs mongodb

# Restart MongoDB
docker compose restart mongodb

# If needed, recreate container
docker compose up -d mongodb
```

### Cannot Connect to DocumentDB (Production)

**Problem:** `manage-documentdb.py` commands fail with connection errors

**Solutions:**

```bash
# 1. Verify you're in the ECS container with activated venv
source .venv/bin/activate

# 2. Check environment variables are set
env | grep DOCUMENTDB

# 3. Test connection with simple list command
python scripts/manage-documentdb.py list

# 4. Check security group allows access from ECS tasks (from local machine)
aws ec2 describe-security-groups --group-ids <docdb-sg-id>

# 5. Verify DocumentDB endpoint (from local machine)
aws docdb describe-db-clusters --db-cluster-identifier mcp-registry-prod
```

### Replica Set Not Initialized (Local)

**Problem:** `rs.status()` shows "not initialized"

**Solutions:**

```bash
# Re-run initialization
docker compose up mongodb-init

# Or manually initialize
docker exec -it mcp-mongodb mongosh --eval 'rs.initiate({_id: "rs0", members: [{_id: 0, host: "mongodb:27017"}]})'
```

### Collections Not Found

**Problem:** `show collections` returns empty

**Solutions:**

```bash
# Verify you're in correct database
# In mongosh:
db.getName()  // Should show "mcp_registry"

# Re-run initialization
docker compose up mongodb-init

# Check if data is in different namespace
db.getCollectionNames()
```

### Slow Queries

**Problem:** Queries taking too long

**Solutions (Local MongoDB CE):**

```javascript
// In mongosh - Check if indexes exist
db.mcp_servers_default.getIndexes()

// Analyze query plan
db.mcp_servers_default.find({ path: "..." }).explain("executionStats")

// Check embeddings indexes
db.mcp_embeddings_1536_default.getIndexes()
```

**Solutions (Production DocumentDB):**

```bash
# Use inspect command to check indexes
python scripts/manage-documentdb.py inspect --collection mcp_servers_default

# Check embeddings collection indexes
python scripts/manage-documentdb.py inspect --collection mcp_embeddings_1536_default
```

---

## Quick Reference

### Connection Strings

**Local MongoDB CE:**
```
mongodb://localhost:27017/mcp_registry
```

**Production DocumentDB:**
```
mongodb://<username>:<password>@<cluster-endpoint>:27017/mcp_registry?tls=true&tlsCAFile=global-bundle.pem&replicaSet=rs0
```

### Common mongosh Commands

| Command | Description |
|---------|-------------|
| `show dbs` | List all databases |
| `use mcp_registry` | Switch to mcp_registry database |
| `show collections` | List all collections |
| `db.stats()` | Database statistics |
| `rs.status()` | Replica set status |
| `db.mcp_servers_default.find()` | List all servers |
| `db.mcp_servers_default.countDocuments()` | Count documents |
| `.pretty()` | Format output nicely |
| `exit` | Exit mongosh |

### Environment Variables Reference

**Local (.env):**
```bash
STORAGE_BACKEND=mongodb-ce
DOCUMENTDB_HOST=mongodb
DOCUMENTDB_PORT=27017
DOCUMENTDB_DATABASE=mcp_registry
DOCUMENTDB_NAMESPACE=default
DOCUMENTDB_USE_TLS=false
```

**Production (ECS Task Definition):**
```bash
STORAGE_BACKEND=documentdb
DOCUMENTDB_HOST=<cluster-endpoint>
DOCUMENTDB_PORT=27017
DOCUMENTDB_DATABASE=mcp_registry
DOCUMENTDB_NAMESPACE=production
DOCUMENTDB_USERNAME=<from-secrets>
DOCUMENTDB_PASSWORD=<from-secrets>
DOCUMENTDB_USE_TLS=true
DOCUMENTDB_TLS_CA_FILE=/app/global-bundle.pem
DOCUMENTDB_REPLICA_SET=rs0
```

---

## See Also

- [Storage Architecture: MongoDB CE & AWS DocumentDB](design/storage-architecture-mongodb-documentdb.md)
- [Datastore Schema Design](database-design.md)
- [Configuration Guide](configuration.md)
- [MongoDB Documentation](https://www.mongodb.com/docs/manual/)
- [AWS DocumentDB Documentation](https://docs.aws.amazon.com/documentdb/)
