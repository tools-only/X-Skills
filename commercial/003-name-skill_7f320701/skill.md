---
name: databricks-enterprise-rbac
description: |
  Configure Databricks enterprise SSO, Unity Catalog RBAC, and organization management.
  Use when implementing SSO integration, configuring role-based permissions,
  or setting up organization-level controls with Unity Catalog.
  Trigger with phrases like "databricks SSO", "databricks RBAC",
  "databricks enterprise", "unity catalog permissions", "databricks SCIM".
allowed-tools: Read, Write, Edit, Bash(databricks:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Databricks Enterprise RBAC

## Overview
Configure enterprise-grade access control using Unity Catalog and SSO.

## Prerequisites
- Databricks Enterprise or Premium tier
- Unity Catalog enabled
- Identity Provider (IdP) with SAML/SCIM support
- Account admin access

## Role Definitions

| Role | Unity Catalog Privileges | Use Case |
|------|-------------------------|----------|
| Data Admin | ALL PRIVILEGES on catalogs | Platform administrators |
| Data Engineer | CREATE, MODIFY, SELECT on schemas | Pipeline development |
| Data Analyst | SELECT on gold tables | Analytics and reporting |
| Data Scientist | SELECT + ML privileges | Model training |
| Viewer | SELECT on specific views | Business stakeholders |

## Instructions

### Step 1: Configure Identity Provider (SCIM)

```python
# scripts/setup_scim.py
from databricks.sdk import AccountClient
from databricks.sdk.service.iam import (
    Group,
    User,
    ServicePrincipal,
)

def setup_account_groups(client: AccountClient) -> dict:
    """
    Set up standard account-level groups.

    These groups are synced from IdP via SCIM.
    """
    standard_groups = [
        {
            "display_name": "data-platform-admins",
            "description": "Full platform administration"
        },
        {
            "display_name": "data-engineers",
            "description": "Data pipeline development"
        },
        {
            "display_name": "data-analysts",
            "description": "Analytics and reporting"
        },
        {
            "display_name": "data-scientists",
            "description": "ML and experimentation"
        },
        {
            "display_name": "data-viewers",
            "description": "Read-only business users"
        },
    ]

    created_groups = {}
    for group_config in standard_groups:
        group = client.groups.create(
            display_name=group_config["display_name"],
        )
        created_groups[group_config["display_name"]] = group.id

    return created_groups

# SCIM Configuration for Azure AD
scim_config = {
    "token_endpoint": "https://accounts.azuredatabricks.net/scim/v2",
    "bearer_token": "<account-level-pat>",
    "attribute_mappings": {
        "userName": "userPrincipalName",
        "displayName": "displayName",
        "emails[type eq \"work\"].value": "mail",
        "groups": "memberOf",
    }
}
```

### Step 2: Unity Catalog Permission Model

```sql
-- Create catalog hierarchy
CREATE CATALOG IF NOT EXISTS enterprise_data;

-- Create standard schemas
CREATE SCHEMA IF NOT EXISTS enterprise_data.bronze
COMMENT 'Raw ingested data';

CREATE SCHEMA IF NOT EXISTS enterprise_data.silver
COMMENT 'Cleansed and conformed data';

CREATE SCHEMA IF NOT EXISTS enterprise_data.gold
COMMENT 'Business-ready aggregations';

CREATE SCHEMA IF NOT EXISTS enterprise_data.ml_features
COMMENT 'ML feature store';

CREATE SCHEMA IF NOT EXISTS enterprise_data.sandbox
COMMENT 'Experimentation sandbox';

-- Grant catalog-level permissions
-- Data Admins: Full control
GRANT ALL PRIVILEGES ON CATALOG enterprise_data TO `data-platform-admins`;

-- Data Engineers: Create/modify in bronze/silver, read gold
GRANT USAGE ON CATALOG enterprise_data TO `data-engineers`;
GRANT CREATE, USAGE ON SCHEMA enterprise_data.bronze TO `data-engineers`;
GRANT CREATE, USAGE ON SCHEMA enterprise_data.silver TO `data-engineers`;
GRANT USAGE, SELECT ON SCHEMA enterprise_data.gold TO `data-engineers`;
GRANT ALL PRIVILEGES ON SCHEMA enterprise_data.sandbox TO `data-engineers`;

-- Data Analysts: Read gold, create in sandbox
GRANT USAGE ON CATALOG enterprise_data TO `data-analysts`;
GRANT USAGE, SELECT ON SCHEMA enterprise_data.gold TO `data-analysts`;
GRANT ALL PRIVILEGES ON SCHEMA enterprise_data.sandbox TO `data-analysts`;

-- Data Scientists: Read silver/gold, manage ML features
GRANT USAGE ON CATALOG enterprise_data TO `data-scientists`;
GRANT USAGE, SELECT ON SCHEMA enterprise_data.silver TO `data-scientists`;
GRANT USAGE, SELECT ON SCHEMA enterprise_data.gold TO `data-scientists`;
GRANT ALL PRIVILEGES ON SCHEMA enterprise_data.ml_features TO `data-scientists`;
GRANT ALL PRIVILEGES ON SCHEMA enterprise_data.sandbox TO `data-scientists`;

-- Viewers: Read-only on gold
GRANT USAGE ON CATALOG enterprise_data TO `data-viewers`;
GRANT USAGE, SELECT ON SCHEMA enterprise_data.gold TO `data-viewers`;
```

### Step 3: Workspace-Level Permissions

```python
# scripts/workspace_permissions.py
from databricks.sdk import WorkspaceClient
from databricks.sdk.service.iam import PermissionLevel

def configure_workspace_permissions(w: WorkspaceClient):
    """Configure workspace-level permissions."""

    # Cluster permissions
    cluster_policies = {
        "data-engineers": "CAN_USE",
        "data-scientists": "CAN_USE",
        "data-analysts": "CAN_ATTACH_TO",
    }

    # Get all cluster policies
    for policy in w.cluster_policies.list():
        for group, permission in cluster_policies.items():
            w.permissions.update(
                object_type="cluster-policies",
                object_id=policy.policy_id,
                access_control_list=[{
                    "group_name": group,
                    "permission_level": permission,
                }]
            )

    # SQL Warehouse permissions
    for warehouse in w.warehouses.list():
        w.permissions.update(
            object_type="sql/warehouses",
            object_id=warehouse.id,
            access_control_list=[
                {"group_name": "data-analysts", "permission_level": "CAN_USE"},
                {"group_name": "data-engineers", "permission_level": "CAN_MANAGE"},
            ]
        )

    # Workspace folder permissions
    folder_permissions = {
        "/Shared/Production": {
            "data-platform-admins": "CAN_MANAGE",
            "data-engineers": "CAN_RUN",
            "data-analysts": "CAN_READ",
        },
        "/Shared/Development": {
            "data-engineers": "CAN_MANAGE",
            "data-scientists": "CAN_MANAGE",
            "data-analysts": "CAN_EDIT",
        },
    }

    for folder, permissions in folder_permissions.items():
        acl = [
            {"group_name": group, "permission_level": level}
            for group, level in permissions.items()
        ]
        w.workspace.set_permissions(path=folder, access_control_list=acl)
```

### Step 4: Service Principal Management

```python
# scripts/service_principals.py
from databricks.sdk import AccountClient, WorkspaceClient
from databricks.sdk.service.iam import ServicePrincipal

def create_service_principal(
    account_client: AccountClient,
    workspace_client: WorkspaceClient,
    name: str,
    purpose: str,
    workspace_permissions: dict,
) -> dict:
    """
    Create and configure a service principal.

    Args:
        account_client: Account-level client
        workspace_client: Workspace-level client
        name: Service principal name
        purpose: Description/purpose
        workspace_permissions: Dict of {object_type: permission_level}

    Returns:
        Service principal details
    """
    # Create at account level
    sp = account_client.service_principals.create(
        display_name=name,
        active=True,
    )

    # Create OAuth secret
    secret = account_client.service_principal_secrets.create(
        service_principal_id=sp.id,
    )

    # Add to workspace groups based on purpose
    group_mappings = {
        "etl-pipeline": "data-engineers",
        "ml-training": "data-scientists",
        "reporting": "data-analysts",
    }

    if purpose in group_mappings:
        group_name = group_mappings[purpose]
        # Add SP to group (via SCIM or API)

    # Set workspace-level permissions
    for obj_type, permission in workspace_permissions.items():
        workspace_client.permissions.update(
            object_type=obj_type,
            object_id="*",  # All objects of type
            access_control_list=[{
                "service_principal_name": sp.application_id,
                "permission_level": permission,
            }]
        )

    return {
        "id": sp.id,
        "application_id": sp.application_id,
        "client_secret": secret.secret,  # Store securely!
    }

# Create standard service principals
service_principals = [
    {
        "name": "etl-pipeline-sp",
        "purpose": "etl-pipeline",
        "permissions": {
            "clusters": "CAN_ATTACH_TO",
            "jobs": "CAN_MANAGE",
        }
    },
    {
        "name": "ml-training-sp",
        "purpose": "ml-training",
        "permissions": {
            "clusters": "CAN_MANAGE",
            "mlflow-experiments": "CAN_MANAGE",
        }
    },
]
```

### Step 5: Audit and Compliance

```sql
-- Query audit logs for access patterns
SELECT
    event_time,
    user_identity.email as user_email,
    service_name,
    action_name,
    request_params,
    response.status_code
FROM system.access.audit
WHERE event_date > current_date() - INTERVAL 7 DAYS
  AND service_name IN ('unityCatalog', 'clusters', 'jobs')
ORDER BY event_time DESC
LIMIT 100;

-- Permission changes audit
SELECT
    event_time,
    user_identity.email as changed_by,
    action_name,
    request_params:full_name_arg as resource,
    request_params:changes as permission_changes
FROM system.access.audit
WHERE action_name LIKE '%Grant%' OR action_name LIKE '%Revoke%'
  AND event_date > current_date() - INTERVAL 30 DAYS
ORDER BY event_time DESC;

-- Privileged access report
SELECT DISTINCT
    grantee,
    privilege,
    object_type,
    object_name
FROM system.information_schema.object_privileges
WHERE privilege IN ('ALL PRIVILEGES', 'CREATE', 'MODIFY')
ORDER BY grantee, object_type;
```

### Step 6: Access Review Automation

```python
# scripts/access_review.py
from databricks.sdk import WorkspaceClient
from datetime import datetime
import pandas as pd

def generate_access_review_report(w: WorkspaceClient, spark) -> pd.DataFrame:
    """
    Generate access review report for compliance.

    Returns DataFrame with all permissions by user/group.
    """
    # Get Unity Catalog permissions
    uc_permissions = spark.sql("""
        SELECT
            grantee,
            grantee_type,
            privilege,
            object_type,
            object_name,
            inherited_from
        FROM system.information_schema.object_privileges
        ORDER BY grantee, object_type
    """).toPandas()

    # Get workspace permissions
    workspace_permissions = []
    for obj_type in ['clusters', 'jobs', 'sql/warehouses']:
        for obj in getattr(w, obj_type.replace('/', '_').replace('sql_', '')).list():
            perms = w.permissions.get(object_type=obj_type, object_id=obj.id)
            for acl in perms.access_control_list or []:
                workspace_permissions.append({
                    "object_type": obj_type,
                    "object_id": obj.id,
                    "principal": acl.user_name or acl.group_name or acl.service_principal_name,
                    "permission": acl.all_permissions[0].permission_level if acl.all_permissions else None,
                })

    ws_df = pd.DataFrame(workspace_permissions)

    # Combine reports
    report = {
        "generated_at": datetime.now().isoformat(),
        "unity_catalog_permissions": uc_permissions,
        "workspace_permissions": ws_df,
    }

    return report

def flag_excessive_permissions(report: dict) -> list[dict]:
    """Flag potentially excessive permissions."""
    flags = []

    uc_df = report["unity_catalog_permissions"]

    # Flag ALL PRIVILEGES
    excessive = uc_df[uc_df["privilege"] == "ALL PRIVILEGES"]
    for _, row in excessive.iterrows():
        if row["grantee_type"] != "GROUP":  # Individual users with ALL PRIVILEGES
            flags.append({
                "type": "EXCESSIVE_PRIVILEGE",
                "grantee": row["grantee"],
                "privilege": row["privilege"],
                "object": row["object_name"],
                "recommendation": "Consider using specific privileges instead of ALL PRIVILEGES",
            })

    # Flag direct user grants (should use groups)
    direct_grants = uc_df[uc_df["grantee_type"] == "USER"]
    unique_users = direct_grants["grantee"].unique()
    for user in unique_users:
        user_grants = len(direct_grants[direct_grants["grantee"] == user])
        if user_grants > 5:
            flags.append({
                "type": "DIRECT_USER_GRANTS",
                "grantee": user,
                "grant_count": user_grants,
                "recommendation": "Consider using groups for permission management",
            })

    return flags
```

## Output
- SCIM integration configured
- Unity Catalog RBAC implemented
- Workspace permissions set
- Service principals created
- Audit logging enabled

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| SSO login fails | Wrong callback URL | Verify IdP config |
| Permission denied | Missing grant | Check inherited permissions |
| SCIM sync failed | Token expired | Regenerate SCIM token |
| Group not found | Case sensitivity | Match exact group name |

## Examples

### Quick Permission Check
```sql
-- Check what a user can access
SELECT * FROM system.information_schema.object_privileges
WHERE grantee = 'user@company.com'
ORDER BY object_type, object_name;
```

### Grant Template
```sql
-- Standard data consumer access
GRANT USAGE ON CATALOG enterprise_data TO `new-team`;
GRANT USAGE ON SCHEMA enterprise_data.gold TO `new-team`;
GRANT SELECT ON SCHEMA enterprise_data.gold TO `new-team`;
```

## Resources
- [Unity Catalog Privileges](https://docs.databricks.com/data-governance/unity-catalog/manage-privileges/privileges.html)
- [SCIM Provisioning](https://docs.databricks.com/administration-guide/users-groups/scim/index.html)
- [Identity Federation](https://docs.databricks.com/administration-guide/users-groups/best-practices.html)

## Next Steps
For major migrations, see `databricks-migration-deep-dive`.
