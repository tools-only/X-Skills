---
name: databricks-multi-env-setup
description: |
  Configure Databricks across development, staging, and production environments.
  Use when setting up multi-environment deployments, configuring per-environment secrets,
  or implementing environment-specific Databricks configurations.
  Trigger with phrases like "databricks environments", "databricks staging",
  "databricks dev prod", "databricks environment setup", "databricks config by env".
allowed-tools: Read, Write, Edit, Bash(databricks:*), Bash(terraform:*)
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Databricks Multi-Environment Setup

## Overview
Configure Databricks across development, staging, and production environments.

## Prerequisites
- Multiple Databricks workspaces (or single workspace with isolation)
- Secret management solution (Databricks Secrets, Azure Key Vault, etc.)
- CI/CD pipeline configured
- Unity Catalog for cross-workspace governance

## Environment Strategy

| Environment | Purpose | Workspace | Data | Permissions |
|-------------|---------|-----------|------|-------------|
| Development | Local dev, experimentation | Shared dev workspace | Sample/synthetic | Broad access |
| Staging | Integration testing, UAT | Dedicated staging | Prod clone (masked) | Team access |
| Production | Live workloads | Dedicated prod | Real data | Restricted |

## Instructions

### Step 1: Asset Bundle Environment Configuration

```yaml
# databricks.yml
bundle:
  name: data-platform

variables:
  catalog:
    description: Unity Catalog name
    default: dev_catalog
  schema_prefix:
    description: Schema prefix for isolation
    default: dev
  warehouse_size:
    description: SQL Warehouse size
    default: "2X-Small"
  cluster_size:
    description: Default cluster workers
    default: 1

include:
  - resources/*.yml

workspace:
  host: ${DATABRICKS_HOST}

targets:
  # Development - personal workspaces
  dev:
    default: true
    mode: development
    variables:
      catalog: dev_catalog
      schema_prefix: "${workspace.current_user.short_name}"
      warehouse_size: "2X-Small"
      cluster_size: 1
    workspace:
      root_path: /Users/${workspace.current_user.userName}/.bundle/${bundle.name}/dev

  # Staging - shared test environment
  staging:
    mode: development
    variables:
      catalog: staging_catalog
      schema_prefix: staging
      warehouse_size: "Small"
      cluster_size: 2
    workspace:
      host: ${DATABRICKS_HOST_STAGING}
      root_path: /Shared/.bundle/${bundle.name}/staging
    run_as:
      service_principal_name: ${STAGING_SERVICE_PRINCIPAL}
    permissions:
      - level: CAN_VIEW
        group_name: developers
      - level: CAN_MANAGE_RUN
        group_name: data-engineers

  # Production - locked down
  prod:
    mode: production
    variables:
      catalog: prod_catalog
      schema_prefix: prod
      warehouse_size: "Medium"
      cluster_size: 4
    workspace:
      host: ${DATABRICKS_HOST_PROD}
      root_path: /Shared/.bundle/${bundle.name}/prod
    run_as:
      service_principal_name: ${PROD_SERVICE_PRINCIPAL}
    permissions:
      - level: CAN_VIEW
        group_name: data-consumers
      - level: CAN_MANAGE_RUN
        group_name: data-engineers
      - level: CAN_MANAGE
        service_principal_name: ${PROD_SERVICE_PRINCIPAL}
```

### Step 2: Environment-Specific Resources

```yaml
# resources/jobs.yml
resources:
  jobs:
    etl_pipeline:
      name: "${bundle.name}-etl-${bundle.target}"

      # Environment-specific schedule
      schedule:
        quartz_cron_expression: >-
          ${bundle.target == "prod" ? "0 0 6 * * ?" :
            bundle.target == "staging" ? "0 0 8 * * ?" : null}
        timezone_id: "America/New_York"
        pause_status: ${bundle.target == "dev" ? "PAUSED" : "UNPAUSED"}

      # Environment-specific notifications
      email_notifications:
        on_failure: ${bundle.target == "prod" ?
          ["oncall@company.com", "pagerduty@company.pagerduty.com"] :
          ["team@company.com"]}

      # Environment-specific cluster sizing
      job_clusters:
        - job_cluster_key: etl_cluster
          new_cluster:
            spark_version: "14.3.x-scala2.12"
            node_type_id: >-
              ${bundle.target == "prod" ? "Standard_DS4_v2" : "Standard_DS3_v2"}
            num_workers: ${var.cluster_size}
            autoscale:
              min_workers: ${bundle.target == "prod" ? 2 : 1}
              max_workers: ${bundle.target == "prod" ? 10 : 4}

            # Spot instances for non-prod
            azure_attributes:
              availability: >-
                ${bundle.target == "prod" ? "ON_DEMAND_AZURE" : "SPOT_AZURE"}
              first_on_demand: 1
```

### Step 3: Unity Catalog Cross-Environment Setup

```sql
-- Create environment-specific catalogs
CREATE CATALOG IF NOT EXISTS dev_catalog;
CREATE CATALOG IF NOT EXISTS staging_catalog;
CREATE CATALOG IF NOT EXISTS prod_catalog;

-- Grant cross-environment read access for data lineage
GRANT USAGE ON CATALOG prod_catalog TO `staging-service-principal`;
GRANT SELECT ON CATALOG prod_catalog TO `staging-service-principal`;

-- Set up data sharing between environments
CREATE SHARE IF NOT EXISTS prod_to_staging;
ALTER SHARE prod_to_staging ADD SCHEMA prod_catalog.reference;

-- Create recipient for staging workspace
CREATE RECIPIENT IF NOT EXISTS staging_workspace
  USING IDENTITY ('staging-workspace-identity');

GRANT SELECT ON SHARE prod_to_staging TO RECIPIENT staging_workspace;
```

### Step 4: Secret Management by Environment

```python
# src/config/secrets.py
from databricks.sdk import WorkspaceClient
import os

class EnvironmentSecrets:
    """Environment-aware secret management."""

    def __init__(self, environment: str = None):
        self.environment = environment or os.getenv("ENVIRONMENT", "dev")
        self.w = WorkspaceClient()
        self._secret_scope = f"{self.environment}-secrets"

    def get_secret(self, key: str) -> str:
        """Get secret for current environment."""
        # In notebooks, use dbutils
        # return dbutils.secrets.get(scope=self._secret_scope, key=key)

        # Via API (for testing)
        return self.w.secrets.get_secret(
            scope=self._secret_scope,
            key=key
        ).value

    def get_database_url(self) -> str:
        """Get environment-specific database URL."""
        return self.get_secret("database_url")

    def get_api_key(self, service: str) -> str:
        """Get API key for service."""
        return self.get_secret(f"{service}_api_key")

# Usage in notebooks
# secrets = EnvironmentSecrets()
# db_url = secrets.get_database_url()
```

### Step 5: Environment Detection

```python
# src/config/environment.py
import os
from dataclasses import dataclass
from typing import Optional

@dataclass
class EnvironmentConfig:
    """Environment-specific configuration."""
    name: str
    catalog: str
    schema_prefix: str
    is_production: bool
    debug_enabled: bool
    max_cluster_size: int

def detect_environment() -> EnvironmentConfig:
    """Detect current environment from context."""
    # From environment variable (set by Asset Bundles)
    env = os.getenv("ENVIRONMENT", "dev")

    # Or from Databricks tags
    # spark.conf.get("spark.databricks.tags.Environment")

    configs = {
        "dev": EnvironmentConfig(
            name="dev",
            catalog="dev_catalog",
            schema_prefix="dev",
            is_production=False,
            debug_enabled=True,
            max_cluster_size=4,
        ),
        "staging": EnvironmentConfig(
            name="staging",
            catalog="staging_catalog",
            schema_prefix="staging",
            is_production=False,
            debug_enabled=True,
            max_cluster_size=8,
        ),
        "prod": EnvironmentConfig(
            name="prod",
            catalog="prod_catalog",
            schema_prefix="prod",
            is_production=True,
            debug_enabled=False,
            max_cluster_size=20,
        ),
    }

    return configs.get(env, configs["dev"])

# Usage
env = detect_environment()
full_table_name = f"{env.catalog}.{env.schema_prefix}_sales.orders"
```

### Step 6: Environment Promotion Pipeline

```python
# scripts/promote_to_prod.py
from databricks.sdk import WorkspaceClient
import subprocess

def promote_to_production(
    version_tag: str,
    dry_run: bool = True,
) -> dict:
    """
    Promote code from staging to production.

    Steps:
    1. Verify staging tests passed
    2. Tag release in git
    3. Deploy to production
    4. Run smoke tests
    5. Enable schedules
    """
    results = {"steps": []}

    # 1. Verify staging tests
    print("Verifying staging tests...")
    staging_result = subprocess.run(
        ["databricks", "bundle", "run", "-t", "staging", "integration-tests"],
        capture_output=True
    )
    if staging_result.returncode != 0:
        raise Exception("Staging tests failed")
    results["steps"].append({"stage": "verify_staging", "status": "passed"})

    # 2. Tag release
    print(f"Tagging release {version_tag}...")
    if not dry_run:
        subprocess.run(["git", "tag", version_tag])
        subprocess.run(["git", "push", "origin", version_tag])
    results["steps"].append({"stage": "tag_release", "status": "done" if not dry_run else "skipped"})

    # 3. Deploy to production
    print("Deploying to production...")
    if not dry_run:
        subprocess.run(["databricks", "bundle", "deploy", "-t", "prod"])
    results["steps"].append({"stage": "deploy_prod", "status": "done" if not dry_run else "skipped"})

    # 4. Run smoke tests
    print("Running smoke tests...")
    if not dry_run:
        subprocess.run(["databricks", "bundle", "run", "-t", "prod", "smoke-tests"])
    results["steps"].append({"stage": "smoke_tests", "status": "done" if not dry_run else "skipped"})

    return results

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    result = promote_to_production(args.version, args.dry_run)
    print(result)
```

## Output
- Multi-environment Asset Bundle configuration
- Environment-specific job settings
- Cross-environment Unity Catalog setup
- Secret management by environment
- Promotion pipeline

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Wrong environment | Missing env var | Check ENVIRONMENT variable |
| Secret not found | Wrong scope | Verify scope name matches environment |
| Permission denied | Missing grants | Add Unity Catalog grants |
| Config mismatch | Target override issue | Check bundle target syntax |

## Examples

### Quick Environment Check
```python
env = detect_environment()
print(f"Running in {env.name} environment")
print(f"Using catalog: {env.catalog}")
print(f"Production mode: {env.is_production}")
```

### Terraform Multi-Workspace
```hcl
# infrastructure/terraform/main.tf
locals {
  environments = {
    dev = {
      workspace_name = "data-platform-dev"
      sku           = "premium"
    }
    staging = {
      workspace_name = "data-platform-staging"
      sku           = "premium"
    }
    prod = {
      workspace_name = "data-platform-prod"
      sku           = "premium"
    }
  }
}

resource "azurerm_databricks_workspace" "workspace" {
  for_each = local.environments

  name                = each.value.workspace_name
  resource_group_name = azurerm_resource_group.main.name
  location            = azurerm_resource_group.main.location
  sku                 = each.value.sku

  tags = {
    environment = each.key
    managed_by  = "terraform"
  }
}
```

## Resources
- [Asset Bundles Targets](https://docs.databricks.com/dev-tools/bundles/settings.html#targets)
- [Unity Catalog Isolation](https://docs.databricks.com/data-governance/unity-catalog/best-practices.html)
- [Workspace Federation](https://docs.databricks.com/administration-guide/workspace/index.html)

## Next Steps
For observability setup, see `databricks-observability`.
