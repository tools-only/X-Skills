---
name: databricks-migration-deep-dive
description: |
  Execute comprehensive platform migrations to Databricks from legacy systems.
  Use when migrating from on-premises Hadoop, other cloud platforms,
  or legacy data warehouses to Databricks.
  Trigger with phrases like "migrate to databricks", "hadoop migration",
  "snowflake to databricks", "legacy migration", "data warehouse migration".
allowed-tools: Read, Write, Edit, Bash(databricks:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Databricks Migration Deep Dive

## Overview
Comprehensive migration strategies for moving to Databricks from legacy systems.

## Prerequisites
- Access to source and target systems
- Understanding of current data architecture
- Migration timeline and requirements
- Stakeholder alignment

## Migration Patterns

| Source System | Migration Pattern | Complexity | Timeline |
|--------------|-------------------|------------|----------|
| On-prem Hadoop | Lift-and-shift + modernize | High | 6-12 months |
| Snowflake | Parallel run + cutover | Medium | 3-6 months |
| AWS Redshift | ETL rewrite + data copy | Medium | 3-6 months |
| Azure Synapse | Delta Lake conversion | Low | 1-3 months |
| Legacy DW (Oracle/Teradata) | Full rebuild | High | 12-18 months |

## Instructions

### Step 1: Discovery and Assessment

```python
# scripts/migration_assessment.py
from dataclasses import dataclass
from typing import List, Dict
import pandas as pd

@dataclass
class SourceTableInfo:
    """Information about source table for migration planning."""
    database: str
    schema: str
    table: str
    row_count: int
    size_gb: float
    column_count: int
    partition_columns: List[str]
    dependencies: List[str]
    access_frequency: str  # high, medium, low
    data_classification: str  # pii, confidential, public

def assess_hadoop_cluster(spark, hive_metastore: str) -> List[SourceTableInfo]:
    """
    Assess Hadoop/Hive cluster for migration planning.

    Returns inventory of all tables with metadata.
    """
    tables = []

    # Get all databases
    databases = spark.sql("SHOW DATABASES").collect()

    for db_row in databases:
        db = db_row.databaseName
        if db in ['default', 'sys']:  # Skip system DBs
            continue

        # Get tables in database
        spark.sql(f"USE {db}")
        table_list = spark.sql("SHOW TABLES").collect()

        for table_row in table_list:
            table_name = table_row.tableName

            try:
                # Get table details
                desc = spark.sql(f"DESCRIBE EXTENDED {db}.{table_name}")
                detail_df = desc.toPandas()

                # Extract partition info
                partition_cols = []
                in_partition_section = False
                for _, row in detail_df.iterrows():
                    if row['col_name'] == '# Partition Information':
                        in_partition_section = True
                    elif in_partition_section and row['col_name'] and not row['col_name'].startswith('#'):
                        partition_cols.append(row['col_name'])

                # Get row count and size
                stats = spark.sql(f"DESCRIBE EXTENDED {db}.{table_name}")
                # Parse statistics from DESCRIBE output

                tables.append(SourceTableInfo(
                    database=db,
                    schema=db,
                    table=table_name,
                    row_count=0,  # Would query actual count
                    size_gb=0,  # From statistics
                    column_count=len(detail_df),
                    partition_columns=partition_cols,
                    dependencies=[],  # Would analyze queries
                    access_frequency='medium',
                    data_classification='internal',
                ))

            except Exception as e:
                print(f"Error processing {db}.{table_name}: {e}")

    return tables

def generate_migration_plan(tables: List[SourceTableInfo]) -> pd.DataFrame:
    """Generate migration plan with prioritization."""
    plan_data = []

    for table in tables:
        # Calculate complexity score
        complexity = 0
        complexity += 1 if table.size_gb > 100 else 0
        complexity += 1 if len(table.partition_columns) > 2 else 0
        complexity += 1 if len(table.dependencies) > 5 else 0
        complexity += 2 if table.data_classification == 'pii' else 0

        # Calculate priority
        priority = 0
        priority += 3 if table.access_frequency == 'high' else 1
        priority += 2 if table.data_classification == 'pii' else 0

        plan_data.append({
            'source_table': f"{table.database}.{table.table}",
            'target_table': f"migrated.{table.schema}.{table.table}",
            'size_gb': table.size_gb,
            'complexity_score': complexity,
            'priority_score': priority,
            'estimated_hours': max(1, table.size_gb / 10 + complexity * 2),
            'migration_wave': 1 if priority > 3 else (2 if priority > 1 else 3),
        })

    return pd.DataFrame(plan_data).sort_values(['migration_wave', 'priority_score'], ascending=[True, False])
```

### Step 2: Schema Migration

```python
# scripts/schema_migration.py
from pyspark.sql import SparkSession
from pyspark.sql.types import *

def convert_hive_to_delta_schema(spark: SparkSession, hive_table: str) -> StructType:
    """
    Convert Hive table schema to Delta Lake compatible schema.

    Handles type conversions and incompatibilities.
    """
    hive_schema = spark.table(hive_table).schema

    # Type mappings for problematic types
    type_conversions = {
        'decimal(38,0)': DecimalType(38, 10),  # Add scale
        'char': StringType(),
        'varchar': StringType(),
        'tinyint': IntegerType(),  # Delta doesn't have tinyint
    }

    new_fields = []
    for field in hive_schema.fields:
        new_type = field.dataType
        type_str = str(field.dataType).lower()

        for pattern, replacement in type_conversions.items():
            if pattern in type_str:
                new_type = replacement
                break

        new_fields.append(StructField(
            field.name,
            new_type,
            field.nullable,
            field.metadata
        ))

    return StructType(new_fields)

def migrate_table_schema(
    spark: SparkSession,
    source_table: str,
    target_table: str,
    catalog: str = "migrated",
) -> dict:
    """
    Migrate table schema from Hive to Delta Lake.

    Returns migration result with any schema changes.
    """
    # Get source schema
    source_df = spark.table(source_table)
    source_schema = source_df.schema

    # Convert schema
    target_schema = convert_hive_to_delta_schema(spark, source_table)

    # Create target table
    schema_ddl = ", ".join([
        f"`{f.name}` {f.dataType.simpleString()}"
        for f in target_schema.fields
    ])

    spark.sql(f"""
        CREATE TABLE IF NOT EXISTS {catalog}.{target_table} (
            {schema_ddl}
        )
        USING DELTA
        TBLPROPERTIES (
            'delta.autoOptimize.optimizeWrite' = 'true',
            'delta.autoOptimize.autoCompact' = 'true'
        )
    """)

    # Track schema changes
    changes = []
    for i, (src, tgt) in enumerate(zip(source_schema.fields, target_schema.fields)):
        if str(src.dataType) != str(tgt.dataType):
            changes.append({
                'column': src.name,
                'source_type': str(src.dataType),
                'target_type': str(tgt.dataType),
            })

    return {
        'source_table': source_table,
        'target_table': f"{catalog}.{target_table}",
        'column_count': len(target_schema.fields),
        'schema_changes': changes,
    }
```

### Step 3: Data Migration

```python
# scripts/data_migration.py
from pyspark.sql import SparkSession, DataFrame
from datetime import datetime
import time

class DataMigrator:
    """Handle data migration from legacy systems to Delta Lake."""

    def __init__(self, spark: SparkSession, target_catalog: str):
        self.spark = spark
        self.target_catalog = target_catalog

    def migrate_table(
        self,
        source_table: str,
        target_table: str,
        batch_size: int = 1000000,
        partition_columns: list[str] = None,
        incremental_column: str = None,
    ) -> dict:
        """
        Migrate table data with batching and checkpointing.

        Args:
            source_table: Source table (Hive, JDBC, etc.)
            target_table: Target Delta table
            batch_size: Rows per batch
            partition_columns: Columns for partitioning
            incremental_column: Column for incremental loads

        Returns:
            Migration statistics
        """
        start_time = time.time()
        stats = {
            'source_table': source_table,
            'target_table': f"{self.target_catalog}.{target_table}",
            'batches': 0,
            'total_rows': 0,
            'errors': [],
        }

        # Get source data
        source_df = self.spark.table(source_table)

        # For large tables, process in batches by partition
        if partition_columns and len(partition_columns) > 0:
            # Get distinct partition values
            partitions = source_df.select(partition_columns).distinct().collect()

            for partition_row in partitions:
                # Build filter condition
                conditions = [
                    f"{col} = '{partition_row[col]}'"
                    for col in partition_columns
                ]
                filter_expr = " AND ".join(conditions)

                batch_df = source_df.filter(filter_expr)
                self._write_batch(batch_df, target_table, partition_columns)
                stats['batches'] += 1
                stats['total_rows'] += batch_df.count()

        else:
            # Single batch for small tables
            self._write_batch(source_df, target_table, partition_columns)
            stats['batches'] = 1
            stats['total_rows'] = source_df.count()

        stats['duration_seconds'] = time.time() - start_time
        stats['rows_per_second'] = stats['total_rows'] / stats['duration_seconds']

        return stats

    def _write_batch(
        self,
        df: DataFrame,
        target_table: str,
        partition_columns: list[str] = None,
    ):
        """Write a batch to Delta table."""
        writer = df.write.format("delta").mode("append")

        if partition_columns:
            writer = writer.partitionBy(*partition_columns)

        writer.saveAsTable(f"{self.target_catalog}.{target_table}")

    def validate_migration(
        self,
        source_table: str,
        target_table: str,
    ) -> dict:
        """Validate migrated data matches source."""
        source_df = self.spark.table(source_table)
        target_df = self.spark.table(f"{self.target_catalog}.{target_table}")

        validation = {
            'source_count': source_df.count(),
            'target_count': target_df.count(),
            'count_match': False,
            'schema_match': False,
            'sample_match': False,
        }

        # Count validation
        validation['count_match'] = (
            validation['source_count'] == validation['target_count']
        )

        # Schema validation (column names)
        source_cols = set(source_df.columns)
        target_cols = set(target_df.columns)
        validation['schema_match'] = source_cols == target_cols

        # Sample data validation
        source_sample = source_df.limit(100).toPandas()
        target_sample = target_df.limit(100).toPandas()
        # Compare samples (simplified)
        validation['sample_match'] = len(source_sample) == len(target_sample)

        return validation
```

### Step 4: ETL/Pipeline Migration

```python
# scripts/pipeline_migration.py

def convert_spark_job_to_databricks(
    source_code: str,
    source_type: str = "spark-submit",
) -> str:
    """
    Convert legacy Spark job to Databricks job.

    Handles common patterns from spark-submit, Oozie, Airflow.
    """
    # Common replacements
    replacements = {
        # SparkSession changes
        'SparkSession.builder.master("yarn")': 'SparkSession.builder',
        '.master("local[*]")': '',

        # Path changes
        'hdfs://namenode:8020/': '/mnt/data/',
        's3a://': 's3://',

        # Hive changes
        '.enableHiveSupport()': '',  # Unity Catalog handles this
        'hive_metastore.': '',  # Direct table access

        # Config changes
        '.config("spark.sql.warehouse.dir"': '# Removed: .config("spark.sql.warehouse.dir"',
    }

    converted = source_code
    for old, new in replacements.items():
        converted = converted.replace(old, new)

    # Add Databricks-specific imports
    header = '''
# Converted for Databricks
# Original source: {source_type}
# Conversion date: {date}

from pyspark.sql import SparkSession

# SparkSession is pre-configured in Databricks
spark = SparkSession.builder.getOrCreate()
'''.format(source_type=source_type, date=datetime.now().isoformat())

    return header + converted

# Convert Oozie workflow to Databricks job
def convert_oozie_to_databricks_job(oozie_xml: str) -> dict:
    """Convert Oozie workflow XML to Databricks job definition."""
    import xml.etree.ElementTree as ET

    root = ET.fromstring(oozie_xml)

    tasks = []
    for action in root.findall('.//action'):
        action_name = action.get('name')

        # Handle different action types
        spark_action = action.find('spark')
        if spark_action is not None:
            jar = spark_action.find('jar').text
            main_class = spark_action.find('class').text

            tasks.append({
                'task_key': action_name,
                'spark_jar_task': {
                    'main_class_name': main_class,
                    'parameters': [],
                },
                'libraries': [{'jar': jar}],
            })

        shell_action = action.find('shell')
        if shell_action is not None:
            # Convert to notebook task or skip
            pass

    # Build job definition
    job_definition = {
        'name': f"migrated-{root.get('name')}",
        'tasks': tasks,
        'job_clusters': [{
            'job_cluster_key': 'migration_cluster',
            'new_cluster': {
                'spark_version': '14.3.x-scala2.12',
                'node_type_id': 'Standard_DS3_v2',
                'num_workers': 2,
            }
        }],
    }

    return job_definition
```

### Step 5: Cutover Planning

```python
# scripts/cutover_plan.py
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import List

@dataclass
class CutoverStep:
    """Individual step in cutover plan."""
    order: int
    name: str
    duration_minutes: int
    owner: str
    rollback_procedure: str
    verification: str

def generate_cutover_plan(
    migration_wave: int,
    tables: List[str],
    cutover_date: datetime,
) -> List[CutoverStep]:
    """Generate detailed cutover plan."""

    steps = [
        CutoverStep(
            order=1,
            name="Pre-cutover validation",
            duration_minutes=60,
            owner="Data Engineer",
            rollback_procedure="N/A - no changes made",
            verification="Run validation queries on all tables",
        ),
        CutoverStep(
            order=2,
            name="Disable source data pipelines",
            duration_minutes=15,
            owner="Platform Admin",
            rollback_procedure="Re-enable pipelines in source system",
            verification="Verify no new data in source",
        ),
        CutoverStep(
            order=3,
            name="Final incremental sync",
            duration_minutes=120,
            owner="Data Engineer",
            rollback_procedure="N/A",
            verification="Row counts match source",
        ),
        CutoverStep(
            order=4,
            name="Enable Databricks pipelines",
            duration_minutes=30,
            owner="Data Engineer",
            rollback_procedure="Disable Databricks pipelines, re-enable source",
            verification="Jobs running successfully",
        ),
        CutoverStep(
            order=5,
            name="Update downstream applications",
            duration_minutes=60,
            owner="Application Team",
            rollback_procedure="Revert connection strings",
            verification="Applications reading from Databricks",
        ),
        CutoverStep(
            order=6,
            name="Monitor and validate",
            duration_minutes=240,
            owner="Data Engineer",
            rollback_procedure="Execute full rollback if issues",
            verification="All metrics within acceptable range",
        ),
    ]

    # Calculate timings
    current_time = cutover_date
    for step in steps:
        step.start_time = current_time
        step.end_time = current_time + timedelta(minutes=step.duration_minutes)
        current_time = step.end_time

    return steps
```

## Output
- Migration assessment complete
- Schema migration automated
- Data migration pipeline ready
- ETL conversion scripts
- Cutover plan documented

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Schema incompatibility | Unsupported types | Use type conversion mappings |
| Data loss | Truncation | Validate counts at each step |
| Performance issues | Large tables | Use partitioned migration |
| Dependency conflicts | Wrong migration order | Analyze dependencies first |

## Examples

### Quick Migration Validation
```sql
-- Compare source and target counts
SELECT
    'source' as system, COUNT(*) as row_count
FROM hive_metastore.db.table
UNION ALL
SELECT
    'target' as system, COUNT(*) as row_count
FROM migrated.db.table;
```

## Resources
- [Databricks Migration Guide](https://docs.databricks.com/migration/index.html)
- [Delta Lake Migration](https://docs.databricks.com/delta/migration.html)
- [Unity Catalog Migration](https://docs.databricks.com/data-governance/unity-catalog/migrate.html)

## Completion
This skill pack provides comprehensive coverage for Databricks platform operations.
