---
name: databricks-core-workflow-a
description: |
  Execute Databricks primary workflow: Delta Lake ETL pipelines.
  Use when building data ingestion pipelines, implementing medallion architecture,
  or creating Delta Lake transformations.
  Trigger with phrases like "databricks ETL", "delta lake pipeline",
  "medallion architecture", "databricks data pipeline", "bronze silver gold".
allowed-tools: Read, Write, Edit, Bash(databricks:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Databricks Core Workflow A: Delta Lake ETL

## Overview
Build production Delta Lake ETL pipelines using medallion architecture.

## Prerequisites
- Completed `databricks-install-auth` setup
- Understanding of Delta Lake concepts
- Unity Catalog configured (recommended)

## Medallion Architecture

```
Raw Sources → Bronze (Raw) → Silver (Cleaned) → Gold (Aggregated)
                  ↓              ↓                  ↓
           Landing Zone    Business Logic    Analytics Ready
```

## Instructions

### Step 1: Bronze Layer - Raw Ingestion
```python
# src/pipelines/bronze.py
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.functions import current_timestamp, input_file_name, lit
from delta.tables import DeltaTable

def ingest_to_bronze(
    spark: SparkSession,
    source_path: str,
    target_table: str,
    source_format: str = "json",
    schema: str = None,
) -> DataFrame:
    """
    Ingest raw data to Bronze layer with metadata.

    Args:
        spark: SparkSession
        source_path: Path to source data
        target_table: Unity Catalog table name (catalog.schema.table)
        source_format: Source file format (json, csv, parquet)
        schema: Optional schema string
    """
    # Read raw data
    reader = spark.read.format(source_format)
    if schema:
        reader = reader.schema(schema)

    df = reader.load(source_path)

    # Add ingestion metadata
    df_with_metadata = (
        df
        .withColumn("_ingested_at", current_timestamp())
        .withColumn("_source_file", input_file_name())
        .withColumn("_source_format", lit(source_format))
    )

    # Write to Delta with merge for idempotency
    df_with_metadata.write \
        .format("delta") \
        .mode("append") \
        .option("mergeSchema", "true") \
        .saveAsTable(target_table)

    return df_with_metadata

# Auto Loader for streaming ingestion
def stream_to_bronze(
    spark: SparkSession,
    source_path: str,
    target_table: str,
    checkpoint_path: str,
    schema_location: str,
) -> None:
    """Stream data to Bronze using Auto Loader."""
    (
        spark.readStream
        .format("cloudFiles")
        .option("cloudFiles.format", "json")
        .option("cloudFiles.schemaLocation", schema_location)
        .option("cloudFiles.inferColumnTypes", "true")
        .load(source_path)
        .withColumn("_ingested_at", current_timestamp())
        .writeStream
        .format("delta")
        .option("checkpointLocation", checkpoint_path)
        .option("mergeSchema", "true")
        .trigger(availableNow=True)
        .toTable(target_table)
    )
```

### Step 2: Silver Layer - Data Cleansing
```python
# src/pipelines/silver.py
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.functions import (
    col, when, trim, lower, to_timestamp,
    regexp_replace, sha2, concat_ws
)
from delta.tables import DeltaTable

def transform_to_silver(
    spark: SparkSession,
    bronze_table: str,
    silver_table: str,
    primary_keys: list[str],
    watermark_column: str = "_ingested_at",
) -> DataFrame:
    """
    Transform Bronze to Silver with cleansing and deduplication.

    Args:
        spark: SparkSession
        bronze_table: Source Bronze table
        silver_table: Target Silver table
        primary_keys: Columns for deduplication/merge
        watermark_column: Column for incremental processing
    """
    # Read Bronze (incremental using CDF)
    bronze_df = spark.readStream \
        .format("delta") \
        .option("readChangeFeed", "true") \
        .table(bronze_table)

    # Apply transformations
    silver_df = (
        bronze_df
        # Standardize strings
        .withColumn("name", trim(col("name")))
        .withColumn("email", lower(trim(col("email"))))
        # Parse timestamps
        .withColumn("created_at", to_timestamp(col("created_at")))
        # Remove PII for analytics
        .withColumn("email_hash", sha2(col("email"), 256))
        # Data quality filters
        .filter(col("email").isNotNull())
        .filter(col("created_at").isNotNull())
        # Generate surrogate key
        .withColumn(
            "_row_key",
            sha2(concat_ws("||", *[col(k) for k in primary_keys]), 256)
        )
    )

    # Merge into Silver (upsert pattern)
    if DeltaTable.isDeltaTable(spark, silver_table):
        delta_table = DeltaTable.forName(spark, silver_table)
        merge_condition = " AND ".join(
            [f"target.{k} = source.{k}" for k in primary_keys]
        )
        (
            delta_table.alias("target")
            .merge(silver_df.alias("source"), merge_condition)
            .whenMatchedUpdateAll()
            .whenNotMatchedInsertAll()
            .execute()
        )
    else:
        silver_df.write \
            .format("delta") \
            .mode("overwrite") \
            .saveAsTable(silver_table)

    return silver_df
```

### Step 3: Gold Layer - Business Aggregations
```python
# src/pipelines/gold.py
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.functions import (
    col, count, sum, avg, max, min,
    date_trunc, window, current_timestamp
)

def aggregate_to_gold(
    spark: SparkSession,
    silver_table: str,
    gold_table: str,
    group_by_columns: list[str],
    aggregations: dict[str, str],
    time_grain: str = "day",
) -> DataFrame:
    """
    Aggregate Silver to Gold for analytics.

    Args:
        spark: SparkSession
        silver_table: Source Silver table
        gold_table: Target Gold table
        group_by_columns: Columns to group by
        aggregations: Dict of {output_col: "agg_func(source_col)"}
        time_grain: Time aggregation grain (hour, day, week, month)
    """
    silver_df = spark.table(silver_table)

    # Build aggregation expressions
    agg_exprs = []
    for output_col, expr in aggregations.items():
        agg_exprs.append(f"{expr} as {output_col}")

    # Apply aggregations
    gold_df = (
        silver_df
        .withColumn("time_period", date_trunc(time_grain, col("created_at")))
        .groupBy(*group_by_columns, "time_period")
        .agg(*[eval(e) for e in agg_exprs])
        .withColumn("_aggregated_at", current_timestamp())
    )

    # Write to Gold (overwrite partition)
    gold_df.write \
        .format("delta") \
        .mode("overwrite") \
        .option("replaceWhere", f"time_period >= '{get_min_date()}'") \
        .saveAsTable(gold_table)

    return gold_df

# Example usage
gold_df = aggregate_to_gold(
    spark=spark,
    silver_table="catalog.silver.events",
    gold_table="catalog.gold.daily_metrics",
    group_by_columns=["region", "product_category"],
    aggregations={
        "total_orders": "count(*)",
        "total_revenue": "sum(amount)",
        "avg_order_value": "avg(amount)",
        "unique_customers": "count(distinct customer_id)",
    },
    time_grain="day"
)
```

### Step 4: Delta Live Tables (DLT) Pipeline
```python
# pipelines/dlt_pipeline.py
import dlt
from pyspark.sql.functions import *

@dlt.table(
    name="bronze_events",
    comment="Raw events from source",
    table_properties={"quality": "bronze"}
)
def bronze_events():
    return (
        spark.readStream
        .format("cloudFiles")
        .option("cloudFiles.format", "json")
        .load("/mnt/landing/events/")
        .withColumn("_ingested_at", current_timestamp())
    )

@dlt.table(
    name="silver_events",
    comment="Cleansed and validated events",
    table_properties={"quality": "silver"}
)
@dlt.expect_or_drop("valid_email", "email IS NOT NULL")
@dlt.expect_or_drop("valid_amount", "amount > 0")
def silver_events():
    return (
        dlt.read_stream("bronze_events")
        .withColumn("email", lower(trim(col("email"))))
        .withColumn("event_time", to_timestamp(col("event_time")))
    )

@dlt.table(
    name="gold_daily_summary",
    comment="Daily aggregated metrics",
    table_properties={"quality": "gold"}
)
def gold_daily_summary():
    return (
        dlt.read("silver_events")
        .groupBy(date_trunc("day", col("event_time")).alias("date"))
        .agg(
            count("*").alias("total_events"),
            sum("amount").alias("total_revenue"),
            countDistinct("customer_id").alias("unique_customers")
        )
    )
```

## Output
- Bronze layer with raw data and metadata
- Silver layer with cleansed, deduplicated data
- Gold layer with business aggregations
- Delta Lake tables with ACID transactions

## Error Handling
| Error | Cause | Solution |
|-------|-------|----------|
| Schema mismatch | Source schema changed | Use `mergeSchema` option |
| Duplicate records | Missing deduplication | Add merge logic with primary keys |
| Null values | Data quality issues | Add expectations/filters |
| Memory errors | Large aggregations | Increase cluster size or partition data |

## Examples

### Complete Pipeline Orchestration
```python
# Run full medallion pipeline
from src.pipelines import bronze, silver, gold

# Bronze ingestion
bronze.ingest_to_bronze(
    spark, "/mnt/landing/orders/", "catalog.bronze.orders"
)

# Silver transformation
silver.transform_to_silver(
    spark, "catalog.bronze.orders", "catalog.silver.orders",
    primary_keys=["order_id"]
)

# Gold aggregation
gold.aggregate_to_gold(
    spark, "catalog.silver.orders", "catalog.gold.order_metrics",
    group_by_columns=["region"], time_grain="day"
)
```

## Resources
- [Delta Lake Guide](https://docs.databricks.com/delta/index.html)
- [Medallion Architecture](https://www.databricks.com/glossary/medallion-architecture)
- [Delta Live Tables](https://docs.databricks.com/delta-live-tables/index.html)
- [Auto Loader](https://docs.databricks.com/ingestion/auto-loader/index.html)

## Next Steps
For ML workflows, see `databricks-core-workflow-b`.
