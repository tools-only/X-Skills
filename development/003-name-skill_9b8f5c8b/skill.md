---
name: guidewire-reference-architecture
description: |
  Enterprise reference architecture for Guidewire InsuranceSuite Cloud deployments.
  Use when designing system architecture, planning integrations,
  or understanding Guidewire cloud infrastructure patterns.
  Trigger with phrases like "guidewire architecture", "system design",
  "integration architecture", "enterprise guidewire", "reference architecture".
allowed-tools: Read, Write, Edit, Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Guidewire Reference Architecture

## Overview

Enterprise reference architecture patterns for Guidewire InsuranceSuite Cloud deployments, including integration patterns, data flows, and scalability considerations.

## Prerequisites

- Understanding of enterprise architecture concepts
- Familiarity with Guidewire InsuranceSuite components
- Knowledge of cloud architecture patterns

## High-Level Architecture

```
                              ┌─────────────────────────────────────────────┐
                              │         External Users & Channels           │
                              │  (Agents, Customers, Partners, Regulators)  │
                              └─────────────────┬───────────────────────────┘
                                                │
                              ┌─────────────────▼───────────────────────────┐
                              │           Digital Experience Layer          │
                              │  ┌─────────┐ ┌─────────┐ ┌─────────────┐   │
                              │  │  Agent  │ │Customer │ │  Partner    │   │
                              │  │ Portal  │ │ Portal  │ │   Portal    │   │
                              │  │ (Jutro) │ │ (Jutro) │ │   (API)     │   │
                              │  └────┬────┘ └────┬────┘ └──────┬──────┘   │
                              └───────┼──────────┼─────────────┼───────────┘
                                      │          │             │
                              ┌───────▼──────────▼─────────────▼───────────┐
                              │              API Gateway                    │
                              │      (Authentication, Rate Limiting)        │
                              └────────────────────┬────────────────────────┘
                                                   │
  ┌────────────────────────────────────────────────┼────────────────────────────────────────────────┐
  │                                    Guidewire Cloud Platform                                     │
  │                                                │                                                │
  │  ┌─────────────────┐  ┌─────────────────┐  ┌──┴──────────────┐  ┌─────────────────┐            │
  │  │   PolicyCenter  │  │   ClaimCenter   │  │  BillingCenter  │  │  Contact        │            │
  │  │                 │  │                 │  │                 │  │  Manager        │            │
  │  │ • Submissions   │  │ • FNOL          │  │ • Invoicing     │  │                 │            │
  │  │ • Quoting       │  │ • Investigation │  │ • Payments      │  │ • Contacts      │            │
  │  │ • Binding       │  │ • Settlement    │  │ • Collections   │  │ • Addresses     │            │
  │  │ • Issuance      │  │ • Payments      │  │ • Commissions   │  │ • Roles         │            │
  │  │ • Endorsements  │  │ • Litigation    │  │                 │  │                 │            │
  │  │ • Renewals      │  │                 │  │                 │  │                 │            │
  │  └────────┬────────┘  └────────┬────────┘  └────────┬────────┘  └────────┬────────┘            │
  │           │                    │                    │                    │                     │
  │           └────────────────────┴────────────────────┴────────────────────┘                     │
  │                                           │                                                    │
  │                              ┌────────────▼────────────┐                                       │
  │                              │     Shared Services     │                                       │
  │                              │ • Document Management   │                                       │
  │                              │ • Workflow Engine       │                                       │
  │                              │ • Rules Engine          │                                       │
  │                              │ • Reporting             │                                       │
  │                              └────────────┬────────────┘                                       │
  │                                           │                                                    │
  │  ┌────────────────────────────────────────┴────────────────────────────────────────┐          │
  │  │                         Integration Layer                                        │          │
  │  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐        │          │
  │  │  │   Cloud API  │  │  App Events  │  │  Integration │  │    Batch     │        │          │
  │  │  │   (REST)     │  │   (Kafka)    │  │   Gateway    │  │   Services   │        │          │
  │  │  └──────────────┘  └──────────────┘  └──────────────┘  └──────────────┘        │          │
  │  └─────────────────────────────────────────────────────────────────────────────────┘          │
  │                                                                                                │
  └────────────────────────────────────────────────────────────────────────────────────────────────┘
                                                │
  ┌─────────────────────────────────────────────┼─────────────────────────────────────────────────┐
  │                              Enterprise Integration Layer                                      │
  │                                             │                                                  │
  │  ┌───────────┐  ┌───────────┐  ┌───────────┴───────────┐  ┌───────────┐  ┌───────────┐       │
  │  │    CRM    │  │  ERP/GL   │  │   Rating Engines      │  │  Document │  │  Legacy   │       │
  │  │ (Salesforce) │ (SAP/Oracle) │  (External/Internal)   │  │   Mgmt    │  │  Systems  │       │
  │  └───────────┘  └───────────┘  └───────────────────────┘  └───────────┘  └───────────┘       │
  │                                                                                                │
  └────────────────────────────────────────────────────────────────────────────────────────────────┘
```

## Integration Patterns

### Pattern 1: Synchronous API Integration

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Client    │────▶│  Cloud API  │────▶│  External   │
│ Application │     │             │     │   Service   │
│             │◀────│             │◀────│             │
└─────────────┘     └─────────────┘     └─────────────┘
       │                   │                   │
       │   HTTP Request    │   HTTP Request    │
       │──────────────────▶│──────────────────▶│
       │                   │                   │
       │   HTTP Response   │   HTTP Response   │
       │◀──────────────────│◀──────────────────│
```

**Use Cases:**
- Real-time policy quoting
- Address validation
- Credit scoring
- Real-time fraud detection

```typescript
// Synchronous integration example
async function getRealTimeQuote(submissionId: string): Promise<Quote> {
  // Call external rating engine
  const ratingResponse = await ratingService.calculatePremium({
    submissionId,
    effectiveDate: submission.effectiveDate,
    coverages: submission.coverages
  });

  // Update Guidewire with results
  return await guidewireClient.updateQuote(submissionId, {
    premium: ratingResponse.premium,
    taxes: ratingResponse.taxes,
    fees: ratingResponse.fees
  });
}
```

### Pattern 2: Asynchronous Event-Driven

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│ InsuranceSuite │──▶│  App Events │──▶│    Kafka    │──▶│  Consumer   │
│             │     │   Service   │     │   Topic     │     │  Service    │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
       │                   │                   │                   │
       │  Business Event   │   Publish Event   │   Consume Event   │
       │──────────────────▶│──────────────────▶│──────────────────▶│
       │                   │                   │                   │
       │                   │                   │   Process Async   │
       │                   │                   │                   │
```

**Use Cases:**
- Policy issued notifications
- Claims status updates
- Billing events
- Data warehouse synchronization

### Pattern 3: Batch Integration

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│    SFTP     │────▶│    Batch    │────▶│  Transform  │────▶│ InsuranceSuite │
│   Server    │     │   Pickup    │     │   & Load    │     │             │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
       │                   │                   │                   │
       │   Drop File       │   Schedule Job    │   Process Data    │
       │──────────────────▶│──────────────────▶│──────────────────▶│
```

**Use Cases:**
- Nightly policy updates
- Bulk claims import
- Premium bordereaux
- Regulatory reporting

## Data Flow Architecture

### Policy Lifecycle Data Flow

```yaml
# Policy data flow through system components
policy_flow:
  1_submission:
    source: Agent Portal / Direct Customer
    target: PolicyCenter
    data:
      - applicant_info
      - coverage_requests
      - risk_data

  2_underwriting:
    source: PolicyCenter
    integrations:
      - external_rating_engine
      - credit_bureau
      - mvr_service
      - loss_history
    data:
      - risk_scores
      - premium_calculations
      - underwriting_decision

  3_binding:
    source: PolicyCenter
    target: BillingCenter
    data:
      - policy_terms
      - premium_schedule
      - payment_plan

  4_document_generation:
    source: PolicyCenter
    target: Document Service
    data:
      - policy_documents
      - dec_pages
      - endorsements

  5_distribution:
    source: Document Service
    targets:
      - customer_email
      - agent_portal
      - document_archive
```

### Claims Data Flow

```yaml
claims_flow:
  1_fnol:
    source: Customer Portal / Call Center
    target: ClaimCenter
    data:
      - loss_details
      - policy_verification
      - initial_reserve

  2_investigation:
    source: ClaimCenter
    integrations:
      - fraud_detection
      - medical_records
      - police_reports
    data:
      - investigation_results
      - liability_assessment

  3_settlement:
    source: ClaimCenter
    target: BillingCenter
    data:
      - payment_authorization
      - vendor_payments
      - subrogation_recovery

  4_reporting:
    source: ClaimCenter
    target: Data Warehouse
    data:
      - claim_metrics
      - loss_ratios
      - regulatory_reports
```

## Security Architecture

```
┌──────────────────────────────────────────────────────────────────────┐
│                        Security Perimeter                            │
│  ┌────────────────────────────────────────────────────────────────┐  │
│  │                         WAF / DDoS                              │  │
│  └────────────────────────────────────────────────────────────────┘  │
│                                   │                                  │
│  ┌────────────────────────────────▼───────────────────────────────┐  │
│  │                      API Gateway                                │  │
│  │  • Rate Limiting    • OAuth2/JWT    • Request Validation       │  │
│  └────────────────────────────────────────────────────────────────┘  │
│                                   │                                  │
│  ┌─────────────────────────────────────────────────────────────────┐ │
│  │                    Guidewire Hub (IdP)                          │ │
│  │  • Identity Federation   • MFA    • Role-Based Access          │ │
│  └─────────────────────────────────────────────────────────────────┘ │
│                                   │                                  │
│  ┌─────────────────────────────────────────────────────────────────┐ │
│  │              InsuranceSuite Applications                        │ │
│  │  • Data Encryption (AES-256)  • PII Masking   • Audit Logging  │ │
│  └─────────────────────────────────────────────────────────────────┘ │
│                                   │                                  │
│  ┌─────────────────────────────────────────────────────────────────┐ │
│  │                     Database Layer                              │ │
│  │  • TDE (Transparent Data Encryption)  • Backup Encryption      │ │
│  └─────────────────────────────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────────────────────┘
```

## Scalability Patterns

### Horizontal Scaling

```yaml
# Auto-scaling configuration
scaling:
  application_tier:
    min_instances: 2
    max_instances: 10
    target_cpu: 70%
    scale_up_cooldown: 300s
    scale_down_cooldown: 600s

  batch_processing:
    strategy: parallel_workers
    worker_count: 4
    queue_threshold: 1000

  database:
    read_replicas: 2
    connection_pool:
      min: 10
      max: 50
```

### Caching Strategy

```
┌─────────────────────────────────────────────────────────────────┐
│                         Cache Tiers                              │
│                                                                  │
│  ┌─────────────┐     ┌─────────────┐     ┌─────────────┐       │
│  │  L1 Cache   │────▶│  L2 Cache   │────▶│  Database   │       │
│  │  (In-Memory)│     │   (Redis)   │     │             │       │
│  │  TTL: 60s   │     │  TTL: 300s  │     │             │       │
│  └─────────────┘     └─────────────┘     └─────────────┘       │
│        │                   │                   │                │
│  • Product Models    • API Responses    • Transactional        │
│  • Typelists         • Session Data     • Master Data          │
│  • User Preferences  • Rate Tables      • Historical           │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Deployment Topology

### Multi-Region Architecture

```
                    ┌─────────────────────────────────────┐
                    │          Global Load Balancer       │
                    │              (CDN/DNS)               │
                    └─────────────────┬───────────────────┘
                                      │
          ┌───────────────────────────┼───────────────────────────┐
          │                           │                           │
          ▼                           ▼                           ▼
┌─────────────────┐       ┌─────────────────┐       ┌─────────────────┐
│   US East       │       │   US West       │       │   EU West       │
│   Region        │       │   Region        │       │   Region        │
│                 │       │                 │       │                 │
│ ┌─────────────┐ │       │ ┌─────────────┐ │       │ ┌─────────────┐ │
│ │ App Cluster │ │       │ │ App Cluster │ │       │ │ App Cluster │ │
│ └─────────────┘ │       │ └─────────────┘ │       │ └─────────────┘ │
│ ┌─────────────┐ │       │ ┌─────────────┐ │       │ ┌─────────────┐ │
│ │  Database   │◀───────▶│ │  Database   │◀───────▶│ │  Database   │ │
│ │  (Primary)  │ │       │ │  (Replica)  │ │       │ │  (Replica)  │ │
│ └─────────────┘ │       │ └─────────────┘ │       │ └─────────────┘ │
└─────────────────┘       └─────────────────┘       └─────────────────┘
```

## Environment Strategy

| Environment | Purpose | Data | Integrations |
|-------------|---------|------|--------------|
| Development | Feature development | Synthetic | Mocked |
| Test/QA | Integration testing | Anonymized | Sandbox endpoints |
| UAT | User acceptance | Anonymized | Sandbox endpoints |
| Staging | Pre-production | Prod subset | Production endpoints |
| Production | Live system | Production | Production endpoints |

## Technology Stack Summary

| Layer | Technology | Purpose |
|-------|------------|---------|
| Frontend | Jutro Digital Platform | React-based portals |
| API Gateway | Guidewire Hub | Auth, routing |
| Core Apps | InsuranceSuite | PC, CC, BC |
| Integration | Integration Gateway | Apache Camel |
| Messaging | Apache Kafka | Event streaming |
| Database | PostgreSQL/Oracle | Relational data |
| Cache | Redis | Session, API cache |
| Search | Elasticsearch | Full-text search |
| Monitoring | Datadog/Splunk | Observability |

## Resources

- [Guidewire Cloud Architecture](https://docs.guidewire.com/cloud/)
- [Integration Framework](https://www.guidewire.com/developers/developer-tools-and-guides/integration-framework)
- [Security Documentation](https://docs.guidewire.com/security/)

## Next Steps

For multi-environment setup, see `guidewire-multi-env-setup`.
