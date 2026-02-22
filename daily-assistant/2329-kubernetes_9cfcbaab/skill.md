---
skill_name: kubernetes
skill_category: devops
description: Kubernetes deployment patterns, resource management, and cluster best practices
allowed_tools: [Read, Edit, Glob, Grep, Write]
token_estimate: 850
version: 1.0
last_updated: 2026-01-13
owner: DevOps Team
status: active

tags: [kubernetes, k8s, containers, anti-pattern, best-practice, deployment]
related_skills: [ci-cd-patterns, docker-patterns, terraform-patterns]

trigger_files: ["**/k8s/*.yaml", "**/k8s/*.yml", "**/kubernetes/*.yaml", "**/kubernetes/*.yml", "**/deployment.yaml", "**/deployment.yml", "Helm*", "**/charts/**/*.yaml"]
trigger_keywords: [kubernetes, k8s, kubectl, deployment, pod, service, ingress, helm, container orchestration]

quality_keywords: [anti-pattern, pattern, validation, best-practice, security, resource-limits]
---

# Kubernetes Patterns

Best practices and anti-patterns for Kubernetes deployments, resource configuration, and cluster management.

## Purpose

- Ensure reliable, scalable Kubernetes deployments
- Prevent resource exhaustion and scheduling failures
- Establish security and observability patterns

---

## Core Patterns

### Pattern 1: Resource Requests and Limits

**When to use:** Every production container deployment.

**Implementation:**
```yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: api-service
spec:
  template:
    spec:
      containers:
        - name: api
          image: api:v1.2.3
          resources:
            requests:
              memory: "256Mi"
              cpu: "250m"
            limits:
              memory: "512Mi"
              cpu: "500m"
```

**Benefits:**
- Predictable scheduling
- Protection against resource exhaustion
- Fair resource sharing across pods

### Pattern 2: Health Probes

**When to use:** All production deployments.

**Implementation:**
```yaml
containers:
  - name: api
    livenessProbe:
      httpGet:
        path: /health
        port: 8080
      initialDelaySeconds: 30
      periodSeconds: 10
      failureThreshold: 3
    readinessProbe:
      httpGet:
        path: /ready
        port: 8080
      initialDelaySeconds: 5
      periodSeconds: 5
      successThreshold: 1
```

**Benefits:**
- Auto-restart unhealthy pods
- Traffic only to ready pods
- Graceful startup handling

### Pattern 3: Pod Disruption Budgets

**When to use:** Production services requiring high availability.

**Implementation:**
```yaml
apiVersion: policy/v1
kind: PodDisruptionBudget
metadata:
  name: api-pdb
spec:
  minAvailable: 2  # or maxUnavailable: 1
  selector:
    matchLabels:
      app: api
```

**Benefits:**
- Controlled disruption during upgrades
- Maintains minimum availability
- Safe cluster maintenance

---

## Anti-Patterns

### Anti-Pattern 1: No Resource Limits

| Aspect | Description |
|--------|-------------|
| **WHY** | Unlimited resources cause node exhaustion, OOM kills, and unpredictable scheduling. One pod can starve others. |
| **DETECTION** | Missing `resources.limits` in container specs. Pods with no resource constraints. |
| **FIX** | Always set requests (guaranteed) and limits (maximum). Start with requests = typical usage. |

**Bad Example:**
```yaml
containers:
  - name: api
    image: api:latest
    # No resource constraints!
```

**Good Example:**
```yaml
containers:
  - name: api
    image: api:v1.2.3
    resources:
      requests:
        memory: "256Mi"
        cpu: "250m"
      limits:
        memory: "512Mi"
        cpu: "500m"
```

### Anti-Pattern 2: Using :latest Tag

| Aspect | Description |
|--------|-------------|
| **WHY** | Non-reproducible deployments. Rollbacks impossible. Different nodes may pull different versions. |
| **DETECTION** | Image tags ending in `:latest` or no tag specified. |
| **FIX** | Use specific semantic version tags. Pin to exact digest for critical workloads. |

**Bad Example:**
```yaml
containers:
  - name: api
    image: api:latest
    imagePullPolicy: Always
```

**Good Example:**
```yaml
containers:
  - name: api
    image: api:v1.2.3
    # Or for maximum reproducibility:
    # image: api@sha256:abc123...
```

### Anti-Pattern 3: Missing Health Probes

| Aspect | Description |
|--------|-------------|
| **WHY** | Kubernetes can't detect unhealthy pods. Traffic goes to broken pods. No auto-recovery. |
| **DETECTION** | Deployments without `livenessProbe` or `readinessProbe`. |
| **FIX** | Add both probes. Liveness for restart, readiness for traffic routing. |

**Bad Example:**
```yaml
containers:
  - name: api
    image: api:v1.2.3
    ports:
      - containerPort: 8080
    # No health probes!
```

**Good Example:**
```yaml
containers:
  - name: api
    image: api:v1.2.3
    ports:
      - containerPort: 8080
    livenessProbe:
      httpGet:
        path: /health
        port: 8080
      initialDelaySeconds: 30
      periodSeconds: 10
    readinessProbe:
      httpGet:
        path: /ready
        port: 8080
      initialDelaySeconds: 5
      periodSeconds: 5
```

### Anti-Pattern 4: Running as Root

| Aspect | Description |
|--------|-------------|
| **WHY** | Container escape vulnerabilities are exploitable as root. Violates least privilege. |
| **DETECTION** | No `securityContext` or `runAsNonRoot: false`. Missing `runAsUser`. |
| **FIX** | Set `runAsNonRoot: true`. Specify non-root user. Drop all capabilities. |

**Bad Example:**
```yaml
containers:
  - name: api
    image: api:v1.2.3
    # Runs as root by default!
```

**Good Example:**
```yaml
containers:
  - name: api
    image: api:v1.2.3
    securityContext:
      runAsNonRoot: true
      runAsUser: 1000
      readOnlyRootFilesystem: true
      capabilities:
        drop:
          - ALL
```

### Anti-Pattern 5: Hardcoded ConfigMaps/Secrets

| Aspect | Description |
|--------|-------------|
| **WHY** | Configuration changes require deployment updates. Secrets exposed in manifests. |
| **DETECTION** | Environment values hardcoded in deployment. Secrets inline in YAML. |
| **FIX** | Use ConfigMaps for config, Secrets for sensitive data. Reference via envFrom. |

**Bad Example:**
```yaml
containers:
  - name: api
    env:
      - name: DATABASE_URL
        value: "postgres://user:password@db:5432/app"
      - name: API_KEY
        value: "sk-secret-key-here"
```

**Good Example:**
```yaml
containers:
  - name: api
    envFrom:
      - configMapRef:
          name: api-config
      - secretRef:
          name: api-secrets
```

### Anti-Pattern 6: Single Replica in Production

| Aspect | Description |
|--------|-------------|
| **WHY** | No fault tolerance. Pod restart = downtime. No rolling update capability. |
| **DETECTION** | `replicas: 1` in production deployments. No PDB defined. |
| **FIX** | Minimum 2 replicas for production. Add PodDisruptionBudget. Use anti-affinity. |

**Bad Example:**
```yaml
spec:
  replicas: 1  # Single point of failure!
```

**Good Example:**
```yaml
spec:
  replicas: 3
  template:
    spec:
      affinity:
        podAntiAffinity:
          preferredDuringSchedulingIgnoredDuringExecution:
            - weight: 100
              podAffinityTerm:
                topologyKey: kubernetes.io/hostname
                labelSelector:
                  matchLabels:
                    app: api
```

---

## Validation Checklist

### Pre-Implementation
- [ ] Review cluster resource capacity
- [ ] Check existing deployments for patterns
- [ ] Understand application resource needs

### Implementation
- [ ] Resource requests and limits set
- [ ] Specific image tags (no :latest)
- [ ] Liveness and readiness probes configured
- [ ] Security context with non-root user
- [ ] ConfigMaps/Secrets for configuration
- [ ] Multiple replicas with anti-affinity
- [ ] PodDisruptionBudget defined

### Post-Implementation
- [ ] Pods schedule successfully
- [ ] Health checks pass consistently
- [ ] Rolling updates work without downtime
- [ ] `kubectl get events` shows no warnings

---

## Quick Reference

| Task | Pattern |
|------|---------|
| Prevent resource exhaustion | Set requests and limits |
| Auto-restart broken pods | Liveness probe |
| Route traffic to ready pods | Readiness probe |
| Reproducible deployments | Specific version tags |
| Secure containers | Non-root, drop capabilities |
| High availability | Multiple replicas + PDB |
| External configuration | ConfigMaps + Secrets |
