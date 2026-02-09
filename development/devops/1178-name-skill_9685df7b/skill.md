---
name: twinmind-deploy-integration
description: |
  Deploy TwinMind integrations to production environments.
  Use when deploying to cloud platforms, configuring production infrastructure,
  or setting up deployment automation.
  Trigger with phrases like "deploy twinmind", "twinmind production deploy",
  "twinmind cloud deployment", "twinmind infrastructure".
allowed-tools: Read, Write, Edit, Bash(npm:*), Bash(docker:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Deploy Integration

## Overview
Deploy TwinMind integrations to production cloud environments.

## Prerequisites
- Completed `twinmind-ci-integration` setup
- Cloud provider account (AWS, GCP, or Azure)
- Docker installed (for containerized deployments)
- Terraform or Pulumi (for IaC)

## Instructions

### Step 1: Docker Configuration

```dockerfile
# Dockerfile
FROM node:20-alpine AS builder

WORKDIR /app

# Install dependencies
COPY package*.json ./
RUN npm ci --only=production

# Copy source
COPY . .

# Build TypeScript
RUN npm run build

# Production image
FROM node:20-alpine AS runner

WORKDIR /app

# Security: run as non-root
RUN addgroup -g 1001 -S nodejs && \
    adduser -S twinmind -u 1001 -G nodejs

# Copy built assets
COPY --from=builder /app/dist ./dist
COPY --from=builder /app/node_modules ./node_modules
COPY --from=builder /app/package.json ./

USER twinmind

ENV NODE_ENV=production
ENV PORT=8080

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
  CMD wget --no-verbose --tries=1 --spider http://localhost:8080/health || exit 1

CMD ["node", "dist/index.js"]
```

```yaml
# docker-compose.yml
version: '3.8'

services:
  twinmind-service:
    build: .
    ports:
      - "8080:8080"
    environment:
      - NODE_ENV=production
      - TWINMIND_API_KEY=${TWINMIND_API_KEY}
      - TWINMIND_WEBHOOK_SECRET=${TWINMIND_WEBHOOK_SECRET}
      - DATABASE_URL=${DATABASE_URL}
      - REDIS_URL=${REDIS_URL}
    healthcheck:
      test: ["CMD", "wget", "-q", "--spider", "http://localhost:8080/health"]
      interval: 30s
      timeout: 10s
      retries: 3
    deploy:
      replicas: 2
      resources:
        limits:
          cpus: '1'
          memory: 512M
    logging:
      driver: "json-file"
      options:
        max-size: "10m"
        max-file: "3"

  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"
    volumes:
      - redis-data:/data

volumes:
  redis-data:
```

### Step 2: AWS Deployment (ECS/Fargate)

```hcl
# terraform/aws/main.tf
terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

provider "aws" {
  region = var.aws_region
}

# Secrets Manager for API keys
resource "aws_secretsmanager_secret" "twinmind" {
  name = "twinmind-api-credentials"
}

resource "aws_secretsmanager_secret_version" "twinmind" {
  secret_id = aws_secretsmanager_secret.twinmind.id
  secret_string = jsonencode({
    api_key        = var.twinmind_api_key
    webhook_secret = var.twinmind_webhook_secret
  })
}

# ECS Cluster
resource "aws_ecs_cluster" "main" {
  name = "twinmind-cluster"

  setting {
    name  = "containerInsights"
    value = "enabled"
  }
}

# Task Definition
resource "aws_ecs_task_definition" "twinmind" {
  family                   = "twinmind-service"
  network_mode             = "awsvpc"
  requires_compatibilities = ["FARGATE"]
  cpu                      = 512
  memory                   = 1024
  execution_role_arn       = aws_iam_role.ecs_execution.arn
  task_role_arn            = aws_iam_role.ecs_task.arn

  container_definitions = jsonencode([
    {
      name  = "twinmind"
      image = "${aws_ecr_repository.twinmind.repository_url}:latest"

      portMappings = [
        {
          containerPort = 8080
          protocol      = "tcp"
        }
      ]

      environment = [
        { name = "NODE_ENV", value = "production" },
        { name = "PORT", value = "8080" }
      ]

      secrets = [
        {
          name      = "TWINMIND_API_KEY"
          valueFrom = "${aws_secretsmanager_secret.twinmind.arn}:api_key::"
        },
        {
          name      = "TWINMIND_WEBHOOK_SECRET"
          valueFrom = "${aws_secretsmanager_secret.twinmind.arn}:webhook_secret::"
        }
      ]

      logConfiguration = {
        logDriver = "awslogs"
        options = {
          awslogs-group         = aws_cloudwatch_log_group.twinmind.name
          awslogs-region        = var.aws_region
          awslogs-stream-prefix = "ecs"
        }
      }

      healthCheck = {
        command     = ["CMD-SHELL", "wget -q --spider http://localhost:8080/health || exit 1"]
        interval    = 30
        timeout     = 5
        retries     = 3
        startPeriod = 60
      }
    }
  ])
}

# ECS Service
resource "aws_ecs_service" "twinmind" {
  name            = "twinmind-service"
  cluster         = aws_ecs_cluster.main.id
  task_definition = aws_ecs_task_definition.twinmind.arn
  desired_count   = 2
  launch_type     = "FARGATE"

  network_configuration {
    subnets          = var.private_subnet_ids
    security_groups  = [aws_security_group.ecs.id]
    assign_public_ip = false
  }

  load_balancer {
    target_group_arn = aws_lb_target_group.twinmind.arn
    container_name   = "twinmind"
    container_port   = 8080
  }

  deployment_circuit_breaker {
    enable   = true
    rollback = true
  }
}

# Auto-scaling
resource "aws_appautoscaling_target" "twinmind" {
  max_capacity       = 10
  min_capacity       = 2
  resource_id        = "service/${aws_ecs_cluster.main.name}/${aws_ecs_service.twinmind.name}"
  scalable_dimension = "ecs:service:DesiredCount"
  service_namespace  = "ecs"
}

resource "aws_appautoscaling_policy" "cpu" {
  name               = "cpu-autoscaling"
  policy_type        = "TargetTrackingScaling"
  resource_id        = aws_appautoscaling_target.twinmind.resource_id
  scalable_dimension = aws_appautoscaling_target.twinmind.scalable_dimension
  service_namespace  = aws_appautoscaling_target.twinmind.service_namespace

  target_tracking_scaling_policy_configuration {
    predefined_metric_specification {
      predefined_metric_type = "ECSServiceAverageCPUUtilization"
    }
    target_value       = 70.0
    scale_in_cooldown  = 300
    scale_out_cooldown = 60
  }
}
```

### Step 3: GCP Deployment (Cloud Run)

```hcl
# terraform/gcp/main.tf
terraform {
  required_providers {
    google = {
      source  = "hashicorp/google"
      version = "~> 5.0"
    }
  }
}

provider "google" {
  project = var.project_id
  region  = var.region
}

# Secret Manager for API keys
resource "google_secret_manager_secret" "twinmind_api_key" {
  secret_id = "twinmind-api-key"

  replication {
    auto {}
  }
}

resource "google_secret_manager_secret_version" "twinmind_api_key" {
  secret      = google_secret_manager_secret.twinmind_api_key.id
  secret_data = var.twinmind_api_key
}

# Cloud Run Service
resource "google_cloud_run_v2_service" "twinmind" {
  name     = "twinmind-service"
  location = var.region

  template {
    scaling {
      min_instance_count = 1
      max_instance_count = 10
    }

    containers {
      image = "${var.region}-docker.pkg.dev/${var.project_id}/twinmind/service:latest"

      ports {
        container_port = 8080
      }

      env {
        name  = "NODE_ENV"
        value = "production"
      }

      env {
        name = "TWINMIND_API_KEY"
        value_source {
          secret_key_ref {
            secret  = google_secret_manager_secret.twinmind_api_key.secret_id
            version = "latest"
          }
        }
      }

      resources {
        limits = {
          cpu    = "1"
          memory = "512Mi"
        }
      }

      startup_probe {
        http_get {
          path = "/health"
          port = 8080
        }
        initial_delay_seconds = 10
        period_seconds        = 10
        failure_threshold     = 3
      }

      liveness_probe {
        http_get {
          path = "/health"
          port = 8080
        }
        period_seconds = 30
      }
    }

    service_account = google_service_account.twinmind.email
  }

  traffic {
    percent = 100
    type    = "TRAFFIC_TARGET_ALLOCATION_TYPE_LATEST"
  }
}

# Service Account
resource "google_service_account" "twinmind" {
  account_id   = "twinmind-service"
  display_name = "TwinMind Service Account"
}

# Allow access to secrets
resource "google_secret_manager_secret_iam_member" "twinmind" {
  secret_id = google_secret_manager_secret.twinmind_api_key.id
  role      = "roles/secretmanager.secretAccessor"
  member    = "serviceAccount:${google_service_account.twinmind.email}"
}
```

### Step 4: Kubernetes Deployment

```yaml
# k8s/deployment.yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: twinmind-service
  labels:
    app: twinmind
spec:
  replicas: 2
  selector:
    matchLabels:
      app: twinmind
  template:
    metadata:
      labels:
        app: twinmind
    spec:
      serviceAccountName: twinmind
      containers:
        - name: twinmind
          image: your-registry/twinmind:latest
          ports:
            - containerPort: 8080
          env:
            - name: NODE_ENV
              value: "production"
            - name: TWINMIND_API_KEY
              valueFrom:
                secretKeyRef:
                  name: twinmind-secrets
                  key: api-key
            - name: TWINMIND_WEBHOOK_SECRET
              valueFrom:
                secretKeyRef:
                  name: twinmind-secrets
                  key: webhook-secret
          resources:
            requests:
              cpu: 250m
              memory: 256Mi
            limits:
              cpu: 1000m
              memory: 512Mi
          readinessProbe:
            httpGet:
              path: /health
              port: 8080
            initialDelaySeconds: 10
            periodSeconds: 5
          livenessProbe:
            httpGet:
              path: /health
              port: 8080
            initialDelaySeconds: 30
            periodSeconds: 10
          securityContext:
            runAsNonRoot: true
            runAsUser: 1001
            allowPrivilegeEscalation: false
---
apiVersion: v1
kind: Service
metadata:
  name: twinmind-service
spec:
  selector:
    app: twinmind
  ports:
    - port: 80
      targetPort: 8080
  type: ClusterIP
---
apiVersion: autoscaling/v2
kind: HorizontalPodAutoscaler
metadata:
  name: twinmind-hpa
spec:
  scaleTargetRef:
    apiVersion: apps/v1
    kind: Deployment
    name: twinmind-service
  minReplicas: 2
  maxReplicas: 10
  metrics:
    - type: Resource
      resource:
        name: cpu
        target:
          type: Utilization
          averageUtilization: 70
```

### Step 5: GitHub Actions Deployment

```yaml
# .github/workflows/deploy.yml
name: Deploy TwinMind Service

on:
  push:
    branches: [main]
    tags: ['v*']

env:
  AWS_REGION: us-east-1
  ECR_REPOSITORY: twinmind

jobs:
  build-and-push:
    runs-on: ubuntu-latest
    outputs:
      image: ${{ steps.build.outputs.image }}

    steps:
      - uses: actions/checkout@v4

      - name: Configure AWS credentials
        uses: aws-actions/configure-aws-credentials@v4
        with:
          aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}
          aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
          aws-region: ${{ env.AWS_REGION }}

      - name: Login to ECR
        id: login-ecr
        uses: aws-actions/amazon-ecr-login@v2

      - name: Build and push image
        id: build
        env:
          ECR_REGISTRY: ${{ steps.login-ecr.outputs.registry }}
          IMAGE_TAG: ${{ github.sha }}
        run: |
          docker build -t $ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG .
          docker push $ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG
          echo "image=$ECR_REGISTRY/$ECR_REPOSITORY:$IMAGE_TAG" >> $GITHUB_OUTPUT

  deploy-staging:
    needs: build-and-push
    runs-on: ubuntu-latest
    environment: staging

    steps:
      - uses: actions/checkout@v4

      - name: Configure AWS credentials
        uses: aws-actions/configure-aws-credentials@v4
        with:
          aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}
          aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
          aws-region: ${{ env.AWS_REGION }}

      - name: Deploy to ECS
        run: |
          aws ecs update-service \
            --cluster twinmind-staging \
            --service twinmind-service \
            --force-new-deployment

      - name: Wait for deployment
        run: |
          aws ecs wait services-stable \
            --cluster twinmind-staging \
            --services twinmind-service

  deploy-production:
    needs: [build-and-push, deploy-staging]
    runs-on: ubuntu-latest
    environment: production

    steps:
      - uses: actions/checkout@v4

      - name: Configure AWS credentials
        uses: aws-actions/configure-aws-credentials@v4
        with:
          aws-access-key-id: ${{ secrets.AWS_ACCESS_KEY_ID }}
          aws-secret-access-key: ${{ secrets.AWS_SECRET_ACCESS_KEY }}
          aws-region: ${{ env.AWS_REGION }}

      - name: Deploy to production
        run: |
          aws ecs update-service \
            --cluster twinmind-production \
            --service twinmind-service \
            --force-new-deployment

      - name: Wait for deployment
        run: |
          aws ecs wait services-stable \
            --cluster twinmind-production \
            --services twinmind-service
```

## Output
- Docker configuration files
- AWS ECS/Fargate Terraform
- GCP Cloud Run Terraform
- Kubernetes manifests
- GitHub Actions deployment workflow

## Error Handling

| Issue | Cause | Solution |
|-------|-------|----------|
| Image pull failed | Registry auth | Verify credentials |
| Health check failed | Service not ready | Increase start period |
| Scaling not working | Metrics missing | Check monitoring config |
| Secrets not found | IAM permissions | Update service account |

## Resources
- [AWS ECS Documentation](https://docs.aws.amazon.com/ecs/)
- [GCP Cloud Run Documentation](https://cloud.google.com/run/docs)
- [Kubernetes Documentation](https://kubernetes.io/docs/)

## Next Steps
For webhook handling, see `twinmind-webhooks-events`.
