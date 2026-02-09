---
name: langfuse-enterprise-rbac
description: |
  Configure Langfuse enterprise organization management and access control.
  Use when implementing team access controls, configuring organization settings,
  or setting up role-based permissions for Langfuse projects.
  Trigger with phrases like "langfuse RBAC", "langfuse teams",
  "langfuse organization", "langfuse access control", "langfuse permissions".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Enterprise RBAC

## Overview
Configure enterprise-grade access control for Langfuse projects and organizations.

## Prerequisites
- Langfuse organization/team account
- Understanding of role-based access patterns
- SSO provider (optional, for enterprise SSO)

## Role Definitions

| Role | Dashboard | API Read | API Write | Settings | Billing |
|------|-----------|----------|-----------|----------|---------|
| Owner | Full | Full | Full | Full | Full |
| Admin | Full | Full | Full | Full | View |
| Member | Full | Full | Full | None | None |
| Viewer | View | Read | None | None | None |
| API Only | None | Full | Full | None | None |

## Instructions

### Step 1: Implement Role-Based Access in Your Application

```typescript
// lib/langfuse/rbac.ts

enum LangfuseRole {
  Owner = "owner",
  Admin = "admin",
  Member = "member",
  Viewer = "viewer",
  ApiOnly = "api_only",
}

interface LangfusePermissions {
  canViewDashboard: boolean;
  canReadTraces: boolean;
  canWriteTraces: boolean;
  canManageSettings: boolean;
  canManageMembers: boolean;
  canManageBilling: boolean;
  canCreateApiKeys: boolean;
  canDeleteData: boolean;
}

const ROLE_PERMISSIONS: Record<LangfuseRole, LangfusePermissions> = {
  [LangfuseRole.Owner]: {
    canViewDashboard: true,
    canReadTraces: true,
    canWriteTraces: true,
    canManageSettings: true,
    canManageMembers: true,
    canManageBilling: true,
    canCreateApiKeys: true,
    canDeleteData: true,
  },
  [LangfuseRole.Admin]: {
    canViewDashboard: true,
    canReadTraces: true,
    canWriteTraces: true,
    canManageSettings: true,
    canManageMembers: true,
    canManageBilling: false,
    canCreateApiKeys: true,
    canDeleteData: true,
  },
  [LangfuseRole.Member]: {
    canViewDashboard: true,
    canReadTraces: true,
    canWriteTraces: true,
    canManageSettings: false,
    canManageMembers: false,
    canManageBilling: false,
    canCreateApiKeys: false,
    canDeleteData: false,
  },
  [LangfuseRole.Viewer]: {
    canViewDashboard: true,
    canReadTraces: true,
    canWriteTraces: false,
    canManageSettings: false,
    canManageMembers: false,
    canManageBilling: false,
    canCreateApiKeys: false,
    canDeleteData: false,
  },
  [LangfuseRole.ApiOnly]: {
    canViewDashboard: false,
    canReadTraces: true,
    canWriteTraces: true,
    canManageSettings: false,
    canManageMembers: false,
    canManageBilling: false,
    canCreateApiKeys: false,
    canDeleteData: false,
  },
};

function hasPermission(
  role: LangfuseRole,
  permission: keyof LangfusePermissions
): boolean {
  return ROLE_PERMISSIONS[role][permission];
}

function checkPermission(
  role: LangfuseRole,
  permission: keyof LangfusePermissions
): void {
  if (!hasPermission(role, permission)) {
    throw new ForbiddenError(
      `Permission denied: ${permission} requires higher role than ${role}`
    );
  }
}
```

### Step 2: Implement API Key Scoping

```typescript
// lib/langfuse/scoped-keys.ts

interface ScopedApiKey {
  id: string;
  publicKey: string;
  secretKey: string;
  name: string;
  scope: ApiKeyScope;
  createdAt: Date;
  createdBy: string;
  expiresAt?: Date;
}

interface ApiKeyScope {
  permissions: ("read" | "write")[];
  projects?: string[];       // Restrict to specific projects
  environments?: string[];   // Restrict to specific environments
  ipAllowlist?: string[];    // Restrict by IP
  rateLimit?: number;        // Custom rate limit
}

class ScopedLangfuseClient {
  private langfuse: Langfuse;
  private scope: ApiKeyScope;

  constructor(apiKey: ScopedApiKey) {
    this.langfuse = new Langfuse({
      publicKey: apiKey.publicKey,
      secretKey: apiKey.secretKey,
    });
    this.scope = apiKey.scope;
  }

  trace(params: Parameters<typeof this.langfuse.trace>[0]) {
    // Check write permission
    if (!this.scope.permissions.includes("write")) {
      throw new ForbiddenError("API key does not have write permission");
    }

    // Check project restriction
    if (this.scope.projects?.length) {
      const projectId = params.metadata?.projectId;
      if (projectId && !this.scope.projects.includes(projectId)) {
        throw new ForbiddenError(`API key not authorized for project ${projectId}`);
      }
    }

    return this.langfuse.trace(params);
  }

  async fetchTraces(params: any) {
    // Check read permission
    if (!this.scope.permissions.includes("read")) {
      throw new ForbiddenError("API key does not have read permission");
    }

    return this.langfuse.fetchTraces(params);
  }
}
```

### Step 3: Implement Project-Based Access Control

```typescript
// lib/langfuse/projects.ts

interface LangfuseProject {
  id: string;
  name: string;
  organizationId: string;
  environment: "development" | "staging" | "production";
  members: ProjectMember[];
}

interface ProjectMember {
  userId: string;
  email: string;
  role: LangfuseRole;
  addedAt: Date;
  addedBy: string;
}

class ProjectAccessController {
  private projects: Map<string, LangfuseProject> = new Map();

  async getUserProjects(userId: string): Promise<LangfuseProject[]> {
    return Array.from(this.projects.values()).filter((project) =>
      project.members.some((m) => m.userId === userId)
    );
  }

  async getUserRole(userId: string, projectId: string): Promise<LangfuseRole | null> {
    const project = this.projects.get(projectId);
    if (!project) return null;

    const member = project.members.find((m) => m.userId === userId);
    return member?.role || null;
  }

  async addMember(
    projectId: string,
    member: Omit<ProjectMember, "addedAt">,
    addedBy: string
  ): Promise<void> {
    // Check if adder has permission
    const adderRole = await this.getUserRole(addedBy, projectId);
    if (!adderRole || !hasPermission(adderRole, "canManageMembers")) {
      throw new ForbiddenError("Cannot add members to this project");
    }

    // Prevent privilege escalation
    if (this.isHigherRole(member.role, adderRole)) {
      throw new ForbiddenError("Cannot add member with higher role than yourself");
    }

    const project = this.projects.get(projectId);
    if (!project) throw new NotFoundError("Project not found");

    project.members.push({
      ...member,
      addedAt: new Date(),
      addedBy,
    });
  }

  private isHigherRole(role1: LangfuseRole, role2: LangfuseRole): boolean {
    const hierarchy = [
      LangfuseRole.Viewer,
      LangfuseRole.ApiOnly,
      LangfuseRole.Member,
      LangfuseRole.Admin,
      LangfuseRole.Owner,
    ];
    return hierarchy.indexOf(role1) > hierarchy.indexOf(role2);
  }
}
```

### Step 4: Implement SSO Integration

```typescript
// lib/langfuse/sso.ts

interface SSOConfig {
  provider: "okta" | "azure" | "google" | "saml";
  domain: string;
  clientId: string;
  clientSecret: string;
  callbackUrl: string;
  groupMapping: Record<string, LangfuseRole>;
}

class LangfuseSSO {
  private config: SSOConfig;

  constructor(config: SSOConfig) {
    this.config = config;
  }

  // Map SSO groups to Langfuse roles
  mapGroupsToRole(groups: string[]): LangfuseRole {
    for (const group of groups) {
      if (this.config.groupMapping[group]) {
        return this.config.groupMapping[group];
      }
    }
    return LangfuseRole.Viewer; // Default role
  }

  // SAML assertion handler
  async handleSAMLAssertion(assertion: any): Promise<{
    userId: string;
    email: string;
    role: LangfuseRole;
    groups: string[];
  }> {
    const email = assertion.user.email;
    const groups = assertion.user.groups || [];

    // Verify domain
    if (!email.endsWith(`@${this.config.domain}`)) {
      throw new UnauthorizedError("Email domain not allowed");
    }

    return {
      userId: assertion.user.id,
      email,
      role: this.mapGroupsToRole(groups),
      groups,
    };
  }
}

// Usage
const ssoConfig: SSOConfig = {
  provider: "okta",
  domain: "company.com",
  clientId: process.env.OKTA_CLIENT_ID!,
  clientSecret: process.env.OKTA_CLIENT_SECRET!,
  callbackUrl: "https://app.company.com/auth/callback",
  groupMapping: {
    "Engineering-Admins": LangfuseRole.Admin,
    "Engineering": LangfuseRole.Member,
    "Data-Team": LangfuseRole.Viewer,
  },
};
```

### Step 5: Implement Audit Logging for RBAC

```typescript
// lib/langfuse/rbac-audit.ts

interface RBACEvent {
  timestamp: Date;
  action:
    | "member_added"
    | "member_removed"
    | "role_changed"
    | "api_key_created"
    | "api_key_revoked"
    | "permission_denied";
  actor: string;
  target: string;
  projectId?: string;
  details: Record<string, any>;
}

class RBACAuditLogger {
  async log(event: Omit<RBACEvent, "timestamp">) {
    const auditEvent: RBACEvent = {
      ...event,
      timestamp: new Date(),
    };

    // Log to audit trail
    console.log("[RBAC Audit]", JSON.stringify(auditEvent));

    // Store in database
    await this.persist(auditEvent);

    // Alert on sensitive actions
    if (event.action === "permission_denied") {
      await this.alertSecurityTeam(auditEvent);
    }
  }

  private async persist(event: RBACEvent) {
    // Store in audit database
  }

  private async alertSecurityTeam(event: RBACEvent) {
    // Send alert for potential security issues
  }

  async query(options: {
    actor?: string;
    action?: string;
    projectId?: string;
    from?: Date;
    to?: Date;
  }): Promise<RBACEvent[]> {
    // Query audit log
    return [];
  }
}

export const rbacAudit = new RBACAuditLogger();
```

## Output
- Role-based permission system
- Scoped API key implementation
- Project-based access control
- SSO integration
- RBAC audit logging

## Access Control Matrix

| Action | Owner | Admin | Member | Viewer | API Only |
|--------|-------|-------|--------|--------|----------|
| View traces | Yes | Yes | Yes | Yes | No |
| Create traces | Yes | Yes | Yes | No | Yes |
| Delete traces | Yes | Yes | No | No | No |
| Manage members | Yes | Yes | No | No | No |
| Create API keys | Yes | Yes | No | No | No |
| Manage billing | Yes | No | No | No | No |

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Permission denied | Insufficient role | Request role upgrade |
| SSO login fails | Wrong group mapping | Update group mapping |
| API key rejected | Expired or scoped | Create new key |
| Project not found | No access | Request project access |

## Resources
- [Langfuse Organizations](https://langfuse.com/docs/teams)
- [Langfuse API Keys](https://langfuse.com/docs/api-reference)
- [SAML 2.0 Specification](https://wiki.oasis-open.org/security/FrontPage)

## Next Steps
For complex migrations, see `langfuse-migration-deep-dive`.
