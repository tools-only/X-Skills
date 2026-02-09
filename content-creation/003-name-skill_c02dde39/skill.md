---
name: twinmind-enterprise-rbac
description: |
  Implement enterprise role-based access control for TwinMind.
  Use when setting up team permissions, configuring SSO/SAML,
  or implementing organization-level access policies.
  Trigger with phrases like "twinmind RBAC", "twinmind permissions",
  "twinmind enterprise access", "twinmind SSO", "twinmind team management".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# TwinMind Enterprise RBAC

## Overview
Enterprise-grade role-based access control for TwinMind with SSO integration, team management, and audit capabilities.

## Prerequisites
- TwinMind Enterprise account
- Admin access to TwinMind organization
- Identity provider (Okta, Azure AD, Google Workspace)
- Understanding of RBAC concepts

## Instructions

### Step 1: Define Roles and Permissions

```typescript
// src/twinmind/rbac/roles.ts
export enum Permission {
  // Transcript permissions
  TRANSCRIPT_CREATE = 'transcript:create',
  TRANSCRIPT_READ = 'transcript:read',
  TRANSCRIPT_READ_ALL = 'transcript:read:all',  // Read any transcript
  TRANSCRIPT_DELETE = 'transcript:delete',
  TRANSCRIPT_DELETE_ALL = 'transcript:delete:all',

  // Summary permissions
  SUMMARY_GENERATE = 'summary:generate',
  SUMMARY_READ = 'summary:read',

  // Action item permissions
  ACTION_ITEM_CREATE = 'action_item:create',
  ACTION_ITEM_ASSIGN = 'action_item:assign',
  ACTION_ITEM_DELETE = 'action_item:delete',

  // Team permissions
  TEAM_VIEW = 'team:view',
  TEAM_MANAGE = 'team:manage',
  TEAM_INVITE = 'team:invite',
  TEAM_REMOVE = 'team:remove',

  // Settings permissions
  SETTINGS_VIEW = 'settings:view',
  SETTINGS_MANAGE = 'settings:manage',
  INTEGRATIONS_MANAGE = 'integrations:manage',

  // Billing permissions
  BILLING_VIEW = 'billing:view',
  BILLING_MANAGE = 'billing:manage',

  // Admin permissions
  ADMIN_FULL = 'admin:*',
  AUDIT_VIEW = 'audit:view',
  COMPLIANCE_MANAGE = 'compliance:manage',
}

export interface Role {
  name: string;
  description: string;
  permissions: Permission[];
  inherits?: string[];  // Inherit from other roles
}

export const predefinedRoles: Record<string, Role> = {
  viewer: {
    name: 'Viewer',
    description: 'Can view own transcripts and summaries',
    permissions: [
      Permission.TRANSCRIPT_READ,
      Permission.SUMMARY_READ,
    ],
  },

  member: {
    name: 'Member',
    description: 'Standard team member with create access',
    permissions: [
      Permission.TRANSCRIPT_CREATE,
      Permission.TRANSCRIPT_READ,
      Permission.TRANSCRIPT_DELETE,
      Permission.SUMMARY_GENERATE,
      Permission.SUMMARY_READ,
      Permission.ACTION_ITEM_CREATE,
    ],
    inherits: ['viewer'],
  },

  manager: {
    name: 'Manager',
    description: 'Can manage team transcripts and members',
    permissions: [
      Permission.TRANSCRIPT_READ_ALL,
      Permission.TRANSCRIPT_DELETE_ALL,
      Permission.ACTION_ITEM_ASSIGN,
      Permission.ACTION_ITEM_DELETE,
      Permission.TEAM_VIEW,
      Permission.TEAM_INVITE,
    ],
    inherits: ['member'],
  },

  admin: {
    name: 'Admin',
    description: 'Full organizational admin access',
    permissions: [
      Permission.TEAM_MANAGE,
      Permission.TEAM_REMOVE,
      Permission.SETTINGS_VIEW,
      Permission.SETTINGS_MANAGE,
      Permission.INTEGRATIONS_MANAGE,
      Permission.BILLING_VIEW,
      Permission.AUDIT_VIEW,
    ],
    inherits: ['manager'],
  },

  owner: {
    name: 'Owner',
    description: 'Organization owner with full access',
    permissions: [
      Permission.ADMIN_FULL,
      Permission.BILLING_MANAGE,
      Permission.COMPLIANCE_MANAGE,
    ],
    inherits: ['admin'],
  },
};

// Resolve all permissions for a role (including inherited)
export function resolvePermissions(roleName: string): Permission[] {
  const role = predefinedRoles[roleName];
  if (!role) return [];

  const permissions = new Set<Permission>(role.permissions);

  for (const inheritedRole of role.inherits || []) {
    const inheritedPermissions = resolvePermissions(inheritedRole);
    inheritedPermissions.forEach(p => permissions.add(p));
  }

  return Array.from(permissions);
}
```

### Step 2: Implement Permission Checking

```typescript
// src/twinmind/rbac/authorization.ts
import { Permission, resolvePermissions } from './roles';

export interface User {
  id: string;
  email: string;
  organizationId: string;
  role: string;
  customPermissions?: Permission[];  // Additional granted permissions
  deniedPermissions?: Permission[];  // Explicitly denied
}

export class AuthorizationService {
  hasPermission(user: User, permission: Permission): boolean {
    // Check for admin wildcard
    const userPermissions = this.getUserPermissions(user);

    if (userPermissions.includes(Permission.ADMIN_FULL)) {
      return true;
    }

    // Check denied list first
    if (user.deniedPermissions?.includes(permission)) {
      return false;
    }

    // Check for exact permission or wildcard
    return userPermissions.includes(permission) ||
           this.matchesWildcard(permission, userPermissions);
  }

  hasAnyPermission(user: User, permissions: Permission[]): boolean {
    return permissions.some(p => this.hasPermission(user, p));
  }

  hasAllPermissions(user: User, permissions: Permission[]): boolean {
    return permissions.every(p => this.hasPermission(user, p));
  }

  private getUserPermissions(user: User): Permission[] {
    const rolePermissions = resolvePermissions(user.role);
    const customPermissions = user.customPermissions || [];

    return [...new Set([...rolePermissions, ...customPermissions])];
  }

  private matchesWildcard(permission: Permission, userPermissions: Permission[]): boolean {
    const [resource] = permission.split(':');
    const wildcardPermission = `${resource}:*` as Permission;
    return userPermissions.includes(wildcardPermission);
  }
}

// Express middleware
export function requirePermission(...permissions: Permission[]) {
  const authService = new AuthorizationService();

  return (req: Request, res: Response, next: NextFunction) => {
    const user = req.user as User;

    if (!user) {
      return res.status(401).json({ error: 'Unauthorized' });
    }

    if (!authService.hasAnyPermission(user, permissions)) {
      return res.status(403).json({
        error: 'Insufficient permissions',
        required: permissions,
        user_role: user.role,
      });
    }

    next();
  };
}

// Resource-level authorization
export function canAccessTranscript(user: User, transcript: Transcript): boolean {
  const authService = new AuthorizationService();

  // Admins can access all
  if (authService.hasPermission(user, Permission.TRANSCRIPT_READ_ALL)) {
    return true;
  }

  // Check if user owns the transcript
  if (transcript.created_by === user.id) {
    return authService.hasPermission(user, Permission.TRANSCRIPT_READ);
  }

  // Check if user was a participant
  if (transcript.participants?.some(p => p.email === user.email)) {
    return authService.hasPermission(user, Permission.TRANSCRIPT_READ);
  }

  return false;
}
```

### Step 3: Configure SSO/SAML Integration

```typescript
// src/twinmind/rbac/sso.ts
import { Strategy as SamlStrategy } from 'passport-saml';
import passport from 'passport';

export interface SSOConfig {
  provider: 'okta' | 'azure' | 'google' | 'onelogin';
  entryPoint: string;          // IdP SSO URL
  issuer: string;              // SP Entity ID
  cert: string;                // IdP Certificate
  callbackUrl: string;         // ACS URL
  signatureAlgorithm?: string;
  identifierFormat?: string;
  attributeMapping: {
    email: string;
    firstName: string;
    lastName: string;
    groups?: string;
  };
  groupRoleMapping?: Record<string, string>;  // IdP group -> TwinMind role
}

export function configureSAML(config: SSOConfig): void {
  const samlStrategy = new SamlStrategy(
    {
      entryPoint: config.entryPoint,
      issuer: config.issuer,
      cert: config.cert,
      callbackUrl: config.callbackUrl,
      signatureAlgorithm: config.signatureAlgorithm || 'sha256',
      identifierFormat: config.identifierFormat || 'urn:oasis:names:tc:SAML:1.1:nameid-format:emailAddress',
    },
    async (profile, done) => {
      try {
        // Extract user attributes
        const email = profile[config.attributeMapping.email];
        const firstName = profile[config.attributeMapping.firstName];
        const lastName = profile[config.attributeMapping.lastName];
        const groups = profile[config.attributeMapping.groups || 'groups'] || [];

        // Map IdP groups to TwinMind roles
        let role = 'member';  // Default role
        if (config.groupRoleMapping) {
          for (const [group, mappedRole] of Object.entries(config.groupRoleMapping)) {
            if (groups.includes(group)) {
              role = mappedRole;
              break;
            }
          }
        }

        // Find or create user
        const user = await findOrCreateUser({
          email,
          firstName,
          lastName,
          role,
          authProvider: 'saml',
          authProviderId: profile.nameID,
        });

        done(null, user);
      } catch (error) {
        done(error as Error);
      }
    }
  );

  passport.use('saml', samlStrategy);
}

// Example Okta configuration
export const oktaConfig: SSOConfig = {
  provider: 'okta',
  entryPoint: process.env.OKTA_SSO_URL!,
  issuer: process.env.OKTA_ISSUER!,
  cert: process.env.OKTA_CERT!,
  callbackUrl: `${process.env.APP_URL}/auth/saml/callback`,
  attributeMapping: {
    email: 'http://schemas.xmlsoap.org/ws/2005/05/identity/claims/emailaddress',
    firstName: 'http://schemas.xmlsoap.org/ws/2005/05/identity/claims/givenname',
    lastName: 'http://schemas.xmlsoap.org/ws/2005/05/identity/claims/surname',
    groups: 'http://schemas.xmlsoap.org/claims/Group',
  },
  groupRoleMapping: {
    'TwinMind-Admins': 'admin',
    'TwinMind-Managers': 'manager',
    'TwinMind-Members': 'member',
    'TwinMind-Viewers': 'viewer',
  },
};
```

### Step 4: Team Management

```typescript
// src/twinmind/rbac/teams.ts
export interface Team {
  id: string;
  name: string;
  organizationId: string;
  members: TeamMember[];
  settings: TeamSettings;
}

export interface TeamMember {
  userId: string;
  email: string;
  role: string;
  joinedAt: Date;
  invitedBy: string;
}

export interface TeamSettings {
  defaultTranscriptVisibility: 'private' | 'team' | 'organization';
  allowExternalSharing: boolean;
  requireApprovalForSharing: boolean;
  autoTranscribe: boolean;
}

export class TeamManager {
  private client = getTwinMindClient();
  private authService = new AuthorizationService();

  async createTeam(name: string, creatorId: string): Promise<Team> {
    const response = await this.client.post('/teams', {
      name,
      settings: {
        defaultTranscriptVisibility: 'team',
        allowExternalSharing: false,
        requireApprovalForSharing: true,
        autoTranscribe: true,
      },
    });

    // Add creator as owner
    await this.addMember(response.data.id, creatorId, 'owner');

    return response.data;
  }

  async addMember(
    teamId: string,
    userId: string,
    role: string,
    invitedBy?: string
  ): Promise<void> {
    await this.client.post(`/teams/${teamId}/members`, {
      user_id: userId,
      role,
      invited_by: invitedBy,
    });

    // Sync to TwinMind organization
    await this.syncToTwinMind(teamId);
  }

  async removeMember(teamId: string, userId: string): Promise<void> {
    await this.client.delete(`/teams/${teamId}/members/${userId}`);
    await this.syncToTwinMind(teamId);
  }

  async updateMemberRole(teamId: string, userId: string, newRole: string): Promise<void> {
    await this.client.patch(`/teams/${teamId}/members/${userId}`, {
      role: newRole,
    });
  }

  async getTeamMembers(teamId: string): Promise<TeamMember[]> {
    const response = await this.client.get(`/teams/${teamId}/members`);
    return response.data;
  }

  async inviteMember(
    teamId: string,
    email: string,
    role: string,
    invitedBy: User
  ): Promise<void> {
    // Check permission
    if (!this.authService.hasPermission(invitedBy, Permission.TEAM_INVITE)) {
      throw new Error('Insufficient permissions to invite team members');
    }

    // Send invitation
    await this.client.post(`/teams/${teamId}/invitations`, {
      email,
      role,
      invited_by: invitedBy.id,
    });

    // Send invitation email
    await sendInvitationEmail(email, teamId, invitedBy);
  }

  private async syncToTwinMind(teamId: string): Promise<void> {
    const members = await this.getTeamMembers(teamId);

    await this.client.put(`/organization/members`, {
      members: members.map(m => ({
        email: m.email,
        role: m.role,
      })),
    });
  }
}
```

### Step 5: Audit Logging for Access

```typescript
// src/twinmind/rbac/audit.ts
export interface AccessAuditEvent {
  timestamp: Date;
  userId: string;
  userEmail: string;
  action: string;
  resource: string;
  resourceId?: string;
  granted: boolean;
  reason?: string;
  ipAddress?: string;
  userAgent?: string;
  metadata?: Record<string, any>;
}

export class AccessAuditLogger {
  async logAccess(event: Omit<AccessAuditEvent, 'timestamp'>): Promise<void> {
    const fullEvent: AccessAuditEvent = {
      ...event,
      timestamp: new Date(),
    };

    // Store in database
    await db.accessAuditLogs.create(fullEvent);

    // Send to SIEM if configured
    if (process.env.SIEM_ENDPOINT) {
      await this.sendToSIEM(fullEvent);
    }

    // Log denied access at warning level
    if (!event.granted) {
      logger.warn('Access denied', {
        user: event.userEmail,
        action: event.action,
        resource: event.resource,
        reason: event.reason,
      });
    }
  }

  async queryAuditLogs(filters: {
    userId?: string;
    action?: string;
    resource?: string;
    granted?: boolean;
    from?: Date;
    to?: Date;
    limit?: number;
  }): Promise<AccessAuditEvent[]> {
    return db.accessAuditLogs.find({
      ...filters,
      timestamp: {
        $gte: filters.from,
        $lte: filters.to,
      },
    }).limit(filters.limit || 100);
  }

  async generateAccessReport(
    userId: string,
    dateRange: { from: Date; to: Date }
  ): Promise<{
    totalAccesses: number;
    deniedAccesses: number;
    resourcesAccessed: string[];
    uniqueResources: number;
  }> {
    const logs = await this.queryAuditLogs({
      userId,
      from: dateRange.from,
      to: dateRange.to,
    });

    const resourcesAccessed = logs
      .filter(l => l.granted)
      .map(l => l.resource);

    return {
      totalAccesses: logs.filter(l => l.granted).length,
      deniedAccesses: logs.filter(l => !l.granted).length,
      resourcesAccessed,
      uniqueResources: new Set(resourcesAccessed).size,
    };
  }

  private async sendToSIEM(event: AccessAuditEvent): Promise<void> {
    await fetch(process.env.SIEM_ENDPOINT!, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        source: 'twinmind-integration',
        event_type: 'access_audit',
        ...event,
      }),
    });
  }
}

// Middleware to audit all accesses
export function auditMiddleware(auditLogger: AccessAuditLogger) {
  return async (req: Request, res: Response, next: NextFunction) => {
    const originalSend = res.send;
    const user = req.user as User;

    res.send = function(body) {
      const granted = res.statusCode < 400;

      auditLogger.logAccess({
        userId: user?.id || 'anonymous',
        userEmail: user?.email || 'anonymous',
        action: `${req.method} ${req.path}`,
        resource: req.path.split('/')[2] || 'unknown',
        resourceId: req.params.id,
        granted,
        reason: !granted ? body?.error : undefined,
        ipAddress: req.ip,
        userAgent: req.headers['user-agent'],
      });

      return originalSend.call(this, body);
    };

    next();
  };
}
```

## Output
- Role and permission definitions
- Permission checking service
- SSO/SAML configuration
- Team management functionality
- Access audit logging

## Role Hierarchy

```
Owner
  └── Admin
        └── Manager
              └── Member
                    └── Viewer
```

## Enterprise Features

| Feature | Description | Configuration |
|---------|-------------|---------------|
| SSO/SAML | Single sign-on integration | IdP configuration |
| SCIM | User provisioning | Automatic sync |
| Custom roles | Define organization roles | Role builder |
| Audit logs | Access tracking | 90-day retention |
| IP allowlisting | Restrict access by IP | Network settings |

## Error Handling

| Issue | Cause | Solution |
|-------|-------|----------|
| SSO login failed | Certificate mismatch | Update IdP cert |
| Role not syncing | SCIM misconfigured | Check provisioning |
| Permission denied | Wrong role assigned | Review role mapping |
| Audit gaps | Middleware not applied | Add to all routes |

## Resources
- [TwinMind Enterprise](https://twinmind.com/enterprise)
- [SAML 2.0 Specification](https://docs.oasis-open.org/security/saml/v2.0/)
- [SCIM Protocol](https://scim.cloud/)

## Next Steps
For migration from other tools, see `twinmind-migration-deep-dive`.
