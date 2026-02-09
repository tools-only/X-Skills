---
name: openevidence-enterprise-rbac
description: |
  Configure OpenEvidence enterprise SSO, role-based access control, and organization management.
  Use when implementing SSO integration, configuring role-based permissions,
  or setting up organization-level controls for clinical AI applications.
  Trigger with phrases like "openevidence SSO", "openevidence RBAC",
  "openevidence enterprise", "openevidence roles", "openevidence permissions".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# OpenEvidence Enterprise RBAC

## Overview
Configure enterprise-grade access control for OpenEvidence clinical AI integrations in healthcare organizations.

## Prerequisites
- OpenEvidence Enterprise tier subscription
- Identity Provider (IdP) with SAML/OIDC support
- Understanding of role-based access patterns
- HIPAA audit logging infrastructure

## Role Definitions

| Role | Permissions | Use Case |
|------|-------------|----------|
| Physician | Full clinical query, DeepConsult | Active patient care |
| Nurse | Clinical query (no DeepConsult) | Nursing support |
| Pharmacist | Drug-focused queries | Medication management |
| Resident | Clinical query (supervised) | Training |
| Admin | Full access, user management | Platform administration |
| Auditor | Read-only audit logs | Compliance review |
| Integration | API access only | System integration |

## Instructions

### Step 1: Role and Permission Definitions
```typescript
// src/rbac/roles.ts

export enum ClinicalRole {
  Physician = 'physician',
  Nurse = 'nurse',
  Pharmacist = 'pharmacist',
  Resident = 'resident',
  Admin = 'admin',
  Auditor = 'auditor',
  Integration = 'integration',
}

export interface ClinicalPermissions {
  clinicalQuery: boolean;
  deepConsult: boolean;
  drugInfo: boolean;
  guidelineAccess: boolean;
  exportResults: boolean;
  viewAuditLogs: boolean;
  manageUsers: boolean;
  manageSettings: boolean;
}

export const ROLE_PERMISSIONS: Record<ClinicalRole, ClinicalPermissions> = {
  [ClinicalRole.Physician]: {
    clinicalQuery: true,
    deepConsult: true,
    drugInfo: true,
    guidelineAccess: true,
    exportResults: true,
    viewAuditLogs: false,
    manageUsers: false,
    manageSettings: false,
  },
  [ClinicalRole.Nurse]: {
    clinicalQuery: true,
    deepConsult: false,
    drugInfo: true,
    guidelineAccess: true,
    exportResults: false,
    viewAuditLogs: false,
    manageUsers: false,
    manageSettings: false,
  },
  [ClinicalRole.Pharmacist]: {
    clinicalQuery: true,
    deepConsult: false,
    drugInfo: true,
    guidelineAccess: true,
    exportResults: true,
    viewAuditLogs: false,
    manageUsers: false,
    manageSettings: false,
  },
  [ClinicalRole.Resident]: {
    clinicalQuery: true,
    deepConsult: false, // Requires attending approval
    drugInfo: true,
    guidelineAccess: true,
    exportResults: false,
    viewAuditLogs: false,
    manageUsers: false,
    manageSettings: false,
  },
  [ClinicalRole.Admin]: {
    clinicalQuery: true,
    deepConsult: true,
    drugInfo: true,
    guidelineAccess: true,
    exportResults: true,
    viewAuditLogs: true,
    manageUsers: true,
    manageSettings: true,
  },
  [ClinicalRole.Auditor]: {
    clinicalQuery: false,
    deepConsult: false,
    drugInfo: false,
    guidelineAccess: false,
    exportResults: false,
    viewAuditLogs: true,
    manageUsers: false,
    manageSettings: false,
  },
  [ClinicalRole.Integration]: {
    clinicalQuery: true,
    deepConsult: true,
    drugInfo: true,
    guidelineAccess: true,
    exportResults: false,
    viewAuditLogs: false,
    manageUsers: false,
    manageSettings: false,
  },
};

export function hasPermission(
  role: ClinicalRole,
  permission: keyof ClinicalPermissions
): boolean {
  return ROLE_PERMISSIONS[role][permission];
}

export function getPermissions(role: ClinicalRole): ClinicalPermissions {
  return ROLE_PERMISSIONS[role];
}
```

### Step 2: SSO Integration (SAML)
```typescript
// src/auth/saml.ts
import { Strategy as SamlStrategy } from 'passport-saml';
import passport from 'passport';

interface SAMLConfig {
  entryPoint: string;
  issuer: string;
  cert: string;
  callbackUrl: string;
  identifierFormat?: string;
}

export function configureSAML(config: SAMLConfig): void {
  passport.use(
    new SamlStrategy(
      {
        entryPoint: config.entryPoint,
        issuer: config.issuer,
        cert: config.cert,
        callbackUrl: config.callbackUrl,
        identifierFormat: config.identifierFormat || 'urn:oasis:names:tc:SAML:1.1:nameid-format:emailAddress',
      },
      async (profile, done) => {
        try {
          // Extract user info from SAML assertion
          const user = await findOrCreateUser({
            email: profile.nameID,
            firstName: profile.firstName,
            lastName: profile.lastName,
            groups: profile.groups || [],
          });

          // Map IdP groups to clinical roles
          const role = mapGroupsToRole(profile.groups);

          return done(null, { ...user, role });
        } catch (error) {
          return done(error);
        }
      }
    )
  );
}

// Map IdP groups to clinical roles
const GROUP_ROLE_MAPPING: Record<string, ClinicalRole> = {
  'Physicians': ClinicalRole.Physician,
  'Attending-Physicians': ClinicalRole.Physician,
  'Nursing': ClinicalRole.Nurse,
  'RN': ClinicalRole.Nurse,
  'Pharmacy': ClinicalRole.Pharmacist,
  'Residents': ClinicalRole.Resident,
  'IT-Admin': ClinicalRole.Admin,
  'Compliance': ClinicalRole.Auditor,
  'Service-Accounts': ClinicalRole.Integration,
};

function mapGroupsToRole(groups: string[]): ClinicalRole {
  // Priority order: Admin > Physician > Pharmacist > Nurse > Resident
  if (groups.some(g => GROUP_ROLE_MAPPING[g] === ClinicalRole.Admin)) {
    return ClinicalRole.Admin;
  }
  if (groups.some(g => GROUP_ROLE_MAPPING[g] === ClinicalRole.Physician)) {
    return ClinicalRole.Physician;
  }
  if (groups.some(g => GROUP_ROLE_MAPPING[g] === ClinicalRole.Pharmacist)) {
    return ClinicalRole.Pharmacist;
  }
  if (groups.some(g => GROUP_ROLE_MAPPING[g] === ClinicalRole.Nurse)) {
    return ClinicalRole.Nurse;
  }
  if (groups.some(g => GROUP_ROLE_MAPPING[g] === ClinicalRole.Resident)) {
    return ClinicalRole.Resident;
  }

  // Default to most restrictive role
  return ClinicalRole.Nurse;
}
```

### Step 3: OAuth2/OIDC Integration
```typescript
// src/auth/oidc.ts
import { Strategy as OpenIDConnectStrategy } from 'passport-openidconnect';
import passport from 'passport';

interface OIDCConfig {
  issuer: string;
  authorizationURL: string;
  tokenURL: string;
  userInfoURL: string;
  clientID: string;
  clientSecret: string;
  callbackURL: string;
  scope: string[];
}

export function configureOIDC(config: OIDCConfig): void {
  passport.use(
    new OpenIDConnectStrategy(
      {
        issuer: config.issuer,
        authorizationURL: config.authorizationURL,
        tokenURL: config.tokenURL,
        userInfoURL: config.userInfoURL,
        clientID: config.clientID,
        clientSecret: config.clientSecret,
        callbackURL: config.callbackURL,
        scope: config.scope,
      },
      async (issuer, profile, done) => {
        try {
          const user = await findOrCreateUser({
            email: profile.emails?.[0]?.value,
            firstName: profile.name?.givenName,
            lastName: profile.name?.familyName,
            providerId: profile.id,
          });

          // Get roles from custom claims or separate lookup
          const role = await getRoleFromClaims(profile._json);

          return done(null, { ...user, role });
        } catch (error) {
          return done(error);
        }
      }
    )
  );
}

async function getRoleFromClaims(claims: any): Promise<ClinicalRole> {
  // Custom claim for role
  if (claims['clinical_role']) {
    return claims['clinical_role'] as ClinicalRole;
  }

  // Fallback to group membership
  if (claims['groups']) {
    return mapGroupsToRole(claims['groups']);
  }

  return ClinicalRole.Nurse; // Most restrictive default
}
```

### Step 4: Permission Middleware
```typescript
// src/middleware/authorization.ts
import { Request, Response, NextFunction } from 'express';
import { ClinicalRole, hasPermission, ClinicalPermissions } from '../rbac/roles';
import { auditLogger } from '../compliance/audit-trail';

interface AuthenticatedRequest extends Request {
  user: {
    id: string;
    email: string;
    role: ClinicalRole;
  };
}

export function requirePermission(permission: keyof ClinicalPermissions) {
  return async (req: AuthenticatedRequest, res: Response, next: NextFunction) => {
    const user = req.user;

    if (!user) {
      return res.status(401).json({ error: 'Authentication required' });
    }

    if (!hasPermission(user.role, permission)) {
      // Audit failed access attempt
      await auditLogger.logAccess({
        userId: user.id,
        userRole: user.role,
        action: 'access_denied',
        resourceType: permission,
        resourceId: req.path,
        ipAddress: req.ip,
        userAgent: req.get('user-agent') || 'unknown',
      });

      return res.status(403).json({
        error: 'Forbidden',
        message: `Permission '${permission}' required for this action`,
        requiredRole: getRolesWithPermission(permission),
      });
    }

    next();
  };
}

function getRolesWithPermission(permission: keyof ClinicalPermissions): ClinicalRole[] {
  return Object.entries(ROLE_PERMISSIONS)
    .filter(([_, perms]) => perms[permission])
    .map(([role]) => role as ClinicalRole);
}

// Usage in routes
app.post('/api/clinical/query',
  requirePermission('clinicalQuery'),
  clinicalQueryHandler
);

app.post('/api/clinical/deepconsult',
  requirePermission('deepConsult'),
  deepConsultHandler
);

app.get('/api/admin/audit-logs',
  requirePermission('viewAuditLogs'),
  auditLogsHandler
);
```

### Step 5: Organization Management
```typescript
// src/rbac/organization.ts

interface Organization {
  id: string;
  name: string;
  openEvidenceOrgId: string;
  ssoEnabled: boolean;
  enforceSso: boolean;
  allowedDomains: string[];
  defaultRole: ClinicalRole;
  settings: OrganizationSettings;
}

interface OrganizationSettings {
  deepConsultEnabled: boolean;
  maxDeepConsultsPerDay: number;
  auditLogRetentionDays: number;
  allowExport: boolean;
  requireMFA: boolean;
}

export class OrganizationManager {
  constructor(private db: Database) {}

  async createOrganization(config: Omit<Organization, 'id'>): Promise<Organization> {
    const org = await this.db.organizations.create({
      data: {
        ...config,
        id: crypto.randomUUID(),
      },
    });

    return org;
  }

  async updateSettings(
    orgId: string,
    settings: Partial<OrganizationSettings>
  ): Promise<Organization> {
    return this.db.organizations.update({
      where: { id: orgId },
      data: { settings },
    });
  }

  async addUser(orgId: string, userId: string, role: ClinicalRole): Promise<void> {
    await this.db.organizationUsers.create({
      data: {
        organizationId: orgId,
        userId,
        role,
        addedAt: new Date(),
      },
    });
  }

  async updateUserRole(orgId: string, userId: string, newRole: ClinicalRole): Promise<void> {
    await this.db.organizationUsers.update({
      where: { organizationId_userId: { organizationId: orgId, userId } },
      data: { role: newRole },
    });
  }

  async removeUser(orgId: string, userId: string): Promise<void> {
    await this.db.organizationUsers.delete({
      where: { organizationId_userId: { organizationId: orgId, userId } },
    });
  }

  async getUserOrganization(userId: string): Promise<Organization | null> {
    const membership = await this.db.organizationUsers.findFirst({
      where: { userId },
      include: { organization: true },
    });

    return membership?.organization || null;
  }
}
```

### Step 6: Session Management
```typescript
// src/auth/session.ts
import session from 'express-session';
import RedisStore from 'connect-redis';

export function configureSession(redis: Redis): session.SessionOptions {
  return {
    store: new RedisStore({ client: redis }),
    secret: process.env.SESSION_SECRET!,
    name: 'clinical.session',
    resave: false,
    saveUninitialized: false,
    cookie: {
      secure: process.env.NODE_ENV === 'production',
      httpOnly: true,
      maxAge: 8 * 60 * 60 * 1000, // 8 hours (typical shift)
      sameSite: 'strict',
    },
    rolling: true, // Extend session on activity
  };
}

// Session timeout warning
export function sessionTimeoutMiddleware(warningMinutes: number = 15) {
  return (req: Request, res: Response, next: NextFunction) => {
    if (req.session?.cookie?.maxAge) {
      const remainingMs = req.session.cookie.maxAge;
      const warningMs = warningMinutes * 60 * 1000;

      if (remainingMs < warningMs) {
        res.set('X-Session-Warning', `Session expires in ${Math.round(remainingMs / 60000)} minutes`);
      }
    }
    next();
  };
}

// Force re-authentication for sensitive operations
export function requireRecentAuth(maxAgeMinutes: number = 15) {
  return async (req: AuthenticatedRequest, res: Response, next: NextFunction) => {
    const lastAuth = req.session?.lastAuthTime;

    if (!lastAuth) {
      return res.status(401).json({ error: 'Re-authentication required' });
    }

    const ageMs = Date.now() - new Date(lastAuth).getTime();
    const maxAgeMs = maxAgeMinutes * 60 * 1000;

    if (ageMs > maxAgeMs) {
      return res.status(401).json({
        error: 'Re-authentication required',
        reason: 'Session too old for sensitive operation',
      });
    }

    next();
  };
}
```

## Output
- Role definitions with clinical permissions
- SAML/OIDC SSO integration
- Permission middleware
- Organization management
- Secure session handling

## RBAC Checklist
- [ ] Roles defined matching clinical workflow
- [ ] IdP group mapping configured
- [ ] SSO integration tested
- [ ] Permission middleware on all routes
- [ ] Audit logging for access denied
- [ ] Session timeout configured
- [ ] MFA enforced for admin roles

## Error Handling
| RBAC Issue | Detection | Resolution |
|------------|-----------|------------|
| SSO login fails | Auth callback error | Check IdP configuration |
| Wrong role assigned | User reports | Review group mappings |
| Permission denied | 403 responses | Check role permissions |
| Session expired | User redirect | Implement session warning |

## Resources
- [SAML 2.0 Specification](https://wiki.oasis-open.org/security/FrontPage)
- [OpenID Connect](https://openid.net/connect/)
- [HIPAA Access Control](https://www.hhs.gov/hipaa/for-professionals/security/guidance/access-control/index.html)

## Next Steps
For EHR integration migrations, see `openevidence-migration-deep-dive`.
