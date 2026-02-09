---
name: lokalise-enterprise-rbac
description: |
  Configure Lokalise enterprise SSO, role-based access control, and team management.
  Use when implementing SSO integration, configuring role-based permissions,
  or setting up organization-level controls for Lokalise.
  Trigger with phrases like "lokalise SSO", "lokalise RBAC",
  "lokalise enterprise", "lokalise roles", "lokalise permissions", "lokalise team".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Lokalise Enterprise RBAC

## Overview
Configure enterprise-grade access control for Lokalise with teams, roles, and SSO.

## Prerequisites
- Lokalise Enterprise or Team plan
- Identity Provider (IdP) for SSO (optional)
- Understanding of role-based access patterns
- Admin access to Lokalise organization

## Lokalise Role Hierarchy

| Role | Scope | Permissions |
|------|-------|-------------|
| Owner | Organization | Full control, billing, delete org |
| Admin | Organization | Manage teams, projects, users |
| Manager | Team/Project | Manage project settings, contributors |
| Developer | Project | Read/write keys, translations |
| Translator | Project | Read/write translations only |
| Reviewer | Project | Review and approve translations |
| Viewer | Project | Read-only access |

## Instructions

### Step 1: Define Role Mappings
```typescript
// Role definitions aligned with Lokalise
enum LokaliseRole {
  Owner = "owner",
  Admin = "admin",
  Manager = "manager",
  Developer = "developer",
  Translator = "translator",
  Reviewer = "reviewer",
  Viewer = "viewer",
}

interface LokalisePermissions {
  manageTeam: boolean;
  manageProject: boolean;
  manageKeys: boolean;
  editTranslations: boolean;
  reviewTranslations: boolean;
  downloadFiles: boolean;
  uploadFiles: boolean;
  viewOnly: boolean;
}

const ROLE_PERMISSIONS: Record<LokaliseRole, LokalisePermissions> = {
  owner: {
    manageTeam: true,
    manageProject: true,
    manageKeys: true,
    editTranslations: true,
    reviewTranslations: true,
    downloadFiles: true,
    uploadFiles: true,
    viewOnly: false,
  },
  admin: {
    manageTeam: true,
    manageProject: true,
    manageKeys: true,
    editTranslations: true,
    reviewTranslations: true,
    downloadFiles: true,
    uploadFiles: true,
    viewOnly: false,
  },
  manager: {
    manageTeam: false,
    manageProject: true,
    manageKeys: true,
    editTranslations: true,
    reviewTranslations: true,
    downloadFiles: true,
    uploadFiles: true,
    viewOnly: false,
  },
  developer: {
    manageTeam: false,
    manageProject: false,
    manageKeys: true,
    editTranslations: true,
    reviewTranslations: false,
    downloadFiles: true,
    uploadFiles: true,
    viewOnly: false,
  },
  translator: {
    manageTeam: false,
    manageProject: false,
    manageKeys: false,
    editTranslations: true,
    reviewTranslations: false,
    downloadFiles: true,
    uploadFiles: false,
    viewOnly: false,
  },
  reviewer: {
    manageTeam: false,
    manageProject: false,
    manageKeys: false,
    editTranslations: false,
    reviewTranslations: true,
    downloadFiles: true,
    uploadFiles: false,
    viewOnly: false,
  },
  viewer: {
    manageTeam: false,
    manageProject: false,
    manageKeys: false,
    editTranslations: false,
    reviewTranslations: false,
    downloadFiles: true,
    uploadFiles: false,
    viewOnly: true,
  },
};

function checkPermission(
  role: LokaliseRole,
  permission: keyof LokalisePermissions
): boolean {
  return ROLE_PERMISSIONS[role][permission];
}
```

### Step 2: Team Management
```typescript
import { LokaliseApi } from "@lokalise/node-api";

const client = new LokaliseApi({
  apiKey: process.env.LOKALISE_API_TOKEN!,
});

// List team users
async function listTeamUsers(teamId: number) {
  const users = await client.teamUsers().list({ team_id: teamId });
  return users.items.map(u => ({
    userId: u.user_id,
    email: u.email,
    fullname: u.fullname,
    role: u.role,
    createdAt: u.created_at,
  }));
}

// Add user to team
async function addTeamUser(teamId: number, email: string, role: LokaliseRole) {
  const user = await client.teamUsers().create({
    team_id: teamId,
    email,
    role,
  });

  console.log(`Added ${email} as ${role} to team ${teamId}`);
  return user;
}

// Update user role
async function updateUserRole(teamId: number, userId: number, newRole: LokaliseRole) {
  const user = await client.teamUsers().update(userId, {
    team_id: teamId,
    role: newRole,
  });

  console.log(`Updated user ${userId} to role ${newRole}`);
  return user;
}

// Remove user from team
async function removeTeamUser(teamId: number, userId: number) {
  await client.teamUsers().delete(userId, { team_id: teamId });
  console.log(`Removed user ${userId} from team ${teamId}`);
}
```

### Step 3: Project-Level Access Control
```typescript
// Add contributor to project with specific role
async function addProjectContributor(
  projectId: string,
  email: string,
  role: LokaliseRole,
  languages?: string[]  // Optional: restrict to specific languages
) {
  const params: any = {
    email,
    is_admin: role === LokaliseRole.Admin || role === LokaliseRole.Manager,
    is_reviewer: role === LokaliseRole.Reviewer,
  };

  // Language-specific access for translators
  if (languages && languages.length > 0) {
    params.languages = languages.map(lang => ({
      lang_iso: lang,
      is_writable: role === LokaliseRole.Translator,
    }));
  }

  const contributor = await client.contributors().create(projectId, params);
  console.log(`Added ${email} to project ${projectId} as ${role}`);
  return contributor;
}

// List project contributors
async function listProjectContributors(projectId: string) {
  const contributors = await client.contributors().list({
    project_id: projectId,
  });

  return contributors.items.map(c => ({
    userId: c.user_id,
    email: c.email,
    fullname: c.fullname,
    isAdmin: c.is_admin,
    isReviewer: c.is_reviewer,
    languages: c.languages,
  }));
}
```

### Step 4: Permission Middleware
```typescript
// Express middleware for permission checks
function requireLokalisePermission(permission: keyof LokalisePermissions) {
  return async (req: Request, res: Response, next: NextFunction) => {
    const user = req.user as { lokaliseRole: LokaliseRole };

    if (!checkPermission(user.lokaliseRole, permission)) {
      return res.status(403).json({
        error: "Forbidden",
        message: `Missing permission: ${permission}`,
        requiredPermission: permission,
        userRole: user.lokaliseRole,
      });
    }

    next();
  };
}

// Usage in routes
app.post("/api/translations/upload",
  requireLokalisePermission("uploadFiles"),
  uploadHandler
);

app.delete("/api/keys/:keyId",
  requireLokalisePermission("manageKeys"),
  deleteKeyHandler
);

app.get("/api/translations",
  requireLokalisePermission("downloadFiles"),
  downloadHandler
);
```

### Step 5: SSO Integration (Enterprise)
```typescript
// Map IdP groups to Lokalise roles
const IDP_GROUP_MAPPING: Record<string, LokaliseRole> = {
  "Engineering": LokaliseRole.Developer,
  "Localization-Admins": LokaliseRole.Admin,
  "Translators-ES": LokaliseRole.Translator,
  "Translators-FR": LokaliseRole.Translator,
  "QA-Team": LokaliseRole.Reviewer,
  "Product": LokaliseRole.Viewer,
};

// SAML callback handler
async function handleSamlCallback(samlResponse: any) {
  const { email, groups } = parseSamlResponse(samlResponse);

  // Determine Lokalise role from IdP groups
  let role = LokaliseRole.Viewer;  // Default
  for (const group of groups) {
    if (IDP_GROUP_MAPPING[group]) {
      const mappedRole = IDP_GROUP_MAPPING[group];
      // Use highest privilege role
      if (roleHierarchy(mappedRole) > roleHierarchy(role)) {
        role = mappedRole;
      }
    }
  }

  // Provision user in Lokalise if needed
  await ensureUserExists(email, role);

  return { email, role };
}

function roleHierarchy(role: LokaliseRole): number {
  const hierarchy: Record<LokaliseRole, number> = {
    owner: 100,
    admin: 90,
    manager: 80,
    developer: 70,
    reviewer: 50,
    translator: 40,
    viewer: 10,
  };
  return hierarchy[role];
}
```

## Output
- Role definitions implemented
- Team management APIs
- Project-level access control
- Permission middleware active

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Permission denied | Wrong role | Check role assignment |
| SSO mismatch | Group mapping wrong | Update IDP_GROUP_MAPPING |
| User not found | Not provisioned | Auto-provision on first login |
| Language access denied | Not in languages array | Update contributor languages |

## Examples

### Quick Permission Check
```typescript
// Before allowing action
if (!checkPermission(user.role, "uploadFiles")) {
  throw new ForbiddenError("You don't have permission to upload files");
}
```

### Audit User Access
```typescript
async function auditProjectAccess(projectId: string) {
  const contributors = await listProjectContributors(projectId);

  console.log(`Project ${projectId} access audit:`);
  contributors.forEach(c => {
    console.log(`  ${c.email}: ${c.isAdmin ? "Admin" : c.isReviewer ? "Reviewer" : "Contributor"}`);
    if (c.languages) {
      console.log(`    Languages: ${c.languages.map(l => l.lang_iso).join(", ")}`);
    }
  });

  return contributors;
}
```

### Bulk Role Update
```typescript
async function updateTeamRoles(
  teamId: number,
  roleUpdates: Array<{ userId: number; newRole: LokaliseRole }>
) {
  for (const update of roleUpdates) {
    await updateUserRole(teamId, update.userId, update.newRole);
    // Respect rate limits
    await new Promise(r => setTimeout(r, 200));
  }
}
```

## Resources
- [Lokalise Team Management](https://docs.lokalise.com/en/articles/1400472-team-management)
- [Lokalise Contributors](https://developers.lokalise.com/reference/list-all-contributors)
- [Lokalise Enterprise](https://lokalise.com/enterprise)

## Next Steps
For major migrations, see `lokalise-migration-deep-dive`.
