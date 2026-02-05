# Actions

| Property | Value |
|----------|-------|
| **Name** | Actions |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/react-router-v7/ACTIONS.md) (⭐ 17) |
| **Original Path** | `skills/react-router-v7/ACTIONS.md` |
| **Category** | communication |
| **Subcategory** | email |
| **Tags** | communication |
| **Created** | 2025-12-21 |
| **Updated** | 2025-12-21 |
| **File Hash** | `18fa59604d903a27...` |

## Description

tsx
{
  path: "/projects/:id",
  action: async ({ request, params }) => {
    const formData = await request.formData();
    const title = formData.get("title");
    await updateProject(params.id, { title });
    return { success: true };
  },
  Component: Project,
}

**Tags:** `communication`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/react-router-v7/ACTIONS.md)*
