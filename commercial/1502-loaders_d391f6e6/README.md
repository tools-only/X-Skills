# Loaders

| Property | Value |
|----------|-------|
| **Name** | Loaders |
| **Repository** | [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/react-router-v7/LOADERS.md) (⭐ 17) |
| **Original Path** | `skills/react-router-v7/LOADERS.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-21 |
| **Updated** | 2025-12-21 |
| **File Hash** | `d391f6e6720be737...` |

## Description

tsx
{
  path: "/teams/:teamId",
  loader: async ({ params, request }) => {
    const url = new URL(request.url);
    const query = url.searchParams.get("q");
    const team = await fetchTeam(params.teamId, query);
    return { team, name: team.name };
  },
  Component: Team,
}

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [existential-birds/beagle](https://raw.githubusercontent.com/existential-birds/beagle/main/skills/react-router-v7/LOADERS.md)*
