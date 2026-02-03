# Anti Patterns

| Property | Value |
|----------|-------|
| **Name** | Anti Patterns |
| **Repository** | [softaworks/agent-toolkit](https://raw.githubusercontent.com/softaworks/agent-toolkit/main/skills/react-useeffect/anti-patterns.md) (⭐ 383) |
| **Original Path** | `skills/react-useeffect/anti-patterns.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2026-01-21 |
| **Updated** | 2026-01-21 |
| **File Hash** | `69146321260d6797...` |

## Description

tsx
// BAD: Extra state + Effect for derived value
function Form() {
  const [firstName, setFirstName] = useState('Taylor');
  const [lastName, setLastName] = useState('Swift');
  const [fullName, setFullName] = useState('');

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [softaworks/agent-toolkit](https://raw.githubusercontent.com/softaworks/agent-toolkit/main/skills/react-useeffect/anti-patterns.md)*
