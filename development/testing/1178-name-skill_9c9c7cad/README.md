# Skill

| Property | Value |
|----------|-------|
| **Name** | Skill |
| **Repository** | [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/langchain-pack/skills/langchain-local-dev-loop/SKILL.md) (⭐ 1.3k) |
| **Original Path** | `plugins/saas-packs/langchain-pack/skills/langchain-local-dev-loop/SKILL.md` |
| **Category** | development |
| **Subcategory** | testing |
| **Tags** | development |
| **Created** | 2026-01-09 |
| **Updated** | 2026-01-09 |
| **File Hash** | `9c9c7cadedad8cfb...` |

## Description

@pytest.fixture
def mock_llm():
    """Mock LLM for unit tests without API calls."""
    mock = MagicMock()
    mock.invoke.return_value = AIMessage(content="Mocked response")
    return mock

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [jeremylongshore/claude-code-plugins-plus-skills](https://raw.githubusercontent.com/jeremylongshore/claude-code-plugins-plus-skills/main/plugins/saas-packs/langchain-pack/skills/langchain-local-dev-loop/SKILL.md)*
