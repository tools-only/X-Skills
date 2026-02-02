# Setup langgraph dev server

| Property | Value |
|----------|-------|
| **Name** | Setup langgraph dev server |
| **Repository** | [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/libs/next_gen_ui_langgraph/LANGRAPH_DEV.md) (⭐ 12) |
| **Original Path** | `libs/next_gen_ui_langgraph/LANGRAPH_DEV.md` |
| **Category** | data-analysis |
| **Subcategory** | visualization |
| **Tags** | data analysis |
| **Created** | 2025-03-13 |
| **Updated** | 2025-03-17 |
| **File Hash** | `f4bcb49560d3656d...` |

## Description

sh
cd libs/next_gen_ui_langgraph
python3 m venv .venv
source .venv/bin/activate
pip install U r requirements.txt
pip install U langgraphcli
pip install U "langgraphcli[inmem]"
pip install forcereinstall ../../dist/next_gen_ui_agent0.0.1py3noneany.whl

**Tags:** `data analysis`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [RedHat-UX/next-gen-ui-agent](https://raw.githubusercontent.com/RedHat-UX/next-gen-ui-agent/main/libs/next_gen_ui_langgraph/LANGRAPH_DEV.md)*
