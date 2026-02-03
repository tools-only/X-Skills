# Best Practices

| Property | Value |
|----------|-------|
| **Name** | Best Practices |
| **Repository** | [K-Dense-AI/claude-scientific-skills](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-skills/main/scientific-skills/pytorch-lightning/references/best_practices.md) (🔥 7.7k) |
| **Original Path** | `scientific-skills/pytorch-lightning/references/best_practices.md` |
| **Category** | productivity |
| **Subcategory** | optimization |
| **Tags** | productivity |
| **Created** | 2025-11-14 |
| **Updated** | 2025-11-14 |
| **File Hash** | `dddba5ae2af7c9ea...` |

## Description

Good:
python
class MyModel(L.LightningModule):
     Research code (what the model does)
    def training_step(self, batch, batch_idx):
        loss = self.compute_loss(batch)
        return loss

**Tags:** `productivity`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [K-Dense-AI/claude-scientific-skills](https://raw.githubusercontent.com/K-Dense-AI/claude-scientific-skills/main/scientific-skills/pytorch-lightning/references/best_practices.md)*
