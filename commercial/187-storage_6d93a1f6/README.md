# fast-ssd

| Property | Value |
|----------|-------|
| **Name** | fast-ssd |
| **Repository** | [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/storage.md) (⭐ 216) |
| **Original Path** | `skills/kubernetes-specialist/references/storage.md` |
| **Category** | commercial |
| **Subcategory** | ecommerce |
| **Tags** | commercial |
| **Created** | 2025-12-15 |
| **Updated** | 2026-01-29 |
| **File Hash** | `6d93a1f63bc7101e...` |

## Description

yaml
apiVersion: storage.k8s.io/v1
kind: StorageClass
metadata:
  name: fastssdgce
provisioner: pd.csi.storage.gke.io
parameters:
  type: pdssd
  replicationtype: regionalpd
volumeBindingMode: WaitForFirstConsumer
allowVolumeExpansion: true
reclaimPolicy: Delete

**Tags:** `commercial`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [Jeffallan/claude-skills](https://raw.githubusercontent.com/Jeffallan/claude-skills/main/skills/kubernetes-specialist/references/storage.md)*
