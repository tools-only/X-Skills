# Qemu Legacy Os

| Property | Value |
|----------|-------|
| **Name** | Qemu Legacy Os |
| **Repository** | [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/benchmarks/trajectory-only/install-windows-3.11/references/qemu_legacy_os.md) (⭐ 44) |
| **Original Path** | `letta/benchmarks/trajectory-only/install-windows-3.11/references/qemu_legacy_os.md` |
| **Category** | development |
| **Subcategory** | coding |
| **Tags** | development |
| **Created** | 2025-12-19 |
| **Updated** | 2025-12-19 |
| **File Hash** | `0c1456fb9d4e2cb9...` |

## Description

bash
qemusystemi386 \
  m 32 \
  hda /path/to/win311.img \
  vnc :1,share=ignoredisconnects \
  qmp unix:/tmp/qmp.sock,server,nowait \
  vga std \
  cpu 486 \
  noacpi \
  nohpet

**Tags:** `development`

---

*This skill is maintained by [SkillFlow](https://github.com/tools-only/SkillFlow)*
*Source: [letta-ai/skills](https://raw.githubusercontent.com/letta-ai/skills/main/letta/benchmarks/trajectory-only/install-windows-3.11/references/qemu_legacy_os.md)*
