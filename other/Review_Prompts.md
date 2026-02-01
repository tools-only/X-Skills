---
name: Review Prompts
source: https://raw.githubusercontent.com/masoncl/review-prompts/main/kernel/btrfs.md
original_path: kernel/btrfs.md
source_repo: masoncl/review-prompts
category: other
subcategory: general
tags: ['other']
collected_at: 2026-02-01T03:15:34.041035
file_hash: 8063ac650a8e16f0386fc9badb39387f9b70e22dfcab3be4a386c86f127d77e2
---

#### BT-001: Extent map field usage correctness

**Risk**: Logic error/Over-read

**Details**: Verify correct field usage:
- em->len: logical file extent size
- em->ram_bytes: uncompressed size of underlying physical extent
- em->disk_num_bytes: compressed size on disk

When expanding operations based on extent size, ensure the expansion respects intended boundaries (file extent vs physical extent)
