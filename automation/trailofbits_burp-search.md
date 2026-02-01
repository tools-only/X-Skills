---
name: trailofbits:burp-search
source: https://raw.githubusercontent.com/trailofbits/skills/main/plugins/burpsuite-project-parser/commands/burp-search.md
original_path: plugins/burpsuite-project-parser/commands/burp-search.md
source_repo: trailofbits/skills
category: automation
subcategory: workflow
tags: ['automation']
collected_at: 2026-02-01T04:15:15.909142
file_hash: 464d94781ae80af4f5717f8530e0c83d173a8a0eb2da6abfe2377c0da07de0e2
---

---
name: trailofbits:burp-search
description: Searches Burp Suite project files for security analysis
argument-hint: "<burp-file> [operation]"
allowed-tools:
  - Bash
  - Read
---

# Search Burp Suite Project Files

**Arguments:** $ARGUMENTS

Parse arguments:
1. **Burp file** (required): Path to .burp project file
2. **Operation** (optional): `auditItems`, `proxyHistory.*`, `responseHeader='...'`, `responseBody='...'`

Invoke the `burpsuite-project-parser` skill with these arguments for the full workflow.
