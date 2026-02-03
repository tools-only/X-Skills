# Pattern: Claude Desktop Server Troubleshooting

**Pattern ID:** `claude-desktop-server-debug-001`
**Category:** System Debugging
**Confidence:** High (tested, successful)
**Last Updated:** 2025-11-01
**Reuse Count:** 0 (newly created)

---

## Pattern Summary

Systematic approach to debugging Claude Desktop crashes or connection errors caused by MCP servers or Claude Extensions.

**Token Cost:** 500-800 tokens when re-explained
**With Pattern:** 100-150 tokens (reference this document)
**Savings:** 400-650 tokens per use

---

## Recognition Triggers

Use this pattern when you encounter:
- ✅ Claude Desktop crashes on startup
- ✅ "Server connection failed" error messages
- ✅ "Transport closed unexpectedly" in logs
- ✅ Slow Claude Desktop startup
- ✅ Unknown/unresponsive MCP servers or extensions

---

## Core Concept

**Key Insight:** Claude Desktop has TWO independent server systems:
1. **MCP Servers** (JSON config)
2. **Claude Extensions** (Registry + directory)

Servers can be installed in BOTH systems simultaneously, requiring cleanup in BOTH locations.

---

## Investigation Checklist

### Phase 1: Symptom Collection (2 min)
```bash
# 1. User reports error/crash
# 2. Note error message exactly
# 3. Determine timing (startup vs runtime)
# 4. List recent config changes
```

### Phase 2: MCP Server Analysis (5 min)
```bash
# Check config
cat "C:/Users/[User]/AppData/Roaming/Claude/claude_desktop_config.json"

# Check logs (parallel)
ls -lh "C:/Users/[User]/AppData/Roaming/Claude/logs/mcp-server-*.log" | sort -k5 -h

# Large logs = problematic servers
tail -100 "C:/Users/[User]/AppData/Roaming/Claude/logs/mcp-server-[name].log"
```

**Red Flags:**
- Log file > 100KB (normal is < 50KB)
- "Server transport closed unexpectedly"
- "Server disconnected"
- DeprecationWarnings
- Repeated connection attempts

### Phase 3: Extensions Analysis (5 min)
```bash
# Check registry
cat "C:/Users/[User]/AppData/Roaming/Claude/extensions-installations.json"

# Check installed extensions
ls "C:/Users/[User]/AppData/Roaming/Claude/Claude Extensions/"

# Check for suspect extension (grep case-insensitive)
grep -i "[server-name]" "C:/Users/[User]/AppData/Roaming/Claude/extensions-installations.json"
```

### Phase 4: Complete Removal (10 min)
```bash
# 1. Backup first (ALWAYS)
cp claude_desktop_config.json "claude_desktop_config.backup.$(date +%Y%m%d_%H%M%S).json"
cp extensions-installations.json "extensions.backup.$(date +%Y%m%d_%H%M%S).json"

# 2. Remove from MCP config (if present)
# Edit claude_desktop_config.json → remove server entry

# 3. Remove from Extensions registry (if present)
# Edit extensions-installations.json → remove extension entry

# 4. Delete extension directory
rm -rf "Claude Extensions/[extension-id]"

# 5. Delete extension settings
rm "Claude Extensions Settings/[extension-id].json"

# 6. Verify cleanup
grep -i "[server-name]" claude_desktop_config.json
grep -i "[server-name]" extensions-installations.json
ls "Claude Extensions/" | grep -i "[server-name]"
```

### Phase 5: Verification (2 min)
```
# Restart Claude Desktop
# Check error message is gone
# Monitor logs for new errors
```

---

## Decision Matrix: Keep vs Remove

| Criterion | Keep | Remove |
|-----------|------|--------|
| **Value** | Solves real problem | Theoretical benefit |
| **Transparency** | Actions visible/loggable | Black box operation |
| **Redundancy** | Unique capability | Duplicates existing feature |
| **Stability** | No errors in logs | Frequent crashes/warnings |
| **Overhead** | Fast startup | Slow/heavy |

**Examples:**
- ✅ **Keep:** everything (filesystem), browser (automation), code-runner
- ❌ **Remove:** memory (black box), aws-api (crashing), sequential-thinking (redundant)

---

## File Locations Reference

```
C:\Users\[User]\AppData\Roaming\Claude\
├── claude_desktop_config.json          ← MCP Servers definition
├── extensions-installations.json       ← Extensions registry
├── Claude Extensions\
│   └── [extension-id]\                 ← Extension files
├── Claude Extensions Settings\
│   └── [extension-id].json             ← Extension configuration
└── logs\
    ├── mcp.log                         ← Main MCP log
    └── mcp-server-[name].log          ← Individual server logs
```

---

## Common Pitfalls

❌ **Mistake:** Only removing from MCP config
✅ **Solution:** Check Extensions too

❌ **Mistake:** Deleting files without backup
✅ **Solution:** Always create timestamped backups

❌ **Mistake:** Incomplete removal (registry but not files)
✅ **Solution:** Remove from ALL 4 locations

❌ **Mistake:** Assuming one config file
✅ **Solution:** Remember the dual system (MCP + Extensions)

---

## Success Indicators

✅ Error message no longer appears
✅ Claude Desktop starts quickly
✅ No crash reports in Crashpad/
✅ Clean logs (no transport errors)
✅ Reduced log file sizes

---

## Estimated Time

- **First Use:** 20-25 minutes
- **With Pattern:** 5-10 minutes
- **Savings:** 10-15 minutes per incident

---

## Related Patterns

- `system-config-optimization` (general config cleanup)
- `tool-selection-framework` (choosing which tools to install)
- `log-analysis-methodology` (interpreting error logs)

---

## Version History

- **v1.0 (2025-11-01):** Initial pattern from AWS MCP Server debugging session
