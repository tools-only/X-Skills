# Quotas, Limits, and Debugging

## Quotas

| Resource | Free Account | Google Workspace |
|----------|-------------|-----------------|
| Script runtime | 6 min / execution | 6 min / execution |
| Time-driven trigger runtime | 30 min | 30 min |
| Triggers total daily runtime | 90 min | 6 hours |
| Triggers total | 20 per user per script | 20 per user per script |
| Email recipients/day | 100 | 1,500 |
| URL Fetch calls/day | 20,000 | 100,000 |
| Properties storage | 500 KB | 500 KB |
| Custom function runtime | 30 seconds | 30 seconds |
| Simultaneous executions | 30 | 30 |
| Calendar events created/day | 5,000 | 10,000 |
| Spreadsheets created/day | 250 | 250 |

## Debugging

- **Logger.log()** / **console.log()** — view in Apps Script editor: View > Execution Log
- **Run manually** — select function in editor dropdown > click Run
- **Executions tab** — in Apps Script editor, shows all recent runs with errors and stack traces
- **Trigger failures** — check at script.google.com > My Projects > select project > Executions
- **Test on a copy** — always test scripts on a copy of the sheet before deploying

## Deployment Checklist

Before deploying to end users:

- [ ] All functions called from HTML dialogs are public (no trailing underscore)
- [ ] `SpreadsheetApp.flush()` called before returning from modifying functions
- [ ] Error handling (try/catch) around external API calls and MailApp
- [ ] Configuration constants at the top of the file
- [ ] Header comment with install instructions
- [ ] Tested on a copy of the sheet
- [ ] Considered multi-user behaviour (different permissions, different active sheet)
- [ ] Long operations use modal progress dialogs
- [ ] No hardcoded sheet names — use configuration constants
- [ ] Checked email quota before batch sends
