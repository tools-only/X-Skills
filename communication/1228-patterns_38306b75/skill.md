# Common Patterns

Code examples for Google Apps Script patterns. Load when a specific pattern is needed.

## Custom Menus

```javascript
function onOpen() {
  const ui = SpreadsheetApp.getUi();
  ui.createMenu('Reports')
    .addItem('Generate Weekly Report', 'generateReport')
    .addItem('Email Report to Team', 'emailReport')
    .addSeparator()
    .addSubMenu(ui.createMenu('Settings')
      .addItem('Configure Recipients', 'configureRecipients'))
    .addToUi();
}
```

- `onOpen()` runs every time the sheet is opened/refreshed
- Menu items reference function names as strings — must be public (no underscore)
- Emojis work in menu titles
- Each menu item requires a unique, named, public function

## Toast Notifications

Quick, non-blocking messages:

```javascript
SpreadsheetApp.getActiveSpreadsheet().toast('Operation complete!', 'Title', 5);
// Arguments: message, title, duration in seconds (-1 = until dismissed)
```

## Alert Dialogs

```javascript
const ui = SpreadsheetApp.getUi();

// Simple alert
ui.alert('Operation complete!');

// Yes/No confirmation
const response = ui.alert('Delete this data?', 'This cannot be undone.',
  ui.ButtonSet.YES_NO);
if (response === ui.Button.YES) {
  // proceed
}

// Prompt for input
const result = ui.prompt('Enter your name:', ui.ButtonSet.OK_CANCEL);
if (result.getSelectedButton() === ui.Button.OK) {
  const name = result.getResponseText();
}
```

## Sidebar Apps

HTML panel on the right side of the sheet. Use `google.script.run` to call server functions.

```javascript
function showSidebar() {
  const html = HtmlService.createHtmlOutput(`
    <h3>Quick Entry</h3>
    <select id="worker">
      <option>Craig</option>
      <option>Steve</option>
    </select>
    <input id="suburb" placeholder="Suburb">
    <button onclick="submit()">Add Job</button>
    <script>
      function submit() {
        const worker = document.getElementById('worker').value;
        const suburb = document.getElementById('suburb').value;
        google.script.run
          .withSuccessHandler(function() { alert('Added!'); })
          .addJob(worker, suburb);
      }
    </script>
  `)
  .setTitle('Job Entry')
  .setWidth(300);

  SpreadsheetApp.getUi().showSidebar(html);
}

// MUST be public (no underscore) for sidebar to call it
function addJob(worker, suburb) {
  const sheet = SpreadsheetApp.getActiveSpreadsheet().getActiveSheet();
  sheet.appendRow([new Date(), worker, suburb]);
}
```

## Triggers

### onEdit (Simple Trigger)

Limited permissions but no auth needed:

```javascript
function onEdit(e) {
  const sheet = e.source.getActiveSheet();
  const range = e.range;

  // Only act on specific sheet
  if (sheet.getName() !== 'Data') return;

  // Only act on specific column (C = 3)
  if (range.getColumn() !== 3) return;

  // Auto-timestamp when column C is edited
  sheet.getRange(range.getRow(), 4).setValue(new Date());
}
```

### Installable Triggers (Full Permissions)

Create via script — run the setup function once manually:

```javascript
function createTriggers() {
  // Time-driven: run every day at 8am
  ScriptApp.newTrigger('dailyReport')
    .timeBased()
    .atHour(8)
    .everyDays(1)
    .create();

  // On edit with full permissions (can send email, fetch URLs)
  ScriptApp.newTrigger('onEditFull')
    .forSpreadsheet(SpreadsheetApp.getActive())
    .onEdit()
    .create();

  // On form submit
  ScriptApp.newTrigger('onFormSubmit')
    .forSpreadsheet(SpreadsheetApp.getActive())
    .onFormSubmit()
    .create();
}
```

## Email from Sheets

```javascript
function emailWeeklySchedule() {
  const ss = SpreadsheetApp.getActiveSpreadsheet();
  const sheet = ss.getActiveSheet();

  // Build HTML email body from sheet data
  const data = sheet.getRange('A2:E10').getDisplayValues();
  let body = '<h2>Weekly Schedule</h2><table border="1" cellpadding="8">';
  body += '<tr><th>Job</th><th>Suburb</th><th>Time</th><th>Price</th></tr>';

  for (const row of data) {
    if (row[0]) {  // skip empty rows
      body += '<tr>' + row.map(cell => '<td>' + cell + '</td>').join('') + '</tr>';
    }
  }
  body += '</table>';

  MailApp.sendEmail({
    to: 'worker@example.com',
    subject: 'Your Schedule - Week ' + sheet.getName(),
    htmlBody: body
  });
}

// Check remaining email quota
function checkQuota() {
  const remaining = MailApp.getRemainingDailyQuota();
  Logger.log('Emails remaining today: ' + remaining);
}
```

## PDF Export

Non-obvious URL construction — the export parameters are undocumented:

```javascript
function exportSheetAsPdf() {
  const ss = SpreadsheetApp.getActiveSpreadsheet();
  const sheet = ss.getActiveSheet();

  const url = ss.getUrl()
    .replace(/\/edit.*$/, '')
    + '/export?exportFormat=pdf'
    + '&format=pdf'
    + '&size=A4'
    + '&portrait=true'
    + '&fitw=true'           // fit to width
    + '&sheetnames=false'
    + '&printtitle=false'
    + '&gridlines=false'
    + '&gid=' + sheet.getSheetId();

  const token = ScriptApp.getOAuthToken();
  const response = UrlFetchApp.fetch(url, {
    headers: { 'Authorization': 'Bearer ' + token }
  });

  // Email as attachment
  MailApp.sendEmail({
    to: 'boss@example.com',
    subject: 'Weekly Report PDF',
    body: 'Please find attached the weekly report.',
    attachments: [response.getBlob().setName('report.pdf')]
  });
}
```

## Row/Column Show/Hide

```javascript
sheet.hideRows(startRow, numRows);
sheet.showRows(startRow, numRows);
const isHidden = sheet.isRowHiddenByUser(rowNumber);

// Same for columns
sheet.hideColumns(startCol, numCols);
sheet.showColumns(startCol, numCols);
```

## Formatting and Colours

```javascript
// Background and text colour
sheet.getRange('A1:Z1').setBackground('#9fc5e8');
sheet.getRange('A1').setFontColor('#ffffff');

// Font styles
sheet.getRange('A1:D1').setFontWeight('bold');
sheet.getRange('A1:D1').setFontSize(12);

// Borders
sheet.getRange('A1:D10').setBorder(true, true, true, true, true, true);

// Number format
sheet.getRange('D2:D100').setNumberFormat('$#,##0.00');

// Alignment and merge
sheet.getRange('A1:D1').setHorizontalAlignment('center');
sheet.getRange('A1:D1').merge();

// Conditional formatting
const rule = SpreadsheetApp.newConditionalFormatRule()
  .whenNumberGreaterThan(1000)
  .setBackground('#b6d7a8')
  .setRanges([sheet.getRange('D2:D100')])
  .build();
sheet.setConditionalFormatRules([rule]);
```

## Data Protection

```javascript
// Protect a range
const protection = sheet.getRange('A1:D10')
  .protect()
  .setDescription('Header row - do not edit');

// Allow specific editors
protection.addEditor('admin@example.com');
protection.removeEditors(protection.getEditors());

// Protect entire sheet, unprotect specific ranges
const sheetProtection = sheet.protect().setDescription('Protected Sheet');
sheetProtection.setUnprotectedRanges([
  sheet.getRange('B3:B50'),
  sheet.getRange('C3:C50')
]);
```

## Data Validation Dropdowns

```javascript
// Dropdown from list
const rule = SpreadsheetApp.newDataValidation()
  .requireValueInList(['Option A', 'Option B', 'Option C'], true)
  .setAllowInvalid(false)
  .setHelpText('Select an option')
  .build();
sheet.getRange('C3:C50').setDataValidation(rule);

// Dropdown from range
const rule2 = SpreadsheetApp.newDataValidation()
  .requireValueInRange(ss.getSheetByName('Lookups').getRange('A1:A100'))
  .build();
sheet.getRange('B3:B50').setDataValidation(rule2);
```

## Properties Service (Persistent Storage)

```javascript
// Script-level (shared by all users)
const scriptProps = PropertiesService.getScriptProperties();
scriptProps.setProperty('lastRun', new Date().toISOString());
const lastRun = scriptProps.getProperty('lastRun');

// User-level (per user)
const userProps = PropertiesService.getUserProperties();

// Document-level (per spreadsheet)
const docProps = PropertiesService.getDocumentProperties();
```

## Working with Multiple Sheets

```javascript
// Get specific sheet by name
const sheet = ss.getSheetByName('Week 01');

// Loop through numbered tabs
for (let i = 1; i <= 52; i++) {
  const tabName = String(i).padStart(2, '0');
  const sheet = ss.getSheetByName(tabName);
  if (sheet) {
    // process each tab
  }
}

// Copy and create sheets
const newSheet = sheet.copyTo(ss);
newSheet.setName('New Sheet Name');
const created = ss.insertSheet('My New Tab');
```

## External API Calls

```javascript
// GET request
function fetchData() {
  const response = UrlFetchApp.fetch('https://api.example.com/data', {
    headers: { 'Authorization': 'Bearer ' + getApiKey() }
  });
  return JSON.parse(response.getContentText());
}

// POST request
function postData(payload) {
  const response = UrlFetchApp.fetch('https://api.example.com/submit', {
    method: 'post',
    contentType: 'application/json',
    payload: JSON.stringify(payload),
    muteHttpExceptions: true  // don't throw on 4xx/5xx
  });

  if (response.getResponseCode() !== 200) {
    throw new Error('API error: ' + response.getContentText());
  }
  return JSON.parse(response.getContentText());
}

// Post to Google Chat webhook
function notifyChat(message) {
  UrlFetchApp.fetch('https://chat.googleapis.com/v1/spaces/XXX/messages?key=YYY', {
    method: 'post',
    contentType: 'application/json',
    payload: JSON.stringify({ text: message })
  });
}
```
