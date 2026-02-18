# Recipes

Complete automation recipes ready to adapt.

## Auto-Archive Completed Rows

Move rows with "Complete" status to an Archive sheet. Processes bottom-up to avoid shifting row indices.

```javascript
function archiveCompleted() {
  const ss = SpreadsheetApp.getActiveSpreadsheet();
  const source = ss.getSheetByName('Active');
  const archive = ss.getSheetByName('Archive');
  const data = source.getDataRange().getValues();
  const statusCol = 4; // column E (0-indexed)

  // Process bottom-up to avoid shifting rows
  for (let i = data.length - 1; i >= 1; i--) {
    if (data[i][statusCol] === 'Complete') {
      archive.appendRow(data[i]);
      source.deleteRow(i + 1); // +1 for 1-indexed rows
    }
  }
  SpreadsheetApp.flush();
}
```

## Duplicate Detection and Highlighting

Scan a column and highlight duplicate values in light red.

```javascript
function highlightDuplicates() {
  const sheet = SpreadsheetApp.getActiveSheet();
  const data = sheet.getRange('A2:A' + sheet.getLastRow()).getValues();
  const seen = {};

  for (let i = 0; i < data.length; i++) {
    const val = data[i][0];
    if (val === '') continue;
    if (seen[val]) {
      sheet.getRange(i + 2, 1).setBackground('#f4cccc'); // light red
      sheet.getRange(seen[val], 1).setBackground('#f4cccc');
    } else {
      seen[val] = i + 2; // store row number
    }
  }
}
```

## Auto-Numbering Rows

Automatically number column A when column B is edited.

```javascript
function onEdit(e) {
  const sheet = e.source.getActiveSheet();
  if (sheet.getName() === 'Data' && e.range.getColumn() === 2) {
    if (e.value && e.value !== '') {
      const row = e.range.getRow();
      const above = sheet.getRange(2, 1, row - 1, 1).getValues()
        .filter(r => r[0] !== '').length;
      sheet.getRange(row, 1).setValue(above + 1);
    }
  }
}
```

## Summary Dashboard Generator

Create a Summary sheet aggregating data from numbered weekly tabs.

```javascript
function generateSummary() {
  const ss = SpreadsheetApp.getActiveSpreadsheet();
  const summary = ss.getSheetByName('Summary') || ss.insertSheet('Summary');
  summary.clear();

  summary.getRange('A1').setValue('Weekly Summary')
    .setFontSize(16).setFontWeight('bold');
  summary.getRange('A3:D3')
    .setValues([['Week', 'Total Jobs', 'Revenue', 'Avg Price']])
    .setFontWeight('bold')
    .setBackground('#4a86c8')
    .setFontColor('#ffffff');

  let row = 4;
  for (let w = 1; w <= 52; w++) {
    const tabName = String(w).padStart(2, '0');
    const sheet = ss.getSheetByName(tabName);
    if (!sheet) continue;

    const jobs = sheet.getRange('D1').getDisplayValue();
    const revenue = sheet.getRange('E55').getDisplayValue();
    const avg = jobs > 0
      ? (parseFloat(revenue.replace(/[^0-9.]/g, '')) / jobs)
      : 0;

    summary.getRange(row, 1, 1, 4).setValues([[
      'Week ' + tabName, jobs, revenue, '$' + avg.toFixed(0)
    ]]);
    row++;
  }

  summary.autoResizeColumns(1, 4);
  SpreadsheetApp.flush();
}
```

## Batch Email Sender

Send personalised emails to a list of recipients from sheet data.

```javascript
function sendBatchEmails() {
  const sheet = SpreadsheetApp.getActiveSpreadsheet().getSheetByName('Recipients');
  const data = sheet.getRange('A2:C' + sheet.getLastRow()).getValues();
  // Columns: A=Email, B=Name, C=Status

  const remaining = MailApp.getRemainingDailyQuota();
  if (remaining < data.length) {
    SpreadsheetApp.getUi().alert(
      'Only ' + remaining + ' emails left today. Need ' + data.length + '.');
    return;
  }

  let sent = 0;
  for (let i = 0; i < data.length; i++) {
    const [email, name, status] = data[i];
    if (!email || status === 'Sent') continue;

    try {
      MailApp.sendEmail({
        to: email,
        subject: 'Your Weekly Update',
        htmlBody: '<p>Hi ' + name + ',</p><p>Here is your update...</p>'
      });
      sheet.getRange(i + 2, 3).setValue('Sent');
      sent++;
    } catch (e) {
      sheet.getRange(i + 2, 3).setValue('Error: ' + e.message);
    }
  }

  SpreadsheetApp.flush();
  SpreadsheetApp.getActiveSpreadsheet().toast('Sent ' + sent + ' emails', 'Done', 5);
}
```
