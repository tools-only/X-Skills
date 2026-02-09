---
name: evernote-hello-world
description: |
  Create a minimal working Evernote example.
  Use when starting a new Evernote integration, testing your setup,
  or learning basic Evernote API patterns.
  Trigger with phrases like "evernote hello world", "evernote example",
  "evernote quick start", "simple evernote code", "create first note".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Evernote Hello World

## Overview

Create your first Evernote note using the Cloud API, demonstrating ENML format and NoteStore operations.

## Prerequisites

- Completed `evernote-install-auth` setup
- Valid access token (OAuth or Developer Token for sandbox)
- Development environment ready

## Instructions

### Step 1: Create Entry File

```javascript
// hello-evernote.js
const Evernote = require('evernote');

// Initialize authenticated client
const client = new Evernote.Client({
  token: process.env.EVERNOTE_ACCESS_TOKEN,
  sandbox: true // Set to false for production
});
```

### Step 2: Understand ENML Format

Evernote uses ENML (Evernote Markup Language), a restricted subset of XHTML:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE en-note SYSTEM "http://xml.evernote.com/pub/enml2.dtd">
<en-note>
  <h1>Note Title</h1>
  <p>This is a paragraph.</p>
  <div>Content goes here</div>
</en-note>
```

**Key ENML Rules:**
- Must include XML declaration and DOCTYPE
- Root element is `<en-note>`, not `<html>` or `<body>`
- All tags must be lowercase and properly closed
- No `<script>`, `<form>`, `<iframe>`, or event handlers
- Inline styles only (no `class` or `id` attributes)

### Step 3: Create Your First Note

```javascript
async function createHelloWorldNote() {
  const noteStore = client.getNoteStore();

  // Define ENML content (required format)
  const content = `<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE en-note SYSTEM "http://xml.evernote.com/pub/enml2.dtd">
<en-note>
  <h1>Hello from Claude Code!</h1>
  <p>This is my first note created via the Evernote API.</p>
  <p>Created at: ${new Date().toISOString()}</p>
  <hr/>
  <div style="color: #666;">
    <en-todo checked="false"/> Learn Evernote API basics<br/>
    <en-todo checked="true"/> Install SDK<br/>
    <en-todo checked="true"/> Configure authentication<br/>
  </div>
</en-note>`;

  // Create note object
  const note = new Evernote.Types.Note();
  note.title = 'Hello World - Evernote API';
  note.content = content;

  // Optional: specify notebook (uses default if not set)
  // note.notebookGuid = 'your-notebook-guid';

  try {
    const createdNote = await noteStore.createNote(note);
    console.log('Note created successfully!');
    console.log('Note GUID:', createdNote.guid);
    console.log('Note Title:', createdNote.title);
    console.log('Created:', new Date(createdNote.created));
    return createdNote;
  } catch (error) {
    console.error('Failed to create note:', error);
    throw error;
  }
}

createHelloWorldNote();
```

### Step 4: Python Version

```python
from evernote.api.client import EvernoteClient
from evernote.edam.type.ttypes import Note
import os

# Initialize client
client = EvernoteClient(
    token=os.environ['EVERNOTE_ACCESS_TOKEN'],
    sandbox=True
)

note_store = client.get_note_store()

# Create note content in ENML format
content = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE en-note SYSTEM "http://xml.evernote.com/pub/enml2.dtd">
<en-note>
  <h1>Hello from Python!</h1>
  <p>This note was created using the Evernote Python SDK.</p>
</en-note>'''

# Create and submit note
note = Note()
note.title = 'Hello World - Python'
note.content = content

created_note = note_store.createNote(note)
print(f'Note created: {created_note.guid}')
```

### Step 5: List Notebooks

```javascript
async function listNotebooks() {
  const noteStore = client.getNoteStore();

  const notebooks = await noteStore.listNotebooks();
  console.log('Your notebooks:');

  notebooks.forEach(notebook => {
    console.log(`- ${notebook.name} (${notebook.guid})`);
    if (notebook.defaultNotebook) {
      console.log('  ^ Default notebook');
    }
  });

  return notebooks;
}
```

### Step 6: Retrieve a Note

```javascript
async function getNote(noteGuid) {
  const noteStore = client.getNoteStore();

  // Options: withContent, withResourcesData, withResourcesRecognition, withResourcesAlternateData
  const note = await noteStore.getNote(noteGuid, true, false, false, false);

  console.log('Title:', note.title);
  console.log('Content:', note.content);
  console.log('Tags:', note.tagGuids);

  return note;
}
```

## Complete Example

```javascript
// complete-hello.js
const Evernote = require('evernote');

async function main() {
  const client = new Evernote.Client({
    token: process.env.EVERNOTE_ACCESS_TOKEN,
    sandbox: true
  });

  const noteStore = client.getNoteStore();
  const userStore = client.getUserStore();

  // 1. Get user info
  const user = await userStore.getUser();
  console.log(`Hello, ${user.username}!`);

  // 2. List notebooks
  const notebooks = await noteStore.listNotebooks();
  console.log(`You have ${notebooks.length} notebooks`);

  // 3. Create a note
  const content = `<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE en-note SYSTEM "http://xml.evernote.com/pub/enml2.dtd">
<en-note>
  <p>Hello, Evernote!</p>
</en-note>`;

  const note = new Evernote.Types.Note();
  note.title = 'My First API Note';
  note.content = content;

  const created = await noteStore.createNote(note);
  console.log(`Created note: ${created.guid}`);

  // 4. Read it back
  const fetched = await noteStore.getNote(created.guid, true, false, false, false);
  console.log('Note content retrieved successfully');
}

main().catch(console.error);
```

## Output

- Working code file with Evernote client initialization
- Successfully created note in your Evernote account
- Console output showing:
```
Hello, yourUsername!
You have 3 notebooks
Created note: 12345678-abcd-1234-efgh-123456789012
Note content retrieved successfully
```

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| `EDAMUserException: BAD_DATA_FORMAT` | Invalid ENML content | Validate against ENML DTD |
| `EDAMNotFoundException` | Note or notebook not found | Check GUID is correct |
| `EDAMSystemException: RATE_LIMIT_REACHED` | Too many requests | Wait for `rateLimitDuration` |
| `Missing DOCTYPE` | ENML missing required header | Add XML declaration and DOCTYPE |

## ENML Quick Reference

```html
<!-- Allowed -->
<p>, <div>, <span>, <br/>, <hr/>
<h1>-<h6>, <b>, <i>, <u>, <strong>, <em>
<ul>, <ol>, <li>, <table>, <tr>, <td>
<a href="...">, <img src="...">
<en-todo checked="false"/>
<en-media type="image/png" hash="..."/>

<!-- NOT Allowed -->
<script>, <form>, <input>, <button>
<iframe>, <object>, <embed>
class="...", id="...", onclick="..."
```

## Resources

- [Creating Notes](https://dev.evernote.com/doc/articles/creating_notes.php)
- [ENML Reference](https://dev.evernote.com/doc/articles/enml.php)
- [Core Concepts](https://dev.evernote.com/doc/articles/core_concepts.php)
- [API Reference](https://dev.evernote.com/doc/reference/)

## Next Steps

Proceed to `evernote-local-dev-loop` for development workflow setup.
