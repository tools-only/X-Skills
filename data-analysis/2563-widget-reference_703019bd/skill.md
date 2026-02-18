# Google Chat Widget Reference

All widget types available in cardsV2 sections.

## textParagraph

Formatted text block. Supports Google Chat formatting (`*bold*`, `_italic_`, `<url|text>`).

```typescript
{
  textParagraph: {
    text: '*Status*: All systems operational\n_Last checked_: 5 minutes ago'
  }
}
```

## decoratedText

Labelled value with optional icons. Most versatile widget for key-value data.

### Basic

```typescript
{
  decoratedText: {
    topLabel: 'Environment',
    text: 'Production',
    bottomLabel: 'Last deployed 2h ago'
  }
}
```

### With Start Icon

```typescript
{
  decoratedText: {
    topLabel: 'Status',
    text: 'Healthy',
    startIcon: { knownIcon: 'STAR' }
  }
}
```

### With Custom Icon URL

```typescript
{
  decoratedText: {
    topLabel: 'GitHub',
    text: 'PR #142 merged',
    startIcon: {
      iconUrl: 'https://github.githubassets.com/favicons/favicon.svg',
      altText: 'GitHub'
    }
  }
}
```

### With Button

```typescript
{
  decoratedText: {
    topLabel: 'Alert',
    text: 'CPU at 95%',
    button: {
      text: 'View',
      onClick: { openLink: { url: 'https://monitoring.example.com' } }
    }
  }
}
```

### Clickable (Whole Widget)

```typescript
{
  decoratedText: {
    text: 'View full report',
    wrapText: true,
    onClick: { openLink: { url: 'https://reports.example.com' } }
  }
}
```

### With Wrap Text

```typescript
{
  decoratedText: {
    topLabel: 'Description',
    text: 'This is a longer description that should wrap to multiple lines instead of being truncated',
    wrapText: true
  }
}
```

## buttonList

One or more action buttons. Buttons open URLs or trigger actions.

### Single Button

```typescript
{
  buttonList: {
    buttons: [{
      text: 'Open Dashboard',
      onClick: { openLink: { url: 'https://dashboard.example.com' } }
    }]
  }
}
```

### Multiple Buttons

```typescript
{
  buttonList: {
    buttons: [
      {
        text: 'Approve',
        onClick: { openLink: { url: 'https://app.example.com/approve/123' } },
        color: { red: 0, green: 0.5, blue: 0, alpha: 1 }
      },
      {
        text: 'Reject',
        onClick: { openLink: { url: 'https://app.example.com/reject/123' } }
      }
    ]
  }
}
```

### Button with Icon

```typescript
{
  buttonList: {
    buttons: [{
      text: 'View on GitHub',
      icon: { knownIcon: 'BOOKMARK' },
      onClick: { openLink: { url: 'https://github.com/org/repo/pull/42' } }
    }]
  }
}
```

## image

Standalone image widget.

```typescript
{
  image: {
    imageUrl: 'https://example.com/chart.png',
    altText: 'Monthly usage chart'
  }
}
```

## divider

Horizontal line separator between widgets.

```typescript
{ divider: {} }
```

## Collapsible Sections

Sections can be collapsed with only the first N widgets visible:

```typescript
{
  header: 'Details',
  collapsible: true,
  uncollapsibleWidgetsCount: 2,  // Show first 2, collapse rest
  widgets: [
    // First 2 always visible
    { decoratedText: { topLabel: 'Status', text: 'Active' } },
    { decoratedText: { topLabel: 'Region', text: 'AU' } },
    // These start collapsed
    { decoratedText: { topLabel: 'Instance', text: 'prod-01' } },
    { decoratedText: { topLabel: 'Memory', text: '2.1 GB' } },
    { decoratedText: { topLabel: 'CPU', text: '45%' } }
  ]
}
```

## Combining Widgets

A typical card section combining multiple widget types:

```typescript
{
  header: 'Deployment Summary',
  widgets: [
    { textParagraph: { text: '*Production deployment complete*' } },
    { divider: {} },
    { decoratedText: { topLabel: 'Version', text: 'v2.3.1', startIcon: { knownIcon: 'STAR' } } },
    { decoratedText: { topLabel: 'Duration', text: '2m 34s', startIcon: { knownIcon: 'CLOCK' } } },
    { decoratedText: { topLabel: 'Commit', text: '`abc1234` Fix auth redirect' } },
    { divider: {} },
    { buttonList: { buttons: [
      { text: 'View Logs', onClick: { openLink: { url: 'https://logs.example.com' } } },
      { text: 'Rollback', onClick: { openLink: { url: 'https://deploy.example.com/rollback' } } }
    ]}}
  ]
}
```
