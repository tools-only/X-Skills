# Elementor Workflows Reference

## Browser Automation: Editing a Page

### Login Flow

```
1. Navigate to https://example.com/wp-admin/
2. Enter username and password
3. Click "Log In"
4. Wait for dashboard to load
```

If already logged in (Chrome MCP with user's session), skip to the editor.

### Open Elementor Editor

```
1. Navigate to: https://example.com/wp-admin/post.php?post={ID}&action=elementor
2. Wait for Elementor loading overlay to disappear (can take 5-10 seconds)
3. Editor is ready when the left sidebar shows widget panels
```

### Edit Text Content

```
1. Click on the text element in the page preview (right panel)
2. The element becomes selected (blue border)
3. The left sidebar shows the element's settings
4. Under "Content" tab, edit the text in the editor field
5. Changes appear live in the preview
6. Click "Update" (green button, bottom left) or Ctrl+S
```

### Edit Heading

```
1. Click the heading in the preview
2. Left sidebar → Content tab → "Title" field
3. Edit the text
4. Optionally adjust: HTML tag (H1-H6), alignment, link
5. Save
```

### Change Image

```
1. Click the image widget in the preview
2. Left sidebar → Content tab → click the image thumbnail
3. Media Library opens
4. Select new image or upload
5. Click "Insert Media"
6. Save
```

### Edit Button

```
1. Click the button in the preview
2. Left sidebar → Content tab:
   - Text: button label
   - Link: URL (click the link icon for options)
   - Icon: optional icon selection
3. Style tab: colours, typography, border, padding
4. Save
```

## Template Management

### Export Template

```
1. Navigate to: https://example.com/wp-admin/edit.php?post_type=elementor_library
2. Hover over the template → "Export Template"
3. Downloads as .json file
```

### Import Template

```
1. Navigate to: https://example.com/wp-admin/edit.php?post_type=elementor_library
2. Click "Import Templates" at the top
3. Choose file → upload .json
4. Template appears in the library
```

### Apply Template via WP-CLI

For duplicating Elementor data between pages:

```bash
# Get source data
SOURCE=$(wp @site post meta get {source_id} _elementor_data)
SETTINGS=$(wp @site post meta get {source_id} _elementor_page_settings)

# Apply to target
wp @site post meta update {target_id} _elementor_data "$SOURCE"
wp @site post meta update {target_id} _elementor_edit_mode "builder"
wp @site post meta update {target_id} _elementor_page_settings "$SETTINGS"

# Clear cache
wp @site elementor flush-css
```

## Safe Text Replacement via WP-CLI

### Pre-flight Checklist

1. Back up the postmeta: `wp @site post meta get {id} _elementor_data > backup.json`
2. Dry run the replacement
3. Verify the dry run matches expectations (correct number of replacements)
4. Execute
5. Flush CSS cache
6. Verify visually

### Replacement Commands

```bash
# Dry run
wp @site search-replace "Old Company Name" "New Company Name" wp_postmeta \
  --include-columns=meta_value --dry-run --precise

# Execute
wp @site search-replace "Old Company Name" "New Company Name" wp_postmeta \
  --include-columns=meta_value --precise

# Flush cache
wp @site elementor flush-css
```

### What's Safe to Replace

| Safe | Risky |
|------|-------|
| Headings text | HTML structure |
| Paragraph text | Widget IDs |
| Button text and URLs | Section/column settings |
| Image URLs (same dimensions) | Layout properties |
| Phone numbers, emails | CSS classes |
| Addresses | Element ordering |

## Elementor Global Widgets

Global widgets are shared across pages. Editing one updates all instances.

```bash
# List global widgets
wp @site post list --post_type=elementor_library --meta_key=_elementor_template_type \
  --meta_value=widget --fields=ID,post_title
```

**Caution**: Replacing text in a global widget's data affects every page that uses it.
