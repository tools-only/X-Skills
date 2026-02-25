# WP-CLI Content Management Reference

## Posts

```bash
# Create
wp @site post create --post_type=post --post_title="Title" --post_content="<p>Content</p>" --post_status=draft

# Create from file
wp @site post create ./content.html --post_title="Title" --post_status=draft

# Update
wp @site post update {id} --post_title="New Title" --post_content="<p>Updated</p>"

# List
wp @site post list --post_type=post --posts_per_page=20 --fields=ID,post_title,post_status,post_date

# Get single post
wp @site post get {id} --fields=ID,post_title,post_status,post_content

# Delete (trash)
wp @site post delete {id}

# Delete permanently
wp @site post delete {id} --force

# Search
wp @site post list --s="search term" --fields=ID,post_title
```

## Pages

Same as posts with `--post_type=page`:

```bash
wp @site post create --post_type=page --post_title="About" --post_status=publish --post_parent=0 --menu_order=5
wp @site post list --post_type=page --fields=ID,post_title,post_parent,menu_order
```

## Media

```bash
# Import from URL
wp @site media import "https://example.com/image.jpg" --title="Photo" --alt="Description"

# Import from local path (on server)
wp @site media import /path/to/image.jpg --title="Upload" --featured_image --post_id={id}

# List media
wp @site post list --post_type=attachment --fields=ID,post_title,guid

# Set featured image
wp @site post meta update {post_id} _thumbnail_id {attachment_id}

# Regenerate thumbnails
wp @site media regenerate --yes
```

## Categories

```bash
# List
wp @site term list category --fields=term_id,name,slug,count

# Create
wp @site term create category "News" --slug=news --description="Updates"

# Create child
wp @site term create category "Tech News" --parent={parent_id}

# Update
wp @site term update category {term_id} --name="Updated Name"

# Assign to post
wp @site post term add {post_id} category {slug_or_id}
```

## Tags

```bash
wp @site term list post_tag --fields=term_id,name,slug,count
wp @site term create post_tag "new-tag"
wp @site post term add {post_id} post_tag tag1 tag2 tag3
```

## Menus

```bash
# List menus
wp @site menu list --fields=term_id,name,slug

# List menu items
wp @site menu item list {menu_name} --fields=db_id,type,title,link,position

# Add page to menu
wp @site menu item add-post {menu_name} {page_id} --title="Page Title"

# Add custom link
wp @site menu item add-custom {menu_name} "External" "https://example.com"

# Add category archive
wp @site menu item add-term {menu_name} category {term_id}

# Update item position
wp @site menu item update {item_id} --position=3

# Delete item
wp @site menu item delete {item_id}
```

## Post Meta / Custom Fields

```bash
# Get all meta
wp @site post meta list {post_id} --fields=meta_key,meta_value

# Get specific meta
wp @site post meta get {post_id} meta_key

# Set meta
wp @site post meta update {post_id} meta_key "meta_value"

# Add meta (allows duplicates)
wp @site post meta add {post_id} meta_key "meta_value"

# Delete meta
wp @site post meta delete {post_id} meta_key
```

## Search and Replace

```bash
# Dry run first
wp @site search-replace "old text" "new text" --dry-run

# Execute
wp @site search-replace "old text" "new text" --precise

# Limit to specific table
wp @site search-replace "old" "new" wp_posts --precise

# Limit to specific column
wp @site search-replace "old" "new" wp_posts post_content --precise
```

## Export / Import

```bash
# Export all content
wp @site export --dir=/tmp/

# Export specific post type
wp @site export --post_type=post --dir=/tmp/

# Import
wp @site import /path/to/file.xml --authors=mapping.csv
```
