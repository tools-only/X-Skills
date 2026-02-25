# WordPress REST API Endpoints

## Authentication

All write operations require authentication via Application Password:

```bash
# Base64 encode credentials
AUTH=$(echo -n "username:xxxx xxxx xxxx xxxx xxxx xxxx" | base64)

# Use in requests
curl -s https://example.com/wp-json/wp/v2/posts \
  -H "Authorization: Basic $AUTH" \
  -H "Content-Type: application/json"
```

## Posts

| Operation | Method | Endpoint |
|-----------|--------|----------|
| List posts | GET | `/wp-json/wp/v2/posts` |
| Get single post | GET | `/wp-json/wp/v2/posts/{id}` |
| Create post | POST | `/wp-json/wp/v2/posts` |
| Update post | PUT | `/wp-json/wp/v2/posts/{id}` |
| Delete post | DELETE | `/wp-json/wp/v2/posts/{id}` |

Create example:

```bash
curl -s https://example.com/wp-json/wp/v2/posts \
  -H "Authorization: Basic $AUTH" \
  -H "Content-Type: application/json" \
  -d '{
    "title": "My New Post",
    "content": "<p>Post content here.</p>",
    "status": "draft",
    "categories": [3, 5],
    "tags": [10, 12],
    "excerpt": "Brief summary",
    "featured_media": 456
  }' | jq '{id, link, status}'
```

## Pages

Same pattern as posts but at `/wp-json/wp/v2/pages`:

```bash
curl -s https://example.com/wp-json/wp/v2/pages \
  -H "Authorization: Basic $AUTH" \
  -H "Content-Type: application/json" \
  -d '{
    "title": "About Us",
    "content": "<h2>Our Story</h2><p>Content...</p>",
    "status": "publish",
    "parent": 0
  }'
```

## Media

```bash
# Upload image
curl -s https://example.com/wp-json/wp/v2/media \
  -H "Authorization: Basic $AUTH" \
  -H "Content-Disposition: attachment; filename=photo.jpg" \
  -H "Content-Type: image/jpeg" \
  --data-binary @photo.jpg | jq '{id, source_url}'
```

## Categories & Tags

| Resource | Endpoint |
|----------|----------|
| Categories | `/wp-json/wp/v2/categories` |
| Tags | `/wp-json/wp/v2/tags` |

```bash
# Create category
curl -s https://example.com/wp-json/wp/v2/categories \
  -H "Authorization: Basic $AUTH" \
  -H "Content-Type: application/json" \
  -d '{"name": "News", "slug": "news", "description": "Company updates"}'
```

## Useful Query Parameters

| Parameter | Purpose | Example |
|-----------|---------|---------|
| `per_page` | Results per page (max 100) | `?per_page=50` |
| `page` | Pagination | `?page=2` |
| `search` | Search term | `?search=keyword` |
| `status` | Filter by status | `?status=draft` |
| `categories` | Filter by category ID | `?categories=3` |
| `orderby` | Sort field | `?orderby=date` |
| `order` | Sort direction | `?order=desc` |
| `_fields` | Limit response fields | `?_fields=id,title,link` |

## Menus (requires plugin or WP 5.9+)

Navigation menus have limited REST API support. The `/wp-json/wp/v2/navigation` endpoint exists for block-based navigation in FSE themes. For classic menus, use WP-CLI or browser automation.
