# Shopify Product CSV Import Format

## Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `Handle` | URL slug (unique per product) | `classic-tshirt` |
| `Title` | Product name (first row per product) | `Classic T-Shirt` |
| `Body (HTML)` | Description in HTML | `<p>Premium cotton</p>` |
| `Vendor` | Brand or manufacturer | `My Brand` |
| `Product Category` | Shopify standard taxonomy | `Apparel & Accessories > Clothing > Shirts & Tops` |
| `Type` | Custom product type | `T-Shirts` |
| `Tags` | Comma-separated tags | `summer, cotton, casual` |
| `Published` | Whether product is visible | `TRUE` or `FALSE` |

## Variant Columns

| Column | Description | Example |
|--------|-------------|---------|
| `Option1 Name` | First option name | `Size` |
| `Option1 Value` | First option value | `Medium` |
| `Option2 Name` | Second option name | `Colour` |
| `Option2 Value` | Second option value | `Black` |
| `Option3 Name` | Third option name | `Material` |
| `Option3 Value` | Third option value | `Cotton` |
| `Variant SKU` | Stock keeping unit | `TSHIRT-M-BLK` |
| `Variant Grams` | Weight in grams | `200` |
| `Variant Inventory Qty` | Stock quantity | `50` |
| `Variant Price` | Variant price | `29.95` |
| `Variant Compare At Price` | Original price (for sales) | `39.95` |
| `Variant Requires Shipping` | Physical product | `TRUE` |
| `Variant Taxable` | Subject to tax | `TRUE` |

## Image Columns

| Column | Description | Example |
|--------|-------------|---------|
| `Image Src` | Image URL | `https://example.com/img.jpg` |
| `Image Position` | Display order (1-based) | `1` |
| `Image Alt Text` | Alt text for accessibility | `Classic T-Shirt front view` |

## SEO Columns

| Column | Description | Example |
|--------|-------------|---------|
| `SEO Title` | Page title tag | `Classic T-Shirt | My Brand` |
| `SEO Description` | Meta description | `Premium cotton tee in 5 colours` |

## Multi-Variant Row Format

The first row has the product title and details. Subsequent rows for the same product have only the `Handle` and variant-specific columns:

```csv
Handle,Title,Body (HTML),Vendor,Type,Tags,Published,Option1 Name,Option1 Value,Variant SKU,Variant Price,Variant Inventory Qty,Image Src
classic-tshirt,Classic T-Shirt,<p>Premium cotton</p>,My Brand,T-Shirts,"summer,cotton",TRUE,Size,Small,TSH-S,29.95,50,https://example.com/tshirt.jpg
classic-tshirt,,,,,,,,Medium,TSH-M,29.95,75,
classic-tshirt,,,,,,,,Large,TSH-L,29.95,60,
```

## Notes

- UTF-8 encoding required
- Maximum 50MB file size
- Handle must be unique per product — duplicate handles update existing products
- Leave variant columns blank on variant rows for fields that don't change
- Images can be on any row — they're associated by Handle
- `Published` = `TRUE` makes the product immediately visible
