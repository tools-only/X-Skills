# Shopify Product GraphQL Mutations

Key mutations for product management via the Admin API.

## Product CRUD

### productCreate

```graphql
mutation productCreate($product: ProductCreateInput!) {
  productCreate(product: $product) {
    product {
      id
      title
      handle
      status
      variants(first: 100) {
        edges {
          node { id title price sku inventoryQuantity }
        }
      }
    }
    userErrors { field message }
  }
}
```

Variables:
```json
{
  "product": {
    "title": "Product Name",
    "descriptionHtml": "<p>Description</p>",
    "vendor": "Brand",
    "productType": "Category",
    "tags": ["tag1", "tag2"],
    "status": "DRAFT",
    "options": ["Size", "Colour"],
    "variants": [
      {
        "optionValues": [
          {"optionName": "Size", "name": "S"},
          {"optionName": "Colour", "name": "Red"}
        ],
        "price": "29.95",
        "sku": "PROD-S-RED",
        "inventoryPolicy": "DENY",
        "inventoryItem": {
          "tracked": true
        }
      }
    ],
    "seo": {
      "title": "SEO Title",
      "description": "Meta description"
    }
  }
}
```

### productUpdate

```graphql
mutation productUpdate($input: ProductInput!) {
  productUpdate(input: $input) {
    product { id title }
    userErrors { field message }
  }
}
```

Variables include `id` (required) plus any fields to update.

### productDelete

```graphql
mutation productDelete($input: ProductDeleteInput!) {
  productDelete(input: $input) {
    deletedProductId
    userErrors { field message }
  }
}
```

## Variant Management

### productVariantsBulkCreate

For adding variants to an existing product:

```graphql
mutation productVariantsBulkCreate($productId: ID!, $variants: [ProductVariantsBulkInput!]!) {
  productVariantsBulkCreate(productId: $productId, variants: $variants) {
    productVariants { id title price }
    userErrors { field message }
  }
}
```

### productVariantsBulkUpdate

```graphql
mutation productVariantsBulkUpdate($productId: ID!, $variants: [ProductVariantsBulkInput!]!) {
  productVariantsBulkUpdate(productId: $productId, variants: $variants) {
    productVariants { id title price }
    userErrors { field message }
  }
}
```

## Image Upload

### stagedUploadsCreate

```graphql
mutation stagedUploadsCreate($input: [StagedUploadInput!]!) {
  stagedUploadsCreate(input: $input) {
    stagedTargets {
      url
      resourceUrl
      parameters { name value }
    }
    userErrors { field message }
  }
}
```

Input per file:
```json
{
  "filename": "image.jpg",
  "mimeType": "image/jpeg",
  "httpMethod": "POST",
  "resource": "IMAGE"
}
```

After uploading to the staged URL, attach with `productCreateMedia`.

### productCreateMedia

```graphql
mutation productCreateMedia($productId: ID!, $media: [CreateMediaInput!]!) {
  productCreateMedia(productId: $productId, media: $media) {
    media { alt status }
    mediaUserErrors { field message }
  }
}
```

## Collection Assignment

### collectionAddProducts

```graphql
mutation collectionAddProducts($id: ID!, $productIds: [ID!]!) {
  collectionAddProducts(id: $id, productIds: $productIds) {
    collection { title productsCount }
    userErrors { field message }
  }
}
```

## Inventory

### inventorySetQuantities

```graphql
mutation inventorySetQuantities($input: InventorySetQuantitiesInput!) {
  inventorySetQuantities(input: $input) {
    inventoryAdjustmentGroup { reason }
    userErrors { field message }
  }
}
```

Input:
```json
{
  "reason": "correction",
  "name": "available",
  "quantities": [{
    "inventoryItemId": "gid://shopify/InventoryItem/123",
    "locationId": "gid://shopify/Location/456",
    "quantity": 50
  }]
}
```

## Useful Queries

### Get all products

```graphql
{
  products(first: 50) {
    edges {
      node {
        id title handle status productType vendor
        variants(first: 10) {
          edges { node { id title price sku inventoryQuantity } }
        }
      }
    }
    pageInfo { hasNextPage endCursor }
  }
}
```

### Get all collections

```graphql
{
  collections(first: 50) {
    edges {
      node { id title handle productsCount }
    }
  }
}
```

### Get inventory locations

```graphql
{
  locations(first: 10) {
    edges {
      node { id name isActive }
    }
  }
}
```
