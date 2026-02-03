---
name: map-generator
description: This skill generates interactive maps using the Leaflet JavaScript library. Use this skill when users need to create geographic visualizations, location-based data displays, or interactive maps for educational textbooks. The skill creates complete MicroSim packages with HTML, CSS, and documentation, optimized for iframe embedding in narrow MkDocs pages with navbar and TOC.
---

# Map Generator Skill

This skill creates interactive Leaflet maps as MicroSims for intelligent textbooks built with MkDocs Material theme.

## Working Templates (REQUIRED REFERENCE)

Before generating any map MicroSim, you MUST reference the working templates in:

```
skills/microsim-generator/assets/templates/map/
├── main-template.html      # HTML shell with Leaflet CDN links
├── style.css               # Complete styling for info panel, legend, responsive
├── script.js               # Core Leaflet logic with choropleth/marker support
├── data-template.json      # Sample data structure for map features
├── index-template.md       # MkDocs documentation page template
├── metadata-template.json  # Dublin Core metadata template
└── choropleth-example.html # Complete working example (US state quality map)
```

**CRITICAL**: Use these templates as your starting point. The `choropleth-example.html` demonstrates the complete working pattern with inline styles for quick reference.

### Template File Structure for Generated MicroSims

Each map MicroSim should produce these 6 files:

```
docs/sims/[map-name]/
├── main.html       # Uses separate CSS/JS files
├── style.css       # Based on templates/map/style.css
├── script.js       # Based on templates/map/script.js
├── data.json       # Map data (features, markers, config)
├── index.md        # MkDocs documentation
└── metadata.json   # Dublin Core metadata
```

## When to Use This Skill

Use this skill when users request:

- Geographic visualizations (regions, countries, cities)
- Location-based data displays (historical sites, landmarks)
- Campus or facility maps
- Travel route visualizations
- Custom maps with markers, layers, or highlighted borders
- Educational geography content

## Workflow

### Step 1: Gather Map Requirements

Ask the user for the following information:

1. **Map Purpose**: What geographic data are we visualizing?
2. **Geographic Region**: Region name, country, city, or coordinates
3. **Markers/Points**: What locations need markers? (names, coordinates, descriptions)
4. **Map Layers**: Do you need multiple layers (satellite, terrain, street map)?
5. **Borders/Regions**: Do any borders or regions need to be highlighted?
6. **Zoom Level**: Initial zoom level (1-18, where 1 is world view, 18 is building level)
7. **Interactive Features**: Popups, custom markers, layer controls?
8. **Educational Context**: Related concepts, Bloom's taxonomy level, target audience

**Example user input**: "Create a map showing major universities in California with markers for each campus"

### Step 2: Create Directory Structure

Create the MicroSim directory:

```
docs/sims/[map-name]/
```

**Naming convention**: Use kebab-case (e.g., `california-universities`, `ancient-rome-map`, `world-capitals`)

### Step 3: Create map-data.json (Optional)

If the map includes markers or GeoJSON data, create a `map-data.json` file:

```json
{
  "center": {
    "lat": 37.7749,
    "lng": -122.4194
  },
  "zoom": 10,
  "title": "Map Title",
  "subtitle": "Map Subtitle",
  "markers": [
    {
      "lat": 37.7749,
      "lng": -122.4194,
      "title": "Location Name",
      "description": "Location description",
      "category": "category-name"
    }
  ]
}
```

**For borders/regions**: Use GeoJSON format for complex geometries.

### Step 4: Create main.html

Create `main.html` based on the template from `assets/templates/map/main-template.html`:

**Key elements**:

- Leaflet CDN links (CSS and JS)
- Map container div with id="map"
- Info panel (top-right overlay for hover details)
- Legend container (below map)
- External style.css and script.js references
- Minimal padding/margins for iframe embedding

**Replace placeholders**:

- `{{TITLE}}` - Map title
- `{{SUBTITLE}}` - Map subtitle
- `{{LEGEND_TITLE}}` - Legend title text
- `{{LEGEND_ITEMS}}` - Legend color swatches HTML

### Step 5: Create style.css

Copy and customize `style.css` from `assets/templates/map/style.css`:

**Critical requirements for iframe embedding**:

- `body { margin: 0; padding: 0; }` - No body margins
- Minimal margins throughout (2px max for headings)
- Fixed height for #map container (default: 420px)
- Responsive breakpoints for mobile
- aliceblue background (repository standard)

**Included styles**:

- `.info` - Floating info panel for hover details
- `.legend-container` - Legend below the map
- `.metric-row` - Data display rows in info panel
- `.better`/`.worse` - Status indicators with colors
- Responsive styles for mobile (max-width: 600px)

**Customization options**:

- Map height (adjust in `#map` and responsive section)
- Info panel max-width and font sizes
- Legend layout and colors
- Marker popup styling

### Step 6: Create script.js

Copy and customize `script.js` from `assets/templates/map/script.js`:

**Core functionality**:

1. Load data from data.json (map data and feature values)
2. Initialize Leaflet map with configurable center and zoom
3. Add tile layer (OpenStreetMap default)
4. Create info control panel with hover updates
5. Load GeoJSON for choropleth maps (optional)
6. Style features based on data values with color scale
7. Handle hover highlight and click-to-zoom interactions
8. Add markers with popups (optional)

**Configuration constants to modify**:

- `MAP_CONFIG` - Center coordinates, zoom levels
- `GEOJSON_URL` - GeoJSON boundary source (null for marker-only maps)
- `DATA_URL` - Local data file path
- `colorScale` - Color definitions for value ranges
- `getColor()` - Threshold values for color assignment

**Common tile layers**:

- OpenStreetMap: `https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png`
- Satellite: `https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}`
- Terrain: `https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png`

### Step 7: Create index.md

Create `index.md` based on the template from `assets/templates/map/index-template.md`:

**Structure**:

- Title and overview
- iframe embed (width="100%", height="560px" for map + legend)
- Link to fullscreen view
- About This MicroSim section
- How to Use section
- Key Concepts Demonstrated
- Data Sources
- Technical Notes
- Related Concepts

**iframe embed format**:

```markdown
<iframe src="main.html" width="100%" height="700" frameborder="0"></iframe>

[View Fullscreen](main.html){:target="_blank"}
```

### Step 8: Create metadata.json

Create `metadata.json` based on the template from `assets/templates/map/metadata-template.json`:

**Dublin Core fields**:

- title, description, subject, creator, date
- version, format, license, language
- audience, educationalLevel

**Map-specific fields**:

- `mapType`: "choropleth", "markers", "route", etc.
- `geoJsonSource`: URL or description of GeoJSON source
- `concepts`: Array of related educational concepts
- `dependencies`: ["leaflet.js"]

### Step 9: Update mkdocs.yml Navigation

Add the new map to the navigation in `mkdocs.yml`:

```yaml
nav:
  - MicroSims:
      - Introduction: sims/index.md
      - [Map Name]: sims/[map-name]/index.md
```

**Naming**: Use Title Case for navigation labels

### Step 10: Test and Validate

Perform these validation steps:

1. **Direct HTML test**: Open `docs/sims/[map-name]/main.html` in browser
   - Verify map loads correctly
   - Test zoom and pan controls
   - Click markers to verify popups
   - Test responsive behavior (resize window)

2. **MkDocs test**: Run `mkdocs serve` and navigate to the map page
   - Verify iframe embedding works
   - Check margins/padding (should be minimal)
   - Test fullscreen link
   - Verify navigation link works

3. **Browser compatibility**: Test in Chrome, Firefox, Safari

4. **Mobile test**: Verify responsive behavior on mobile devices

## Default Layout: Map with Info Panel

The default map layout includes:

1. **Title** (top center) - 18px font, 8px top margin
2. **Subtitle** (below title) - 12px font, gray color
3. **Map container** - 420px height, full width
4. **Info panel** (top-right overlay) - Shows details on hover
5. **Legend** (below map) - Horizontal color scale with labels
6. **Hover instruction** (bottom) - Italic helper text

### Info Panel Features

The info panel automatically:

- Appears on hover over map features
- Displays feature name with styled header
- Shows metric rows with value and status indicators
- Supports grade badges with color backgrounds
- Hides when mouse leaves feature

### Choropleth Color Scale

Default color scale for value-based maps:

| Value Range | Color | Hex Code |
|-------------|-------|----------|
| > 80 | Dark Green (Best) | #1a5e1a |
| 60-80 | Green | #4CAF50 |
| 40-60 | Yellow (Average) | #FFC107 |
| 20-40 | Orange | #FF9800 |
| < 20 | Dark Red (Worst) | #8B0000 |

Customize the `getColor()` function in script.js for different threshold values.

## Common Map Patterns

### Simple Marker Map

Basic map with multiple markers:

```javascript
const markers = [
  { lat: 40.7128, lng: -74.0060, title: "New York", description: "The Big Apple" },
  { lat: 34.0522, lng: -118.2437, title: "Los Angeles", description: "City of Angels" }
];

markers.forEach(marker => {
  L.marker([marker.lat, marker.lng])
    .bindPopup(`<b>${marker.title}</b><br>${marker.description}`)
    .addTo(map);
});
```

### Highlighted Region (GeoJSON)

Show a highlighted border or region:

```javascript
const region = {
  "type": "Feature",
  "geometry": {
    "type": "Polygon",
    "coordinates": [[...]]
  }
};

L.geoJSON(region, {
  style: { color: 'red', weight: 3, fillOpacity: 0.2 }
}).addTo(map);
```

### Layer Controls

Toggle between map types:

```javascript
const street = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png');
const satellite = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}');

const baseMaps = {
  "Street": street,
  "Satellite": satellite
};

L.control.layers(baseMaps).addTo(map);
street.addTo(map); // Default layer
```

### Custom Marker Icons

Use custom icons for different categories:

```javascript
const universityIcon = L.icon({
  iconUrl: 'university-icon.png',
  iconSize: [32, 32],
  iconAnchor: [16, 32],
  popupAnchor: [0, -32]
});

L.marker([lat, lng], { icon: universityIcon }).addTo(map);
```

## Educational Considerations

### Bloom's Taxonomy Alignment

- **Remember**: Identify locations on a map
- **Understand**: Explain geographic relationships
- **Apply**: Use maps for problem-solving (routes, distances)
- **Analyze**: Compare geographic patterns
- **Evaluate**: Assess geographic data quality
- **Create**: Design custom maps for specific purposes

### Accessibility

- Ensure marker popups have descriptive text
- Provide text alternatives for visual information
- Use high-contrast colors for highlighted regions
- Include keyboard navigation support

### Performance

- Limit markers to <100 for optimal performance
- Use marker clustering for large datasets:
  ```javascript
  const markers = L.markerClusterGroup();
  markers.addLayer(L.marker([lat, lng]));
  map.addLayer(markers);
  ```

## Troubleshooting

### Map not displaying

- Check that Leaflet CDN links are correct
- Verify `#map` div has a fixed height in CSS
- Ensure coordinates are in decimal format (not DMS)

### Markers not appearing

- Verify latitude/longitude order (lat first, lng second)
- Check that coordinates are within valid ranges (-90 to 90 lat, -180 to 180 lng)
- Ensure markers are added after map initialization

### Iframe not fitting properly

- Verify `body { margin: 0; padding: 0; }` in CSS
- Check iframe height in index.md (adjust if needed)
- Ensure container has minimal margins

## References

- [Leaflet Documentation](https://leafletjs.com/reference.html)
- [Leaflet Tutorials](https://leafletjs.com/examples.html)
- [OpenStreetMap Tile Usage Policy](https://operations.osmfoundation.org/policies/tiles/)
- [GeoJSON Format Specification](https://geojson.org/)

## Version History

- v1.0 (2025-01-16): Initial release with basic marker maps, layer controls, and GeoJSON support
