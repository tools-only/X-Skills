# Causal Loop MicroSim Generator

## Overview

The causal-loop-microsim-generator skill creates interactive Causal Loop Diagram (CLD) MicroSims using the vis-network JavaScript library for systems thinking education. Each MicroSim visualizes cause-and-effect relationships, feedback loops, and system dynamics through an interactive node-edge diagram that can be embedded via iframe in educational content.

## Purpose

This skill transforms complex system dynamics into interactive, explorable visualizations that enable students to understand feedback loops, reinforcing and balancing dynamics, and cause-and-effect relationships. CLDs are fundamental tools in systems thinking, helping learners recognize patterns in complex systems and identify intervention points.

## Key Features

- **Interactive Causal Diagrams**: Node-link diagrams showing causal relationships with polarity indicators
- **Feedback Loop Visualization**: Reinforcing (R) and Balancing (B) loop indicators with automatic placement
- **Polarity Indicators**: Green (+) for positive relationships, Red (-) for negative relationships
- **Click-to-Detail**: Interactive panels showing node, edge, and loop descriptions
- **Systems Archetypes**: Support for common patterns like "limits to growth", "fixes that fail", "shifting the burden"
- **Educational Content**: Built-in discussion questions, key insights, and learning objectives
- **MicroSim Architecture**: Standardized patterns for iframe embedding at 500px height
- **Width Responsive**: Adapts to container width with re-centering

## When to Use

Use this skill when users need to:

- Create causal loop diagrams for systems thinking courses
- Visualize feedback loops in business, ecology, or social systems
- Explain reinforcing and balancing dynamics
- Demonstrate systems archetypes (limits to growth, tragedy of the commons, etc.)
- Build interactive cause-and-effect visualizations
- Teach system dynamics concepts

## Common Trigger Phrases

- "Create a CLD showing..."
- "Visualize the feedback loops in..."
- "Build a causal loop diagram for..."
- "Show the reinforcing loop between..."
- "Create a systems thinking diagram of..."
- "Generate a causal diagram showing..."

## MicroSim Architecture

Educational MicroSims occupy the intersection of:

1. **Simplicity**: Focused scope, transparent code
2. **Accessibility**: Browser-native, universal embedding
3. **AI Generation**: Standardized patterns, prompt-compatible design

## Folder Structure

Each causal loop MicroSim contains:

```
/docs/sims/$MICROSIM_NAME/
├── index.md               # Main documentation with iframe
├── main.html              # HTML container with vis-network CDN
├── $MICROSIM_NAME.js      # JavaScript code for CLD rendering
├── data.json              # Node, edge, and loop definitions
└── style.css              # Custom CSS for layout and legend
```

## Educational Requirements Specification

Before generating, the skill gathers:

1. **MicroSim Name**: kebab-case identifier (e.g., `ai-flywheel`, `climate-feedback`)
2. **Title**: Display title for the diagram
3. **Description**: What system is being modeled
4. **Nodes**: Variables in the system with descriptions
5. **Edges**: Causal relationships with polarity (positive/negative)
6. **Loops**: Feedback loop identification (reinforcing R or balancing B)

## CLD Concepts

### Causal Relationships (Edges)

- **Positive Polarity (+)**: "More of A leads to more of B" (same direction change)
- **Negative Polarity (-)**: "More of A leads to less of B" (opposite direction change)

### Feedback Loop Types

| Type | Symbol | Behavior | Example |
|------|--------|----------|---------|
| Reinforcing | R | Exponential growth or collapse | Viral adoption, compound interest |
| Balancing | B | Goal-seeking, stability | Thermostat, predator-prey equilibrium |

### Loop Determination Rules

- **Reinforcing Loop**: All edges same polarity OR even number of negative edges
- **Balancing Loop**: Odd number of negative edges in the loop

## Systems Thinking Archetypes

The skill recognizes and can generate common systems archetypes:

| Archetype | Loops | Pattern |
|-----------|-------|---------|
| `limits-to-growth` | R + B | Growth eventually limited by constraint |
| `fixes-that-fail` | B + R | Quick fix creates side effects |
| `shifting-the-burden` | 2B + R | Symptomatic vs fundamental solution |
| `success-to-the-successful` | 2R | Winner-take-all dynamics |
| `tragedy-of-the-commons` | Multiple R + B | Shared resource depletion |
| `escalation` | 2R | Arms race pattern |
| `drifting-goals` | B + R | Standards erosion over time |

## JSON Data Schema

The `data.json` file follows a comprehensive schema:

```json
{
  "metadata": {
    "id": "example-cld",
    "title": "Example System",
    "archetype": "limits-to-growth",
    "description": "Description of the system",
    "learning_objectives": ["..."],
    "version": "1.0.0"
  },
  "nodes": [
    {
      "id": "node_1",
      "label": "Variable Name",
      "position": {"x": 300, "y": 150},
      "type": "variable",
      "description": "What this variable represents"
    }
  ],
  "edges": [
    {
      "id": "node_1_to_node_2",
      "source": "node_1",
      "target": "node_2",
      "polarity": "positive",
      "description": "How these relate"
    }
  ],
  "loops": [
    {
      "id": "main_loop",
      "type": "reinforcing",
      "path": ["node_1", "node_2", "node_1"],
      "label": "Growth Cycle",
      "position": {"x": 300, "y": 300}
    }
  ],
  "educational_content": {
    "discussion_questions": ["..."],
    "key_insights": ["..."]
  }
}
```

### Node Types

| Type | Description | Use Case |
|------|-------------|----------|
| `stock` | Accumulation variable | Resources, inventory, population |
| `variable` | Flow or rate variable | Rates, decisions, activities |
| `parameter` | External constant | Policies, fixed factors |

## Node Positioning Guidelines

The skill uses these positioning patterns for clean layouts:

### Canvas Dimensions
- Standard canvas: 600x600 pixels
- Center point: (300, 300)
- Margin: 100px from edges

### Layout Patterns

**3-Node Triangle:**
```
Node 1: (300, 100)  - Top
Node 2: (450, 400)  - Bottom Right
Node 3: (150, 400)  - Bottom Left
```

**4-Node Diamond:**
```
Node 1: (300, 150)  - Top
Node 2: (450, 300)  - Right
Node 3: (300, 450)  - Bottom
Node 4: (150, 300)  - Left
```

**5-Node Pentagon:**
```
Node 1: (300, 100)  - Top
Node 2: (470, 220)  - Top Right
Node 3: (410, 420)  - Bottom Right
Node 4: (190, 420)  - Bottom Left
Node 5: (130, 220)  - Top Left
```

## Visual Configuration

### Polarity Colors
- **Positive (+)**: Green `#28a745`
- **Negative (-)**: Red `#dc3545`

### Loop Indicators
- **Reinforcing (R)**: Red ellipse background `#dc3545`
- **Balancing (B)**: Green ellipse background `#28a745`

### Edge Curves
Curves prevent edge overlap:

```javascript
{
  "curve": {
    "type": "curvedCW",    // or "curvedCCW"
    "roundness": 0.4       // 0.1 to 0.5
  }
}
```

## Interactive Features

### Click Events
- **Node Click**: Shows variable details (type, description, examples)
- **Edge Click**: Shows relationship details (polarity, strength, delay)
- **Loop Click**: Shows loop information (type, path, behavior pattern)

### URL Parameters
- `menu=true`: Show header and details panel (default hidden for iframe)
- `file=filename`: Load a specific JSON file

## Educational Content Features

The skill generates:

### Learning Objectives (Bloom's Taxonomy)
- **Remember**: Identify loop types and polarity
- **Understand**: Explain causal relationships
- **Apply**: Predict system behavior
- **Analyze**: Identify leverage points
- **Evaluate**: Assess intervention strategies
- **Create**: Design system modifications

### Discussion Questions
Generated questions explore "what if" scenarios, challenge assumptions, and connect to real-world examples.

### Key Insights
Highlights non-obvious relationships, counterintuitive behaviors, delay effects, and unintended consequences.

## Best Practices

### Diagram Design
1. **Start Simple**: Begin with core feedback loop, add complexity gradually
2. **Clear Labels**: Use concise, descriptive variable names (max 20 characters)
3. **Meaningful Relationships**: Include only significant causal links
4. **Loop Labels**: Name loops descriptively (e.g., "Growth Engine", "Quality Control")

### Educational Value
1. **Context**: Explain what real-world system the CLD represents
2. **Guided Exploration**: Provide questions to focus investigation
3. **Intervention Points**: Identify where to apply leverage
4. **Assessment**: Include comprehension questions

### Accessibility
1. **Color and Symbol**: Use both color and +/- symbols for polarity
2. **Click Details**: Provide text descriptions for all elements
3. **Legend**: Include legend explaining visual conventions
4. **Alternative Text**: Document provides text-based alternative

## Output Files

1. **index.md**: Documentation with iframe embed, learning objectives, and explanation
2. **main.html**: HTML container with vis-network CDN and layout structure
3. **[sim-name].js**: JavaScript code for loading data and rendering the CLD
4. **data.json**: Complete CLD data with nodes, edges, loops, and educational content
5. **style.css**: CSS for layout, legend, and details panel styling

## Post-Generation Actions

After generating files, the skill reminds users to:

1. Take a screenshot and save as `[sim-name].png` in the MicroSim directory
2. Update `mkdocs.yml` navigation (entries should be alphabetically ordered)
3. Test the iframe embed in the documentation page

## Integration with Other Skills

**Primary Integrations:**

- **microsim-p5**: For more complex dynamic simulations with animation
- **vis-network**: Base library for network visualization
- **learning-graph-generator**: Can visualize concept dependencies as CLD

**Content Integrations:**

- **chapter-content-generator**: Embed CLDs in systems thinking chapters
- **quiz-generator**: Create questions about system dynamics
- **glossary-generator**: Link CLD terms to glossary definitions

## Technical Requirements

- **vis-network.js**: Loaded from CDN (unpkg.com)
- **Modern Browser**: Chrome, Firefox, Safari, Edge
- **No Server**: Runs entirely client-side
- **Responsive**: Container-based sizing

## Example Use Cases

1. **AI Adoption Flywheel**: Reinforcing loop of AI usage → data → accuracy → adoption
2. **Climate Feedback**: Multiple loops showing warming → ice melt → albedo effects
3. **Technical Debt**: Balancing pressures of speed vs. quality in software development
4. **Market Dynamics**: Supply and demand with price equilibration
5. **Learning Curve**: Skill → productivity → practice time feedback

## Troubleshooting

### Issue: Loops not displaying correctly
**Solution**: Verify loop paths form complete cycles, check polarity count for loop type

### Issue: Edges overlapping
**Solution**: Use opposite curve directions (curvedCW/curvedCCW) for parallel edges

### Issue: Nodes too close together
**Solution**: Adjust position coordinates, maintain 150px minimum spacing

### Issue: Details panel not updating
**Solution**: Verify click handlers in JavaScript, check element IDs match

## References

- **vis-network Documentation**: https://visjs.github.io/vis-network/docs/network/
- **Systems Thinking Archetypes**: Senge, P. (1990) *The Fifth Discipline*
- **Causal Loop Diagrams**: Sterman, J. (2000) *Business Dynamics*
- **System Dynamics Society**: https://systemdynamics.org/
