# MicroSim Type Matcher Skill

Run the skill-generator skill to create a new skill called `microsim-type-matcher`.
It will take as input a single learning objective document that is created by the chapter-content-generator.  These descriptions are currently stored in the `details` XML element withing chapter content.

Given this input document, this skill will create a ranked-list of the microsims generators in the the intelligent-textbook microsim generator skill list.  It will also give the pros and cons of each microsim generator for this learning objective.  It will prompt the user with this list, and the user can tell the skill to proceed to build the microsim with a given microsim generator.

Currently the list of MicroSims includes the following:

1. bubble-chart-generator - 2D bubble chart using chart.js
2. infographic-generator-p5 - allows users to hover over regions an click through
3. microsim-p5 - baseline microsim using p5.js - very general for interactive simulations - ideal for 2D bouncing balls, geometry, physics
4. mermaid-generator - ideal for workflow diagrams
5. causal-loop-generator - ideal for generating casual loop drawings in systems thinking
6. timeline-generator - perfect for timelines with hover tooltips and click-throughs
7. map-movement-generator - ideal for showing geographic movement such as migrations of populations
8. circuit-simulator - ideal for showing current moving through wires in a circuit
9. venn-diagram-generator - good for venn diagrams
10. complementary-products - good for displaying two products on a yin-yang diagram
11. healthcare-graph - network graphs with healthcare icons
12. it-management-graph - network graph with IT icons

The systematic process for matching learning objectives with appropriate MicroSim types. This is a great pedagogical design challenge that can significantly improve learning outcomes.

## Framework for Matching Learning Objectives to MicroSim Types

### Step 1: Analyze the Learning Objective

First, classify the learning objective along these dimensions:

**Cognitive Level (Bloom's Taxonomy):**

- **Remember/Understand**: Static visualizations work well
- **Apply/Analyze**: Interactive simulations needed  
- **Evaluate/Create**: Complex, manipulable systems required

**Content Type:**

- **Temporal**: Events over time → Timeline, Animation
- **Spatial**: Geographic or physical relationships → Maps, Drawings
- **Process**: Sequential steps → Workflow, Animation
- **System**: Interconnected components → Graph-network, Causal-loop
- **Quantitative**: Data relationships → Charts, Graphs
- **Conceptual**: Abstract relationships → Infographic, Diagram

### Step 2: Match Patterns to MicroSim Types

Here's a decision matrix based on common learning objective patterns:

**"Understand how X changes over time"**

- Primary: Timeline (discrete events)
- Alternative: Animation (continuous change)

**"Analyze the relationship between X and Y"**

- Primary: Chart (quantitative)
- Alternative: Causal-loop diagram (systemic)

**"Explain the process of X"**

- Primary: Workflow (linear process)
- Alternative: Animation (dynamic process)

**"Compare and contrast X and Y"**

- Primary: Infographic (side-by-side)
- Alternative: Chart (quantitative comparison)

**"Identify the causes and effects of X"**

- Primary: Causal-loop diagram
- Alternative: Graph-network (complex causation)

**"Demonstrate how X works"**

- Primary: Simulation (interactive)
- Alternative: Animation (passive viewing)

**"Map the spread/movement of X"**

- Primary: Map with motion arrows
- Alternative: Timeline (if geographic component is secondary)

### Step 3: Decision Algorithm

Here's a practical decision flow:

1. **Extract key verbs** from the learning objective (analyze, compare, demonstrate, etc.)

2. **Identify the subject matter**:
   - Is it about change? (temporal dimension)
   - Is it about relationships? (connections)
   - Is it about location? (spatial)
   - Is it about quantity? (numerical)

3. **Determine interactivity needs**:
   - Does the learner need to manipulate variables? → Simulation
   - Does the learner need to observe dynamics? → Animation
   - Is static representation sufficient? → Infographic/Chart

### Step 4: Validation Questions

Before finalizing the MicroSim type, ask:

- **Complexity**: Can this type handle the complexity of the concept?
- **Engagement**: Will this type maintain learner attention?
- **Accessibility**: Is this type appropriate for the target audience?
- **Technical feasibility**: Can this be reasonably implemented?

## Example Mappings

**Learning Objective**: "Students will understand the water cycle"

- Analysis: Process-based, cyclical, involves transformation
- Recommendation: Animation (primary) or Causal-loop diagram (alternative)

**Learning Objective**: "Students will analyze population growth trends"

- Analysis: Quantitative, temporal, comparative
- Recommendation: Chart (line graph) with time series

**Learning Objective**: "Students will trace the spread of the Renaissance"

- Analysis: Spatial, temporal, movement-based
- Recommendation: Map with motion arrows

## Implementation Suggestion

Consider creating a scoring rubric where each MicroSim type gets points based on:

- Content type fit (0-3 points)
- Cognitive level match (0-3 points)
- Interactivity appropriateness (0-3 points)
- Implementation complexity (0-3 points, inverse scoring)

The highest scoring MicroSim type becomes your primary recommendation, with the second highest as an alternative.

Would you like me to help you develop this into a more detailed rubric or create a specific tool for automating these recommendations?