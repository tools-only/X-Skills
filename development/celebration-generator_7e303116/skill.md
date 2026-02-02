# Celebration Animation Generator

## Overview

The **Celebration Animation Generator** skill creates self-contained p5.js celebration animations for educational MicroSims. These animations provide visual feedback and rewards when students complete tasks correctly, enhancing engagement and motivation in interactive learning environments.

Each animation is a single JavaScript file that can be copied into any MicroSim folder to provide visual celebration feedback. The animations are designed to work with the p5.js library and follow a consistent API pattern.

## When to Use This Skill

Use this skill when you need to create:

- **Particle effects** for student rewards and achievements
- **Celebration animations** (confetti, fireworks, sparkles)
- **Visual feedback** for correct answers in educational games
- **Motion-based animations** (floating, bursting, falling objects)
- **Custom themed celebrations** (baseballs, butterflies, stars, hearts)

## Key Features

- **Self-contained modules**: Each animation is a single JavaScript file
- **Consistent API**: Four standard functions for all animations
- **Speed control**: Adjustable animation speed (slow, medium, fast)
- **Rainbow color palette**: Built-in colorful defaults
- **Multiple motion patterns**: Burst, float, fall, explode, zoom, bounce
- **p5.js integration**: Works seamlessly with p5.js MicroSims

## Animation API Pattern

All celebration animations follow this consistent API:

```javascript
// Initialize and trigger the animation
create[AnimationName](centerX, startY, speedMultiplier)

// Update physics and render particles (call in draw loop)
updateAndDraw[AnimationName]()

// Check if animation is still playing
is[AnimationName]Active()

// Stop animation immediately
clear[AnimationName]()
```

### Speed Multiplier Values

| Value | Description | Use Case |
|-------|-------------|----------|
| 0.5 | Slow | Younger students, accessibility |
| 1.0 | Medium | Default speed |
| 1.8 | Fast | Quick feedback, older students |

## Available Motion Patterns

| Pattern | Description | Example Use |
|---------|-------------|-------------|
| **Burst Up** | Objects shoot upward from a point with gravity | Book Burst, Alphabet Fireworks |
| **Float Up** | Objects gently float upward | Yellow Stars, Balloons |
| **Fall Down** | Objects fall from top | Happy Star Sprinkle, Confetti, Spark Shower |
| **Explode Out** | Objects radiate outward from center | Rainbow Sparkle Burst, Magic Book Bloom |
| **Zoom Across** | Objects move horizontally | Reading Rocket Zoom |
| **Pop/Bounce** | Objects appear, bounce, then pop | Giggle Glitter Pop |

## Standard Color Palette

All animations use the rainbow color palette for visual variety:

```javascript
const rainbowColors = [
  '#FF6B6B', // red
  '#FF8E53', // orange
  '#FFD93D', // yellow
  '#6BCB77', // green
  '#4D96FF', // blue
  '#9B59B6', // purple
  '#FF6B9D'  // pink
];
```

## Common Particle Properties

Each particle in an animation typically has these properties:

```javascript
{
  x, y,           // Position
  vx, vy,         // Velocity
  size,           // Size/radius
  alpha,          // Transparency (0-255)
  fadeRate,       // How fast alpha decreases
  rotation,       // Current angle
  rotationSpeed,  // Angular velocity
  color,          // p5.js color object
  gravity,        // Downward acceleration (burst patterns)
  wobble,         // Oscillation factor (floating patterns)
  trail: []       // Past positions (trail effects)
}
```

## Workflow

### Step 1: Parse the Animation Request

Extract from the user's description:

1. **Object/Shape**: What is being animated (baseballs, hearts, butterflies)
2. **Motion Pattern**: How it moves (exploding, floating, falling, zooming)
3. **Origin Point**: Where it starts (bottom center, top, sides, center)
4. **Suggested Name**: Derive a kebab-case filename (e.g., `baseball-explosion.js`)

### Step 2: Generate the Animation File

Create a new JavaScript file following requirements:

- Use a **unique particle array name** (e.g., `baseballExplosionParticles`)
- Use a **unique suffix for helper functions** to avoid conflicts
- Include all four standard API functions
- Support `speedMultiplier` parameter
- Use the standard rainbow color palette

### Step 3: Test with Animation Library Tester

Integrate with the animation-lib-tester MicroSim:

1. Add script import to main.html
2. Add animation to the `animationTypes` array
3. Add draw and clear function calls
4. Add trigger case in switch statement

### Step 4: Update Documentation

Add the new animation to the README:

- Add row to "Available Animations" table
- Add API documentation with function signatures

## File Naming Convention

- **Filename**: kebab-case (`baseball-explosion.js`)
- **Function names**: PascalCase (`createBaseballExplosion()`)
- **Particle array**: camelCase (`baseballExplosionParticles`)
- **Helper functions**: Add unique suffix (`drawBaseballBE()`)

## Example: Creating "Baseball Explosion"

For request: "Baseballs exploding from the bottom middle of the screen"

1. **Filename**: `baseball-explosion.js`
2. **Motion**: Burst Up pattern (like Book Burst)
3. **Object**: Baseball with red stitching
4. **Functions**:
   - `createBaseballExplosion(centerX, startY, speedMultiplier)`
   - `updateAndDrawBaseballExplosion()`
   - `isBaseballExplosionActive()`
   - `clearBaseballExplosion()`

## Integration with MicroSims

To use a celebration animation in your MicroSim:

```html
<!-- In main.html, add the animation script -->
<script src="../shared/animations/baseball-explosion.js"></script>
```

```javascript
// In your MicroSim's sketch.js
function checkAnswer() {
  if (isCorrect) {
    // Trigger celebration from bottom center
    createBaseballExplosion(width/2, height - 50, 1.0);
  }
}

function draw() {
  // ... your drawing code ...

  // Update and render the celebration
  updateAndDrawBaseballExplosion();
}
```

## Educational Benefits

Celebration animations support learning by:

- **Positive reinforcement**: Immediate visual reward for correct answers
- **Engagement**: Colorful effects maintain student interest
- **Motivation**: Celebrations encourage continued participation
- **Dopamine response**: Visual rewards trigger positive associations
- **Accessibility**: Speed control accommodates different needs

## Technical Details

- **Library**: p5.js (loaded via CDN)
- **File Location**: `/docs/sims/shared/animations/`
- **Testing Tool**: `/docs/sims/animation-lib-tester/`
- **Particle Count**: Typically 20-50 particles per animation
- **Duration**: Most animations complete in 2-4 seconds
- **Performance**: Optimized for smooth 60fps rendering

## Best Practices

1. **Keep particle counts reasonable** (20-50) for performance
2. **Use fade rates** so animations end naturally
3. **Provide clear API documentation** for each animation
4. **Test at different speeds** to ensure good experience
5. **Match animation theme** to educational content when possible
6. **Ensure animations don't obscure** important UI elements

## Related Skills

- **MicroSim P5 Generator** - For creating the base MicroSim
- **MicroSim Standardization** - For quality validation
- **MicroSim Screen Capture** - For capturing animation screenshots

## References

- [p5.js Reference](https://p5js.org/reference/)
- [Particle Systems Tutorial](https://p5js.org/examples/hello-p5-flocking.html)
- Animation Library Tester: `/docs/sims/animation-lib-tester/`
- Shared Animations: `/docs/sims/shared/animations/`
