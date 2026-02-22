# Editing & Composition Guide

Load this reference when the user provides input images for editing or multi-reference composition. These principles do not apply to text-to-image generation.

---

## Core Principle

Connect the dots, don't describe. Reference images provide the visual context. The prompt's job is to tell the model what goes where — not re-describe what the model can already see. Over-describing creates competing signals between text and image that degrade output quality.

```
Bad:  "Replace the front figure with a woman who has light fair skin, delicate oval
       face, dark sunglasses, dark brown headscarf, sleeveless beige linen dress..."

Good: "Replace the front figure with the woman from image 3"
```

### The Pointing vs Naming Boundary

Point to reference images for elements that will be used whole (a face, a background, a scene). Explicitly name elements when extracting or transferring parts from a reference to a different context, because generic references like "outfit from image 2" don't tell the model which visual elements to isolate from surrounding context.

```
Bad:  "Replace the clothing with the outfit from image 2"
Good: "Replace the clothing with the white linen button-up shirt and camel
       brown wide-leg trousers"
```

---

## The Edit Grammar

### Reference Block

Start multi-image edit prompts with a reference block that explicitly labels what each image represents. This disambiguates image roles before the directive.

```
Image 3: [role/description]
Image 2: [role/description]
Image 1: [role/description]

[Main directive]
```

Example:
```
Image 2: Reference outfit - white linen shirt and camel trousers
Image 1: Base scene - woman on street to preserve

Replace only the clothing with the white linen shirt and camel trousers.
Keep her exact pose, facial features, and the street background unchanged.
```

### Three-Sentence Directive Structure

Each sentence has one job. Compounding directives into a single sentence dilutes each one's weight.

1. **Replace** (scope) — what changes
2. **Match** (constraint) — how new elements behave
3. **Keep** (preservation) — what stays

### Language Patterns

Use directive verbs ("replace", "match", "keep") not passive phrasing ("should be", "could be"). Directives produce stronger adherence.

"Replace" anchors to the base scene (model edits in-place). "Change" allows full recomposition (model may discard the scene). Always prefer "replace" for in-place edits.

"Only" after the verb constrains scope: "replace only the front figure" is tighter than "replace the front figure."

Sentence order affects spatial placement. The element mentioned last in a spatial assignment tends to get placed in the more prominent position.

---

## Image Ordering & Numbering

Base image (the canvas being edited) goes last in the `--input` list. Gemini numbers images in reverse from input order:

```
--input ref_a.jpg    -> image 3 in prompt
--input ref_b.jpg    -> image 2 in prompt
--input base.jpg     -> image 1 in prompt
```

---

## Multi-Image Composition

Reference images provide visual context — the prompt connects them. Point to images by number, assign elements to positions, and describe only what the model cannot infer from the images themselves. Supports up to 14 reference images.

---

## Character Consistency

- Use follow-up edits for multiple views of the same character
- Reference distinctive features explicitly in follow-ups
- Include "exact same character" or "maintain all design details"
- Save successful designs as reference for future prompts

---

## Semantic Masking

No manual masking needed. Language creates the edit boundary — name the element to define the mask, specify the replacement, constrain the scope:

```
"Using the provided image of a living room, change only the blue sofa
to be a vintage, brown leather chesterfield sofa. Keep the rest of the room,
including the pillows on the sofa and the lighting, unchanged."
```

"Only" defines scope. The element name ("the blue sofa") defines the mask region. The preservation clause ("keep the rest...unchanged") protects everything else.
