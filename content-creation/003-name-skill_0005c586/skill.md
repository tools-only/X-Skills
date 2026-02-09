---
name: video
description: >
  Generate AI-cinematic videos using fal.ai. Primary: Kling 3.0 V3 Pro
  (cinematic quality) and Kling 3.0 Omni (character locking via Elements).
  Multi-shot generation is default — one prompt produces connected scenes
  with character lock, automatic transitions, and native audio. Supports
  Grok Imagine (fast prototyping), DreamActor M2.0 (character driving),
  Vidu Q3 (dialogue), Hailuo 2.3 (VFX). Assembly in DaVinci Resolve /
  CapCut for production; FFmpeg for headless/programmatic use only.
  Claude acts as Creative Director — directing scenes, not engineering prompts.
---

# AI Cinema Video Generator (/video)

Autonomous skill that generates cinematic product/promo videos using AI
multi-shot video generation. Claude directs — writing shot lists, engineering
the blueprint image, and running quality gates before spending on animation.

## Input

$ARGUMENTS

## Architecture

```
/video "promo for X/Twitter showing the stop hook catching a lie"
    |
    +-- Phase 0: Activation (state file, dependency check, platform detection)
    +-- Phase 1: Creative Direction (heavy analysis — 3 agents)
    |   +-- Visual Strategist: storyboard, shot list, visual identity
    |   +-- Reference Researcher: SOTA visual language, platform conventions
    |   +-- Contrarian Critic: what will look like AI slop, failure modes
    |
    +-- Phase 2: Storyboard Design (synthesize agents, define shot plan)
    +-- Phase 3: Directing (director-level prompts per shot — SHORT > LONG)
    +-- Phase 4: Blueprint Image (one reference → entire film)
    +-- Phase 4.5: Quality Gate (Gemini Vision evaluation)
    +-- Phase 5: Multi-Shot Generation (Kling 3.0, native audio, 6 shots/call)
    |   +-- Multi-shot sequences → Kling V3 Pro (DEFAULT — narrative, transitions)
    |   +-- Character-locked sequences → Kling Omni (Elements, 4-image ref)
    |   +-- VFX/physics shots → Hailuo 2.3
    |   +-- Fast prototyping → Grok Imagine (#1 ranked, cheapest)
    |   +-- Character driving → DreamActor M2.0 (single image → performance)
    |
    +-- Phase 5.5: Audio (native Kling audio, ElevenLabs VO, ElevenLabs Music)
    +-- Phase 6: Assembly (DaVinci Resolve / CapCut; FFmpeg for headless)
    +-- Phase 7: Video Evaluation (7-dimension final rubric)
    +-- Phase 8: Iterate or Complete
```

## Triggers

- `/video <description>`
- "create a video", "generate a promo video"
- "make a product video", "hero video"
- "cinematic promo", "AI-generated video"

---

## Phase 0: Activation

### State File

Create `.claude/autonomous-state.json` with `"mode": "video"` at activation.

### Dependency Check

```bash
python3 -c "import fal_client" 2>/dev/null || pip3 install fal-client
python3 -c "from google import genai" 2>/dev/null || pip3 install google-genai
# FFmpeg optional — only needed for headless assembly
which ffmpeg >/dev/null 2>&1 || echo "NOTE: ffmpeg not found — headless assembly unavailable"
```

### API Keys

Required:
- `FAL_API_KEY` or `FAL_KEY` — fal.ai API key
- `GEMINI_API_KEY` — Google Gemini (quality gate + video evaluation)

Optional:
- `ELEVENLABS_API_KEY` — ElevenLabs (voiceover + music, if needed)

If missing, ask the user ONCE at start. Then proceed autonomously.

### Platform Detection

| Signal | Platform | Aspect Ratio | Duration | Notes |
|--------|----------|-------------|----------|-------|
| "X", "Twitter", "tweet" | X/Twitter | 1:1 square | 45-60s | Muted autoplay, burned-in captions mandatory |
| "GitHub", "README", "hero" | GitHub README | 16:9 | 10-20s | Autoplay muted, looping, compressed |
| "LinkedIn", "professional" | LinkedIn | 1:1 or 16:9 | 30-60s | Professional tone |
| "YouTube", "shorts" | YouTube Shorts | 9:16 vertical | 30-60s | Vertical, hook in first 2s |
| "website", "landing page" | Web embed | 16:9 | 15-30s | Looping, no audio assumed |
| No clear signal | Default | 16:9 | 15-30s | Versatile |

**X/Twitter specifics**:
- 1:1 square gets 30-35% more views than 16:9
- Muted autoplay is default — burned-in captions are MANDATORY
- No external links in main tweet (kills reach) — put links in reply
- Replies weighted 75x over likes by X algorithm
- 45-60s sweet spot for engagement

---

## Phase 1: Creative Direction (Heavy Analysis — 3 Agents)

**CRITICAL**: Before writing a single prompt, run creative direction with 3 parallel agents.
Do NOT skip this. The difference between generic AI slop and tasteful cinema is in the thinking
that happens before generation.

Launch ALL 3 agents in a SINGLE message with multiple Task tool calls:

### Agent 1: Visual Strategist

```
Task(
  subagent_type="general-purpose",
  description="Visual Strategist: storyboard design",
  model="opus",
  prompt="""You are the video's visual strategist. Design the storyboard and
visual identity for a cinematic AI-generated product video.

BRIEF: [PASTE USER REQUEST + PRODUCT CONTEXT]
PLATFORM: [DETECTED PLATFORM]
ASPECT RATIO: [RATIO]
TARGET DURATION: [DURATION]

## Your Mission

1. **Story Concept**: Design a narrative arc. Proven structures:
   - "The Lie": Something claims to be done → catches the lie → fixes → verified
   - "The Graveyard Shift": Show unglamorous reality → tool makes it elegant
   - "Before/After": With emotional weight, not feature comparison
   - "The Unexpected": Start somewhere the viewer doesn't expect

2. **Scene Breakdown**: For each scene:
   - Duration (vary for rhythm — never all equal)
   - Narrative beat (what happens emotionally)
   - Key visual element
   - Camera movement (director language: dolly, crane, FPV, rack focus)
   - Audio direction (dialogue, ambient, music, silence)

3. **Visual Identity**: Color language, lighting, typography, surface

4. **The Landing**: Last scene earns its stillness. Tagline + visual treatment.

5. **Model Strategy**: For each scene, recommend:
   - Multi-shot sequence (Kling V3 Pro) — narrative, dialogue, transitions
   - Character-locked (Kling Omni) — Elements with 4-image reference
   - VFX/physics (Hailuo 2.3) — destruction, water, particles
   - Fast prototype (Grok Imagine) — cheapest iteration

6. **Blueprint Image**: Design the ONE reference image that anchors the film.
   Every multi-shot sequence traces back to this single source of truth.

## Output
- Story Concept (title + 2-sentence pitch)
- Scene Breakdown (table: scene, duration, beat, visual, camera, audio, model)
- Visual Identity (colors, lighting, typography, surface)
- Blueprint Image description
- The Landing (tagline + visual treatment)
"""
)
```

### Agent 2: Reference Researcher

```
Task(
  subagent_type="general-purpose",
  description="Reference Researcher: SOTA visual language",
  model="sonnet",
  prompt="""Research current SOTA in AI video generation. Find concrete references.

BRIEF: [PASTE USER REQUEST]
PLATFORM: [DETECTED PLATFORM]

## Research (use mcp__x-ai__x_search_tweets)

1. Product video aesthetics:
   - Search: "AI video" product ad 2026
   - Search: "kling 3.0" multi-shot workflow
   What visual language do the best product videos use?

2. Key practitioners to learn from:
   - @CharaspowerAI (multi-shot prompts, FPV action)
   - @Diesol (Blueprint Image technique, short films)
   - @Ror_Fly (product videos, multi-tool pipeline)
   - @Jonnyvandel (content factory, batch production)
   - @ailker (fal.ai workflows, Super Bowl ad in one day)
   - @mrlnonai (12-hour satire film pipeline)
   What techniques do they use?

3. Emerging models:
   - Search: "Seedance 2.0" video
   - Search: "DreamActor" character driving
   Anything that changes the game?

## Output
- Reference Analysis (3-5 examples, what works and why)
- Practitioner Techniques (what the best creators actually do)
- Model Selection Guide (which model for which shot, from experience)
- What to Avoid (patterns that signal "AI-generated")
"""
)
```

### Agent 3: Contrarian Critic

```
Task(
  subagent_type="general-purpose",
  description="Contrarian Critic: AI slop failure modes",
  model="opus",
  prompt="""You are the video's quality immune system. Identify every way
this AI-generated video could look like slop, and how to prevent it.

BRIEF: [PASTE USER REQUEST]
TARGET AUDIENCE: [WHO WILL SEE THIS]

## The AI Video Slop Detector

**Visual Tells**: Uniform soft focus, plastic textures, algorithmic symmetry,
purposeless particles, concept art aesthetic, galaxy backgrounds, AI glow/HDR.

**Motion Tells**: Uniform speed, morphing, unnatural cloth, text warping,
camera ignoring physics, face distortion, object drift.

**Audio Tells (NEW)**: Nonsense dialogue (Kling invents speech if not directed),
generic ambient that doesn't match scene, robotic voice cadence.

**Compositional Tells**: Everything centered, no hierarchy, background as
busy as foreground, too many things competing.

## Your Mission

1. **Failure Mode Catalog**: For THIS video, list 10 most likely failures
2. **Authenticity Test**: What signals "designed" vs "generated"?
3. **Prompt Guardrails**: Words/phrases to NEVER use in generation prompts
4. **Audio Guardrails**: How to prevent nonsense dialogue and generic audio
5. **The Nuclear Option**: If AI cinema can't serve this audience, what instead?

## Output
- Failure Mode Catalog (10 items with prevention)
- Authenticity Checklist
- Prompt Blacklist (words to avoid)
- Audio Direction Rules
- Alternative Assessment
"""
)
```

---

## Phase 2: Storyboard Design

After all 3 agents return, synthesize into a concrete storyboard.

### Synthesis Checklist

1. **Resolve agent tensions**: Critic's failure modes are HARD CONSTRAINTS.
   Strategist's vision operates WITHIN those constraints.

2. **Lock the story**: Write in one sentence:
   "The video tells the story of [WHAT HAPPENS] to show that [PRODUCT VALUE]."

3. **Design the Blueprint Image**: ONE reference frame that defines:
   - The character/subject in their hero pose
   - The lighting setup and color palette
   - The "neutral" state of the story
   This image anchors every multi-shot generation.

4. **Write the Global Style prefix**: One line that forces visual consistency:
   ```
   Global Style: Cinematic, shot on 35mm KODAK film stock, warm amber tones,
   shallow depth of field, natural grain, single-source lighting.
   ```

5. **Assign models per shot**:
   - Narrative sequences with character continuity → Kling V3 Multi-Shot
   - Character-locked sequences needing Elements → Kling Omni
   - VFX-heavy shots → Hailuo 2.3
   - Budget iteration → Grok Imagine
   - Character driving from single image → DreamActor M2.0

6. **Group shots into multi-shot batches**: Consecutive narrative scenes that
   share characters/setting → ONE multi-shot API call (up to 6 shots).
   Individual shots ONLY when scenes are visually distinct.

7. **Audio direction per shot**: Every shot MUST specify one of:
   - `Dialogue: "exact words"` (with accent/tone if needed)
   - `Audio: ambient description, no dialogue`
   - `Audio: none` (for silent shots / title cards)
   NEVER leave audio unspecified — Kling invents nonsense.

### Storyboard Rules

- **Variable duration**: Never all scenes the same length. Rhythm from variation.
- **Hard cuts only**: CROSSFADE_DURATION = 0.0. Dissolves feel like slideshow.
- **Text budget**: Max 2 lines of text per frame. Rendered BY the AI model.
- **Element budget**: One subject, one light source, one surface per shot.
- **Hands below frame**: AI cannot render hands in motion (typing, gesturing).
  Use over-the-shoulder, medium shots, or frame hands out. Never close-up hands.
- **Short prompts**: 1-3 sentences per shot. Direct like a filmmaker.
- **The Landing earns stillness**: Last scene nearly static. Time to absorb.
- **First iteration without audio**: Use `--no-audio` for iteration 1 to assess
  visuals independently. Add audio on the final pass — it's cheaper and you
  avoid regenerating expensive clips just to fix audio direction.

---

## Phase 3: Directing

**THE PARADIGM SHIFT.** You are not "prompt engineering." You are DIRECTING.

A prompt engineer describes what they want to see.
A director tells the camera what to do and the actors how to perform.

> "The trick: don't write 'girl in apartment with gummies.' That gives you slop.
> Write it like you're behind the camera." — @heyrobinai

> "Short prompts work better with Kling 3.0. That's incredible!" — @BarrakAli

> "Finally, an AI video tool that feels like a director's tool and not just
> a prompt lottery." — @Camus194878

### The Core Rule

**Short, precise direction > verbose description.** On Kling 3.0, practitioners
report that concise 1-3 line prompts per shot consistently outperform
paragraph-length descriptions. Every word must earn its place.

### Multi-Shot Prompt Format (PRIMARY)

This is the format practitioners use for Kling 3.0 multi-shot:

```
Global Style: [Visual treatment — film stock, palette, lighting, mood]

SHOT1: [Camera type + angle], [subject] [action]. Camera: [movement].
[One atmospheric detail]. Audio: [direction].

SHOT2: [Camera type + angle], [new perspective]. Camera: [movement].
[What changes from SHOT1]. Audio: [direction].

SHOT3: [Shot type], [resolution/payoff]. Camera: [movement].
[Final emotional state]. Audio: [direction].
```

**Real example (@CharaspowerAI — F1 race, 174 likes, 8.3K views):**
```
SHOT1 Wide tracking shot, a Formula 1 car rockets down a long straight
under golden-hour light. Heat distortion shimmers off the asphalt.
Grandstands blur past, the engine screams.

SHOT2 Low-angle front view, the car's nose fills the frame as it brakes
hard into a hairpin. Tire smoke erupts. Camera locked to apex.

SHOT3 Drone shot from above, the car threads through the chicane.
Long shadows stretch across the circuit.
```

**Real example (@ai_artworkgen — Global Style technique):**
```
Global Style: Dark Cinematic Period Drama, shot on 35mm film with
high contrast and flickering torchlight.

Shot 1 (0:00-0:05): Slow, tense Dolly In moves toward the figure.
Shot 2 (0:05-0:10): Close-up, firelight plays across weathered face.
Shot 3 (0:10-0:15): Wide reveals the full scene. Torches gutter.
```

### Blueprint Image Prompt (Nano Banana Pro)

For the ONE reference image that anchors the film:

```
[Purpose frame: what kind of image]. [Subject: physical description,
materials, proportions, pose]. [Background: physical surface, not absence
— "painted studio cyclorama" not "void"]. [Lighting: equipment, position].
[Composition: subject size, position]. [Format: aspect ratio].
[Negatives: No lens flares. No particles. No fog. Clean studio photography.]
```

Keep this prompt detailed and physically-grounded. The blueprint image is
the exception to the "short > long" rule — it needs precision because
it defines the visual identity of the entire film.

### FPV Dynamic Shot Formula

For action/establishing shots (@CharaspowerAI's signature technique):

```
Dynamic hyperspeed FPV shot [verb] across [environment] during [condition],
[verb] through [obstacle] and [verb] into [new space]. [Light description].
[Film stock reference].
```

Example: "Dynamic hyperspeed FPV shot darting across a frostbitten tundra
during a whiteout, jumps through a crack in an iceberg and enters a flooded
cave system. Bioluminescent glow. Shot on 35mm film."

### Camera Terms That Work

| Term | Effect | Best For |
|------|--------|----------|
| `FPV shot` / `First-person POV` | Immersive, through-the-action | Action, chase |
| `Wide tracking shot` | Full scene with movement | Establishing shots |
| `Medium close-up` | Character emotion | Dialogue, reaction |
| `Drone shot` | Aerial, sweeping | Landscapes, reveals |
| `Low-angle tracking` | Power, drama, speed | Vehicles, heroes |
| `Side tracking shot` | Parallel movement | Speed, proximity |
| `Macro photography` | Extreme detail | Textures, product |
| `Nearly static, faint handheld drift` | Grounded realism | Character moments |

**Movement phrases that work**: "Camera races above and slightly behind",
"Camera matches speed down the cliff", "Camera trembles with every stomp",
"Slow push-in follows through", "Quick tilt-up to reveal",
"Imperceptible pull-back over 6 seconds".

### Film Stock References (They Work)

Practitioners confirm these reliably influence aesthetic:
- `shot on 35mm film` — warm grain, organic
- `KODAK film stock` — golden warmth
- `shot on Arri Alexa` — digital cinema, clean
- `high contrast and flickering torchlight` — period drama
- `low saturation, film-grain` — anti-AI-sheen, organic

### Negative Prompt Checklist

From practitioners (@DaveSkorupski, @rezkhere). Include as negative prompt
or append to end of main prompt:

```
Quality: blurry, noisy, low resolution, compressed, worst quality
Anatomy: deformed, extra limbs, missing fingers, disproportionate
Motion: jittery, frozen pose, robotic motion, stuttering
Artifacts: watermark, text overlay, split screen, duplicate elements
Lighting: overexposed, flat lighting, uniform glow, HDR haze
Style: AI glow, plastic skin, concept art look, stock footage sheen
```

### Audio Direction (CRITICAL)

**WARNING**: Kling 3.0 INVENTS NONSENSE DIALOGUE if you do not specify audio.
Every shot MUST include audio direction. Three patterns:

```
# Pattern 1: Dialogue
Audio: He says "alright guys, look at this. what a beautiful machine."
       American accent, casual tone.

# Pattern 2: Ambient
Audio: ambient office hum, keyboard clicks, no dialogue.

# Pattern 3: Silent
Audio: none.
```

Multiple accents work — Indian, African, posh English, cockney all confirmed
by @Uncanny_Harry. Performance direction works — "whispered", "shouted",
"with quiet determination" all produce nuanced vocal performance.

### Screen Content in Prompts

For shots that include monitors, terminals, or screens — describe the content
as COLOR and MOOD, not readable text. The viewer at phone distance reads
the light on the character's face, not the screen text.

```
# Describe screen as ambient light source:
"Green terminal glow on developer's face — tests passing."
"Screen shifts to harsh red — failures detected."
"Cool blue code scrolling, rapid changes."
```

### Words to NEVER Use in Prompts

These trigger AI defaults and produce generic output:

- beautiful, stunning, cinematic (ironically), elegant, sleek, gorgeous
- epic, dramatic, atmospheric, moody, ethereal
- high-quality, professional-looking, award-winning
- cyberpunk, steampunk, sci-fi, futuristic
- minimalist, Jony Ive, Apple-designed
- Blade Runner, Tron, Matrix

Instead: describe the physical scene. What camera, what light, what action.

---

## Phase 4: Blueprint Image Generation

The Blueprint Image is ONE reference frame that defines your entire film.
Generate it with Nano Banana Pro, then feed it as the reference for every
multi-shot generation.

> "Every sequence in this short film starting with this same image in Kling 3.0.
> I call this my 'blueprint image'" — @Diesol (322K views)

**Model**: `fal-ai/nano-banana-pro`
**Cost**: ~$0.08-0.15 per image
**Resolution**: 2K

```python
import fal_client

result = fal_client.subscribe(
    "fal-ai/nano-banana-pro",
    arguments={
        "prompt": blueprint_prompt,
        "aspect_ratio": "1:1",
        "resolution": "2K",
        "num_images": 1,
        "output_format": "png",
    },
)
blueprint_url = result["images"][0]["url"]
```

### Character Consistency Strategy

**For Kling V3 Multi-Shot**: Upload the blueprint image as `image_url`.
The model extracts visual traits and maintains them across all shots in
the multi-shot sequence. ONE strong reference + multi-shot = consistency.

> "Character consistency from a single frame combined with Kling 3.0's
> multishot system is just insane." — @CharaspowerAI

**For Kling Omni (when Elements are needed)**: Upload 4 reference images:
1. Front face
2. 3/4 angle
3. Side profile
4. Full body

> "Kling 3.0 character references are fantastic. With Omni 3.0 you can
> upload 4 images to get stronger consistency." — @jerrod_lew

**CRITICAL**: Character Elements work in **Omni only**. In base V3, Elements
often BREAK generation. Multiple practitioners confirm this.

> "With base Kling 3.0 when u choose the character from the elements to
> 'lock consistency' the gen breaks and simply does not follow orders."
> — @Mr_impossible5

### When to Generate Blueprint Images

- **Always** for the first multi-shot batch (establishes character + setting)
- **Reuse batch 1's blueprint** for subsequent batches — this is the primary
  cross-batch character consistency mechanism. ONE blueprint → entire film.
- **Always** for title cards — but hold static via FFmpeg, don't animate with Kling
- **Never** generate separate blueprints per batch unless the visual world changes

### Static Title Cards (FFmpeg Hold, Not Kling Animation)

Title cards with typography should be generated as a blueprint image (Nano Banana
Pro renders text well) then held as a static video via FFmpeg — NOT animated with
Kling. Kling warps text during animation, and you're paying ~$0.42 for 6 seconds
of motion that actively degrades the result.

```python
# In pipeline code, use "static" model type for title cards:
ffmpeg -y -loop 1 -i title_card.png -t 6 \
    -vf "scale=1080:1080,noise=alls=2:allf=t,setsar=1" \
    -c:v libx264 -preset medium -crf 18 -pix_fmt yuv420p \
    title_card.mp4
```

---

## Phase 4.5: Blueprint Quality Gate

**DO NOT SKIP.** The quality gate prevents spending on animating bad references.

### 6-Dimension Evaluation (Gemini Vision)

Upload each blueprint to Gemini 2.5 Flash:

| Dimension | Weight | What It Measures |
|-----------|--------|-----------------|
| aesthetic_authenticity | 30% | Does it look DESIGNED or GENERATED? |
| composition_restraint | 15% | Breathing room? Or cluttered AI maximalism? |
| color_fidelity | 15% | Colors match prompt? Background actually dark? |
| prompt_adherence | 15% | All prompted elements present and correct? |
| lighting_atmosphere | 15% | Light from specific source with falloff? |
| animation_readiness | 10% | Will this animate well? Clean shapes? |

### Hard Ceilings

| Condition | Ceiling |
|-----------|---------|
| Instant "AI" recognition | aesthetic_authenticity capped at 5 |
| >3 AI tells present | aesthetic_authenticity capped at 6 |
| >2 unprompted decorative elements | composition_restraint capped at 6 |
| Background texture when plain surface prompted | composition_restraint capped at 6 |
| Uniform bloom/glow on everything | lighting capped at 4 |
| Wrong dominant color | color_fidelity capped at 4 |

### Pass Criteria

- Overall weighted score >= 7
- No dimension below 5
- aesthetic_authenticity >= 6

### Retry Logic

If a blueprint fails: read `regeneration_guidance`, append to prompt, regenerate.
Max 2 retries. After 3 failures, proceed with warning.

---

## Phase 5: Multi-Shot Generation

Multi-shot is the DEFAULT. Individual shots are the exception.

> "This isn't clips stitched together anymore. It's cinematic sequencing."
> — @rahulnanda86

### Model Registry (February 2026)

| Model | Endpoint (fal.ai) | Cost | Max Duration | Best For |
|-------|-------------------|------|-------------|----------|
| **Kling V3 Pro** | `fal-ai/kling-video/v3/pro/image-to-video` | ~$0.07/s | 15s | Cinematic quality, single shots AND multi-shot (via `multi_prompt` param) |
| ~~Kling V3 Multi-Shot~~ | ~~`fal-ai/kling-video/v3/pro/multi-shot`~~ | — | — | **DOES NOT EXIST.** Multi-shot uses the V3 Pro endpoint with `multi_prompt` parameter. |
| **Kling V3 Standard** | `fal-ai/kling-video/v3/standard/image-to-video` | ~$0.05/s | 15s | Budget single shots |
| **Kling Omni Pro** | `fal-ai/kling-video/o3/pro/image-to-video` | ~$0.224/s | 15s | Character Elements (4-image ref), controllability |
| **Kling Omni Std** | `fal-ai/kling-video/o3/standard/image-to-video` | ~$0.168/s | 15s | Budget character work |
| **Grok Imagine** | Via X API | ~$0.07/s | 15s | #1 ranked, fastest, cheapest prototyping |
| **Hailuo 2.3** | `fal-ai/hailuo/video` | varies | 6s | VFX, physics, water, explosions |
| **DreamActor M2.0** | `fal-ai/dreamactor-m2` | TBD | varies | Single-image character driving |
| **Vidu Q3** | `fal-ai/vidu/q3` | TBD | 1-16s | Native dialogue/sound, up to 1080p |

**NOT on fal.ai** (do not use fal endpoints for these):
- Veo 3.1 — Google Vertex AI or VideoSOS only
- Seedance 2.0 — ByteDance, CapCut/Jimeng only. Monitor for fal integration.
  Being called "the DeepSeek of video." If it lands on fal, priority-test it.

### V3 vs Omni: The Critical Distinction

```
V3 Pro = Best raw cinematic quality. Use for final production shots.
         Character Elements BREAK in V3. Do NOT use V3 for Element locking.

Omni (O3) = Best controllability. Character Elements WORK here.
            Use when you need 4-image character reference.
            Use when you need fine-grained control over motion/expression.

V3 Standard = Budget V3. Same quality tier, lower cost. Good for iteration.
```

> "V3 = high quality cinematic outputs. O3 = smart model, video editing,
> controllability. Same idea as Kling 2.6 vs O1." — @gokayfem (fal Academy)

### Multi-Shot API Call (DEFAULT)

**CRITICAL**: There is NO separate multi-shot endpoint. Multi-shot uses the
standard V3 Pro `image-to-video` endpoint with the `multi_prompt` parameter
instead of `prompt`. The endpoint `v3/pro/multi-shot` returns 404.

```python
# Multi-shot: use multi_prompt (array), NOT prompt (string)
result = fal_client.subscribe(
    "fal-ai/kling-video/v3/pro/image-to-video",   # Same endpoint as single-shot
    arguments={
        "multi_prompt": [                           # Array of per-shot prompts
            {"prompt": "Global Style: ...\n\nSHOT1: ...", "duration": "5"},
            {"prompt": "Global Style: ...\n\nSHOT2: ...", "duration": "5"},
            {"prompt": "Global Style: ...\n\nSHOT3: ...", "duration": "5"},
        ],
        "start_image_url": blueprint_image_url,     # NOT "image_url"
        "generate_audio": True,
        "shot_type": "customize",                   # Required for multi_prompt
    },
)
video_url = result["video"]["url"]
```

**Key API details:**
- `multi_prompt` and `prompt` are **mutually exclusive** — provide one or the other
- Each element: `{"prompt": str, "duration": str}` (duration: "3"-"15")
- `shot_type: "customize"` is **required** when using `multi_prompt`
- Use `start_image_url` (NOT `image_url`) for the blueprint reference image
- Prepend the Global Style to each shot's prompt for consistency
- Per-shot durations must sum to total desired duration

**Multi-shot eliminates manual stitching.** One call → multiple connected
scenes with natural transitions, consistent characters, synchronized audio.

**When to use multi-shot vs individual shots:**
- Multi-shot (DEFAULT): 2-6 consecutive narrative scenes, shared characters
- Individual: Visually distinct scenes (different locations, styles)
- Individual: Need per-shot model selection (Hailuo for VFX, Omni for char)
- Individual + static: Title cards (generate blueprint image, hold via FFmpeg)

### Kling Omni — Character Reference

For character-driven videos that need Elements:

```python
result = fal_client.subscribe(
    "fal-ai/kling-video/o3/pro/image-to-video",
    arguments={
        "prompt": director_prompt,
        "image_url": character_front_url,
        "reference_image_urls": [
            character_34_url,
            character_side_url,
            character_full_url,
        ],
        "duration": "12",
        "generate_audio": True,
    },
)
```

### Grok Imagine — Fast Prototyping

#1 on Artificial Analysis Arena. Fastest generation. Use for:
- Quick drafts before committing to Kling Pro
- Budget iteration ($4.20/min, highest quality-per-dollar)
- Product ads (full ad for under $1 — @jamiebxne)

### Shot Duration Strategy (from practitioners)

- **3 seconds/shot** — Fast-paced action montages (@mxvdxn, @Ror_Fly)
- **5 seconds/shot** — Balanced narrative (@heyrobinai, @PJaccetturo)
- **15 seconds single** — Long takes with internal camera moves (@LudovicCreator)

---

## Phase 5.5: Audio Integration

Audio separates professional video from demo reels. The 2026 shift:
models now generate audio WITH the video.

### Option A: Kling Native Audio (DEFAULT for dialogue)

Kling 3.0 generates synchronized audio natively — dialogue with lip-sync,
ambient sounds, vocal performance. Set `generate_audio: True`.

**Best for**: Character dialogue, ambient scenes, narration-over-action.

> "Kling 3.0 adding native synchronized audio to AI video is a real workflow
> shift. Dialogue, ambience and SFX arriving with the clip makes drafts feel
> like scenes." — @BlueLightningTV

> "The lip sync is ridiculous." — @SebJefferies (stress-tested extensively)

Multiple accents work. Performance direction (whispered, shouted, etc.) works.

### Option B: ElevenLabs (Narration + Music)

For voiceover that doesn't need lip-sync, or background music:
- ElevenLabs voice cloning for narration
- ElevenLabs Music for score generation (used by @mrlnonai for full film scores)
- Add in DaVinci Resolve during assembly

### Option C: CassetteAI / Suno (Sound Effects + Music)

- CassetteAI for additional SFX (@ailker uses this)
- Suno for instrumental soundtrack generation

### Audio Strategy Per Video Type

| Video Type | Audio Approach |
|-----------|---------------|
| X/Twitter promo | Muted autoplay — burned-in captions, native audio as bonus |
| Narrative short | Kling native for dialogue + ElevenLabs Music for score |
| Product demo | ElevenLabs narration + light ambient |
| GitHub README | No audio (autoplay muted, looping) |
| UGC ad | Kling native lip-sync (killed the talking head workflow — @SebJefferies) |

---

## Phase 6: Assembly

### For Production: DaVinci Resolve or CapCut (Recommended)

Multi-shot output is already a single coherent video. Import into your NLE for:
- Color grading (CRITICAL — different generations diverge in color/contrast)
- Burned-in captions/subtitles
- Audio mixing (VO + music + native audio)
- Transitions between multi-shot segments
- Export with platform-specific encoding

> "Editing is where the movie is made. AI generates clips — editing makes it
> a movie. Without editing, your AI video looks amateur." — @Riodefirst

> "Another real hidden problem in AI video — when combining shots from different
> generations, color and contrast can diverge. What remains is boring color
> grading as manual matching." — @MateuszStys1

DaVinci Resolve is free. CapCut is free. Both are professional-grade.

### For Headless/Programmatic Use: FFmpeg

When running from a script without NLE access (CI/CD, automated pipelines,
Claude Code autonomous mode), FFmpeg assembly is acceptable:

```bash
# Normalize resolution between clips
for clip in clips/*.mp4; do
    ffmpeg -y -i "$clip" -vf "scale=1080:1080,setsar=1" -c:v libx264 "norm/$(basename $clip)"
done

# Concat with hard cuts
ffmpeg -y -f concat -safe 0 -i concat_list.txt \
    -vf "noise=alls=2:allf=t" \
    -c:v libx264 -preset medium -crf 18 \
    -pix_fmt yuv420p -movflags +faststart \
    output/final.mp4
```

**FFmpeg is ONLY for**: Resolution normalization, format conversion,
concatenation, subtitle burn-in, film grain.

**FFmpeg is NOT for**: Text rendering, visual compositing, color grading,
screen content. AI models handle all visual content.

---

## Phase 7: Video Evaluation (7-Dimension Rubric)

After assembly, evaluate the complete video using Gemini Vision.

| Dimension | Weight | What It Measures |
|-----------|--------|-----------------|
| hook_and_narrative | 20% | Story arc works? Clear problem-solution-proof? Scroll-stop opening? |
| motion_craft | 20% | Beyond static? Purposeful camera movement? Or slideshow? |
| character_consistency | 15% | Same face/body across shots? Same wardrobe, hair, build? |
| pacing_and_rhythm | 15% | Video breathes? Tension → release? Variable shot duration? |
| audio_integration | 10% | Audio matches action? No nonsense dialogue? Emotional sound? |
| emotional_impact | 10% | Viewer FEELS something? Frustration → satisfaction? |
| visual_identity | 10% | Coherent color system? Typography? Professional? |

### Hard Ceilings

| Condition | Ceiling |
|-----------|---------|
| Slideshow feel (static frames with fades) | motion_craft capped at 4 |
| No clear narrative arc | hook_and_narrative capped at 5 |
| Generic "dark mode dev tool" aesthetic | visual_identity capped at 5 |
| Face distortion or object morphing | motion_craft capped at 3 |
| AI glow / HDR haze on surfaces | visual_identity capped at 4 |
| Character inconsistency between shots | character_consistency capped at 3 |
| Different face across 2+ shots | character_consistency capped at 4 |
| Nonsense dialogue audible | audio_integration capped at 2 |
| Generic ambient that doesn't match scene | audio_integration capped at 5 |
| Verbose/unnatural dialogue | audio_integration capped at 5 |

### Pass: overall >= 7. Target: 8+.

If < 7, identify weakest dimension and iterate:
- motion_craft low → switch to multi-shot for better continuity
- hook_and_narrative low → restructure storyboard
- character_consistency low → use Omni with 4-image reference
- audio_integration low → rewrite audio direction, regenerate
- visual_identity low → refine blueprint image, regenerate

---

## Phase 8: Iterate or Complete

### Iteration Loop

1. Identify the weakest dimension
2. Determine fix layer:
   - **Prompt fix** (cheapest): Rewrite prompts, regenerate specific shots
   - **Model switch**: Omni for character work, Hailuo for VFX
   - **Multi-shot consolidation**: Merge individual shots into one multi-shot
   - **Blueprint fix**: Regenerate the reference image
   - **Structure fix** (most expensive): Restructure storyboard
3. Only regenerate what changed
4. Re-evaluate
5. Max 3 full iterations before asking user for direction

### Cost Awareness

| Video Type | Blueprint | Animation | Audio | Total |
|-----------|-----------|-----------|-------|-------|
| GitHub README (15s, 1 multi-shot) | $0.15 | $1.05 | $0 | ~$1.20 |
| X/Twitter promo (45s, 2 multi-shot) | $0.30 | $3.15 | $0 | ~$3.50 |
| Narrative with dialogue (60s, 3 batches) | $0.45 | $4.20 | $0 | ~$4.65 |
| With VO + music (60s) | $0.45 | $4.20 | $2.00 | ~$6.65 |
| With 2 iteration cycles | x1.5-1.8 | | | |

Multi-shot is significantly cheaper than individual shots — fewer API calls,
built-in transitions, no manual stitching cost.

### Completion Checkpoint

```json
{
  "self_report": {
    "is_job_complete": true,
    "code_changes_made": true,
    "linters_pass": true,
    "category": "pattern"
  },
  "reflection": {
    "what_was_done": "Generated 45s X/Twitter promo with 2 multi-shot batches, scored 8/10",
    "what_remains": "none",
    "key_insight": "Reusable lesson about video generation approach (>50 chars)",
    "search_terms": ["video", "promo", "ai-cinema", "multi-shot"]
  },
  "evidence": {
    "video_path": "video/output/hero_video.mp4",
    "evaluation_score": 8,
    "cost_usd": 3.50,
    "duration_seconds": 45,
    "multishot_batches": 2,
    "iterations": 1,
    "models_used": ["nano-banana-pro", "kling-v3-multi-shot"]
  }
}
```

---

## Anti-Slop Visual Patterns

### AI Tells to Eliminate

| AI Tell | What It Looks Like | Prevention |
|---------|-------------------|------------|
| Uniform soft focus | Everything equally blurry | Specify depth of field |
| Plastic textures | Surfaces look waxy | Describe material: "frosted matte glass" |
| Algorithmic symmetry | Everything perfectly centered | Asymmetric compositions |
| Purposeless particles | Floating sparkles | Explicitly exclude particles |
| Concept art aesthetic | ArtStation render look | Describe "product photography" |
| Symmetric radial glow | Light radiating equally | Single directional source |
| Galaxy/nebula backgrounds | Cosmic backdrops | "Painted studio cyclorama" |
| Circuit board textures | PCB for "tech" | Flat uninterrupted dark surface |
| AI glow / HDR haze | Uniform bloom | Negative: "No AI glow, no HDR" |
| Face distortion | Warped features during motion | Kling V3 (best face stability) |
| Hands in motion | Fingers merge, extra digits, claw shapes | Frame hands out of shot. Over-the-shoulder or medium shots. |
| Object drift | Subject changes shape | Reference images + multi-shot |
| Stock footage sheen | Generic clean look | Material specifics, film grain |
| Nonsense dialogue | Kling invents speech | Always specify audio direction |
| Verbose prompt artifacts | Scene looks "described" | Short director-style prompts |
| Character Element breaks | Elements in V3 crash gen | Use Omni for Elements |

### The Authenticity Test

Before approving any frame: "Would a senior designer at Linear or Stripe
approve this?" If "it looks AI-generated," reject it.

**Designed** = intentional restraint, motivated lighting, earned negative
space, typography with hierarchy, color that means something.

**Generated** = everything the model can do at once, uniform treatment,
decorative complexity, impressive but purposeless.

---

## Pipeline Script

The pipeline lives at `video/generate_hero_video.py`.

### Required Structure

```python
I2V_MODELS = {
    "kling_v3_pro": "fal-ai/kling-video/v3/pro/image-to-video",
    "kling_v3_multi": "fal-ai/kling-video/v3/pro/image-to-video",  # Same endpoint, multi_prompt param
    "kling_v3_std": "fal-ai/kling-video/v3/standard/image-to-video",
    "kling_omni_pro": "fal-ai/kling-video/o3/pro/image-to-video",
    "kling_omni_std": "fal-ai/kling-video/o3/standard/image-to-video",
    "hailuo": "fal-ai/hailuo/video",
    "dreamactor": "fal-ai/dreamactor-m2",
    "vidu_q3": "fal-ai/vidu/q3",
}

T2V_MODELS = {
    "kling_v3_t2v": "fal-ai/kling-video/v3/pro/text-to-video",
    "kling_omni_t2v": "fal-ai/kling-video/o3/pro/text-to-video",
}

KEYFRAME_MODEL = "fal-ai/nano-banana-pro"

STORYBOARD = [
    {
        "id": "batch_1",
        "name": "Setup and Turn",
        "prompt": "Global Style: ...\n\nSHOT1: ...\nSHOT2: ...",
        "blueprint_prompt": "...",  # Nano Banana Pro prompt
        "model": "kling_v3_multi",
        "duration": 15,
        "generate_audio": True,
        "is_multishot": True,
    },
]
```

### CLI Flags

```
--dry-run          Preview plan + cost without calling APIs
--skip-eval        Skip blueprint quality gate
--skip-assembly    Generate assets only, skip assembly
--regenerate X     Force regeneration of specific shot or 'all'
--model M          Override model for all shots
--no-audio         Skip audio generation
--blueprint-only   Generate blueprint image only, stop before animation
```

---

## Troubleshooting

### "Generic-looking frames"
Prompts describe aesthetic CONCEPTS instead of physical OBJECTS.
Fix: Rewrite as director instructions — camera, action, light.

### "Kling generates nonsense dialogue"
Audio is enabled but no dialogue/audio direction in prompt.
Fix: ALWAYS specify `Audio: [description]` or `Dialogue: "[words]"` per shot.

### "Character breaks between shots"
Using V3 base for character Elements.
Fix: Use Omni (O3) for character Elements. V3 Elements are broken.

### "Background filling (galaxies, circuits)"
Abstract descriptions like "void" or "dark space."
Fix: "Painted studio cyclorama, seamless, no texture."

### Text hallucination
Complex text compositions with multiple elements.
Fix: Max 2 lines per frame. Large text. Isolated composition.

### Face distortion during camera moves
Model struggles with face geometry during motion.
Fix: Kling V3 Pro has best face stability. Negative: "No face distortion."

### Object morphing / drift
Subject slowly changes shape over the clip.
Fix: Use reference images + multi-shot. Negative: "No morphing."

### Color mismatch between clips
Different generations produce different color/contrast.
Fix: Color grade in DaVinci Resolve. Apply consistent LUT across project.

### "Veo 3.1 endpoint not found on fal.ai"
Veo 3.1 is NOT on fal.ai. Only available via Google Vertex AI or VideoSOS.
Fix: Use Kling V3 Pro for photorealistic shots, or access Veo through Vertex.

### Animation too expensive for iteration
Fix: Use Grok Imagine (~$4.20/min) or Kling V3 Standard for draft iterations.
Switch to V3 Pro only for final pass.

### "Path /v3/pro/multi-shot not found" (404)
There is NO separate multi-shot endpoint. Multi-shot uses `multi_prompt`
parameter on the standard `image-to-video` endpoint.
Fix: Use `fal-ai/kling-video/v3/pro/image-to-video` with `multi_prompt` array
instead of `prompt` string. Add `shot_type: "customize"`.

### "image_url is not a valid parameter"
The V3 API uses `start_image_url`, not `image_url`.
Fix: Change `image_url` to `start_image_url` in your API arguments.

### Hands look distorted or have extra fingers
AI models cannot reliably render hands in motion.
Fix: Frame hands below/outside shot. Use over-the-shoulder or medium shots.
Never write "hands typing" or "fingers on keyboard" in prompts.

### Title card text warps during animation
Kling distorts typography when animating static text.
Fix: Generate the title card as a blueprint image (Nano Banana Pro renders text
well), then hold it as a static video via FFmpeg. Do not animate title cards.

### Burned-in captions for X/Twitter
FFmpeg `drawtext` filter works for caption burn-in. Use `enable='between(t,start,end)'`
for timing. Single quotes in the filter value protect commas from being parsed
as filter separators. Place text at `y=h-th-100` for lower-third positioning.

---

## Model Quick Reference

| Need | Model | Why |
|------|-------|-----|
| Multi-shot narrative (DEFAULT) | Kling V3 Pro + `multi_prompt` | Character lock, transitions, native audio |
| Single cinematic shot | Kling V3 Pro + `prompt` | Best raw quality |
| Static title card | Blueprint image + FFmpeg hold | No animation cost, no text warping |
| Character Elements + locking | Kling Omni (O3) | Only model where Elements work reliably |
| Fast cheap prototyping | Grok Imagine | #1 ranked, fastest, $4.20/min |
| VFX / physics | Hailuo 2.3 | Best physics simulation |
| Character driving | DreamActor M2.0 | Single image → performance |
| Dialogue-heavy | Vidu Q3 | Native multi-speaker dialogue |
| Blueprint images | Nano Banana Pro | Best character consistency in images |
| Voiceover + music | ElevenLabs | Voice cloning + score generation |
| Upscaling (if needed) | Topaz Video AI | Industry standard (Kling 3.0 outputs 4K natively) |

---

## Key Practitioners (Follow for Ongoing SOTA)

| Creator | Handle | Specialty |
|---------|--------|-----------|
| CharaspowerAI | @CharaspowerAI | Multi-shot prompts, FPV action, prompt sharing |
| Diesol | @Diesol | Blueprint Image technique, short films |
| Ror_Fly | @Ror_Fly | Product videos, multi-tool pipelines |
| Uncanny_Harry | @Uncanny_Harry | Performance direction, character acting |
| PJ Accetturo | @PJaccetturo | Cinematic short films (4.2M views) |
| Jonnyvandel | @Jonnyvandel | Content factory, batch production |
| ailker | @ailker | fal.ai workflows, rapid production |
| mrlnonai | @mrlnonai | Complete film pipelines |
| gokayfem | @gokayfem | fal Academy, model comparisons |
| jerrod_lew | @jerrod_lew | Tutorials, Omni vs V3 |
| LudovicCreator | @LudovicCreator | YAML prompting, film-grain aesthetic |

---

## Philosophy

This skill encodes three hard-won lessons:

### 1. Taste comes from constraint, not capability

The model can generate anything. That's the problem. Without constraint,
it generates the statistical average of its training data — "AI-generated."

Every layer of the skill imposes constraints: physically-grounded blueprints,
director-level language, model selection per shot, quality gates, hard cuts,
platform rules, anti-slop patterns.

### 2. AI handles everything visual. Code handles assembly.

AI video models render scenes, characters, text, screen content, lighting,
camera movement, perspective, reflections, audio, dialogue. They handle
these naturally because that's what they're trained on.

Code (FFmpeg / NLE) handles assembly — concatenation, normalization,
encoding, format conversion. These are deterministic operations.

**The anti-pattern**: Using FFmpeg drawtext to composite text onto scenes,
or Remotion to render overlays. This forces manual pixel estimation,
fights camera movement, and produces artificial overlays. The AI model
renders all visual content naturally if you describe it in the prompt.

**Rule**: If it requires understanding the visual scene → AI model does it.
If it's a mathematical transform → code does it.

### 3. Direct, don't describe

The shift from "prompt engineering" to "directing" is real. Short, precise
direction outperforms verbose description on every modern video model.
Write like you're behind the camera giving instructions to actors and
crew, not like you're writing alt text for a screen reader.

Claude is the Creative Director. The APIs generate pixels. The NLE
assembles them. The constraints produce quality.
