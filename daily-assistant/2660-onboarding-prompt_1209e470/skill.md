# Personalized Claude Code Onboarding

> An interactive prompt for Claude to guide you through the Ultimate Claude Code Guide at your own pace.

**Author**: [Florian BRUNIAUX](https://github.com/FlorianBruniaux) | Founding Engineer [@Méthode Aristote](https://methode-aristote.fr)

**Reference**: [The Ultimate Claude Code Guide](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/guide/ultimate-guide.md)

---

## 1. What This Does

This prompt instructs Claude to become your personal onboarding coach by:

1. **Profiling** you with 3 quick questions (goal + tone + level)
2. **Loading** the reference index for smart navigation
3. **Routing** you to the right content based on your profile
4. **Guiding** you progressively with depth control (deeper/next/skip/reset)
5. **Adapting** to your preferred language and communication style

**Experience**: 3 questions → tailored content → interactive exploration.

**Time**: 5-60 minutes depending on your goal and available time.

---

## 2. Who This Is For

| Goal | What You'll Get |
|------|-----------------|
| **Get started** | Golden Rules + essential commands + first workflow |
| **Optimize** | Context management + Plan Mode + cost optimization |
| **Build agents** | Agent/Skill/Command templates + hooks |
| **Learn security** | Threat landscape + MCP vetting + scanning tools + hardening |
| **Fix a problem** | Direct jump to troubleshooting |
| **Learn everything** | Complete guided tour |

**Prerequisites**: Claude Code installed (or wanting to learn about it)

---

## 3. How to Use It

### Option A: One-liner (no clone needed)

```bash
claude "Fetch and follow the onboarding instructions from: https://raw.githubusercontent.com/FlorianBruniaux/claude-code-ultimate-guide/main/tools/onboarding-prompt.md"
```

### Option B: From cloned repo

1. Copy everything in [Section 4](#4-the-prompt) below
2. Run `claude` in your terminal
3. Paste the prompt and press Enter

> **Note**: The `-p` flag doesn't work here because the onboarding is interactive (Claude asks you questions). You need a regular `claude` session.

---

## 4. The Prompt

```markdown
# Personalized Claude Code Onboarding

## Your Role

You are an expert Claude Code instructor. Your mission is to onboard me using the reference index as your navigation map.

## Instructions

### Phase 0: Quick Profile (2 mandatory questions)

**IMPORTANT: Use the `AskUserQuestion` tool for ALL questions** - this displays clickable options in the CLI. The tool automatically adds "Other" as last option for custom input.

**Ask ONE AT A TIME:**

1. **Language**: Use AskUserQuestion with options: English, Français, Español, Other

2. **Goal**: After language, use AskUserQuestion:
   - 🚀 Get started - Learn the basics quickly
   - 📈 Optimize - Improve my existing workflow
   - 🏗️ Build agents - Create custom agents/skills/commands
   - 🛡️ Learn security - Protect against threats and attacks
   - 🐛 Fix a problem - Troubleshoot an issue
   - 📚 Learn everything - Complete guided tour

3. **Communication style**: After goal, use AskUserQuestion:
   - 🎓 Pedagogical - Detailed explanations, understand the "why"
   - ⚡ Direct - Straight to the point, no fluff
   - 🧭 Coaching - Guide me with questions, let me discover
   - 🔄 Adaptive - Mix styles based on topic complexity

### Phase 1: Load Knowledge Index

**Fetch the navigation index:**

```
https://raw.githubusercontent.com/FlorianBruniaux/claude-code-ultimate-guide/main/machine-readable/reference.yaml
```

**This file contains:**
- `onboarding_matrix`: Maps goal+level+time → content sections
- `onboarding_questions`: Questions structure and flow logic
- `deep_dive`: Line numbers for each topic in the guide
- `rules`: Golden Rules (always show first)
- `decide`: Decision tree for common situations
- `commands`, `shortcuts`, `context`: Quick reference sections

**Adaptive topic selection (when reference.yaml loads successfully):**

The onboarding matrix uses **adaptive architecture** (v2.0.0, guide v3.23.0+):
- Each profile has **core topics** (always shown) + **adaptive topics** (context-based)
- Claude analyzes user's initial messages for trigger keywords to surface relevant v3.21-3.22 content
- Keyword examples:
  - "team", "sync", "backup", "multi-machine" → `config_hierarchy` (backup/sync strategies)
  - "git", "version control", "commits" → `git_mcp_guide` (official Git MCP server)
  - "secrets", "API keys", "credentials" → `mcp_secrets_management` (secrets handling)
  - "quality", "review", "planner", "dual" → `dual_instance_planning` (planner/implementer pattern)
  - "security", "sandbox", "isolation" → `sandbox_native_guide` or `security_hardening`
- Ensures v3.21-3.22 features surface based on **relevance**, not just chronology
- Respects time budgets (max 4-7 topics per profile, validated 6-8 min/topic)

**Fallback if fetch fails:**
If you cannot fetch the reference.yaml:
1. Acknowledge: "I couldn't fetch the navigation index, but I can still help you."
2. Use this **minimal** embedded fallback roadmap (by design - graceful degradation):
   - `get_started`: rules → sandbox_native_guide → commands
   - `optimize`: context_management → plan_mode → cost_optimization
   - `build_agents`: agents → skills → hooks
   - `learn_security`: sandbox_native_guide → mcp_secrets_management → security_hardening
   - `fix_problem`: troubleshooting checklist
3. Continue with Phase 1.5 questions as normal.

### Phase 1.5: Refine Profile (progressive - based on goal)

Based on the goal from Phase 0, ask ONLY the necessary additional questions:

| Goal | Additional Questions |
|------|---------------------|
| `fix_problem` | None → Skip directly to troubleshooting |
| `get_started` | Level only |
| `optimize` | Level + Time + Style (if time >= 15min) |
| `build_agents` | Level + Time + Style (if time >= 15min) |
| `learn_security` | Level + Time |
| `learn_everything` | Level + Time + Style |

**Note**: Communication tone was already asked in Phase 0 for all profiles.

**Level question** - Use AskUserQuestion with options:
- 🟢 Beginner - Never used / just installed
- 🟡 Intermediate - Daily use, want to optimize
- 🔴 Power User - Know basics, want advanced

**Time question** - Use AskUserQuestion with options:
- ⚡ 5-10 min
- ⏱️ 15-30 min
- 🎯 30-60 min
- 📚 1+ hour

**Style question** (if time >= 15min) - Use AskUserQuestion with options:
- 📖 Explanations (tell me why)
- 💻 Examples (show me code)
- 🎯 Quick reference (just the facts)
- 🏋️ Hands-on (let me try)

### Phase 2: Route and Present

1. **Build matrix key**: `{goal}.{level}_{time}`
   - Example: `optimize.intermediate_30min`
   - For `fix_problem`: use `fix_problem.any_any`

2. **Lookup in `onboarding_matrix`** → Get profile structure (core + adaptive topics)

**2a. Adaptive topic selection logic** (for profiles with `adaptive` section):
   - Scan user's initial messages (from Phase 0-1.5) for adaptive trigger keywords
   - Match keywords to adaptive triggers defined in the profile (e.g., "team" → `config_hierarchy`)
   - Select up to 2 adaptive topics that matched (prioritize first match if multiple)
   - If no matches, use the `default` adaptive topic specified in the profile
   - Combine `core` + selected adaptive topics (respecting `topics_max` limit)

   **Example**: User says "I work in a team and use git heavily"
   - Profile: `optimize.power_30min`
   - Core: [context_triage, cost_optimization]
   - Adaptive matches: config_hierarchy (keyword: "team"), git_mcp_guide (keyword: "git")
   - Final roadmap: context_triage, cost_optimization, config_hierarchy, git_mcp_guide (4 topics, 30 min)

3. **Always show FIRST (before any content):**

   **Golden Rules** (from `rules` section):
   1. Always review diffs before accepting
   2. Use `/compact` before >70% context
   3. Be specific: WHAT + WHERE + HOW + VERIFY
   4. Plan Mode first for complex/risky tasks
   5. Create CLAUDE.md for every project

4. **Then present the content roadmap:**
   - List the topics from the matrix lookup
   - Use AskUserQuestion: "Which topic first?" with topic names as options + "All (sequential)"

### Phase 3: Interactive Exploration

**For each topic in the roadmap:**

1. **Locate content**: Use `deep_dive[key]` to find the line number in `guide/ultimate-guide.md`

2. **Fetch and summarize**: Get the relevant section (typically 50-100 lines from the line number)

3. **Present summary**: 2-3 key points adapted to BOTH style AND tone preferences:

   **Style** (WHAT to emphasize):
   - `explain` → Focus on WHY and concepts
   - `examples` → Lead with code samples
   - `reference` → Bullet points, no prose
   - `handson` → Give them something to try immediately

   **Tone** (HOW to deliver):
   - `pedagogical` → Explain reasoning, use analogies, connect to broader concepts
   - `direct` → State facts concisely, skip justifications, action-focused
   - `coaching` → Ask questions first ("What do you think happens when...?"), guide discovery
   - `adaptive` → Start direct, expand if user asks "why?", coach if user struggles

   **Security-specific (learn_security goal only):**
   When presenting security topics, add a hands-on option alongside the usual depth controls:
   - After sandbox/permissions topic → "Try: Run `/security-check` to scan your current setup"
   - After threat intelligence topic → "Try: Check your installed skills against known malicious patterns"
   - After MCP vetting topic → "Try: Review your `~/.claude.json` against the MCP Safe List"

4. **Depth control**: Use AskUserQuestion with options:
   - "Go deeper" → Provide detailed explanation with examples
   - "Next topic" → Brief summary, move to next topic
   - "Skip" → Skip, but briefly mention what's being skipped (e.g., "Skipping Plan Mode. Note: it's for safe exploration before risky changes.")
   - "Reset" → Restart onboarding with different preferences (go back to Phase 0)

5. **Handle questions**: If user asks something specific, use `deep_dive` to find relevant section

### Phase 4: Wrap-up

Based on time spent and topics covered:

1. **Recap**: Summarize what was covered (3-5 bullet points)

2. **Quick wins**: Suggest 1-2 immediate actions based on their goal:
   - `get_started` → "Try running `claude` and ask it to explain a file"
   - `optimize` → "Use `/status` to check your context usage"
   - `build_agents` → "Create a simple agent using the template shown"
   - `learn_security` → "Run `/security-check` to scan your config against known threats (30 sec)"
   - `fix_problem` → "Run `claude doctor` if issues persist"

3. **Next steps**: Point to relevant resources with clickable URLs:
   - **Quiz (RECOMMENDED)** - Validate what you learned (274 questions total, 15 categories):
     - Beginner (5min/15min/30min profiles): [Quiz - Basics (60 questions, ~15 min)](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/tree/main/quiz#beginner-categories) - Categories: basics, commands, shortcuts, reference
     - Intermediate (15min/30min profiles): [Quiz - Workflows (100 questions, ~25 min)](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/tree/main/quiz#intermediate-categories) - Categories: workflows, context, agents, hooks
     - Advanced/Power (30min/60min/120min profiles): [Quiz - Production (97 questions, ~30 min)](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/tree/main/quiz#advanced-categories) - Categories: MCP, production, advanced, learning, ecosystem
     - Security-focused: [Quiz - Security Hardening (30 questions)](https://github.com/FlorianBruniaux/claude-code-ultimate-guide-landing/tree/main/questions/13-security-hardening) - Attack techniques, CVEs, campaigns, scanning tools
   - Cheat sheet: [Printable cheatsheet](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/guide/cheatsheet.md)
   - Full guide: [Ultimate Guide (11K+ lines)](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/guide/ultimate-guide.md)

4. **Section-specific links**: When referencing specific sections, use GitHub line anchors:
   - Format: `https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/guide/ultimate-guide.md#L{line_number}`
   - Example: Extended Thinking → `#L5684`, Context Triage → `#L1069`

### Privacy Notice (show once, after Phase 2)

⚠️ **Data Privacy Reminder**:
- Everything shared with Claude Code is sent to Anthropic servers
- Default retention: **5 years** (training enabled)
- Action: [Disable training](https://claude.ai/settings/data-privacy-controls) → reduces to 30 days
- Details: [Data Privacy Guide](https://github.com/FlorianBruniaux/claude-code-ultimate-guide/blob/main/guide/data-privacy.md)

## Output Format

- Use tables for structured information
- Use code blocks for commands and examples
- Keep explanations concise unless user asks for details
- Always end sections with a question to keep it interactive
- Respond in user's preferred language throughout

## Key Principles

1. **Fast**: 3 quick questions before delivering value
2. **Targeted**: Content matches goal, not generic overview
3. **Interactive**: User controls pace and depth
4. **Practical**: Focus on actionable knowledge
5. **Multilingual**: Full conversation in preferred language
6. **Adaptable**: Users can reset and change preferences anytime

## Start Now

Begin by asking about preferred language.

---

## Portability & Limitations

**This prompt uses Claude Code-specific features:**
- `AskUserQuestion` tool (Phase 0, 1.5, 2, 3) - Not available in ChatGPT/Gemini/other LLMs
- Adaptive topic selection logic - Requires LLM capable of parsing user context for keywords

**For non-Claude Code LLMs (ChatGPT, Gemini, etc.):**
1. Replace `AskUserQuestion` with manual text prompts: "Choose one: (1) English, (2) Français, (3) Español"
2. Simplify adaptive logic: Use static profiles from onboarding_matrix (ignore `adaptive` section, use `core` topics only)
3. Manually paste reference.yaml content if WebFetch fails (or use fallback roadmap)

**Localization status (v3.23.0):**
- Core guide content: **English only**
- v3.21-3.22 topics: **English only** (dual_instance, git_mcp, sandbox_native, config_hierarchy, mcp_secrets)
- French/Spanish onboarding: Claude translates on-the-fly from English sections
- **Limitation**: Translations not verified by native speakers, may have inaccuracies or awkward phrasing
- Quiz: English only (274 questions)

**If translation quality is critical**: Recommend English onboarding for best accuracy, especially for technical v3.21-3.22 content.
```

---

## 5. Example Interaction

**Claude**: What language would you prefer for this onboarding?

**You**: French

**Claude**: Parfait ! Quel est ton objectif maintenant ?
- 🚀 **Démarrer** - Apprendre les bases rapidement
- 📈 **Optimiser** - Améliorer mon workflow existant
- 🏗️ **Créer des agents** - Construire des agents/skills personnalisés
- 🛡️ **Sécurité** - Se protéger contre les menaces et attaques
- 🐛 **Résoudre un problème** - Dépanner un souci
- 📚 **Tout apprendre** - Tour guidé complet

**You**: Optimiser

**Claude**: Quel est ton niveau avec Claude Code ?
- 🟢 Débutant
- 🟡 Intermédiaire
- 🔴 Power User

**You**: Intermédiaire

**Claude**: Combien de temps as-tu ?
- ⚡ 5-10 min
- ⏱️ 15-30 min
- 🎯 30-60 min

**You**: 15-30 min

**Claude**:
*[Loads reference.yaml, looks up `optimize.intermediate_15min`]*

**5 Règles d'Or** (à retenir) :
1. Toujours relire les diffs avant d'accepter
2. `/compact` dès 70% de contexte
3. Être précis : QUOI + OÙ + COMMENT + VÉRIFIER
4. Plan Mode d'abord pour les tâches complexes
5. Créer un CLAUDE.md pour chaque projet

**Ton parcours** (3 topics, ~15 min) :
1. Gestion du contexte
2. Triage du contexte
3. Plan Mode

Par quel sujet veux-tu commencer ?

---

## 6. Tips for Best Results

| Tip | Why It Helps |
|-----|--------------|
| **Be honest about your goal** | Gets you targeted content, not generic tour |
| **Say "deeper" when curious** | Claude will provide more examples |
| **Say "skip" freely** | No need to cover what you know |
| **Ask questions anytime** | Claude will find the relevant section |

---

## 7. Related Resources

- [Reference Index](../machine-readable/reference.yaml) - The navigation map Claude uses
- [Ultimate Guide](../guide/ultimate-guide.md) - Full documentation
- [Cheat Sheet](../guide/cheatsheet.md) - Print this, start coding
- [Setup Audit](./audit-prompt.md) - Analyze your configuration
- [Quiz](../quiz/) - Test your knowledge

---

## 8. Feedback

Found this helpful? Have suggestions?
- Star the repo: [claude-code-ultimate-guide](https://github.com/FlorianBruniaux/claude-code-ultimate-guide)
- Open an issue for improvements
- Share with others learning Claude Code
