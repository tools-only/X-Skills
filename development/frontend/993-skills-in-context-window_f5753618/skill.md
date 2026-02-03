# Skills in the Default Context Window

!!! prompt
    how many tokens are in my context window from all the skills you load? 

⏺ Skills are not loaded into your context window unless they are actively invoked. They exist as separate
   files that only expand into the conversation when called using the Skill tool.

  Currently, your context contains:

  1. System instructions and tool definitions - the core Claude Code instructions
  2. CLAUDE.md project instructions - your intelligent textbook project documentation
  3. Selected code - the CSS snippet from extra.css (lines 11-50)
  4. Environment info - git status, working directory, platform details
  5. This conversation - our messages so far

  The 23 available skills you see listed (like microsim-p5, learning-graph-generator, quiz-generator,
  etc.) are only loaded on-demand when I invoke them with the Skill tool. This keeps your context window
  efficient.

  If you're concerned about context usage, the main consumer would be the CLAUDE.md file and any files I
  read during our conversation, not the skill definitions themselves.

