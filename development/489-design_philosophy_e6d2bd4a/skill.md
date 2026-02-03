# User Interface Design Philosophy

## The NeXT Step in Uniformity

The Tunacode TUI (Terminal User Interface) is built upon a strong design philosophy heavily inspired by the classic **NeXTSTEP** operating system. This is not merely an aesthetic choice, but a functional one aimed at maximizing user clarity and control.

### Core Principles

1.  **Uniformity**: Every interaction should feel consistent. Whether you are configuring a model, reviewing a code diff, or reading an error log, the interface behaves in a predictable manner. The "widget" styling, keybindings, and layout logic are shared across the entire application.

2.  **Transparency (Keep the User Informed)**:
    - **No Magic**: In the era of AI, "magic" is often a synonym for "unpredictable behavior". Tunacode strives to show you exactly what is happening.
    - **State Visibility**: The user should always know the agent's current state (Thinking, Coding, Waiting, Error).
    - **Visual Feedback**: Every action—from a file write to a network request—should have a corresponding visual indicator in the UI.

3.  **Aesthetic Functionality**:
    - The retro-modern look pays homage to the clean, object-oriented design of NeXTSTEP.
    - High contrast and distinct borders help separate information density in a terminal environment.
    - The use of the `rich` and `textual` libraries ensures that despite the retro feel, the terminal capabilities are modern (true color, mouse support, resizing).

## Implementation Details

The UI is implemented using [Textual](https://textual.textualize.io/), a Python framework for rapid application development in the terminal.

- **Widgets**: Custom widgets are designed to mimic the "chunky" and distinct look of early GUI systems.
- **Layout**: We prioritize information hierarchy. The most important information (the code or the conversation) takes center stage, while status and tools are peripheral but always visible.
