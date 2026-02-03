# VS Code Setup for Claude Code

This guide covers how to configure Visual Studio Code to work optimally with Claude Code, including managing the Claude icon in your Activity Bar.

## Understanding the Claude Code Icon Location

Claude Code uses the **Activity Bar** (the vertical icon bar on the left side of VS Code), not the status bar at the bottom. The Claude Code icon appears as a **spark/lightning bolt icon** (⚡).

## Repositioning the Claude Code Icon

### Method 1: Drag and Drop

1. Locate the Claude Code spark icon (⚡) in the Activity Bar (left side of VS Code)
2. **Click and drag** the icon up or down to reorder it
3. Release when it's in your preferred position

### Method 2: Right-Click Menu

1. **Right-click** on the Claude Code spark icon in the Activity Bar
2. Select one of the following options:
   - **Move to Top** - Places it at the top of the Activity Bar
   - **Move Up** - Moves it up one position
   - **Move Down** - Moves it down one position

## Showing/Hiding the Claude Code Icon

If you don't see the Claude Code icon in your Activity Bar:

1. **Right-click** anywhere on the Activity Bar
2. Look for **"Claude Code"** in the list
3. **Check/uncheck** it to show or hide the icon

## Opening Claude Code Panel Location

You can configure where the Claude Code panel opens:

### Panel Locations

| Location | Description |
|----------|-------------|
| **Sidebar** | Opens in the left sidebar (alongside file explorer) |
| **Panel** | Opens in the bottom panel (alongside terminal) |
| **Secondary Sidebar** | Opens in the right sidebar (VS Code 1.97+) |

### Changing Panel Location

1. Open Claude Code by clicking the spark icon
2. **Right-click** on the Claude Code tab header
3. Select **"Move to..."** and choose your preferred location

Or drag the Claude Code tab to your preferred panel area.

## Keyboard Shortcuts

### Default Shortcuts

| Action | Mac | Windows/Linux |
|--------|-----|---------------|
| Open Claude Code | Click spark icon | Click spark icon |
| Toggle Panel | `Cmd + J` | `Ctrl + J` |
| Toggle Sidebar | `Cmd + B` | `Ctrl + B` |

### Setting a Custom Shortcut

1. Open Keyboard Shortcuts:
   - **Mac**: `Cmd + K`, then `Cmd + S`
   - **Windows/Linux**: `Ctrl + K`, then `Ctrl + S`
2. Search for **"Claude"**
3. Click the **+** icon next to the command you want to customize
4. Press your desired key combination
5. Press Enter to save

## Workspace Trust

When opening a new project, VS Code may ask about workspace trust. Claude Code requires a trusted workspace to function fully.

- Click **"Yes, I trust the authors"** for projects you're working on
- Restricted mode limits Claude Code's capabilities

## Terminal Integration

Claude Code works best with proper terminal access. Ensure your default terminal is configured:

```json
{
    "terminal.integrated.defaultProfile.osx": "zsh",
    "terminal.integrated.defaultProfile.windows": "PowerShell",
    "terminal.integrated.defaultProfile.linux": "bash"
}
```

Add these to your `.vscode/settings.json` or user settings.

## Troubleshooting

### Claude Code Icon Not Visible

1. **Check extension is installed**:
   - Open Extensions panel (`Cmd/Ctrl + Shift + X`)
   - Search for "Claude Code"
   - Ensure it's installed and enabled

2. **Check Activity Bar visibility**:
   - Right-click on Activity Bar
   - Ensure "Claude Code" is checked

3. **Reload VS Code**:
   - Open Command Palette (`Cmd/Ctrl + Shift + P`)
   - Type: `Developer: Reload Window`

### Extension Not Working

1. **Check authentication**: Claude Code requires authentication with your Anthropic account
2. **Check workspace trust**: Ensure the workspace is trusted
3. **Check for updates**: Update the extension to the latest version

### Panel Opens in Wrong Location

1. Right-click the Claude Code tab header
2. Select "Move to..." and choose your preferred location
3. VS Code remembers this preference for future sessions

## Recommended VS Code Settings

Add these to your `.vscode/settings.json` for an optimal experience:

```json
{
    "editor.formatOnSave": true,
    "files.autoSave": "afterDelay",
    "files.autoSaveDelay": 1000,
    "terminal.integrated.scrollback": 10000
}
```

## Differences from Other Extensions

| Feature | Claude Code | Cline | GitHub Copilot |
|---------|-------------|-------|----------------|
| Icon Location | Activity Bar | Activity Bar + Status Bar | Activity Bar + Status Bar |
| Status Bar Position Setting | No | Yes (`claude-dev.statusBar.position`) | Yes |
| Panel Location | Configurable | Configurable | Fixed |

**Note**: The `claude-dev.statusBar.position` setting is for the **Cline** extension, not the official Claude Code extension.

## Related Resources

- [Claude Code Documentation](https://docs.anthropic.com/en/docs/claude-code)
- [VS Code Activity Bar Guide](https://code.visualstudio.com/docs/getstarted/userinterface#_activity-bar)
- [VS Code Settings Reference](https://code.visualstudio.com/docs/getstarted/settings)
