#!/bin/bash
# X-Skills One-Line Installer
# Install with: curl -fsSL https://raw.githubusercontent.com/tools-only/X-Skills/main/install.sh | bash

set -e

echo "🚀 Installing X-Skills..."

# Detect Python
PYTHON_CMD="python3"
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is required"
    exit 1
fi

# Install via pip
echo "📦 Installing xskills package..."
pip3 install --user git+https://github.com/tools-only/X-Skills.git@main 2>/dev/null || \
    pip3 install --user -e .

# Add to PATH
BASH_RC="$HOME/.bashrc"
ZSH_RC="$HOME/.zshrc"

if [ -n "$ZSH_VERSION" ]; then
    RC_FILE="$ZSH_RC"
else
    RC_FILE="$BASH_RC"
fi

# Check if PATH already includes ~/.local/bin
if ! grep -q '~/.local/bin' "$RC_FILE" 2>/dev/null; then
    echo ""
    echo "export PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$RC_FILE"
    echo "✓ Added to PATH in $RC_FILE"
fi

echo ""
echo "✅ Installation complete!"
echo ""
echo "Start a new terminal or run: source $RC_FILE"
echo ""
echo "Usage:"
echo "  xsk scenarios              # List scenarios"
echo "  xsk install @research      # Install research scenario"
echo "  xsk list                  # List installed"
