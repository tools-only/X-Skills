"""X-Skills CLI - Main entry point

Install with: pip install xskills
Or: pip install git+https://github.com/tools-only/X-Skills.git
"""

import os
import sys
import subprocess
from pathlib import Path


def get_skillflow_dir():
    """Get or install SkillFlow directory"""
    skillflow_dir = Path.home() / ".skillflow"

    if not (skillflow_dir / "src" / "xskills_cli_new.py").exists():
        # Auto-install SkillFlow
        print("📦 Installing SkillFlow...")
        try:
            subprocess.run(
                ["git", "clone", "--depth", "1",
                 "https://github.com/tools-only/SkillFlow.git",
                 str(skillflow_dir)],
                check=True
            )
            print("✓ SkillFlow installed")
        except subprocess.CalledProcessError as e:
            print(f"⚠️  Installation failed: {e}")
            print("Please install manually:")
            print("   git clone https://github.com/tools-only/SkillFlow.git ~/.skillflow")
            sys.exit(1)
        except FileNotFoundError:
            print("⚠️  git is not installed")
            print("Please install git first, then run:")
            print("   git clone https://github.com/tools-only/SkillFlow.git ~/.skillflow")
            sys.exit(1)

    return skillflow_dir


def main():
    """Main CLI entry point"""
    skillflow_dir = get_skillflow_dir()

    # Set up Python path
    sys.path.insert(0, str(skillflow_dir))
    os.chdir(skillflow_dir)

    # Import and run
    from src.xskills_cli_new import xskills_cli
    xskills_cli()


if __name__ == "__main__":
    main()
