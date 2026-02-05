"""X-Skills Plugin System - Setup Script"""

from setuptools import setup, find_packages
import os

# Read README
with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="xskills",
    version="1.0.0",
    description="X-Skills Plugin System for Claude Code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="X-Skills Team",
    url="https://github.com/tools-only/X-Skills",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "click>=8.1.0",
        "rich>=13.0",
        "pyyaml>=6.0",
        "fuzzywuzzy>=0.18.0",
        "python-Levenshtein>=0.21.0",
    ],
    entry_points={
        "console_scripts": [
            "xsk=xskills.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.8",
)
