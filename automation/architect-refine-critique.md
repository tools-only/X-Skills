---
name: architect-refine-critique
source: https://raw.githubusercontent.com/NTCoding/claude-skillz/main/architect-refine-critique/SKILL.md
original_path: architect-refine-critique/SKILL.md
source_repo: NTCoding/claude-skillz
category: automation
subcategory: scripting
tags: ['automation']
collected_at: 2026-01-31T19:32:40.648712
file_hash: d97ecba8f003911af471a2745dcdc592142775ce81d2dc566b4409b8cebb9aaf
---

---
name: architect-refine-critique
description: "Three-phase design review. Chain architect → refiner → critique subagents. Triggers on: 'design review', 'architecture review', '/arc', system design proposals, significant refactoring decisions, new service or module design."
version: 1.3.6
---

# Architect-Refine-Critique

Chain three subagents sequentially from the main conversation.

## Usage

`/arc [name] [target]`

## Execution

Use the architect subagent to create the initial design for [target] with name=[name],
then use the refiner subagent to improve the design with name=[name],
then use the critique subagent to challenge the design with name=[name].

After all three complete, tell the user: "Run `/arc-review [name]` to discuss findings."

Do not analyze. Do not summarize. Do not add value. Just chain and hand off.
