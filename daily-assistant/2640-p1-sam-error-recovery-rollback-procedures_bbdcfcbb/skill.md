---
name: 'SAM: Error Recovery / Rollback Procedures'
description: Define explicit procedure when a task fails irrecoverably. How to undo artifact changes? How to restore artifact plane to consistent state after failure?
metadata:
  topic: sam-error-recovery-rollback-procedures
  source: Gap analysis of SAM framework
  added: '2026-02-01'
  priority: completed
  type: Feature
  status: done
  issue: '#85'
  plan: ''
---

**Suggested location**: [`stateless-software-engineering-framework.md`](https://github.com/bitflight-devops/stateless-agent-methodology/blob/main/stateless-software-engineering-framework.md) (new Appendix or Part 6 addition)

**Research first**: How do GSD, BMAD-METHOD, AutoGPT, and traditional CI/CD handle rollback? What patterns exist for transactional artifact updates?