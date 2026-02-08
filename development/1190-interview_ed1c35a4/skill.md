---
description: Interview user to clarify plan details and resolve ambiguities
---

Use ultrathink to thoroughly consider all aspects of this task.

## Overview

Before implementing a feature or making significant changes, use this command to ensure complete understanding of requirements. This prevents wasted effort from misunderstood requirements and surfaces edge cases early.

## When to Use

Invoke `/interview` when:
- Starting a new feature with multiple possible approaches
- The task description is vague or leaves room for interpretation
- You're unsure about edge cases or error handling expectations
- Integration points with existing systems are unclear
- The feature has user-facing components where UX matters

## Interview Process

### Round 1: Core Requirements

Start by clarifying the fundamental aspects:

**Functionality**
- What is the primary goal of this feature?
- What inputs does it accept?
- What outputs should it produce?
- Are there any existing examples or mockups to reference?

**Scope**
- What is explicitly in scope?
- What is explicitly out of scope?
- Are there related features that should wait for later?

### Round 2: Edge Cases & Error Handling

Probe for scenarios that may not have been considered:

**Invalid Input**
- What happens with empty input?
- What happens with invalid data types?
- What happens with extremely large inputs?
- What happens with special characters or unicode?

**Error States**
- How should errors be displayed to users?
- Should errors be logged? Where?
- Are there retry mechanisms needed?
- What happens if external services are unavailable?

**Concurrency**
- What if multiple users do this simultaneously?
- What if the same user triggers it twice quickly?
- Are there race conditions to consider?

### Round 3: Integration & Dependencies

Understand how this fits with existing systems:

**Data Flow**
- Where does input data come from?
- Where does output data go?
- What existing APIs or services are involved?
- Are there database schema changes needed?

**Authentication & Authorization**
- Who is allowed to use this feature?
- Are there different permission levels?
- How is access controlled?

**Third-party Services**
- Are external APIs involved?
- What are the rate limits?
- How are API keys managed?

### Round 4: Performance & Scale

Clarify non-functional requirements:

**Performance**
- What response time is acceptable?
- Are there specific SLA requirements?
- Should results be cached?

**Scale**
- How many users will use this?
- What's the expected data volume?
- Are there peak usage times to consider?

**Resource Constraints**
- Are there memory limits?
- Are there storage limits?
- Are there cost constraints?

### Round 5: User Experience

For user-facing features:

**Interaction Design**
- What does the UI look like?
- What feedback should users see during processing?
- How should success/failure be communicated?

**Accessibility**
- Are there accessibility requirements?
- Should it work on mobile?
- Are there internationalization needs?

### Round 6: Testing & Validation

Define what "done" looks like:

**Test Criteria**
- What test cases are essential?
- What does a successful test look like?
- Are there existing test patterns to follow?

**Acceptance Criteria**
- Who will verify this works correctly?
- What environments should it be tested in?
- Are there specific browsers or devices to support?

### Round 7: Deployment & Rollout

Plan the release:

**Deployment**
- Is this a feature flag rollout?
- Are database migrations needed?
- Is there a rollback plan?

**Documentation**
- Does user documentation need updating?
- Are there API docs to write?
- Should there be release notes?

## Question Patterns

Use the AskUserQuestion tool with well-structured questions:

**For choices with clear options:**
```
header: "Error handling"
question: "How should validation errors be displayed to users?"
options:
  - Inline below each field (Recommended)
  - Summary at top of form
  - Toast notification
```

**For open-ended requirements:**
```
header: "Performance"
question: "What response time is acceptable for this operation?"
options:
  - Under 100ms (fast, may need caching)
  - Under 1 second (standard)
  - Under 5 seconds (complex operation)
```

**For edge cases:**
```
header: "Edge cases"
question: "What should happen if the user submits while offline?"
options:
  - Queue for later submission
  - Show error immediately
  - Disable submit button when offline
```

## Best Practices

1. **Ask in batches** - Group related questions to avoid question fatigue
2. **Offer recommendations** - Put the best default first with "(Recommended)"
3. **Be specific** - Vague questions get vague answers
4. **Use multiSelect** - When options aren't mutually exclusive
5. **Build on answers** - Use previous answers to inform follow-up questions
6. **Know when to stop** - Don't over-interview; move to implementation when clear

## Output

After completing the interview, summarize:
1. Confirmed requirements
2. Decisions made
3. Remaining unknowns (if any)
4. Recommended approach based on answers

Do not proceed with implementation until the plan is fully clarified and the user confirms understanding.
