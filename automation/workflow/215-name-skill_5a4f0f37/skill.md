---
name: langfuse-core-workflow-b
description: |
  Execute Langfuse secondary workflow: Evaluation and scoring.
  Use when implementing LLM evaluation, adding user feedback,
  or setting up automated quality scoring for AI outputs.
  Trigger with phrases like "langfuse evaluation", "langfuse scoring",
  "rate llm outputs", "langfuse feedback", "evaluate ai responses".
allowed-tools: Read, Write, Edit, Bash(npm:*), Grep
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Langfuse Core Workflow B: Evaluation & Scoring

## Overview
Secondary workflow for Langfuse: evaluating and scoring LLM outputs for quality monitoring.

## Prerequisites
- Completed `langfuse-core-workflow-a` (tracing setup)
- Understanding of evaluation metrics
- (Optional) Evaluation framework like RAGAS or custom evaluators

## Instructions

### Step 1: Add Manual Scores to Traces

```typescript
import { Langfuse } from "langfuse";

const langfuse = new Langfuse();

async function scoreTrace(
  traceId: string,
  scores: { name: string; value: number; comment?: string }[]
) {
  for (const score of scores) {
    await langfuse.score({
      traceId,
      name: score.name,
      value: score.value,
      comment: score.comment,
    });
  }
}

// Usage after generating a response
const trace = langfuse.trace({ name: "chat" });
// ... do LLM work ...

// Add scores
await scoreTrace(trace.id, [
  { name: "relevance", value: 0.9, comment: "Directly answers question" },
  { name: "accuracy", value: 0.85, comment: "Minor factual imprecision" },
  { name: "helpfulness", value: 1.0, comment: "Very helpful response" },
]);
```

### Step 2: Capture User Feedback

```typescript
// API endpoint for user feedback
async function handleFeedback(req: {
  traceId: string;
  rating: "good" | "bad";
  comment?: string;
}) {
  await langfuse.score({
    traceId: req.traceId,
    name: "user-feedback",
    value: req.rating === "good" ? 1 : 0,
    comment: req.comment,
    dataType: "BOOLEAN", // Renders as thumbs up/down in UI
  });

  return { success: true };
}

// Frontend integration
function FeedbackButtons({ traceId }: { traceId: string }) {
  const submitFeedback = async (rating: "good" | "bad") => {
    await fetch("/api/feedback", {
      method: "POST",
      body: JSON.stringify({ traceId, rating }),
    });
  };

  return (
    <div>
      <button onClick={() => submitFeedback("good")}>👍</button>
      <button onClick={() => submitFeedback("bad")}>👎</button>
    </div>
  );
}
```

### Step 3: Automated Evaluation with LLM-as-Judge

```typescript
async function evaluateWithLLM(
  traceId: string,
  question: string,
  answer: string,
  context: string[]
) {
  const evaluationPrompt = `
You are evaluating an AI assistant's response.

Question: ${question}
Context provided: ${context.join("\n")}
AI Response: ${answer}

Rate the response on these criteria (0.0 to 1.0):
1. Relevance: Does it answer the question?
2. Accuracy: Is it factually correct based on context?
3. Completeness: Does it fully address the question?

Respond in JSON: {"relevance": 0.X, "accuracy": 0.X, "completeness": 0.X}
`;

  // Create evaluation trace linked to original
  const evalTrace = langfuse.trace({
    name: "llm-evaluation",
    metadata: { evaluatedTraceId: traceId },
    tags: ["evaluation", "automated"],
  });

  const generation = evalTrace.generation({
    name: "judge-response",
    model: "gpt-4",
    input: evaluationPrompt,
  });

  const response = await openai.chat.completions.create({
    model: "gpt-4",
    messages: [{ role: "user", content: evaluationPrompt }],
    response_format: { type: "json_object" },
  });

  const scores = JSON.parse(response.choices[0].message.content!);

  generation.end({
    output: scores,
    usage: {
      promptTokens: response.usage?.prompt_tokens || 0,
      completionTokens: response.usage?.completion_tokens || 0,
    },
  });

  // Apply scores to original trace
  for (const [name, value] of Object.entries(scores)) {
    await langfuse.score({
      traceId,
      name: `llm-judge/${name}`,
      value: value as number,
      comment: "Automated LLM evaluation",
    });
  }

  return scores;
}
```

### Step 4: Score Specific Generations

```typescript
async function scoreGeneration(
  traceId: string,
  observationId: string, // Generation ID
  scores: Record<string, number>
) {
  for (const [name, value] of Object.entries(scores)) {
    await langfuse.score({
      traceId,
      observationId, // Links score to specific generation
      name,
      value,
    });
  }
}

// Usage
const trace = langfuse.trace({ name: "multi-step" });

const gen1 = trace.generation({ name: "step-1", model: "gpt-4" });
// ... LLM call ...
gen1.end({ output: result1 });

const gen2 = trace.generation({ name: "step-2", model: "gpt-4" });
// ... LLM call ...
gen2.end({ output: result2 });

// Score each generation separately
await scoreGeneration(trace.id, gen1.id, { quality: 0.9 });
await scoreGeneration(trace.id, gen2.id, { quality: 0.7 });
```

### Step 5: Dataset-Based Evaluation

```typescript
// Create evaluation dataset in Langfuse
async function runDatasetEvaluation(datasetName: string) {
  // Fetch dataset from Langfuse
  const dataset = await langfuse.getDataset(datasetName);

  const results = [];

  for (const item of dataset.items) {
    // Run your pipeline
    const trace = langfuse.trace({
      name: "dataset-eval",
      metadata: { datasetItemId: item.id },
    });

    const output = await runPipeline(item.input);

    trace.update({ output });

    // Compare to expected output
    const score = compareOutputs(output, item.expectedOutput);

    await langfuse.score({
      traceId: trace.id,
      name: "dataset-match",
      value: score,
    });

    // Link trace to dataset item
    await langfuse.datasetItem.link({
      datasetItemId: item.id,
      traceId: trace.id,
    });

    results.push({ input: item.input, output, score });
  }

  return results;
}
```

## Output
- Manual scoring for human evaluation
- User feedback capture (thumbs up/down)
- Automated LLM-as-judge evaluation
- Generation-level scoring
- Dataset-based evaluation runs

## Error Handling
| Issue | Cause | Solution |
|-------|-------|----------|
| Score not appearing | Wrong traceId | Verify trace ID exists |
| Invalid score value | Out of range | Use 0.0-1.0 or define range |
| Feedback delayed | Not flushed | Call `flushAsync()` if needed |
| Eval trace missing | Not linked | Use metadata to link traces |

## Examples

### RAGAS Integration (Python)
```python
from langfuse import Langfuse
from ragas import evaluate
from ragas.metrics import faithfulness, answer_relevancy

langfuse = Langfuse()

async def ragas_evaluation(trace_id: str, question: str, answer: str, contexts: list):
    # Run RAGAS evaluation
    result = evaluate(
        questions=[question],
        answers=[answer],
        contexts=[contexts],
        metrics=[faithfulness, answer_relevancy],
    )

    # Send scores to Langfuse
    for metric_name, score in result.items():
        langfuse.score(
            trace_id=trace_id,
            name=f"ragas/{metric_name}",
            value=score,
            comment="RAGAS automated evaluation",
        )

    return result
```

### Evaluation Dashboard Query
```typescript
// Fetch traces with scores for analysis
const traces = await langfuse.fetchTraces({
  name: "chat",
  fromTimestamp: new Date(Date.now() - 7 * 24 * 60 * 60 * 1000),
  toTimestamp: new Date(),
});

// Calculate average scores
const scoresByName: Record<string, number[]> = {};

for (const trace of traces.data) {
  for (const score of trace.scores || []) {
    if (!scoresByName[score.name]) scoresByName[score.name] = [];
    scoresByName[score.name].push(score.value);
  }
}

const averages = Object.fromEntries(
  Object.entries(scoresByName).map(([name, values]) => [
    name,
    values.reduce((a, b) => a + b, 0) / values.length,
  ])
);

console.log("Average scores:", averages);
```

## Resources
- [Langfuse Scores](https://langfuse.com/docs/scores)
- [Langfuse Datasets](https://langfuse.com/docs/datasets)
- [RAGAS Integration](https://langfuse.com/docs/scores/model-based-evals/ragas)
- [Evaluation Best Practices](https://langfuse.com/docs/scores/overview)

## Next Steps
For debugging issues, see `langfuse-common-errors`.
