---
name: speak-hello-world
description: |
  Create a minimal working Speak language learning example.
  Use when starting a new Speak integration, testing your setup,
  or learning basic Speak API patterns for language tutoring.
  Trigger with phrases like "speak hello world", "speak example",
  "speak quick start", "simple speak lesson".
allowed-tools: Read, Write, Edit
version: 1.0.0
license: MIT
author: Jeremy Longshore <jeremy@intentsolutions.io>
---

# Speak Hello World

## Overview
Minimal working example demonstrating core Speak functionality for AI-powered language learning.

## Prerequisites
- Completed `speak-install-auth` setup
- Valid API credentials configured
- Development environment ready
- Microphone access for speech input (optional for testing)

## Instructions

### Step 1: Create Entry File
Create a new file for your hello world lesson.

### Step 2: Import and Initialize Client
```typescript
// src/speak-hello-world.ts
import { SpeakClient, LessonSession, AITutor } from '@speak/language-sdk';

const client = new SpeakClient({
  apiKey: process.env.SPEAK_API_KEY!,
  appId: process.env.SPEAK_APP_ID!,
  language: 'es', // Learning Spanish
});
```

### Step 3: Create Your First Lesson Session
```typescript
async function startFirstLesson() {
  // Initialize AI tutor
  const tutor = new AITutor(client, {
    targetLanguage: 'es',
    nativeLanguage: 'en',
    proficiencyLevel: 'beginner',
  });

  // Start a conversation practice session
  const session = await tutor.startSession({
    topic: 'greetings',
    duration: 5, // 5 minutes
  });

  // Get AI tutor's opening prompt
  const opening = await session.getPrompt();
  console.log('AI Tutor says:', opening.text);
  console.log('Audio URL:', opening.audioUrl);

  // Simulate user response (in production, use speech recognition)
  const feedback = await session.submitResponse({
    text: 'Hola, me llamo Juan',
    audioData: null, // Optional: audio buffer for pronunciation scoring
  });

  console.log('Feedback:', feedback.message);
  console.log('Pronunciation:', feedback.pronunciationScore);
  console.log('Grammar:', feedback.grammarCorrections);

  return session;
}

startFirstLesson().catch(console.error);
```

### Step 4: Run Your First Lesson
```bash
npx tsx src/speak-hello-world.ts
```

## Output
- Working code file with Speak client initialization
- Successful lesson session creation
- AI tutor prompt and feedback displayed
- Console output showing:
```
AI Tutor says: Hola! Welcome to your Spanish lesson. Let's practice greetings. Can you say "Hello, my name is..." in Spanish?
Audio URL: https://speak.com/audio/abc123.mp3
Feedback: Great job! Your pronunciation was clear.
Pronunciation: 85
Grammar: [No corrections needed]
```

## Error Handling
| Error | Cause | Solution |
|-------|-------|----------|
| Import Error | SDK not installed | Verify with `npm list @speak/language-sdk` |
| Auth Error | Invalid credentials | Check environment variables are set |
| Session Timeout | Exceeded lesson duration | Extend session or start new one |
| Rate Limit | Too many requests | Wait and retry with exponential backoff |
| Language Not Supported | Invalid language code | Use supported language codes |

## Examples

### TypeScript: Pronunciation Practice
```typescript
import { SpeakClient, PronunciationPractice } from '@speak/language-sdk';

const client = new SpeakClient({
  apiKey: process.env.SPEAK_API_KEY!,
  appId: process.env.SPEAK_APP_ID!,
  language: 'ko',
});

async function practicePronunciation() {
  const practice = new PronunciationPractice(client);

  // Get a phrase to practice
  const phrase = await practice.getPhrase({
    difficulty: 'easy',
    category: 'daily_conversation',
  });

  console.log('Practice this phrase:');
  console.log('Korean:', phrase.text);
  console.log('Romanization:', phrase.romanization);
  console.log('English:', phrase.translation);
  console.log('Listen:', phrase.audioUrl);
}

practicePronunciation();
```

### Python: Basic Lesson
```python
from speak_sdk import SpeakClient, AITutor

client = SpeakClient(
    api_key=os.environ.get('SPEAK_API_KEY'),
    app_id=os.environ.get('SPEAK_APP_ID'),
    language='ja'
)

async def first_lesson():
    tutor = AITutor(client, target_language='ja', native_language='en')
    session = await tutor.start_session(topic='self_introduction')

    prompt = await session.get_prompt()
    print(f"Tutor: {prompt.text}")

    # Submit response
    feedback = await session.submit_response(text="Hajimemashite, watashi wa John desu")
    print(f"Feedback: {feedback.message}")

asyncio.run(first_lesson())
```

### Vocabulary Quiz Example
```typescript
async function vocabularyQuiz() {
  const quiz = await client.vocabulary.createQuiz({
    language: 'es',
    topic: 'food',
    questionCount: 5,
    type: 'multiple_choice',
  });

  for (const question of quiz.questions) {
    console.log(`Q: ${question.prompt}`);
    console.log(`Options: ${question.options.join(', ')}`);

    // In production, get user input
    const answer = await quiz.submitAnswer(question.id, 'answer_a');
    console.log(`Correct: ${answer.correct}`);
  }

  console.log(`Final Score: ${quiz.score}/${quiz.total}`);
}
```

## Resources
- [Speak Getting Started](https://developer.speak.com/docs/getting-started)
- [Speak API Reference](https://developer.speak.com/api)
- [Speak Lesson Types](https://developer.speak.com/docs/lesson-types)
- [Speech Recognition Guide](https://developer.speak.com/docs/speech-recognition)

## Next Steps
Proceed to `speak-local-dev-loop` for development workflow setup.
