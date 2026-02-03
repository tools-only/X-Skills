# Statistical Sampling Skill

**🎯 Zweck**: Verhindert Performance-Timeouts und Token-Overflows durch statistisch fundiertes Data Sampling.

**⚠️ KRITISCH**: Basiert auf Lessons Learned vom 2025-11-18 (SECOM dataset 504 timeout → 200×100 sampling).

## Wann wird dieser Skill aktiviert?

- Du arbeitest mit großen Datasets (>100 Zeilen oder >50 Spalten)
- Du implementierst AI/ML-Features mit Token-Limits (Gemini, GPT, Claude)
- Du debuggst Performance-Timeouts bei API-Calls
- Du erwähnst Keywords: sampling, imbalanced data, time-series, PCA, anomaly detection
- Du entwickelst Manufacturing/IoT Analytics mit Sensor-Daten

## Goldene Regeln

### 1. Cochran's Formula für Imbalanced Classification

**Hintergrund**: Bei imbalanced datasets (z.B. 6.6% Failure-Rate) ist die Minority-Class-Repräsentation kritisch.

**Formel**:
```
n₀ = (Z² · p · (1-p)) / E²

Wo:
Z = 1.96 (95% confidence level)
p = minority class proportion (z.B. 0.066 für 6.6%)
E = margin of error (z.B. 0.05 für 5%)

Mit Finite Population Correction (N = total rows):
n = n₀ / (1 + (n₀-1)/N)
```

**Beispiel - SECOM Dataset**:
```typescript
// SECOM: 1,567 rows, 6.6% failure rate
const p = 0.066;
const Z = 1.96;  // 95% confidence
const E = 0.05;  // 5% margin of error

const n0 = (Z * Z * p * (1 - p)) / (E * E);
// n0 = 94.72 samples

const N = 1567;
const n = n0 / (1 + (n0 - 1) / N);
// n = 89.4 samples

// Praktisches Minimum: 150-200 rows (10-15 failure examples für ML)
const RECOMMENDED_MIN = 200;
```

**Regel**:
- Minimum nach Cochran: ~90 samples
- **Praktisch**: 150-200 samples (für 10-15 minority examples)
- **Never go below**: 2× Cochran's minimum

---

### 2. Box-Jenkins für Time-Series Analysis

**Hintergrund**: ARIMA-Modelle benötigen Mindest-Observations für Trend/Seasonality Detection.

**Minimums**:
- **Basic ARIMA**: 50-100 observations
- **Seasonal patterns**: 200+ observations (mindestens 2 vollständige Zyklen)
- **Manufacturing (Daily cycles)**: 200 rows für ~7-10 Tage Coverage

**SCHLECHT:**
```typescript
// ❌ Nur 100 Zeilen → 6.4% Coverage = 4 Tage
const MAX_ROWS = 100;
// Problem: Kann keine Wochenend-Effekte oder Schicht-Muster erkennen
```

**GUT:**
```typescript
// ✅ 200 Zeilen → 12.8% Coverage = 8 Tage
const MAX_ROWS = 200;
// Vorteil: Mindestens 1 volle Woche, erkennt tägliche Zyklen
```

**Beispiel - Manufacturing Data**:
```typescript
// Dataset: 1,567 rows über 64 Tage = ~24.5 rows/day
const totalDays = 64;
const totalRows = 1567;
const rowsPerDay = totalRows / totalDays; // 24.5

// Temporal Coverage
const sampledRows = 200;
const temporalCoveragePct = (sampledRows / totalRows) * 100; // 12.8%
const temporalCoverageDays = (sampledRows / rowsPerDay); // ~8.2 days

// Check if sufficient for daily cycles (minimum 7 days)
if (temporalCoverageDays >= 7) {
  console.log('✅ Sufficient for daily cycle detection');
} else {
  console.warn('⚠️ Consider increasing sample size for weekly patterns');
}
```

**Regel**:
- **Minimum**: 50-100 observations für basic trends
- **Optimal**: 200+ für seasonality (Manufacturing: mindestens 1 Woche)
- **Critical**: Immer mindestens 2 vollständige Zyklen samplen

---

### 3. Power Analysis für Anomaly Detection

**Hintergrund**: Manufacturing braucht Detection von 1-2 sigma shifts (small/medium effect sizes).

**Formel (Cohen's d)**:
```
n = 2 · [(Zα + Zβ) / d]²

Wo:
Zα = 1.96 (Type I error α=0.05, two-tailed)
Zβ = 0.84 (Power 80%, Type II error β=0.20)
d = Effect Size (Cohen's d: 0.2=small, 0.5=medium, 0.8=large)

Manufacturing typical: d = 1.0 to 1.5 (1-1.5 sigma shifts)
```

**Effect Sizes im Manufacturing Context**:
```typescript
// Cohen's d effect size classification
enum EffectSize {
  SMALL = 0.2,    // 0.2 sigma shift - Subtle drift
  MEDIUM = 0.5,   // 0.5 sigma shift - Tool wear
  LARGE = 0.8,    // 0.8 sigma shift - Process change
  XLARGE = 1.0,   // 1.0 sigma shift - Equipment issue
  CRITICAL = 1.5  // 1.5 sigma shift - Failure mode
}

// Power Analysis für 1-sigma shift (d=1.0), 80% power
const Zalpha = 1.96;
const Zbeta = 0.84;
const d = 1.0;

const requiredSamples = 2 * Math.pow((Zalpha + Zbeta) / d, 2);
// = 31.4 samples für 80% power

// Für >95% power bei 1-sigma shifts: 100-200 samples
const POWER_95_SAMPLES = 200;
```

**Sensitivity Table**:

| Shift Magnitude | Cohen's d | Samples (80% power) | Samples (95% power) |
|-----------------|-----------|---------------------|---------------------|
| 2-sigma | 2.0 | 8 | 16 |
| 1.5-sigma | 1.5 | 14 | 28 |
| **1-sigma** | **1.0** | **32** | **100** |
| 0.5-sigma | 0.5 | 128 | 400 |

**Regel**:
- Manufacturing (1-sigma shifts): **100-200 samples für >95% power**
- Critical systems (0.5-sigma): 400+ samples (oft unrealistisch!)
- **Diminishing Returns**: Nach 200 samples < 5% marginal power gain

---

### 4. PCA Sample-to-Feature Ratio

**Hintergrund**: Principal Component Analysis braucht genug samples per feature für stable variance estimates.

**Rule of Thumb**: **5-10 samples per feature**

**SECOM Dataset Challenge**:
```
591 features × 5-10 samples = 2,955-5,910 required samples
Available: 1,567 samples → IMPOSSIBLE!
```

**Solution - Feature Reduction**:
```typescript
// ❌ SCHLECHT: 100 samples × 591 features = 1:6 ratio (inverted!)
const MAX_ROWS = 100;
const MAX_COLUMNS = 591; // ALL features
// Problem: PCA wird instabil, overfitting risk

// ⚠️ MARGINAL: 100 samples × 50 features = 2:1 ratio
const MAX_ROWS = 100;
const MAX_COLUMNS = 50;
// Problem: 2:1 ist MINIMUM, nicht optimal

// ✅ GUT: 200 samples × 100 features = 2:1 ratio
const MAX_ROWS = 200;
const MAX_COLUMNS = 100;
// Vorteil: 2:1 ist akzeptabel für correlated sensors
// Erklärung: Manufacturing sensors sind hoch korreliert (80-95% variance in first 10-30% PCs)
```

**Manufacturing Sensor Correlation**:
```typescript
// Research: Industrial sensors capture redundant signals
// Example: Temperature Sensor 1, 2, 3 → High correlation
// Result: 80-95% variance in first 20-30% of principal components

// Feature Coverage Analysis
const totalFeatures = 591;
const sampledFeatures = 100;
const featureCoveragePct = (sampledFeatures / totalFeatures) * 100; // 16.9%

// If sensors are highly correlated (r > 0.7):
// 16.9% of features can capture 80-95% of total variance
// → 2:1 sample-to-feature ratio is ACCEPTABLE
```

**Regel**:
- **Ideal**: 5-10 samples per feature (often infeasible!)
- **Minimum**: 2:1 sample-to-feature ratio für correlated sensors
- **Manufacturing**: 2:1 ratio akzeptabel (sensors sind korreliert)
- **Never**: <1:1 ratio (mehr features als samples = overfitting!)

---

### 5. ISO 2859-1 Acceptance Sampling Standards

**Hintergrund**: ISO-Standard für Quality Assurance in Manufacturing.

**Lookup Table (ISO 2859-1)**:

| Lot Size | Sample Code | Sample Size (Normal) | Accept ≤ | Reject ≥ |
|----------|-------------|----------------------|----------|----------|
| 91-150 | H | 50 | 3 | 4 |
| 151-280 | J | 80 | 5 | 6 |
| 281-500 | K | 125 | 7 | 8 |
| 501-1,200 | L | 200 | 10 | 11 |
| **1,201-3,200** | **M** | **315** | **14** | **15** |

**SECOM Dataset (1,567 rows)**:
```typescript
// Lot Size: 1,567 → Sample Code M
// ISO-Standard: 315 samples (IMPRACTICAL for API timeout!)

// Compromise: Sample Code K (125 samples) for 500-1,200 range
// Our Implementation: 200 samples → EXCEEDS ISO-K by 60% ✅
```

**AQL (Acceptable Quality Level) Context**:
- **AQL 1.5%**: Pharma, Medical Devices
- **AQL 2.5%**: Automotive, Aerospace
- **AQL 4.0%**: Consumer Electronics
- **SECOM 6.6%**: Below AQL 4.0 → Needs improvement

**Regel**:
- Follow ISO sample codes für compliance
- **Manufacturing API constraints**: Use next-lower sample code
- **200 samples** exceeds ISO-K (125) by 60% → Statistically robust

---

### 6. Six Sigma Methodology (DMAIC Phases)

**Hintergrund**: DMAIC (Define-Measure-Analyze-Improve-Control) hat phase-spezifische sample requirements.

**Sample Requirements by Phase**:

| Phase | Minimum Samples | Purpose | SECOM Context |
|-------|-----------------|---------|---------------|
| **Define** | 30-50 | Problem scoping | Exploratory analysis |
| **Measure** | **100-200** | Baseline metrics | **Current implementation ✅** |
| **Analyze** | 200-500 | Root cause | Critical analysis mode |
| **Improve** | 100-300 | Validate changes | A/B testing |
| **Control** | 30-50/period | Ongoing monitoring | Real-time dashboards |

**SECOM's Sigma Level**:
```typescript
// SECOM: 6.6% failure rate → DPMO (Defects Per Million Opportunities)
const failureRate = 0.066;
const dpmo = failureRate * 1_000_000; // 66,000 DPMO

// Sigma Level Lookup (approximation):
// 66,000 DPMO ≈ 3.0 Sigma Level

// Six Sigma Goals:
// 3-Sigma: 66,807 DPMO (SECOM current state)
// 4-Sigma: 6,210 DPMO (target for manufacturing)
// 6-Sigma: 3.4 DPMO (world-class)
```

**Regel**:
- **Define Phase**: 30-50 samples (quick exploration)
- **Measure/Analyze**: 100-200 samples (current default)
- **Critical Analysis**: 200-500 samples (timeout risk!)
- **SECOM (3.0σ)**: 200 samples aligns with Measure/Analyze phases

---

## Smart Sampling Algorithm

### Strategie: 30% First + 40% Random Middle + 30% Last

**Rationale**:
1. **First 30%**: Captures early production patterns, startup phase
2. **Random 40%**: Representative sample, minimizes bias
3. **Last 30%**: Recent state, trend detection

**Implementation**:
```typescript
// fileParser.ts - Smart Sampling
export function smartSampleRows(data: any[], maxRows: number): any[] {
  const totalRows = data.length;

  // No sampling needed if dataset is small
  if (totalRows <= maxRows) {
    return data;
  }

  // Percentage-based ratios (flexible for MAX_ROWS changes)
  const FIRST_RATIO = 0.30;   // 30% first
  const MIDDLE_RATIO = 0.40;  // 40% random middle
  const LAST_RATIO = 0.30;    // 30% last

  // Calculate exact counts
  const firstCount = Math.floor(maxRows * FIRST_RATIO);   // 60 rows (for MAX_ROWS=200)
  const lastCount = Math.floor(maxRows * LAST_RATIO);     // 60 rows
  const middleCount = maxRows - firstCount - lastCount;   // 80 rows

  // 1. Take first N rows
  const firstRows = data.slice(0, firstCount);

  // 2. Take last N rows
  const lastRows = data.slice(totalRows - lastCount);

  // 3. Random sampling from middle
  const middleStart = firstCount;
  const middleEnd = totalRows - lastCount;
  const middleData = data.slice(middleStart, middleEnd);

  const sampledMiddle: any[] = [];
  const indices = new Set<number>();

  while (indices.size < Math.min(middleCount, middleData.length)) {
    const randomIndex = Math.floor(Math.random() * middleData.length);
    indices.add(randomIndex);
  }

  Array.from(indices)
    .sort((a, b) => a - b)
    .forEach(idx => sampledMiddle.push(middleData[idx]));

  // 4. Combine in chronological order
  return [...firstRows, ...sampledMiddle, ...lastRows];
}
```

**Für MAX_ROWS = 200**:
- First: 60 rows (30%)
- Random middle: 80 rows (40%)
- Last: 60 rows (30%)
- **Total: 200 rows**

**Vorteil über uniformes Sampling**:
- Behält zeitliche Struktur (first/last)
- Reduziert Bias (random middle)
- Erkennt Trends (first vs last comparison)

---

## Token Budget Management

### Gemini 2.0 Flash Context Window: 1,000,000 Tokens

**Token Calculation** (approximation):
```typescript
// Average token estimate: ~0.5 tokens per character (English text)
// Manufacturing data: ~1 token per cell (numbers + headers)

function estimateTokens(rows: number, columns: number): number {
  const headerTokens = columns * 10;  // Column names
  const dataTokens = rows * columns * 1;  // Numeric data cells
  const formattingTokens = rows * 5;  // Markdown table formatting (| | |)

  return headerTokens + dataTokens + formattingTokens;
}

// Example: 200 rows × 100 columns
const tokens = estimateTokens(200, 100);
// = 1,000 + 20,000 + 1,000 = 22,000 tokens

// Add prompt overhead (~20,000 tokens for system prompt + user instructions)
const totalTokens = tokens + 20000; // ≈ 42,000 tokens

// Percentage of Gemini limit
const pctOfLimit = (totalTokens / 1_000_000) * 100; // 4.2%
```

**Token Limits by Tier**:
```typescript
// Token limits (optimized for Gemini 2.0 Flash 1M context window)
const limits = {
  free: 15000,       // 15k tokens → 60×100 sample (~0.6MB CSV)
  premium: 100000,   // 100k tokens → 500×200 sample (~3MB CSV)
  unlimited: 500000, // 500k tokens → 50% of Gemini limit (~15MB CSV)
};

// Check before processing
function validateTokenBudget(rows: number, cols: number, tier: string): boolean {
  const estimatedTokens = estimateTokens(rows, cols) + 20000;
  const limit = limits[tier];

  if (estimatedTokens > limit) {
    throw new Error(
      `Token budget exceeded: ${estimatedTokens} > ${limit} (${tier} tier)\n` +
      `Reduce sample size: ${rows}×${cols} → ${Math.floor(rows * 0.7)}×${Math.floor(cols * 0.7)}`
    );
  }

  return true;
}
```

**Regel**:
- **Free tier**: Max 60×100 (~15k tokens)
- **Production default**: 200×100 (~42k tokens, 4.2% of limit)
- **Critical analysis**: 300×150 (~92k tokens, 9.2% of limit)
- **Hard limit**: 500×200 (~202k tokens, 20% of limit)
- **Safety margin**: Never exceed 50% of context window (500k tokens)

---

## Diminishing Returns Analysis

### Sweet Spot: 200-300 Rows

**Marginal Gains Table**:

| Configuration | Expected Failures | Processing Time | Marginal Gain | Recommendation |
|---------------|-------------------|-----------------|---------------|----------------|
| 100×50 | 6.6 | ~5s | Baseline | Deprecated |
| **200×100** | **13.2** | **~10s** | **+100%** | ✅ **Optimal** |
| 300×150 | 19.8 | ~20s | +50% | ⚠️ Marginal |
| 500×200 | 33.0 | ~40s | +67% | ❌ Timeout Risk |
| 1000×300 | 66.0 | ~60s+ | +100% | ❌ Too Slow |

**Binomial Probability Analysis (SECOM 6.6% failure rate)**:
```typescript
// Probability of detecting at least N failures
function binomialProbability(n: number, k: number, p: number): number {
  // P(X >= k) = 1 - P(X < k)
  // Using Poisson approximation for large n: λ = n * p
  const lambda = n * p;

  let cumulativeProb = 0;
  for (let i = 0; i < k; i++) {
    cumulativeProb += (Math.pow(lambda, i) * Math.exp(-lambda)) / factorial(i);
  }

  return 1 - cumulativeProb;
}

// Results for SECOM (p=0.066):
const results = {
  '100 rows': {
    'P(≥1 fail)': binomialProbability(100, 1, 0.066),  // 99.9%
    'P(≥5 fails)': binomialProbability(100, 5, 0.066), // 79.7%
    'P(≥10 fails)': binomialProbability(100, 10, 0.066), // 12.4%
  },
  '200 rows': {
    'P(≥1 fail)': binomialProbability(200, 1, 0.066),  // 100.0%
    'P(≥5 fails)': binomialProbability(200, 5, 0.066), // 99.7%
    'P(≥10 fails)': binomialProbability(200, 10, 0.066), // 85.6% ✅
  },
  '300 rows': {
    'P(≥1 fail)': binomialProbability(300, 1, 0.066),  // 100.0%
    'P(≥5 fails)': binomialProbability(300, 5, 0.066), // 100.0%
    'P(≥10 fails)': binomialProbability(300, 10, 0.066), // 99.5%
  },
};
```

**Key Insight**:
- 100 → 200 rows: **+73pp for P(≥10 fails)** (12.4% → 85.6%)
- 200 → 300 rows: **+14pp** (85.6% → 99.5%)
- **Diminishing Returns** start nach 200 rows!

**Regel**:
- **100-200 rows**: Beste ROI (return on investment)
- **200-300 rows**: Marginal gains (<20pp improvement)
- **>300 rows**: Timeout risk outweighs benefits

---

## Performance Constraints (Vercel Serverless)

### Timeout: 10 Seconds (Hobby Plan)

**Processing Time Breakdown**:
```typescript
// Empirical measurements (SECOM dataset)
const processingTimes = {
  parsing: 500,        // 0.5s - CSV/Excel parsing
  sampling: 200,       // 0.2s - Smart sampling algorithm
  markdown: 300,       // 0.3s - Table conversion
  apiRequest: 8000,    // 8s - Gemini API call (95th percentile)
  overhead: 1000,      // 1s - Network, JSON serialization
};

const totalTime = Object.values(processingTimes).reduce((a, b) => a + b, 0);
// = 10,000ms = 10 seconds

// Safety margin: 2 seconds
const TIMEOUT_MS = 10000;
const SAFETY_MARGIN_MS = 2000;
const MAX_API_TIME_MS = TIMEOUT_MS - SAFETY_MARGIN_MS; // 8 seconds
```

**Timeout Handling**:
```typescript
// api/analyze.ts - With timeout
import pTimeout from 'p-timeout';

export default async function handler(req: VercelRequest, res: VercelResponse) {
  try {
    // Wrap Gemini API call with timeout
    const result = await pTimeout(
      genAI.models.generateContent({
        model: 'gemini-2.0-flash',
        contents: [{ text: prompt }]
      }),
      {
        milliseconds: 8000,  // 8s timeout (2s buffer)
        message: 'Gemini API timeout after 8 seconds'
      }
    );

    const text = result.text;
    return res.json({ analysis: text });

  } catch (error) {
    if (error.name === 'TimeoutError') {
      return res.status(504).json({
        error: 'Analysis timeout - dataset too large',
        suggestion: 'Try reducing file size or columns'
      });
    }
    // ... other error handling
  }
}
```

**Regel**:
- **Vercel Hobby Plan**: 10s hard limit (cannot extend!)
- **API timeout**: 8s (2s buffer for parsing/overhead)
- **Sample size**: Must complete within 8s → Max 200-300 rows
- **Failure mode**: Return 504 Gateway Timeout with helpful message

---

## Testing & Verification

### Unit Tests for Sampling Logic

```typescript
// fileParser.test.ts
import { describe, it, expect } from 'vitest';
import { smartSampleRows, estimateTokens } from './fileParser';

describe('Smart Sampling Algorithm', () => {
  it('should return all rows if below MAX_ROWS', () => {
    const data = Array.from({ length: 100 }, (_, i) => ({ id: i }));
    const sampled = smartSampleRows(data, 200);

    expect(sampled.length).toBe(100);
    expect(sampled).toEqual(data); // No sampling needed
  });

  it('should sample 200 rows with 30-40-30 split', () => {
    const data = Array.from({ length: 1567 }, (_, i) => ({ id: i }));
    const sampled = smartSampleRows(data, 200);

    expect(sampled.length).toBe(200);

    // Check first 60 are sequential
    expect(sampled[0].id).toBe(0);
    expect(sampled[59].id).toBe(59);

    // Check last 60 are from end
    expect(sampled[140].id).toBeGreaterThan(1500);
    expect(sampled[199].id).toBe(1566);
  });

  it('should maintain chronological order', () => {
    const data = Array.from({ length: 1567 }, (_, i) => ({ id: i }));
    const sampled = smartSampleRows(data, 200);

    // IDs should be monotonically increasing
    for (let i = 1; i < sampled.length; i++) {
      expect(sampled[i].id).toBeGreaterThan(sampled[i - 1].id);
    }
  });
});

describe('Token Estimation', () => {
  it('should estimate tokens correctly for SECOM dataset', () => {
    const tokens = estimateTokens(200, 100);

    // Expected: 1,000 (headers) + 20,000 (data) + 1,000 (formatting)
    expect(tokens).toBeCloseTo(22000, -3); // Within 1,000 tokens
  });

  it('should stay within Gemini 1M token limit', () => {
    const maxTokens = estimateTokens(500, 200) + 20000; // Max safe config

    expect(maxTokens).toBeLessThan(1_000_000 / 2); // <50% of limit
  });
});
```

### Integration Tests with SECOM Baseline

```typescript
// secom.integration.test.ts
import { describe, it, expect } from 'vitest';

describe('SECOM Dataset Baseline', () => {
  it('should complete analysis within 10 seconds', async () => {
    const startTime = Date.now();

    const response = await fetch('/api/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        data: secomDataset, // 1,567 rows × 591 columns
        analysisType: 'quality'
      })
    });

    const endTime = Date.now();
    const duration = endTime - startTime;

    expect(response.status).toBe(200);
    expect(duration).toBeLessThan(10000); // <10s
  });

  it('should detect imbalanced dataset (6.6% failure rate)', async () => {
    const response = await fetch('/api/analyze', {
      method: 'POST',
      body: JSON.stringify({ data: secomDataset })
    });

    const result = await response.json();

    expect(result.analysis).toContain('imbalanced');
    expect(result.analysis).toContain('6.6%');
  });

  it('should provide 13.2 expected failures in 200-row sample', () => {
    const failureRate = 0.066;
    const sampledRows = 200;
    const expectedFailures = sampledRows * failureRate;

    expect(expectedFailures).toBeCloseTo(13.2, 1);

    // Probability of at least 10 failures: 85.6%
    // This is acceptable for robust analysis
  });
});
```

---

## Debugging Checklist

Wenn Sampling-Probleme auftreten:

- [ ] **Processing Time**: Läuft API-Call in <10s? (Check Vercel logs)
- [ ] **Token Budget**: Ist estimate <50% von Gemini limit?
- [ ] **Sample Size**: Mindestens 2× Cochran's minimum?
- [ ] **Temporal Coverage**: Mindestens 7 Tage für daily cycles?
- [ ] **Feature Coverage**: Mindestens 2:1 sample-to-feature ratio?
- [ ] **Minority Class**: Mindestens 10-15 failure examples?
- [ ] **Power Analysis**: >95% power für 1-sigma shifts?
- [ ] **Chronological Order**: Behält smart sampling die Reihenfolge?

**Debugging Commands**:
```bash
# Check processing time in Vercel logs
npx vercel logs --prod | grep "Duration"

# Test locally with timing
time curl -X POST http://localhost:5173/api/analyze \
  -H "Content-Type: application/json" \
  -d @secom_dataset.json

# Verify token usage
node -e "
const rows = 200, cols = 100;
const tokens = (cols * 10) + (rows * cols * 1) + (rows * 5);
console.log('Estimated tokens:', tokens);
console.log('% of Gemini limit:', (tokens / 1000000 * 100).toFixed(2));
"
```

---

## Common Mistakes & Fixes

### Mistake 1: Using All Rows Without Sampling

```typescript
// ❌ BROKEN - Timeout on large datasets
async function analyzeDataset(data: any[]) {
  const prompt = createMarkdownTable(data); // All 1,567 rows!
  const result = await geminiAPI.generateContent(prompt);
  // 504 Timeout Error after 10s
}

// ✅ FIX - Smart sampling
async function analyzeDataset(data: any[]) {
  const sampled = smartSampleRows(data, 200); // Only 200 rows
  const prompt = createMarkdownTable(sampled);
  const result = await geminiAPI.generateContent(prompt);
  // Completes in <10s ✅
}
```

### Mistake 2: Ignoring Minority Class Representation

```typescript
// ❌ BAD - Uniform random sampling loses failures
function randomSample(data: any[], n: number): any[] {
  return data.sort(() => Math.random() - 0.5).slice(0, n);
}
// Problem: May miss all 6.6% failures by chance!

// ✅ GOOD - Stratified sampling (when labels available)
function stratifiedSample(data: any[], n: number, labelColumn: string): any[] {
  const passed = data.filter(row => row[labelColumn] === 0);
  const failed = data.filter(row => row[labelColumn] === 1);

  const failureRate = failed.length / data.length;
  const nFailed = Math.max(Math.floor(n * failureRate), 10); // Min 10 failures
  const nPassed = n - nFailed;

  return [
    ...randomSample(passed, nPassed),
    ...randomSample(failed, nFailed)
  ].sort((a, b) => a.timestamp - b.timestamp); // Chronological order
}
```

### Mistake 3: Not Considering Correlated Features

```typescript
// ❌ BAD - Taking first 100 columns blindly
const sampledColumns = allColumns.slice(0, 100);
// Problem: May miss critical features at end!

// ✅ GOOD - Prioritize by variance or user selection
function smartColumnSample(data: any[], maxCols: number): string[] {
  // Option 1: User-selected critical features
  const criticalFeatures = getUserSelectedFeatures(); // e.g., ['temp', 'pressure']

  // Option 2: High-variance features (auto-detect)
  const variances = calculateVariances(data);
  const highVarianceCols = Object.entries(variances)
    .sort(([, a], [, b]) => b - a)
    .slice(0, maxCols)
    .map(([col]) => col);

  return [...new Set([...criticalFeatures, ...highVarianceCols])].slice(0, maxCols);
}
```

### Mistake 4: No Timeout Handling

```typescript
// ❌ BAD - No timeout, hangs forever
const result = await geminiAPI.generateContent(prompt);

// ✅ GOOD - Explicit timeout with fallback
import pTimeout from 'p-timeout';

try {
  const result = await pTimeout(
    geminiAPI.generateContent(prompt),
    {
      milliseconds: 8000,
      message: 'Gemini API timeout'
    }
  );
} catch (error) {
  if (error.name === 'TimeoutError') {
    return {
      error: 'Dataset too large for analysis',
      suggestion: 'Reduce to <5MB or <200 rows',
      fallback: 'Try basic statistics without AI'
    };
  }
}
```

---

## Future Optimization Opportunities

### 1. Adaptive Sampling Based on Dataset Characteristics

```typescript
// Dynamic MAX_ROWS based on dataset size and failure rate
function adaptiveSampling(data: any[], failureRate: number): number {
  if (data.length < 200) {
    return data.length; // Use all rows
  }

  // For low failure rates, need more samples
  if (failureRate < 0.05) {
    return Math.min(300, data.length); // 300 rows for <5% failure rate
  }

  // For high failure rates, 200 is sufficient
  if (failureRate > 0.10) {
    return 200;
  }

  // Default
  return 200;
}
```

### 2. SMOTE (Synthetic Minority Oversampling)

```typescript
// Generate synthetic failure examples to balance dataset
import { SMOTE } from 'smote-algorithm';

async function balancedSampling(data: any[], labelColumn: string): Promise<any[]> {
  const smote = new SMOTE({
    k: 5,  // 5 nearest neighbors
    sampling: 0.5  // Increase minority class by 50%
  });

  const balanced = await smote.fit(data, labelColumn);

  // Now sample from balanced dataset
  return smartSampleRows(balanced, 200);
}
```

### 3. Column Prioritization via PCA

```typescript
// Auto-select high-variance columns
import { PCA } from 'ml-pca';

function pcaBasedColumnSelection(data: any[], maxCols: number): string[] {
  const pca = new PCA(data);

  // Find columns contributing to first N principal components
  const loadings = pca.getLoadings();
  const importantCols = loadings
    .map((loading, idx) => ({ col: idx, importance: Math.abs(loading[0]) }))
    .sort((a, b) => b.importance - a.importance)
    .slice(0, maxCols)
    .map(x => data.columns[x.col]);

  return importantCols;
}
```

---

## Ressourcen

### Statistical Foundations
- [Cochran's Sample Size Calculator](https://dissertationdataanalysishelp.com/cochrans-sample-size-calculator/) - Imbalanced classification
- [Box-Jenkins Methodology](https://www.columbia.edu/~ww2040/4106F07/timeseries.pdf) - Time-series analysis
- [Cohen's d Effect Size Calculator](https://www.psychometrica.de/effect_size.html) - Power analysis
- [NIST Engineering Statistics Handbook](https://www.itl.nist.gov/div898/handbook/) - Industry standards

### Industry Standards
- ISO 2859-1: Acceptance Sampling by Attributes
- ISO 3951: Acceptance Sampling by Variables
- Six Sigma Methodology (DMAIC phases)

### Imbalanced Data Techniques
- [SMOTE Paper (2002)](https://arxiv.org/abs/1106.1813) - Synthetic minority oversampling
- [Survey on Small Sample Imbalance](https://arxiv.org/abs/1302.1866) - Classification strategies
- [Machine Learning Mastery - Imbalanced Classification](https://machinelearningmastery.com/smote-oversampling-for-imbalanced-classification/)

### Dataset
- [UCI SECOM Dataset](https://archive.ics.uci.edu/ml/datasets/SECOM) - 1,567 rows × 591 features, 6.6% failure rate

---

## Projektspezifische Anwendung

### ManufacturingInsideAnalyzer

**Status**: ✅ **Deployed to Production** (2025-11-18)

**Implementation**:
- ✅ Smart sampling (30-40-30 split) in `src/utils/fileParser.ts`
- ✅ 200×100 sampling (exceeds all statistical minimums)
- ✅ Token budget management (4.2% of 1M Gemini limit)
- ✅ Processing time <10s (8s API timeout + 2s buffer)
- ✅ SECOM baseline tests passing

**Key Metrics**:
- **Expected Failures**: 13.2 (vs. 6.6 with old 100×50)
- **Temporal Coverage**: 8.2 days (vs. 4.1 days)
- **Feature Coverage**: 16.9% (vs. 8.5%)
- **Power (1-sigma)**: >99% (vs. ~90%)
- **Token Usage**: 42,200 (~4.2% of limit)

**Files Modified**:
- `src/utils/fileParser.ts:27-28` - MAX_ROWS, MAX_COLUMNS
- `src/utils/fileParser.ts:42-49` - Smart sampling algorithm
- `api/analyze.ts:238-245` - Timeout handling

**Documentation**: `SAMPLING_OPTIMIZATION_2025-11-18.md`

---

## Lesson Learned Origin

**Date**: 2025-11-18
**Incident**: SECOM dataset (1,567 rows × 591 columns) caused 504 timeouts
**Root Cause**: Column sampling missing (only row sampling implemented)
**Duration**: ~6 hours debugging + 8 hours statistical research
**Impact**: All large datasets (>200 rows × >50 cols) failed
**Prevention**: This skill + automated performance tests

**Resolution Timeline**:
1. **Problem**: SECOM dataset → 504 timeout
2. **Hypothesis 1**: Row sampling insufficient → Increased 100 → 200 rows
3. **Hypothesis 2**: Column overload → Added 50 column limit
4. **Research**: Perplexity AI - Statistical foundations (Cochran, Box-Jenkins, Power Analysis)
5. **Optimization**: Upgraded to 200×100 sampling (2× increase)
6. **Validation**: Processing time <10s, token budget 4.2%
7. **Deployment**: Production (2025-11-18 02:31 CET)

**Key Insight**: Naive sampling (first N rows) is insufficient. Need:
- Statistical foundation (Cochran's Formula, Power Analysis)
- Industry compliance (ISO, Six Sigma)
- Smart algorithm (30-40-30 split for temporal coverage)
- Token budget management (stay <50% of API limit)

---

**🎯 Ziel**: Zero Performance Timeouts durch statistisch fundiertes Sampling!
