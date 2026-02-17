# Claude Trading Skills - TODO List

## Completed ✅

- [x] **Institutional Flow Tracker** (2025-11-09)
  - 13F filings analysis with tier-based quality framework
  - Multi-quarter trend identification
  - FMP API integration
  - Comprehensive reference guides

## In Progress 🔄

None currently.

## Pending Skills 📋

### 1. Relative Strength Momentum Scanner (High Priority)

**Purpose:** Identify stocks with strong momentum using IBD-style RS Rating

**Why This Skill:**
- Fills gap in current skill set (no momentum/growth stock focus)
- Complements value-oriented skills (Value Dividend Screener, Institutional Flow Tracker)
- IBD RS Rating is proven methodology used by successful growth investors
- Useful for swing traders and growth stock investors

**Key Features:**
- RS Rating calculation (3/6/12-month price performance vs market)
- Sector-relative strength ranking
- Momentum acceleration/deceleration detection
- Volume-confirmed breakout candidates
- Stage analysis (Accumulation → Markup → Distribution → Decline)

**Data Requirements:**
- FMP API (price history, volume data)
- Free tier sufficient (250 calls/day)

**Implementation Plan:**
1. Create directory structure: `relative-strength-momentum-scanner/`
2. Write SKILL.md with RS Rating methodology
3. Create references:
   - `rs_rating_methodology.md` - IBD-style RS calculation
   - `momentum_patterns.md` - Acceleration/deceleration patterns
   - `breakout_checklist.md` - Volume confirmation criteria
4. Write scripts:
   - `calculate_rs_rating.py` - RS Rating for single stock
   - `screen_momentum_stocks.py` - Sector-wide screening
   - `detect_breakouts.py` - Volume-confirmed breakout detection
5. Write README.md with usage examples
6. Update main README.md and CLAUDE.md
7. Commit and push

**Estimated Effort:** 4-6 hours

**Integration Opportunities:**
- Combine with Institutional Flow Tracker (momentum + smart money confirmation)
- Combine with Technical Analyst (momentum + chart pattern confirmation)
- Use with Sector Analyst to identify leading sectors

---

### 2. Macro Indicator Dashboard (Medium Priority)

**Purpose:** Automated macro economic indicator tracking and analysis

**Why This Skill:**
- Current Market Environment Analysis requires manual data entry
- FRED API provides free access to thousands of economic indicators
- Automates weekly/monthly macro reviews
- Helps with sector rotation and risk management decisions

**Key Features:**
- Automated data fetching from FRED API (free, unlimited)
- Key indicators: Fed Funds Rate, CPI, PPI, NFP, GDP, ISM PMI, Yield Curve
- Trend visualization (MoM, YoY, moving averages)
- Economic cycle detection (Early/Mid/Late Cycle, Recession)
- Market implication generation (sector recommendations, risk posture)
- Alert system for threshold breaches

**Data Requirements:**
- FRED API (free, no rate limits for reasonable use)
- FMP API for complementary market data (optional)

**Implementation Plan:**
1. Create directory structure: `macro-indicator-dashboard/`
2. Write SKILL.md with economic indicator framework
3. Create references:
   - `fred_api_guide.md` - FRED API usage and key series IDs
   - `economic_cycle_indicators.md` - Leading/lagging indicators
   - `sector_cycle_mapping.md` - Sector performance by cycle phase
   - `interpretation_guide.md` - How to read macro data
4. Write scripts:
   - `fetch_macro_data.py` - Pull data from FRED API
   - `analyze_economic_cycle.py` - Determine current cycle phase
   - `generate_dashboard.py` - Create comprehensive macro report
5. Write README.md with indicator catalog
6. Update main README.md and CLAUDE.md
7. Commit and push

**Estimated Effort:** 6-8 hours

**Integration Opportunities:**
- Combine with Sector Analyst for cycle-based rotation
- Combine with US Market Bubble Detector for risk assessment
- Feed into Portfolio Manager for macro-driven rebalancing
- Support Stanley Druckenmiller Investment Advisor with macro themes

---

## Future Enhancements (Backlog) 💡

### 3. Short Interest Analyzer (Low Priority)
- Track short interest changes
- Identify short squeeze opportunities
- FMP API integration
- Estimated effort: 3-4 hours

### 4. Insider Trading Monitor (Low Priority)
- Form 4 filing analysis
- Cluster detection (multiple insiders buying/selling)
- FMP API integration
- Estimated effort: 4-5 hours

### 5. Swing Trade Setup Scanner (Low Priority)
- Price action pattern recognition
- Volume profile analysis
- Specific entry/stop/target levels
- Estimated effort: 6-8 hours

### 6. Correlation Matrix Analyzer (Low Priority)
- Portfolio correlation analysis
- Hidden concentration risk detection
- Diversification scoring
- Estimated effort: 3-4 hours

---

## Repository Maintenance Tasks

### Documentation Updates Needed
- [ ] Create Japanese README.ja.md for Institutional Flow Tracker
- [ ] Generate ZIP package for Institutional Flow Tracker (for web app users)
- [ ] Update skill catalog count in main README.md

### Testing Needed
- [ ] Test Institutional Flow Tracker scripts with real FMP API
- [ ] Validate reference guide accuracy against current SEC rules
- [ ] Test integration with other skills (Value Dividend Screener, etc.)

### Code Quality
- [ ] Add error handling improvements to all FMP API scripts
- [ ] Add retry logic with exponential backoff for rate limits
- [ ] Standardize output format across all skills

---

## Priority Order for Next Implementation

1. **Relative Strength Momentum Scanner** (completes momentum/growth coverage)
2. **Macro Indicator Dashboard** (automates manual workflow, high value-add)
3. Short Interest Analyzer (nice-to-have, lower priority)
4. Insider Trading Monitor (complementary to Institutional Flow)
5. Others as needed

---

## Notes

- All new skills should follow existing directory structure:
  ```
  <skill-name>/
  ├── SKILL.md
  ├── README.md
  ├── references/
  ├── scripts/
  └── assets/
  ```

- API keys should support both environment variables and CLI arguments
- All scripts should include help text and usage examples
- Reference guides should be comprehensive with examples
- Integration with existing skills should be documented

**Last Updated:** 2025-11-09
