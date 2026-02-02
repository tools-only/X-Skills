# Equation Counting Bug Fix

## Problem

The `count_equations_in_file()` function in `book-metrics.py` was significantly over-counting LaTeX equations. For example, a signal processing textbook showed 2,421 equations when the actual count was much lower.

## Root Causes

The original implementation had two major bugs:

### 1. Double Counting Display Math

The original regex patterns:
```python
inline = len(re.findall(r'\$[^$]+\$', content))      # Inline math
display = len(re.findall(r'\$\$[^$]+\$\$', content))  # Display math
```

When a display math equation like `$$x = a + b$$` appeared in the content:
- The display pattern correctly matched `$$x = a + b$$` (counted as 1)
- The inline pattern then matched `$x = a + b$` (the content between the first and second `$`)
- The inline pattern also matched text between consecutive equations like `$\n\n### Title\n$`

This caused display math to be counted multiple times.

### 2. False Positives from Dollar Amounts

The inline pattern `\$[^$]+\$` matched dollar amounts like:
- `$500`
- `$1,000`
- `$25 to $75`

These are not LaTeX equations but were being counted as such.

## Solution

The fixed implementation:

```python
def count_equations_in_file(self, markdown_file: Path) -> int:
    """Count LaTeX equations in a single markdown file."""
    try:
        with open(markdown_file, 'r', encoding='utf-8') as f:
            content = f.read()

            # Count display math: $$...$$ (must come first)
            display_matches = re.findall(r'\$\$[^$]+?\$\$', content, re.DOTALL)
            display = len(display_matches)

            # Remove all display math blocks to avoid double-counting
            content_no_display = re.sub(r'\$\$[^$]+?\$\$', '', content, flags=re.DOTALL)

            # Count inline math: $...$
            # Negative lookahead (?!\d) ensures we don't match dollar amounts like $500
            inline_matches = re.findall(r'\$(?!\d)([^\$]+?)\$', content_no_display)
            inline = len(inline_matches)

            return inline + display
    except Exception as e:
        print(f"Warning: Could not read {markdown_file}: {e}")
        return 0
```

### Key Improvements

1. **Prevent Double Counting**: First count display math `$$...$$`, then remove all display math blocks from the content before counting inline math
2. **Exclude Dollar Amounts**: Use negative lookahead `(?!\d)` to ensure we don't match expressions like `$500` where a digit immediately follows the opening `$`
3. **Non-greedy Matching**: Use `+?` instead of `+` for better handling of adjacent equations

## Test Results

Created test file `equation-count-test.md` with:
- 5 display math equations
- 7 inline math equations
- Multiple dollar amounts (false positives)
- Total: **12 actual equations**

### Original Function Results
- **Counted**: 24 equations
- **Error**: +12 (100% over-count)
- **Status**: ❌ FAILED

### Fixed Function Results
- **Counted**: 12 equations
- **Error**: 0
- **Status**: ✅ PASSED

## Impact

For a typical intelligent textbook:
- **Before**: Counts could be 2-3x higher than actual
- **After**: Accurate counts that reflect actual LaTeX equation usage

This fix ensures the book metrics accurately represent the mathematical content in educational textbooks.

## Files Modified

- `src/book-metrics/book-metrics.py` - Fixed `count_equations_in_file()` method
- `src/book-metrics/equation-count-test.md` - Test file with 12 known equations
- `src/book-metrics/test-equation-count.py` - Original test script
- `src/book-metrics/test-equation-count-fixed.py` - Test comparing original vs fixed versions

## Testing

To verify the fix works:

```bash
cd src/book-metrics
python test-equation-count-fixed.py
```

Expected output: `✅ Fixed version PASSED!`
