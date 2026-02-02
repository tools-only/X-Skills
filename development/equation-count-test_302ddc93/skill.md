# Equation Count Test File

This test file contains exactly 10 LaTeX equations for testing the equation counter.

## Display Math Equations (5 equations)

### 1. Quadratic Formula
$$x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$$

### 2. Pythagorean Theorem
$$a^2 + b^2 = c^2$$

### 3. Euler's Identity
$$e^{i\pi} + 1 = 0$$

### 4. Normal Distribution
$$f(x) = \frac{1}{\sigma\sqrt{2\pi}} e^{-\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2}$$

### 5. Fourier Transform
$$F(\omega) = \int_{-\infty}^{\infty} f(t) e^{-i\omega t} dt$$

## Inline Math Equations (7 equations)

The area of a circle is $A = \pi r^2$ where $r$ is the radius.

Einstein's famous equation is $E = mc^2$ which relates energy and mass.

The derivative of $f(x) = x^2$ is $f'(x) = 2x$ by the power rule.

A sentence with math $x = 5$ and more text $y = 10$ has two equations.

## Potential False Positives (should NOT be counted)

This project costs $500 and will save us $1,000 per year.

The price increased from $25 to $75, a change of $50.

## Notes

The inline count is 7 because:
1. `A = πr²` - area of circle
2. `r` - variable reference in "where r is"
3. `E = mc²` - Einstein's equation
4. `f(x) = x²` - function
5. `f'(x) = 2x` - derivative
6. `x = 5` - assignment
7. `y = 10` - assignment

**Expected Count: 12 equations total (5 display + 7 inline)**
