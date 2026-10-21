
# Logarithmic-Hamming-Numbers
## $$H = 2^i \cdot 3^j \cdot 5^k$$
---
A geometric &amp; logarithmic solution to solving for the nth Hamming Number, allowing for faster calculation of very large Hamming Numbers compared to Dijkstra’s three-pointer method

[Desmos Visualization](https://www.desmos.com/3d/d7d7piepsc)

[Hamming Numbers Wikipedia](https://en.wikipedia.org/wiki/Regular_number)

[Codewars Kata](https://www.codewars.com/kata/526d84b98f428f14a60008da)

---

<img width="1920" height="1080" alt="desmos" src="https://github.com/user-attachments/assets/06f7712d-846e-4a0e-92fb-5a673e78f24d" />

A Hamming Number is a positive integer whose only prime divisors are 2, 3, and 5. Expressed as $$H = 2^i \cdot 3^j \cdot 5^k$$. Standard approaches (like Dijkstra's algorithm) sequentially store every number in memory. This makes it impossible to calculate very high Hamming numbers, but isn't a problem for, say, $n<5000$. Calculating $n=1,000,000,000$ requires a different approach, which is outlined here.

We translate the problem into a 3D geometry problem using logarithmic space. This allows us to compute massive values like $n=1,000,000,000$ in milliseconds with essentially no memory overhead.
By taking the natural logarithm of the standard Hamming equation we get an additive linear inequality: $$\ln(H) = i \ln(2) + j \ln(3) + k \ln(5)$$. Finding the $n$-th Hamming number is equivalent to finding a limit B such that exactly $n$ integer coordinate points $(i, j, k)$ sit at or below B. 

We then count how many valid $(i, j, k)$ combinations exist below a proposed logarithmic bound B. Because the boundary forms a 3D tetrahedron (see [Desmos](https://www.desmos.com/3d/d7d7piepsc)), the **count(B)** function iterates through $k$ and $j$, quickly calculating the remaining valid values for $i$ in constant time using division. 

### Finding the boundary
---
#### Equations
**Geometric Volume Estimate (initial guess):** Calculates the initial limit $B_0$ by approximating $n$ points as the volume of a 3D tetrahedron, adjusted by half the sum of the logarithmic bases to account for the discrete nature of the integer grid.

 $$B_0 = \sqrt[3]{6n \ln(2) \ln(3) \ln(5)} - \frac{\ln(2) + \ln(3) + \ln(5)}{2}$$

**Continuous Density Derivative (used as slope in the Newton-Raphson equation):** Represents the theoretical rate of change (the surface area) of the tetrahedron's volume at a given bound $B$, serving as the slope for the Newton-Raphson step.

 $$f'(B) = \frac{B^2}{2 \ln(2) \ln(3) \ln(5)}$$

**Newton-Raphson Update:** Refines the initial geometric guess by calculating the difference between the actual evaluated count $C_0$ and the target index $n$, scaled by the theoretical slope $f'(B_0)$.

 $$B_1 = B_0 - \frac{C_0 - n}{f'(B_0)}$$

**Secant Method Empirical Slope:** Calculates the literal slope $m$ of the step function by evaluating the change in actual point counts ($C_0$ to $C_1$) over the change between the two most recent logarithmic boundary guesses ($B_0$ to $B_1$).

 $$m = \frac{C_1 - C_0}{B_1 - B_0}$$

**Secant Method Iterative Update:** Computes the next boundary guess $B_{\text{next}}$ using the empirical slope $m$ to rapidly close the gap between the most recent count $C_1$ and the target count $n$.

 $$B_{\text{next}} = B_1 - \frac{C_1 - n}{m}$$

**Constant-Time Lattice Point Count:** Determines the exact number of valid integer points along the $i$-axis for a specific $j$ and $k$ coordinate, bypassing an iterative loop through direct algebraic division and flooring.

 $$i_{\text{count}} = \left\lfloor \frac{B - j \ln(3) - k \ln(5)}{\ln(2)} \right\rfloor + 1$$

Instead of guessing $B$ sequentially, the **hamming(n)** function uses a geometric continuous volume estimate, refined with a Newton-Raphson step, and finally the Secant Method to rapidly converge on the precise boundary $B$ where the number of lattice points exactly equals $n$.

Because floating-point math can introduce slight imprecision at these scales, the algorithm establishes a microscopic window ```[lo, hi]``` around the final bound $B$. Instead of generating all $n$ Hamming numbers, it only computes the tiny fraction of candidates that fall inside this specific band, sorts them, and returns the exact $n$-th value.
