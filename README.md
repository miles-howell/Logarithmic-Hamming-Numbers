# Logarithmic-Hamming-Numbers
## $$H = 2^i \cdot 3^j \cdot 5^k$$
---
A geometric &amp; logarithmic solution to solving for the nth Hamming Number, allowing for faster calculation of very large Hamming Numbers compared to Dijkstra’s three-pointer method

[Desmos Visualization](https://www.desmos.com/3d/d7d7piepsc)
[Hamming Numbers Wikipedia](https://en.wikipedia.org/wiki/Regular_number)
[Codewars Kata](https://www.codewars.com/kata/526d84b98f428f14a60008da)

---

A Hamming Number is a positive integer whose only prime divisors are 2, 3, and 5. Expressed as $$H = 2^i \cdot 3^j \cdot 5^k$$. Standard approaches (like Dijkstra's algorithm) sequentially store every number in memory. This makes it impossible to calculate very high Hamming numbers, but isn't a problem for, say, *n<5000*. Calculating *n=1,000,000,000* requires a different approach, which is outlined here.

We translate the problem into a 3D geometry problem using logarithmic space. This allows us to compute massive values like *n=1,000,000,000* in milliseconds with essentially no memory overhead.
By taking the natural logarithm of the standard Hamming equation we get an additive linear inequality: $$\ln(H) = i \ln(2) + j \ln(3) + k \ln(5)$$. Finding the *n*-th Hamming number is equivalent to finding a limit B such that exactly *n* integer coordinate points *(i, j, k)* sit at or below B. 

We then count how many valid *(i, j, k)* combinations exist below a proposed logarithmic bound B. Because the boundary forms a 3D tetrahedron (see [Desmos](https://www.desmos.com/3d/d7d7piepsc)), the **count(B)** function iterates through *k* and *j*, quickly calculating the remaining valid values for *i* in constant time using division. 

### Finding the boundary
---
#### Equations
#### **Geometric Volume Estimate (initial guess):** $$B_0 = \sqrt[3]{6n \ln(2) \ln(3) \ln(5)} - \frac{\ln(2) + \ln(3) + \ln(5)}{2}$$
#### **Continuous Density Derivative (used as slope in the Newton-Raphson equation):** $$f'(B) = \frac{B^2}{2 \ln(2) \ln(3) \ln(5)}$$
#### **Newton-Raphson Update:** $$B_1 = B_0 - \frac{C_0 - n}{f'(B_0)}$$
#### **Secant Method Iterative Update:** $$B_{\text{next}} = B_1 - \frac{C_1 - n}{m}$$
#### **Constant-Time Lattice Point Count:** $$i_{\text{count}} = \left\lfloor \frac{B - j \ln(3) - k \ln(5)}{\ln(2)} \right\rfloor + 1$$

Instead of guessing ***B*** sequentially, the **hamming(n)** function uses a geometric continuous volume estimate, refined with a Newton-Raphson step, and finally the Secant Method to rapidly converge on the precise boundary ***B*** where the number of lattice points exactly equals *n*.

Because floating-point math can introduce slight imprecision at these scales, the algorithm establishes a microscopic window '''[lo, hi]''' around the final bound ***B***. Instead of generating all *n* Hamming numbers, it only computes the tiny fraction of candidates that fall inside this specific band, sorts them, and returns the exact *n*-th value.
