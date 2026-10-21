# Logarithmic-Hamming-Numbers

## $$H = 2^i \cdot 3^j \cdot 5^k$$

---

> A geometric & logarithmic solution to solving for the nth Hamming Number, allowing for faster calculation of very large Hamming Numbers compared to Dijkstra’s three-pointer method

[Desmos Visualization](https://www.desmos.com/3d/d7d7piepsc) · [Hamming Numbers Wikipedia](https://en.wikipedia.org/wiki/Regular_number) · [Codewars Kata](https://www.codewars.com/kata/526d84b98f428f14a60008da)

---

<p align="center">
  <img width="800" alt="desmos" src="https://github.com/user-attachments/assets/06f7712d-846e-4a0e-92fb-5a673e78f24d" />
</p>

---

A Hamming Number is a positive integer whose only prime divisors are 2, 3, and 5. Expressed as 

$$H = 2^i \cdot 3^j \cdot 5^k$$

Standard approaches (like Dijkstra's algorithm) sequentially store every number in memory. This makes it impossible to calculate very high Hamming numbers, but isn't a problem for, say, $n<5000$. Calculating $n=1,000,000,000$ requires a different approach, which is outlined here.

We translate the problem into a 3D geometry problem using logarithmic space. This allows us to compute massive values like $n=1,000,000,000$ much more efficiently, see the table below.

By taking the natural logarithm of the standard Hamming equation we get an additive linear inequality:

$$i \ln(2) + j \ln(3) + k \ln(5) \le B$$

where $B$ is some boundary beneath which is exactly $n$ Hamming numbers. Finding the $n$-th Hamming number is equivalent to finding a limit B such that exactly $n$ integer coordinate points $(i, j, k)$ sit at or below B.

We then count how many valid $(i, j, k)$ combinations exist below a proposed logarithmic bound B. Because the space below the boundary within the first octant forms a 3D tetrahedron (see [Desmos](https://www.desmos.com/3d/d7d7piepsc)), the **count(B)** function iterates through $k$ and $j$, quickly calculating the remaining valid values for $i$ in constant time using division.

---

# Hamming Numbers Algorithm Comparison

| Metric / Feature | Dijkstra's Three-Pointer Algorithm | Geometric / Logarithmic Approach |
| :--- | :--- | :--- |
| **Time Complexity** | $\mathcal{O}(n)$ | $\mathcal{O}(n^{2/3})$ |
| **Space Complexity** | $\mathcal{O}(n)$ | $\mathcal{O}(1)$ or $\mathcal{O}(n^{2/3})$ auxiliary storage |
| **Core Mechanism** | Sequential generation using 3 pointers/queues ($2 \cdot H_i, 3 \cdot H_j, 5 \cdot H_k$) | Counting and searching lattice points in 3D logarithmic space ($\log 2 \cdot i + \log 3 \cdot j + \log 5 \cdot k \le C$) |
| **Arithmetic Precision** | Exact integer arithmetic (requires big integer types for high $n$) | Floating-point logarithmic bounds (requires arbitrary-precision floats or exact rational bounds for large $n$) |
| **Sequential Dependency** | Strictly sequential; requires calculating $H_1, \dots, H_{n-1}$ to find $H_n$ | Direct range evaluation; enables binary search directly for the $n$-th value |
| **Parallelization** | Hard to parallelize due to iterative pointer state dependency | Highly parallelizable (independent range counts across lattice slices) |
| **Best Use Case** | Small to moderate $n$ where exact integer multiplication is fast | Extremely large $n$ (e.g., $n > 10^{12}$) where linear time scaling becomes prohibitive |

---

# Finding the boundary

#### Equations

### Geometric Volume Estimate (initial guess)

Calculates the initial limit $B_0$ by approximating $n$ points as the volume of a 3D tetrahedron, adjusted by half the sum of the logarithmic bases to account for the discrete nature of the integer grid.

$$B_0 = \sqrt[3]{6n \ln(2) \ln(3) \ln(5)} - \frac{\ln(2) + \ln(3) + \ln(5)}{2}$$

### Continuous Density Derivative

Represents the theoretical rate of change of the tetrahedron's volume at a given bound $B$, serving as the "slope" for the Newton-Raphson step.

$$f'(B) = \frac{B^2}{2 \ln(2) \ln(3) \ln(5)}$$

### Newton-Raphson Update

Refines the initial geometric guess by calculating the difference between the actual evaluated count $C_0$ and the target index $n$, scaled by the theoretical slope $f'(B_0)$.

$$B_1 = B_0 - \frac{C_0 - n}{f'(B_0)}$$

### Secant Method Empirical Slope

Calculates the literal slope $m$ of the step function by evaluating the change in actual point counts ($C_0$ to $C_1$) over the change between the two most recent logarithmic boundary guesses ($B_0$ to $B_1$).

$$m = \frac{C_1 - C_0}{B_1 - B_0}$$

### Secant Method Iterative Update

Computes the next boundary guess $B_{\text{next}}$ using the empirical slope $m$ to rapidly close the gap between the most recent count $C_1$ and the target count $n$.

$$B_{\text{next}} = B_1 - \frac{C_1 - n}{m}$$

### Constant-Time Lattice Point Count

Determines the exact number of valid integer points along the $i$-axis for a specific $j$ and $k$ coordinate, bypassing an iterative loop through direct algebraic division and flooring.

$$i_{\text{count}} = \left\lfloor \frac{B - j \ln(3) - k \ln(5)}{\ln(2)} \right\rfloor + 1$$

---

Instead of guessing $B$ sequentially, the **hamming(n)** function uses a geometric continuous volume estimate, refined with a Newton-Raphson step, and finally the Secant Method to rapidly converge on the precise boundary $B$ where the number of lattice points exactly equals $n$.

Because count is a step function, the secant iteration lands somewhere in the interval $[\ln H_n, \ln H_{n+1})$ rather than exactly on $\ln H_n$. The algorithm therefore builds a bracket $[lo, hi]$ starting at $B \pm 10^{-8}$ and widens it until count(lo) < n ≤ count(hi). It then enumerates only the lattice points inside $(lo, hi]$ (a handful), sorts them, and returns the element at index $n - \text{count}(lo)$.

Each count call costs $O(n^{2/3})$ time, since it loops over $(j,k)$ pairs, and the memory used is only the band. Floating-point error limits accuracy for extremely large $n$.
