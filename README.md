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

We translate the problem into a 3D geometry problem using logarithmic space. This allows us to compute massive values like $n=1,000,000,000$ much more efficiently: with the C extension the $10^9$-th Hamming number takes about 11 ms, and the $10^{16}$-th (a 181,949-digit number) takes under two minutes, with every result verified exactly (see [Benchmarks](#benchmarks)).

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
| **Arithmetic Precision** | Exact integer arithmetic (requires big integer types for high $n$) | Floating-point logarithmic bounds: 64-bit `double` is reliable to about $n = 10^{12}$, 80-bit `long double` to at least $10^{16}$ (checked exactly). The final answer is built with exact integers either way |
| **Sequential Dependency** | Strictly sequential; requires calculating $H_1, \dots, H_{n-1}$ to find $H_n$ | Direct range evaluation; enables binary search directly for the $n$-th value |
| **Parallelization** | Hard to parallelize due to iterative pointer state dependency | Highly parallelizable (independent range counts across lattice slices) |
| **Best Use Case** | Small to moderate $n$ where exact integer multiplication is fast | Extremely large $n$ (e.g., $n > 10^{12}$) where linear time scaling becomes prohibitive |

---

# Benchmarks

Measured with [`benchmark.py`](benchmark.py) on a Runpod CPU pod: AMD EPYC 9575F (up to 5.0 GHz), 8 vCPU, Python 3.13.18, using the C extension (see [Building the C extension](#building-the-c-extension)). Times are the median of 5 runs for $n \le 10^7$, 3 runs for $10^8$, and a single run from $10^9$ up. Peak memory is the Python heap measured with `tracemalloc` (skipped from $10^{13}$, where the run is not repeated). Dijkstra was not run past $10^7$, where it already needs almost 1 GB. Sizes below $10^6$ each finish in under 150 µs.

| $n$ | Digits | Logarithmic | `count()` calls | Peak memory | Exact rank check | Dijkstra | Dijkstra peak memory |
| ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: |
| 1,000,000 | 84 | 140 µs | 7 | 404 B | ✓ | 208 ms | 63 MB |
| 10,000,000 | 182 | 504 µs | 5 | 388 B | ✓ | 2.34 s | 934 MB |
| 100,000,000 | 392 | 2.2 ms | 5 | 892 B | ✓ | — | — |
| 1,000,000,000 | 844 | 11.4 ms | 7 | 9 KB | ✓ | — | — |
| 10,000,000,000 | 1,819 | 17.3 ms | 5 | 65 KB | ✓ | — | — |
| 100,000,000,000 | 3,920 | 71.0 ms | 10 | 2 KB | ✓ | — | — |
| 1,000,000,000,000 | 8,445 | 241 ms | 6 | 14 KB | ✓ | — | — |
| 10,000,000,000,000 | 18,195 | 1.12 s | 6 | — | ✓ | — | — |
| 100,000,000,000,000 | 39,199 | 6.12 s | 8 | — | ✓ | — | — |
| 1,000,000,000,000,000 | 84,453 | 24.3 s | 6 | — | ✓ | — | — |
| 10,000,000,000,000,000 | 181,949 | 107 s | 5 | — | ✓ | — | — |

### How the results are checked

Every row is verified two independent ways in the same script:

- `hamming(n)` matches Dijkstra for every $n$ from 1 to 20,000, and returns the known values $H_{1691} = 2125764000$, $H_{10^6} = 2^{55} \cdot 3^{47} \cdot 5^{64}$ and $H_{10^9} = 2^{1334} \cdot 3^{335} \cdot 5^{404}$.
- The **Exact rank check** column runs `_hamming_exact.c`, which takes the exponents $(i, j, k)$ of the result and counts every $(i, j, k)$ with $2^i 3^j 5^k \le$ the result. It does not reuse the floating-point code it is auditing: it compares exponent-weighted $\log_2$ values in 80-bit `long double`, redoes any cell that lands near an integer boundary in 128-bit quad precision, and raises an error instead of guessing if a cell still cannot be decided. A ✓ means the rank is exactly $n$. For $n = 10^9$ to $10^{11}$ the checker's verdict was also confirmed by an independent pure-integer count, with no floating point.

This check is what caught a real bug. An earlier version computed `count()` in 64-bit `double` with repeated subtraction. From about $n = 10^{13}$ neighbouring Hamming numbers are closer together in log space (about $2.45 / B^2$) than double rounding error, and the returned value was wrong: its rank was $n + 6$ at $10^{13}$ and $n + 210$ at $10^{14}$. The C code now uses `long double` and computes each remainder directly as $B - k \ln 5 - j \ln 3$ instead of subtracting repeatedly. The pure-Python fallback still uses `double`, so it is only reliable up to about $n = 10^{12}$.

Takeaways:

- The logarithmic method is ahead at every size tested. At $n = 10^7$ it is about 4,600× faster than Dijkstra and uses about 2.4 million times less memory (388 B against 934 MB).
- Time grows about 4–5.5× per 10× increase in $n$ from $10^{12}$ to $10^{16}$ (241 ms, 1.12 s, 6.12 s, 24.3 s, 107 s), close to the 4.6× that the $n^{2/3}$ cost of a count predicts.
- From $10^6$ up a search needs only 5–10 counts: a geometric first guess, one Newton step, a few secant steps, then a bracket of about $\pm 10^{-8}$ around the result. The bracket widens in doubling steps, so small $n$, where neighbouring Hamming numbers are far apart in log space, needs a few more counts.
- Memory stays tiny. Apart from the answer itself, the only thing held is the small band of candidate numbers around it.

Run it yourself with `python benchmark.py` (the Dijkstra reference at $10^7$ dominates the small sizes, and the $10^{16}$ row takes about three minutes on 8 cores with the exact check) or `python benchmark.py --quick` (up to $10^7$ for the logarithmic method, up to $10^6$ for Dijkstra).

## Building the C extension

```bash
python setup.py build_ext --inplace
```

This needs a C compiler and the Python headers. It builds two modules:

- `_hamming_count`: `count()` and the band walk, in 80-bit `long double` and parallelized with OpenMP (set `OMP_NUM_THREADS` to limit threads). x86-64 and GCC are assumed.
- `_hamming_exact`: the exact rank checker `exact_count(i, j, k)`. It needs OpenMP and `libquadmath`, and it is optional.

`hamming.py` falls back to pure Python if the extension is missing (reliable only to about $n = 10^{12}$), and `benchmark.py` prints which backend ran.

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

Because count is a step function, the secant iteration lands somewhere in the interval $[\ln H_n, \ln H_{n+1})$ rather than exactly on $\ln H_n$. The algorithm therefore builds a bracket $[lo, hi]$ starting at $B \pm 10^{-8}$ and widens it, doubling the step each time, until count(lo) < n ≤ count(hi). It then enumerates only the lattice points inside $(lo, hi]$ (a handful), sorts them, and returns the element at index $n - \text{count}(lo)$.

Each count call costs $O(n^{2/3})$ time, since it loops over $(j,k)$ pairs, and the memory used is only the band. Floating-point error limits accuracy for extremely large $n$.
