import math

# Pre-calculated natural logarithms and their reciprocals
a, b, c = math.log(2), math.log(3), math.log(5)
inv_a, inv_b, inv_c = 1.0 / a, 1.0 / b, 1.0 / c
density_factor = 2.0 * a * b * c  # 2 * a * b * c for Newton density fallback

def _count_py(B):
    """Inlined logarithmic point counter with floating-point precision protection."""
    if B < 0:
        return 0
    total = 0
    rem_k = B
    while rem_k >= 0:
        rem_j = rem_k
        while rem_j >= 0:
            # 1e-11 guards against precision loss (e.g. 1.99999999999 -> 2)
            total += int(rem_j * inv_a + 1e-11) + 1
            rem_j -= b
        rem_k -= c
    return total

def _band_py(lo, hi):
    """Hamming numbers with lo < ln(value) <= hi, unsorted."""
    band = []
    p5 = 1
    rem_k_hi, rem_k_lo = hi, lo

    while rem_k_hi >= 0:
        p35 = p5
        rem_j_hi, rem_j_lo = rem_k_hi, rem_k_lo

        while rem_j_hi >= 0:
            top_hi = int(rem_j_hi * inv_a + 1e-11)
            top_lo = int(rem_j_lo * inv_a + 1e-11) if rem_j_lo >= 0 else -1

            for i in range(top_lo + 1, top_hi + 1):
                band.append(p35 << i)  # Fast power of 2 via bitwise shift

            p35 *= 3
            rem_j_hi -= b
            rem_j_lo -= b

        p5 *= 5
        rem_k_hi -= c
        rem_k_lo -= c
    return band

# Prefer the C extension (python setup.py build_ext --inplace); count() and the band walk return identical results.
try:
    from _hamming_count import count, band_terms

    def _band(lo, hi):
        return [(3**j * 5**k) << i for j, k, i in band_terms(lo, hi)]

    backend = "C extension"
except ImportError:
    count = _count_py
    _band = _band_py
    backend = "pure Python"

def hamming(n):
    """Calculates the 1-indexed n-th Hamming number using the Secant Method."""
    if n <= 0:
        raise ValueError("n must be a positive integer.")

    # 1. Initial geometric estimate (B0)
    B0 = (6 * n * a * b * c) ** (1/3) - (a + b + c) / 2
    C0 = count(B0)

    # 2. Warm start: Use 1 Newton step to find a high-quality second point (B1)
    f_prime = (B0 * B0) / density_factor if B0 > 0 else 1.0
    B1 = B0 - (C0 - n) / f_prime
    if abs(B1 - B0) < 1e-5:
        B1 = B0 + 1e-5
    C1 = count(B1)

    # 3. Secant Method loop (typically converges in 2-4 iterations)
    for _ in range(10):
        if C1 == n or abs(B1 - B0) < 1e-12:
            break

        # Calculate empirical slope between the last two evaluation points
        slope = (C1 - C0) / (B1 - B0)
        
        # Fallback to theoretical density if local slope is zero (e.g., flat step)
        if slope == 0:
            slope = (B1 * B1) / density_factor if B1 > 0 else 1.0

        # Compute next bound update
        B_next = B1 - (C1 - n) / slope

        # Shift variables
        B0, C0 = B1, C1
        B1, C1 = B_next, count(B_next)

    # 4. Guarantee that n is trapped inside the band [lo, hi].
    # C1 = count(B1) is already known, so it supplies one side of the band for free.
    # Steps double so sparse (small n) bands are reached in a few counts.
    if C1 >= n:
        hi = B1
        lo = B1 - 1e-8
        step = 1e-5
        below = count(lo)
        while below >= n:
            lo -= step
            step *= 2
            below = count(lo)
    else:
        lo, below = B1, C1
        hi = B1 + 1e-8
        step = 1e-5
        while count(hi) < n:
            hi += step
            step *= 2

    # 5. Extract and sort only the candidate numbers inside the narrow band
    band = _band(lo, hi)
    band.sort()
    return band[n - below - 1]
