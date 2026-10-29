"""Benchmarks hamming.hamming(n) against Dijkstra's three-pointer method.

Usage: python benchmark.py [--quick]
"""
import platform
import statistics
import sys
import time
import tracemalloc

import hamming


def dijkstra_sequence(n):
    """Reference O(n) three-pointer generator; returns the first n Hamming numbers."""
    h = [1] * n
    i2 = i3 = i5 = 0
    x2, x3, x5 = 2, 3, 5
    for idx in range(1, n):
        nxt = min(x2, x3, x5)
        h[idx] = nxt
        if nxt == x2:
            i2 += 1
            x2 = 2 * h[i2]
        if nxt == x3:
            i3 += 1
            x3 = 3 * h[i3]
        if nxt == x5:
            i5 += 1
            x5 = 5 * h[i5]
    return h


def dijkstra(n):
    return dijkstra_sequence(n)[-1]


def timed(fn, n, repeats):
    """Returns (result, median seconds) over `repeats` runs."""
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        result = fn(n)
        times.append(time.perf_counter() - start)
    return result, statistics.median(times)


def peak_memory(fn, n):
    """Peak Python heap allocation in bytes during one call."""
    tracemalloc.start()
    fn(n)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


def count_calls(n):
    """Number of count(B) evaluations hamming(n) makes."""
    original = hamming.count
    calls = 0

    def wrapper(B):
        nonlocal calls
        calls += 1
        return original(B)

    hamming.count = wrapper
    try:
        hamming.hamming(n)
    finally:
        hamming.count = original
    return calls


def fmt_time(s):
    if s < 1e-3:
        return f"{s * 1e6:.0f} µs"
    if s < 1:
        return f"{s * 1e3:.1f} ms"
    return f"{s:.2f} s"


def fmt_bytes(b):
    for unit in ("B", "KB", "MB", "GB"):
        if b < 1024:
            return f"{b:.0f} {unit}"
        b /= 1024
    return f"{b:.1f} TB"


# Exponents (i, j, k) of known Hamming numbers, from Rosetta Code / OEIS A051037
KNOWN = {
    1691: (5, 12, 3),             # 2125764000
    10**6: (55, 47, 64),
    10**9: (1334, 335, 404),
}


def main():
    quick = "--quick" in sys.argv
    print(f"Python {platform.python_version()} on {platform.processor() or platform.machine()}")

    # Correctness: every n up to N against Dijkstra, plus known large values
    N = 5000 if quick else 20000
    h = dijkstra_sequence(N)
    assert all(hamming.hamming(n) == h[n - 1] for n in range(1, N + 1))
    print(f"hamming(n) matches Dijkstra for every n in 1..{N}")
    for n, (i, j, k) in KNOWN.items():
        if quick and n > 10**6:
            continue
        assert hamming.hamming(n) == 2**i * 3**j * 5**k, n
        print(f"hamming({n:,}) = 2^{i} * 3^{j} * 5^{k}  ✓")

    sizes = [10**3, 10**4, 10**5, 10**6, 10**7, 10**8, 10**9]
    if quick:
        sizes = sizes[:5]
    dijkstra_max = 10**6 if quick else 10**7

    print()
    print("| n | Digits | Logarithmic | count() calls | Peak memory | Dijkstra | Dijkstra peak memory |")
    print("| ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for n in sizes:
        repeats = 5 if n <= 10**7 else 3
        result, t = timed(hamming.hamming, n, repeats)
        mem = peak_memory(hamming.hamming, n)
        calls = count_calls(n)
        if n <= dijkstra_max:
            d_result, d_t = timed(dijkstra, n, 3 if n <= 10**6 else 1)
            assert d_result == result
            d_mem = peak_memory(dijkstra, n)
            d_cols = f"{fmt_time(d_t)} | {fmt_bytes(d_mem)}"
        else:
            d_cols = "— | —"
        print(f"| {n:,} | {len(str(result)):,} | {fmt_time(t)} | {calls} | {fmt_bytes(mem)} | {d_cols} |",
              flush=True)


if __name__ == "__main__":
    main()
