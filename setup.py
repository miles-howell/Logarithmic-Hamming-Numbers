import platform

from setuptools import Extension, setup

# SSE3 makes long double -> integer conversion a single fisttp instead of a control-word dance.
X86 = platform.machine().lower() in ("x86_64", "amd64")

setup(
    name="hamming-count",
    version="0.1",
    ext_modules=[
        Extension(
            "_hamming_count",
            ["_hamming_count.c"],
            # 80-bit long double arithmetic (x86-64 only); OpenMP for the cell loops.
            extra_compile_args=["-O2", "-ffp-contract=off", "-fopenmp"] + (["-msse3"] if X86 else []),
            extra_link_args=["-fopenmp"],
            libraries=["m"],
        ),
        # Independent exact rank checker (x86-64, needs OpenMP + libquadmath); optional.
        Extension(
            "_hamming_exact",
            ["_hamming_exact.c"],
            extra_compile_args=["-O2", "-fopenmp", "-ffp-contract=off"],
            extra_link_args=["-fopenmp"],
            libraries=["quadmath", "m"],
            optional=True,
        ),
    ],
)
