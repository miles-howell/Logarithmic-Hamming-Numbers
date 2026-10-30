from setuptools import Extension, setup

setup(
    name="hamming-count",
    version="0.1",
    ext_modules=[
        Extension(
            "_hamming_count",
            ["_hamming_count.c"],
            # Keep IEEE rounding identical to the pure-Python count().
            extra_compile_args=["-O2", "-ffp-contract=off"],
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
