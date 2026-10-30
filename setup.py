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
        )
    ],
)
