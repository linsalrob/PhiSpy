from setuptools import setup, Extension

# C++ extension module
ext_modules = [
    Extension(
        "PhiSpyRepeatFinder",
        sources=["src/repeatFinder.cpp"],
        language='c++',
        define_macros=[("PHISPY_PYTHON_MODULE", "1")]
    )
]

setup(ext_modules=ext_modules)

