from setuptools import setup, Extension

# C++ extension module
ext_modules = [
    Extension(
        "PhiSpyRepeatFinder",
        sources=["src/repeatFinder.cpp"],
        language='c++'
    )
]

setup(ext_modules=ext_modules)

