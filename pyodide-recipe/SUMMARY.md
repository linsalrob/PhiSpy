# PhiSpy Pyodide Recipe - Summary

This directory contains everything needed to build and distribute PhiSpy for Pyodide, enabling PhiSpy to run in web browsers and Node.js environments via WebAssembly.

## What's Included

### Recipe Files (for pyodide-recipes repository)

1. **meta.yaml** - The core recipe file
   - Package metadata (name, version, URLs)
   - Build configuration for C++ extension compilation
   - Dependency specifications
   - Import tests

2. **test_phispy.py** - Automated tests
   - Validates all imports work correctly
   - Tests dependencies are available
   - Checks basic functionality

### Documentation

3. **README.md** - Quick start guide
   - Overview of the recipe
   - Dependencies and requirements
   - Basic build instructions
   - Usage overview

4. **SUBMISSION_GUIDE.md** - Detailed submission instructions
   - Step-by-step process for submitting to pyodide-recipes
   - Testing and validation procedures
   - PR submission guidelines
   - Troubleshooting common issues

5. **USAGE.md** - Comprehensive usage documentation
   - Installation in Pyodide environments
   - Code examples for browser use
   - Working with GenBank files
   - Performance considerations
   - Advanced web interface examples

## Quick Start for Recipe Submission

1. Fork https://github.com/pyodide/pyodide-recipes
2. Copy `meta.yaml` and `test_phispy.py` to `packages/phispy/`
3. Build and test: `pyodide build-recipes phispy`
4. Submit a pull request

See **SUBMISSION_GUIDE.md** for complete details.

## Quick Start for Using PhiSpy in Pyodide

```python
import micropip
await micropip.install('phispy')

from PhiSpyModules import main
# ... use PhiSpy
```

See **USAGE.md** for complete examples.

## Technical Details

- **PhiSpy Version:** 5.0.1
- **License:** MIT
- **Homepage:** https://github.com/linsalrob/PhiSpy
- **PyPI:** https://pypi.org/project/PhiSpy/

### Dependencies

| Package | Status in Pyodide | Type |
|---------|-------------------|------|
| biopython | ✅ Built-in | C extension |
| numpy | ✅ Built-in | C extension |
| scikit-learn | ✅ Built-in | C extension |
| bcbio-gff | ✅ Included in recipe | Pure Python |

### C++ Extension

PhiSpy includes a C++ extension module (`PhiSpyRepeatFinder`) for efficient repeat sequence finding. The recipe configures this to compile to WebAssembly with:
- `-fwasm-exceptions` flag for exception support
- Emscripten toolchain for browser compatibility

## Architecture

```
PhiSpy Package
├── PhiSpyModules/ (Python modules)
│   ├── main.py
│   ├── classification.py
│   ├── helper_functions.py
│   └── ... (other modules)
└── PhiSpyRepeatFinder (C++ extension → WebAssembly)
```

## Files in This Directory

```
pyodide-recipe/
├── meta.yaml              # Main recipe file (required)
├── test_phispy.py        # Test script (optional but recommended)
├── README.md             # Quick reference
├── SUBMISSION_GUIDE.md   # Detailed submission instructions
├── USAGE.md              # Usage examples and documentation
└── SUMMARY.md            # This file
```

## Building Locally

```bash
# Clone pyodide-recipes
git clone --recurse-submodules https://github.com/pyodide/pyodide-recipes.git
cd pyodide-recipes

# Copy recipe
cp -r /path/to/PhiSpy/pyodide-recipe packages/phispy

# Install build tools
pip install ./pyodide-build

# Build
pyodide build-recipes phispy

# Output will be in dist/
```

## Testing

After building, the package can be tested:

```bash
cd dist
python -m http.server 8000
# Open browser to localhost:8000
# In browser console:
# await loadPyodide()
# await pyodide.loadPackage('phispy')
# await pyodide.runPython('import PhiSpyModules')
```

## Integration with PhiSpy Repository

This recipe is maintained in the main PhiSpy repository under `/pyodide-recipe/`. 

For the pyodide-recipes repository submission, only `meta.yaml` and optionally `test_phispy.py` are needed. The other documentation files are for reference and can be included in the PhiSpy repository documentation.

## Maintenance

When updating PhiSpy:
1. Update version number in `meta.yaml`
2. Update source URL and SHA256 hash
3. Test build with new version
4. Update recipe in pyodide-recipes repository

## Support and Issues

- PhiSpy issues: https://github.com/linsalrob/PhiSpy/issues
- Pyodide issues: https://github.com/pyodide/pyodide/issues
- Recipe issues: https://github.com/pyodide/pyodide-recipes/issues

## References

- [PhiSpy GitHub](https://github.com/linsalrob/PhiSpy)
- [PhiSpy on PyPI](https://pypi.org/project/PhiSpy/)
- [Pyodide Documentation](https://pyodide.org/en/stable/)
- [Pyodide Recipes Repository](https://github.com/pyodide/pyodide-recipes)
- [Pyodide Recipe Format](https://pyodide.org/en/stable/development/meta-yaml.html)

## Citation

If you use PhiSpy in your research, please cite:

```
Akhter, S., Aziz, R. K., & Edwards, R. A. (2012). 
PhiSpy: a novel algorithm for finding prophages in bacterial genomes.
Nucleic Acids Research, 40(16), e126.
```

---

**Status:** ✅ Ready for submission to pyodide-recipes

**Last Updated:** February 2026
