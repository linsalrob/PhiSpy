# PhiSpy Pyodide Recipe

This directory contains the Pyodide recipe for building PhiSpy for WebAssembly.

## About

PhiSpy is a prophage finder that uses multiple metrics to identify prophages in bacterial and archaeal genomes. This recipe enables PhiSpy to run in web browsers and Node.js environments via Pyodide.

## Recipe Contents

- `meta.yaml` - The main recipe file defining how to build PhiSpy for Pyodide

## Dependencies

PhiSpy requires the following packages, all of which are available in Pyodide:
- **biopython** - For biological sequence handling
- **numpy** - For numerical computations
- **scikit-learn** - For machine learning functionality

Note: PhiSpy also depends on `bcbio-gff`, which is a pure Python package and can be installed via micropip if needed.

## C++ Extension

PhiSpy includes a C++ extension module (`PhiSpyRepeatFinder`) for efficient repeat finding. The recipe is configured to compile this extension using Emscripten with WebAssembly exception support.

## Building

To build this recipe for Pyodide:

1. Clone the pyodide-recipes repository:
   ```bash
   git clone --recurse-submodules https://github.com/pyodide/pyodide-recipes.git
   cd pyodide-recipes
   ```

2. Copy this recipe to the packages directory:
   ```bash
   cp -r /path/to/PhiSpy/pyodide-recipe packages/phispy
   ```

3. Build the package:
   ```bash
   pip install ./pyodide-build
   pyodide build-recipes phispy
   ```

## Usage in Pyodide

Once built and published, PhiSpy can be used in Pyodide:

```python
import micropip
await micropip.install('phispy')
# If bcbio-gff is needed and not pre-installed:
await micropip.install('bcbio-gff')

# Use PhiSpy
from PhiSpyModules import main
# ... use PhiSpy functionality
```

## Submitting to pyodide-recipes

To submit this recipe to the official pyodide-recipes repository:

1. Fork the [pyodide-recipes repository](https://github.com/pyodide/pyodide-recipes)
2. Copy the contents of this directory to `packages/phispy/` in your fork
3. Test the build locally
4. Submit a pull request

## License

PhiSpy is licensed under the MIT License. See the main repository LICENSE file for details.

## References

- [PhiSpy GitHub](https://github.com/linsalrob/PhiSpy)
- [Pyodide Documentation](https://pyodide.org/en/stable/)
- [Pyodide Recipe Guide](https://pyodide.org/en/stable/development/building-packages-using-recipe.html)
