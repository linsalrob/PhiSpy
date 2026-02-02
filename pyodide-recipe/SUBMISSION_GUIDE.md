# Submitting PhiSpy to the Pyodide Recipes Repository

This guide explains how to submit the PhiSpy recipe to the official [pyodide-recipes repository](https://github.com/pyodide/pyodide-recipes).

## Prerequisites

Before submitting, ensure you have:
- Python 3.13 or higher installed
- Git installed with submodule support
- A GitHub account
- Basic knowledge of Git and GitHub workflows

## Step-by-Step Submission Process

### 1. Fork the pyodide-recipes Repository

1. Go to https://github.com/pyodide/pyodide-recipes
2. Click the "Fork" button in the top-right corner
3. This creates a copy of the repository in your GitHub account

### 2. Clone Your Fork

```bash
git clone --recurse-submodules https://github.com/YOUR_USERNAME/pyodide-recipes.git
cd pyodide-recipes
```

**Important:** Use the `--recurse-submodules` flag to ensure all submodules are initialized.

If you already cloned without this flag, run:
```bash
git submodule update --init --recursive
```

### 3. Copy the PhiSpy Recipe

Copy the recipe files from this repository to the pyodide-recipes packages directory:

```bash
# From the PhiSpy repository root
cp -r pyodide-recipe /path/to/pyodide-recipes/packages/phispy
```

Or manually create the `packages/phispy` directory and copy these files:
- `meta.yaml`
- `test_phispy.py` (optional but recommended)

### 4. Install pyodide-build

```bash
cd /path/to/pyodide-recipes
pip install ./pyodide-build
```

### 5. Test the Build Locally

Before submitting, test that the recipe builds successfully:

```bash
pyodide build-recipes phispy
```

This will:
- Download the PhiSpy source package
- Compile the C++ extension for WebAssembly
- Build the wheel package
- Run the tests (if provided)

If the build fails, review error messages and adjust the `meta.yaml` accordingly.

### 6. Verify the Package

After building, test the package:

```bash
# The built package will be in the dist directory
cd dist
python -m http.server 8000
```

Then open your browser to test in Pyodide's console.

### 7. Create a Pull Request

Once the build succeeds:

```bash
cd /path/to/pyodide-recipes
git checkout -b add-phispy-recipe
git add packages/phispy/
git commit -m "Add PhiSpy package recipe

PhiSpy is a prophage finder that identifies prophages in bacterial
and archaeal genomes using multiple metrics.

Package includes:
- C++ extension module (PhiSpyRepeatFinder) for efficient repeat finding
- Dependencies: biopython, numpy, scikit-learn (all available in Pyodide)

Tested with PhiSpy 5.0.1 from PyPI.
"
git push origin add-phispy-recipe
```

### 8. Submit the Pull Request on GitHub

1. Go to your fork on GitHub
2. Click "Compare & pull request"
3. Fill in the PR description with:
   - Package name and version
   - Brief description of what PhiSpy does
   - List of dependencies
   - Any special build requirements
   - Test results

Example PR description:
```markdown
## Add PhiSpy package

**Package:** phispy 5.0.1
**Homepage:** https://github.com/linsalrob/PhiSpy
**PyPI:** https://pypi.org/project/PhiSpy/

### Description
PhiSpy identifies prophages in Bacterial and Archaeal genomes using 
multiple metrics including genome sequence composition, gene orientation, 
and protein function predictions.

### Dependencies
- biopython (available in Pyodide)
- numpy (available in Pyodide)
- scikit-learn (available in Pyodide)
- bcbio-gff (pure Python, included in recipe)

### Build Notes
- Includes C++ extension module (PhiSpyRepeatFinder)
- Configured with `-fwasm-exceptions` for WebAssembly compatibility
- All tests pass successfully

### Testing
Tested imports and basic functionality successfully in Pyodide environment.
```

4. Submit the pull request

### 9. Respond to Review Feedback

Maintainers may request changes. Common requests include:
- Adjusting build flags
- Adding or removing dependencies
- Improving test coverage
- Documentation updates

Make requested changes and push to your branch:
```bash
# Make changes to packages/phispy/meta.yaml
git add packages/phispy/
git commit -m "Address review feedback: [describe changes]"
git push origin add-phispy-recipe
```

## Tips for Success

1. **Test thoroughly** - Ensure the build completes without errors
2. **Minimal recipe** - Keep the recipe as simple as possible
3. **Follow conventions** - Look at similar packages in pyodide-recipes for examples
4. **Document** - Add comments in meta.yaml for non-obvious settings
5. **Be responsive** - Reply to reviewer comments promptly

## Troubleshooting

### Build Fails with Compiler Errors

If the C++ extension fails to compile:
- Check the `cxxflags` and `ldflags` in meta.yaml
- Review Emscripten compatibility issues
- May need to add patches for WebAssembly compatibility

### Import Errors After Build

If imports fail during testing:
- Verify `top-level` in meta.yaml lists all importable modules
- Check that dependencies are correctly listed
- Ensure the package structure is preserved

### Dependency Not Available

If a dependency isn't in Pyodide:
- For pure Python packages: note in PR that micropip can install it
- For packages with C extensions: may need to submit a separate recipe first

## Additional Resources

- [Pyodide Building Packages Guide](https://pyodide.org/en/stable/development/building-packages-using-recipe.html)
- [meta.yaml Specification](https://pyodide.org/en/stable/development/meta-yaml.html)
- [pyodide-recipes Repository](https://github.com/pyodide/pyodide-recipes)
- [Pyodide Packages List](https://pyodide.org/en/stable/usage/packages-in-pyodide.html)

## Support

If you encounter issues:
- Check existing issues in the pyodide-recipes repository
- Ask questions in the Pyodide community channels
- Review similar package recipes for examples
