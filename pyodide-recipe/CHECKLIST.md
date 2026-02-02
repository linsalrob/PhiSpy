# Pre-Submission Checklist for pyodide-recipes

Use this checklist before submitting the PhiSpy recipe to the pyodide-recipes repository.

## Prerequisites

- [ ] Python 3.13+ installed
- [ ] Git installed with submodule support
- [ ] GitHub account created
- [ ] Forked pyodide-recipes repository

## Recipe Preparation

- [ ] `meta.yaml` file is complete and correct
  - [ ] Package name matches: `phispy`
  - [ ] Version is current: `5.0.1`
  - [ ] Source URL is correct
  - [ ] SHA256 hash is verified
  - [ ] Dependencies are listed correctly
  - [ ] Build flags are appropriate (`-fwasm-exceptions`)

- [ ] `test_phispy.py` is included and tests critical imports
  - [ ] Tests PhiSpyModules import
  - [ ] Tests PhiSpyRepeatFinder (C++ extension) import
  - [ ] Tests dependency availability

## Local Testing

- [ ] Cloned pyodide-recipes with `--recurse-submodules`
- [ ] Copied recipe files to `packages/phispy/`
- [ ] Installed pyodide-build: `pip install ./pyodide-build`
- [ ] Built package successfully: `pyodide build-recipes phispy`
- [ ] No build errors or warnings
- [ ] Import tests pass

## Build Verification

- [ ] C++ extension compiles without errors
- [ ] All dependencies resolve correctly
- [ ] Package size is reasonable
- [ ] Import test succeeds:
  ```python
  import PhiSpyModules
  import PhiSpyRepeatFinder
  ```

## Documentation Review

- [ ] Recipe follows pyodide-recipes conventions
- [ ] Comments in `meta.yaml` explain non-obvious choices
- [ ] Test file includes docstrings

## Git Preparation

- [ ] Created feature branch: `add-phispy-recipe`
- [ ] Added recipe files to git
- [ ] Committed with descriptive message
- [ ] Pushed to your fork

## Pull Request Content

### PR Title
```
Add PhiSpy package recipe
```

### PR Description Template
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
- bcbio-gff (pure Python, installable via micropip)

### Build Notes
- Includes C++ extension module (PhiSpyRepeatFinder)
- Configured with `-fwasm-exceptions` for WebAssembly compatibility
- All tests pass successfully

### Testing
- [x] Built successfully with `pyodide build-recipes phispy`
- [x] Import tests pass
- [x] C++ extension compiles and loads correctly
- [x] All dependencies available

### Maintainer Checklist
- [x] Follows recipe conventions
- [x] Tests included
- [x] Build tested locally
```

- [ ] PR description is complete
- [ ] All checklist items marked

## Submission

- [ ] Opened pull request on pyodide-recipes
- [ ] Assigned to appropriate reviewers (if known)
- [ ] Monitoring for review comments
- [ ] Ready to respond to feedback

## Post-Submission

- [ ] Respond to reviewer comments within 48 hours
- [ ] Make requested changes promptly
- [ ] Test changes after each iteration
- [ ] Wait for approval and merge

## After Merge

- [ ] Update PhiSpy documentation about Pyodide availability
- [ ] Announce Pyodide support to users
- [ ] Monitor for issues
- [ ] Plan for future version updates

## Notes

- Keep this checklist for future PhiSpy version updates
- Each new version will need a recipe update
- Document any special issues or solutions discovered

---

**Checklist Version:** 1.0
**Date:** February 2026
