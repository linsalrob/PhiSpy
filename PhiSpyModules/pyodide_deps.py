"""
Module for handling Pyodide-specific dependencies.

This module provides utilities for lazily installing Python packages
via micropip when running under Pyodide (WebAssembly environment).
"""

import sys
import asyncio


def is_pyodide() -> bool:
    """
    Check if running under Pyodide.
    
    Returns:
        bool: True if running under Pyodide, False otherwise.
    """
    return sys.platform == "emscripten"


async def ensure_bcbio_gff_async() -> None:
    """
    Asynchronously ensure bcbio-gff is available in Pyodide.
    
    If running under Pyodide and bcbio-gff is not installed,
    this function will install it via micropip.
    
    This function should be called with await in async contexts:
        await ensure_bcbio_gff_async()
    
    Raises:
        ImportError: If micropip is not available in Pyodide environment.
    """
    if not is_pyodide():
        return
    
    # Check if already installed
    try:
        import BCBio.GFF  # noqa: F401
        return  # Already installed
    except ImportError:
        pass
    
    # Install via micropip
    try:
        import micropip  # type: ignore
        await micropip.install("bcbio-gff")
    except ImportError as e:
        raise ImportError(
            "micropip is required for installing packages in Pyodide but is not available"
        ) from e


def ensure_bcbio_gff_sync() -> None:
    """
    Synchronously ensure bcbio-gff is available in Pyodide.
    
    This is a fallback for sync contexts. It will:
    - Do nothing if not running under Pyodide
    - Return immediately if bcbio-gff is already installed
    - Attempt to install via micropip if needed
    
    In Jupyter/Pyodide notebooks with an already-running event loop,
    this may raise a RuntimeError. In such cases, use the async version:
        await ensure_pyodide_deps()
    
    Raises:
        RuntimeError: If running in Pyodide with an active event loop that
                     prevents synchronous installation. Use ensure_pyodide_deps() instead.
    """
    if not is_pyodide():
        return
    
    # Check if already installed
    try:
        import BCBio.GFF  # noqa: F401
        return  # Already installed
    except ImportError:
        pass
    
    # Try to check if there's a running event loop
    try:
        loop = asyncio.get_running_loop()
        # There's a running loop, we can't use asyncio.run()
        raise RuntimeError(
            "Cannot install bcbio-gff synchronously in Pyodide with an active event loop. "
            "Please use: await phispy.ensure_pyodide_deps() before using GFF output features."
        )
    except RuntimeError as e:
        # No running loop (expected for most CLI usage)
        if "no running event loop" in str(e).lower() or "There is no current event loop" in str(e):
            # Safe to create a new event loop
            asyncio.run(ensure_bcbio_gff_async())
        else:
            # Re-raise if it's our custom error message
            raise


async def ensure_pyodide_deps() -> None:
    """
    Public async helper to ensure all Pyodide dependencies are installed.
    
    In Pyodide environments (e.g., JupyterLite notebooks), call this once
    before using features that require optional dependencies:
    
        import phispy
        await phispy.ensure_pyodide_deps()
    
    This is the recommended way to handle dependencies in async contexts
    like Jupyter notebooks running in Pyodide.
    """
    await ensure_bcbio_gff_async()


def _import_bcbio_gff():
    """
    Import and return BCBio.GFF module.
    
    This function performs the actual import of BCBio.GFF after ensuring
    the dependency is available. It should be called after ensure_bcbio_gff_sync()
    or in contexts where bcbio-gff is already installed.
    
    Returns:
        module: The BCBio.GFF module.
        
    Raises:
        ImportError: If bcbio-gff is not installed.
    """
    from BCBio import GFF
    return GFF
