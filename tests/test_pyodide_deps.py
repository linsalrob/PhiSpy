"""
Test suite for pyodide_deps.py
"""

import unittest
import sys
from unittest.mock import patch, MagicMock, AsyncMock
import asyncio


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class PyodideDepsTest(unittest.TestCase):
    
    def test_is_pyodide_returns_false_on_cpython(self):
        """Test that is_pyodide() returns False on standard CPython"""
        from PhiSpyModules.pyodide_deps import is_pyodide
        
        # On standard CPython, sys.platform should not be "emscripten"
        self.assertFalse(is_pyodide())
        self.assertNotEqual(sys.platform, "emscripten")
    
    @patch('sys.platform', 'emscripten')
    def test_is_pyodide_returns_true_on_pyodide(self):
        """Test that is_pyodide() returns True when sys.platform is emscripten"""
        from PhiSpyModules.pyodide_deps import is_pyodide
        
        self.assertTrue(is_pyodide())
    
    def test_ensure_bcbio_gff_sync_noop_on_cpython(self):
        """Test that ensure_bcbio_gff_sync() is a no-op on CPython"""
        from PhiSpyModules.pyodide_deps import ensure_bcbio_gff_sync
        
        # Should not raise any errors
        ensure_bcbio_gff_sync()
    
    def test_ensure_pyodide_deps_async_noop_on_cpython(self):
        """Test that ensure_pyodide_deps() is a no-op on CPython"""
        from PhiSpyModules.pyodide_deps import ensure_pyodide_deps
        
        # Should complete without errors
        asyncio.run(ensure_pyodide_deps())
    
    @patch('sys.platform', 'emscripten')
    def test_ensure_bcbio_gff_sync_returns_if_already_installed(self):
        """Test that ensure_bcbio_gff_sync returns early if BCBio is importable"""
        from PhiSpyModules.pyodide_deps import ensure_bcbio_gff_sync
        
        # Mock BCBio.GFF as already installed
        mock_bcbio = MagicMock()
        with patch.dict('sys.modules', {'BCBio': mock_bcbio, 'BCBio.GFF': mock_bcbio.GFF}):
            # Should return without attempting to install
            ensure_bcbio_gff_sync()
    
    @patch('sys.platform', 'emscripten')
    @patch('asyncio.get_running_loop')
    def test_ensure_bcbio_gff_sync_raises_with_running_loop(self, mock_get_loop):
        """Test that ensure_bcbio_gff_sync raises error when event loop is running"""
        from PhiSpyModules.pyodide_deps import ensure_bcbio_gff_sync
        
        # Mock BCBio.GFF import to raise ImportError (simulating not installed)
        def mock_import(name, *args, **kwargs):
            if 'BCBio' in name:
                raise ImportError(f"No module named '{name}'")
            # Use actual __import__ for everything else
            return __import__(name, *args, **kwargs)
        
        with patch('builtins.__import__', side_effect=mock_import):
            # Mock a running event loop
            mock_loop = MagicMock()
            mock_get_loop.return_value = mock_loop
            
            # Should raise RuntimeError with helpful message
            with self.assertRaises(RuntimeError) as context:
                ensure_bcbio_gff_sync()
            
            self.assertIn("Cannot install bcbio-gff synchronously", str(context.exception))
            self.assertIn("await phispy.ensure_pyodide_deps()", str(context.exception))
    
    @patch('sys.platform', 'emscripten')
    @patch('asyncio.get_running_loop')
    @patch('asyncio.run')
    def test_ensure_bcbio_gff_sync_uses_asyncio_run_when_no_loop(
        self, mock_asyncio_run, mock_get_loop
    ):
        """Test that ensure_bcbio_gff_sync uses asyncio.run when no loop is running"""
        from PhiSpyModules.pyodide_deps import ensure_bcbio_gff_sync
        
        # Simulate BCBio not being installed
        def mock_import(name, *args, **kwargs):
            if 'BCBio' in name:
                raise ImportError(f"No module named '{name}'")
            return __import__(name, *args, **kwargs)
        
        with patch('builtins.__import__', side_effect=mock_import):
            # Simulate no running event loop
            mock_get_loop.side_effect = RuntimeError("no running event loop")
            
            # Call the function
            ensure_bcbio_gff_sync()
            
            # Verify asyncio.run was called
            mock_asyncio_run.assert_called_once()
    
    @patch('sys.platform', 'emscripten')
    async def test_ensure_bcbio_gff_async_returns_if_already_installed(self):
        """Test that ensure_bcbio_gff_async returns early if BCBio is importable"""
        from PhiSpyModules.pyodide_deps import ensure_bcbio_gff_async
        
        # Mock BCBio.GFF as already installed
        mock_bcbio = MagicMock()
        with patch.dict('sys.modules', {'BCBio': mock_bcbio, 'BCBio.GFF': mock_bcbio.GFF}):
            # Should return without attempting to install
            await ensure_bcbio_gff_async()
    
    @patch('sys.platform', 'emscripten')
    async def test_ensure_bcbio_gff_async_installs_via_micropip(self):
        """Test that ensure_bcbio_gff_async installs via micropip when needed"""
        from PhiSpyModules.pyodide_deps import ensure_bcbio_gff_async
        
        # Mock micropip
        mock_micropip = MagicMock()
        mock_micropip.install = AsyncMock()
        
        # Simulate BCBio not being installed initially
        def mock_import(name, *args, **kwargs):
            if 'BCBio' in name:
                raise ImportError(f"No module named '{name}'")
            if name == 'micropip':
                return mock_micropip
            return __import__(name, *args, **kwargs)
        
        with patch('builtins.__import__', side_effect=mock_import):
            await ensure_bcbio_gff_async()
            
            # Verify micropip.install was called with correct package
            mock_micropip.install.assert_called_once_with("bcbio-gff")
    
    @patch('sys.platform', 'emscripten')
    async def test_ensure_bcbio_gff_async_raises_if_no_micropip(self):
        """Test that ensure_bcbio_gff_async raises ImportError if micropip is not available"""
        from PhiSpyModules.pyodide_deps import ensure_bcbio_gff_async
        
        # Simulate both BCBio and micropip not being available
        def mock_import(name, *args, **kwargs):
            if 'BCBio' in name or name == 'micropip':
                raise ImportError(f"No module named '{name}'")
            return __import__(name, *args, **kwargs)
        
        with patch('builtins.__import__', side_effect=mock_import):
            with self.assertRaises(ImportError) as context:
                await ensure_bcbio_gff_async()
            
            self.assertIn("micropip is required", str(context.exception))
    
    def test_import_bcbio_gff_returns_gff_module(self):
        """Test that _import_bcbio_gff returns the GFF module"""
        from PhiSpyModules.pyodide_deps import _import_bcbio_gff
        
        # This test assumes bcbio-gff is installed (which it should be in the test environment)
        GFF = _import_bcbio_gff()
        
        # Verify we got a module with expected attributes
        self.assertTrue(hasattr(GFF, 'write') or hasattr(GFF, 'parse'))
    
    def test_import_bcbio_gff_raises_if_not_installed(self):
        """Test that _import_bcbio_gff raises ImportError if bcbio-gff is not installed"""
        from PhiSpyModules.pyodide_deps import _import_bcbio_gff
        
        # Mock BCBio as not installed
        def mock_import(name, *args, **kwargs):
            if 'BCBio' in name:
                raise ImportError(f"No module named '{name}'")
            return __import__(name, *args, **kwargs)
        
        with patch('builtins.__import__', side_effect=mock_import):
            with self.assertRaises(ImportError):
                _import_bcbio_gff()
    
    def test_public_api_exports(self):
        """Test that the public API functions are properly exported"""
        import PhiSpyModules
        
        # Check that the new functions are accessible
        self.assertTrue(hasattr(PhiSpyModules, 'is_pyodide'))
        self.assertTrue(hasattr(PhiSpyModules, 'ensure_pyodide_deps'))
        
        # Verify they are the correct functions
        from PhiSpyModules.pyodide_deps import is_pyodide, ensure_pyodide_deps
        self.assertEqual(PhiSpyModules.is_pyodide, is_pyodide)
        self.assertEqual(PhiSpyModules.ensure_pyodide_deps, ensure_pyodide_deps)


if __name__ == '__main__':
    unittest.main()
