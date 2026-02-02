"""
Test script for PhiSpy in Pyodide environment.
This script tests basic imports and functionality of PhiSpy.
"""

def test_imports():
    """Test that all essential PhiSpy modules can be imported."""
    print("Testing PhiSpy imports...")
    
    # Test main module imports
    import PhiSpyModules
    print("✓ PhiSpyModules imported successfully")
    
    # Test C++ extension import
    import PhiSpyRepeatFinder
    print("✓ PhiSpyRepeatFinder (C++ extension) imported successfully")
    
    # Test submodules
    from PhiSpyModules import main
    print("✓ PhiSpyModules.main imported successfully")
    
    from PhiSpyModules import helper_functions
    print("✓ PhiSpyModules.helper_functions imported successfully")
    
    from PhiSpyModules import classification
    print("✓ PhiSpyModules.classification imported successfully")
    
    print("\nAll imports successful! PhiSpy is ready to use in Pyodide.")
    return True


def test_dependencies():
    """Test that PhiSpy dependencies are available."""
    print("\nTesting PhiSpy dependencies...")
    
    import numpy
    print(f"✓ numpy {numpy.__version__}")
    
    import Bio
    print(f"✓ biopython {Bio.__version__}")
    
    import sklearn
    print(f"✓ scikit-learn {sklearn.__version__}")
    
    print("\nAll dependencies available!")
    return True


def test_basic_functionality():
    """Test basic PhiSpy functionality."""
    print("\nTesting basic PhiSpy functionality...")
    
    # Test that we can access version info
    from PhiSpyModules import version
    print(f"✓ PhiSpy version: {version.__version__}")
    
    # Test repeat finder is callable
    import PhiSpyRepeatFinder
    print("✓ PhiSpyRepeatFinder module is accessible")
    
    print("\nBasic functionality tests passed!")
    return True


if __name__ == "__main__":
    try:
        test_imports()
        test_dependencies()
        test_basic_functionality()
        print("\n" + "="*50)
        print("All tests passed! PhiSpy is working in Pyodide.")
        print("="*50)
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        raise
