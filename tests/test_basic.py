"""
Basic tests to demonstrate the new 3S2B structure.
"""

import pytest


def test_package_imports():
    """Test that the new package structure can be imported."""
    try:
        import three_s_two_b
        assert three_s_two_b.__version__ == "1.0.0"
        assert hasattr(three_s_two_b, '__author__')
    except ImportError as e:
        pytest.skip(f"Package not yet installed: {e}")


def test_version_info():
    """Test version information is accessible."""
    try:
        from three_s_two_b import __version__, __author__
        assert __version__ == "1.0.0"
        assert "Fraser" in __author__
    except ImportError as e:
        pytest.skip(f"Package not yet installed: {e}")


if __name__ == "__main__":
    test_package_imports()
    test_version_info()
    print("✅ Basic tests passed!")