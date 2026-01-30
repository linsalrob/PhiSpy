"""
Test suite for formatting.py
"""

import unittest
from PhiSpyModules.formatting import Colors
from PhiSpyModules.errors import ColorNotFoundError


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class FormattingTest(unittest.TestCase):
    
    def test_colors_class_has_required_colors(self):
        """Test that Colors class contains all expected color definitions"""
        colors = Colors()
        required_colors = [
            'HEADER', 'OKBLUE', 'OKGREEN', 'WARNING', 'FAIL', 
            'ENDC', 'BOLD', 'UNDERLINE', 'PINK', 'BLUE', 
            'GREEN', 'YELLOW', 'RED', 'WHITE'
        ]
        for color_name in required_colors:
            with self.subTest(color=color_name):
                self.assertIn(color_name, colors.color, 
                            f"Color '{color_name}' not found in Colors.color dictionary")
    
    def test_colors_ansi_codes(self):
        """Test that color codes start with ANSI escape sequence"""
        colors = Colors()
        for color_name, color_code in colors.color.items():
            with self.subTest(color=color_name):
                self.assertTrue(color_code.startswith('\033['), 
                              f"Color code for '{color_name}' doesn't start with ANSI escape sequence")


if __name__ == '__main__':
    unittest.main()
