"""Dependency-free contract tests; these do not validate model inference."""
import subprocess
import sys
import unittest
from pathlib import Path
from smoke_fold import require_finite, SEQUENCE


class ContractTests(unittest.TestCase):
    def test_public_smoke_sequence(self):
        self.assertEqual(len(SEQUENCE), 76)
        self.assertTrue(set(SEQUENCE) <= set("ACDEFGHIKLMNPQRSTVWY"))

    def test_finite(self):
        require_finite([0, 1, 0.5])

    def test_bad_values(self):
        for values in [[], [float("nan")], [float("inf")]]:
            with self.assertRaises(ValueError):
                require_finite(values)

    def test_help_no_heavy_dependencies(self):
        result = subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_fold.py")), "--help"], capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
