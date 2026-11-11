"""Offline input/security tests. Inputs: none; outputs: unittest status.
Example: python scripts/test_alphagenome2.py --help or python scripts/test_alphagenome2.py
"""
from __future__ import annotations
import argparse
import contextlib
import importlib.util
import io
from pathlib import Path
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('ag2_helper', Path(__file__).with_name('alphagenome2.py'))
assert SPEC and SPEC.loader
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


class Inputs(unittest.TestCase):
    """Reject unsafe/invalid requests before accessing the network."""
    def test_defaults(self) -> None:
        self.assertEqual(HELPER.parse_args(['conditions', '--output', 'out']).organism, 'human')

    def test_invalid_inputs(self) -> None:
        for extra in [['--timeout', '0'], ['--timeout', '-1'], ['--timeout', '301'],
                      ['--timeout', 'nan'], ['--organism', 'fish']]:
            with self.subTest(extra=extra), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                HELPER.parse_args(['conditions', '--output', 'out'] + extra)

    def test_mouse_smoke(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            HELPER.parse_args(['smoke', '--organism', 'mouse', '--output', 'out'])

    def test_no_consent(self) -> None:
        with self.assertRaisesRegex(SystemExit, 'Confirm noncommercial'):
            HELPER.run(HELPER.parse_args(['conditions', '--output', 'out']))

    def test_missing_key(self) -> None:
        with patch.dict('os.environ', {}, clear=True), self.assertRaisesRegex(SystemExit, 'Set ALPHAGENOME_API_KEY'):
            HELPER.run(HELPER.parse_args(['conditions', '--output', 'out', '--confirm-noncommercial']))


def main() -> None:
    """Run offline tests without third-party dependencies or a credential."""
    argparse.ArgumentParser(description=__doc__).parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Inputs))
    if not result.wasSuccessful():
        raise SystemExit(1)
    print(f'✓ Offline tests completed: {result.testsRun} tests')


if __name__ == '__main__':
    main()
