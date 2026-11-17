"""The pencil verifier must reject missing or repeated local receipts."""
import shutil
import subprocess
import sys
import tempfile
import unittest
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ReceiptCoverageTest(unittest.TestCase):
    def test_complete_and_incomplete_local_receipts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / 'scripts' / 'pencil_selmer'
            shutil.copytree(ROOT / 'scripts' / 'pencil_selmer', destination,
                            ignore=shutil.ignore_patterns('__pycache__'))
            verifier = destination / 'verify.py'

            def run():
                return subprocess.run([sys.executable, str(verifier)], cwd=root,
                                      capture_output=True, text=True)

            valid = run()
            self.assertEqual(valid.returncode, 0, valid.stderr)
            for name in ('local_5_7.txt', 'local_odd.txt', 'local_remaining.txt'):
                path = destination / name
                original = path.read_text()
                blocks = original.split('PRIME ')
                variants = {
                    'empty': '',
                    'missing': blocks[0] + ''.join('PRIME ' + b for b in blocks[1:-1]),
                    'duplicate': original + 'PRIME ' + blocks[1],
                    'missing row': re.sub(r'\[[^\]]+\]', '', original, count=1),
                    'invalid bit': re.sub(r'\[\s*[01]', '[2', original, count=1),
                }
                for case, text in variants.items():
                    with self.subTest(file=name, case=case):
                        path.write_text(text)
                        result = run()
                        self.assertNotEqual(result.returncode, 0, result.stdout)
                        self.assertIn('AssertionError', result.stderr)
                path.write_text(original)


if __name__ == '__main__':
    unittest.main()
