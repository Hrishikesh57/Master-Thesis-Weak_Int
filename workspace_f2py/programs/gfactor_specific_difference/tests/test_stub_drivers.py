from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

from _support import PROJECT_ROOT


class StubDriverTest(unittest.TestCase):
    def test_placeholder_scripts_run(self) -> None:
        script_dir = PROJECT_ROOT / "scripts"
        scripts = [
            "run_hlike_pb.py",
            "run_hlike_u.py",
            "run_lilike_pb.py",
            "run_lilike_u.py",
            "run_specific_difference.py",
        ]
        for script_name in scripts:
            with self.subTest(script=script_name):
                completed = subprocess.run(
                    [sys.executable, str(script_dir / script_name), "--placeholder", "--json"],
                    cwd=PROJECT_ROOT,
                    text=True,
                    capture_output=True,
                    check=True,
                )
                self.assertIn("placeholder", completed.stdout)


if __name__ == "__main__":
    unittest.main()
