import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class DockLauncherTest(unittest.TestCase):
    def launch(self, build_status):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, content in {
                'make': '#!/bin/sh\nprintf "%s\\n" "$@" > "$LAUNCH_LOG"\nexit "$BUILD_STATUS"\n',
                'python3': '#!/bin/sh\nprintf "%s\\n" "$@" > "$RUN_LOG"\n',
            }.items():
                path = root / name
                path.write_text(content)
                path.chmod(0o755)
            env = dict(os.environ, PATH=str(root)+os.pathsep+os.environ['PATH'],
                       LAUNCH_LOG=str(root/'build.log'), RUN_LOG=str(root/'run.log'),
                       BUILD_STATUS=str(build_status))
            result = subprocess.run(['bash', str(ROOT/'dock.sh'), 'OR5K1', '(1R,2S,5R)-(-)-menthol'],
                                    env=env, cwd=ROOT, capture_output=True, text=True, timeout=15)
            self.assertIn('bin/ringflip', (root/'build.log').read_text().splitlines())
            if build_status:
                self.assertEqual(result.returncode, build_status)
                self.assertFalse((root/'run.log').exists())
            else:
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual((root/'run.log').read_text().splitlines(),
                                 ['run.py', 'OR5K1', '(1R,2S,5R)-(-)-menthol'])

    def test_builds_the_writer_before_starting_docking(self):
        self.launch(0)

    def test_failed_build_does_not_launch_stale_programs(self):
        self.launch(7)


if __name__ == '__main__':
    unittest.main()
