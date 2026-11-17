"""Short and long modeling options select identical parser settings."""
import contextlib
import io
from pathlib import Path
import re
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import reec_backend
import modeling_backend


class ModelAliases(unittest.TestCase):
    def test_every_long_option_has_a_unique_documented_short_alias(self):
        for repair in (False, True):
            parser = reec_backend.parser(repair)
            help_text = parser.format_help()
            used = set()
            for action in parser._actions:
                longs = [flag for flag in action.option_strings if flag.startswith('--')]
                if not longs:
                    continue
                shorts = [flag for flag in action.option_strings if re.fullmatch('-[A-Za-z]{1,2}', flag)]
                self.assertEqual(len(shorts), 1, action.dest)
                self.assertNotIn(shorts[0], used)
                used.add(shorts[0])
                self.assertIn(shorts[0], help_text)
                self.assertIn(longs[0], help_text)
                self.assertTrue(action.help, action.dest)
                if action.dest == 'help':
                    continue
                value = ([] if action.nargs == 0 else ['208:ND2', 'O10', '2.66'] if action.nargs == 3
                         else [str(action.choices[0])] if action.choices
                         else ['7'] if action.type in (int, float) else ['example'])
                long = parser.parse_args(['OR52E8', longs[0], *value])
                short = parser.parse_args(['OR52E8', shorts[0], *value])
                self.assertEqual(vars(long), vars(short))

    def test_backend_alias_selects_the_engine(self):
        for flag in (['-k', 'modeller'], ['-kmodeller'], ['-k=modeller']):
            self.assertFalse(reec_backend.requested(flag, default='reec'))
            self.assertEqual(reec_backend.legacy_backend_args(flag), [])
            self.assertEqual(modeling_backend.options(['OR52E8', *flag], repair=True).backend, 'modeller')

    def test_no_arguments_show_complete_help(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), self.assertRaises(SystemExit):
            modeling_backend.options([], repair=True)
        for action in reec_backend.parser(True)._actions:
            for flag in action.option_strings:
                self.assertIn(flag, output.getvalue())

    def test_fixfail_help_runs_without_a_modeling_job(self):
        for flag in ('--help', '-h'):
            result = subprocess.run([sys.executable, str(ROOT / 'hm/fixfail.py'), flag],
                                    cwd=ROOT, capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('-lc', result.stdout)
            self.assertIn('--ligand-contact', result.stdout)


if __name__ == '__main__':
    unittest.main()
