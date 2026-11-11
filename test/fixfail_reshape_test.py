import contextlib
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'hm'))
import fixfail


class RepairReshapeTest(unittest.TestCase):
    def test_ligand_alias_and_spaces_use_canonical_docking_filename(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'out/OR1/OR1A1/OR1A1~Example_ligand.active.model1.pdb'
            source.parent.mkdir(parents=True)
            source.write_text('END\n')
            odor = {'full_name': 'Example ligand', 'oid': 'example'}
            with patch.object(fixfail, 'root_dir', str(root)), \
                    patch('data.odorutils.find_odorant', return_value=odor):
                for name in ('Example ligand', 'example', 'alias'):
                    candidate, resolved = fixfail.repair_input(name, 'OR1A1', 'OR1', 'active')
                    self.assertEqual(candidate, source)
                    self.assertEqual(resolved, odor)

    def test_explicit_pdb_can_resolve_without_an_odorant_name(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / 'protein.pdb'
            source.write_text('END\n')
            with patch('data.odorutils.find_odorant', return_value=False):
                candidate, odor = fixfail.repair_input(str(source), 'OR1A1', 'OR1', 'active')
            self.assertEqual(candidate, source)
            self.assertFalse(odor)

    def test_family_scripts_and_missing_script_fallback(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = root / 'data'
            data.mkdir()
            with patch.object(fixfail, 'root_dir', str(root)):
                message = io.StringIO()
                with contextlib.redirect_stdout(message):
                    self.assertIsNone(fixfail.repair_reshape_script('TAAR5', 'active'))
                self.assertIn('repairing the supplied conformation directly', message.getvalue())
                for receptor, stem in [('OR1A1', 'OR_ClassII'), ('OR51E2', 'OR_ClassI'),
                                       ('TAAR5', 'TAAR'), ('VN1R1', 'VN1R')]:
                    for mode, state in [('active', 'a'), ('inactive', 'i')]:
                        with self.subTest(receptor=receptor, mode=mode):
                            script = data / f'{stem}_{state}.ic'
                            script.write_text('')
                            self.assertEqual(fixfail.repair_reshape_script(receptor, mode), script)
                            self.assertIsNone(fixfail.repair_reshape_script(receptor, mode, skip=True))

    def test_taar5_entry_point_reaches_template_copy_without_reshape(self):
        class TemplateReached(Exception):
            pass

        previous = os.getcwd()
        argv = sys.argv[:]
        try:
            with tempfile.TemporaryDirectory() as temporary, \
                    patch.object(fixfail, 'prepare'), \
                    patch.object(fixfail.subprocess, 'run'), \
                    patch('data.globals.wait_cool_cpu'), \
                    patch('data.protutils.ensure_bw_metadata'), \
                    patch.object(fixfail.shutil, 'copyfile', side_effect=TemplateReached) as copy, \
                    contextlib.redirect_stdout(io.StringIO()):
                source = Path(temporary) / 'TAAR5~N-methylpiperidine.active.model1.pdb'
                source.write_bytes((ROOT / 'testdata/reec-small.pdb').read_bytes())
                with self.assertRaises(TemplateReached):
                    fixfail.main(['TAAR5', str(source), '--output', str(Path(temporary) / 'result.pdb')])
                self.assertEqual(Path(copy.call_args.args[0]), source)
        finally:
            os.chdir(previous)
            sys.argv = argv


if __name__ == '__main__':
    unittest.main()
