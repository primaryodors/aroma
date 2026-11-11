"""Completed diagnostics leave receptor storage without losing validation evidence."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hm'))
from modeling_backend import publication_guard
from dohm_recovery import copy_pair
from reec_reports import archive_report, report_path


class ReportStorageTest(unittest.TestCase):
    def test_publication_removes_sidecar_and_copy_retains_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            model = Path(temporary) / 'pdbs/OR1/OR1A1.active.pdb'
            model.parent.mkdir(parents=True)
            sidecar = Path(str(model)+'.json')
            with publication_guard(model):
                model.write_text('accepted coordinates')
                sidecar.write_text('{"valid":true,"contacts":[]}')
                self.assertTrue(sidecar.is_file())
            self.assertFalse(sidecar.exists())
            self.assertTrue(json.loads(report_path(model).read_text())['valid'])
            self.assertNotIn('pdbs', report_path(model).relative_to(temporary).parts)
            copied = model.with_name('copy.active.pdb')
            with publication_guard(copied):
                copy_pair(model, copied)
            self.assertFalse(Path(str(copied)+'.json').exists())
            self.assertEqual(report_path(copied).read_bytes(), report_path(model).read_bytes())

    def test_changed_coordinates_cannot_reuse_archived_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            model = Path(temporary) / 'pdbs/OR1/OR1A1.active.pdb'
            model.parent.mkdir(parents=True)
            model.write_text('first coordinates')
            Path(str(model)+'.json').write_text('{"valid":true}')
            archive_report(model)
            model.write_text('different coordinates')
            self.assertFalse(report_path(model).exists())

    def test_failed_publication_restores_coordinates_and_archived_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            model = Path(temporary) / 'pdbs/OR1/OR1A1.active.pdb'
            model.parent.mkdir(parents=True)
            model.write_text('previous coordinates')
            sidecar = Path(str(model)+'.json')
            sidecar.write_text('{"valid":true}')
            archive_report(model)
            with self.assertRaises(ValueError):
                with publication_guard(model):
                    model.write_text('rejected coordinates')
                    sidecar.write_text('{"valid":false}')
                    raise ValueError('Rejected')
            self.assertEqual(model.read_text(), 'previous coordinates')
            self.assertFalse(sidecar.exists())
            self.assertTrue(json.loads(report_path(model).read_text())['valid'])


if __name__ == '__main__':
    unittest.main()
