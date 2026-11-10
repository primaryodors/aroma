from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from data import protutils as pu


class BwNumberingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pu.load_prots()

    def test_or52e8_deletion_and_residue_filter(self):
        self.assertEqual(pu.bw_insdel('OR52E8', 5, 42), 1)
        self.assertEqual(pu.bw_insdel('OR52E8', '5', '42'), 1)
        self.assertEqual(pu.resno_from_bw('OR52E8', 'N5.42'), 208)
        self.assertEqual(pu.resno_from_bw('OR52E8', '5.43'), 209)
        self.assertEqual(pu.resno_from_bw('OR52E8', 'STNQ5.43'), 0)
        self.assertEqual(pu.resno_from_bw('OR52E8', '5.49'), 0)
        self.assertEqual(pu.resno_from_bw('OR52E8', '5.50'), 215)

    def test_repeated_insertions_at_a_deleted_position(self):
        receptor = {'sequence': 'A'*100, 'bw': {'7.50': '50'},
                    'deletion': ['7.36'], 'insertion': ['7.36', '7.36']}
        with patch.dict(pu.prots, {'numbering_fixture': receptor}):
            self.assertEqual(pu.resno_from_bw('numbering_fixture', '7.35'), 34)
            self.assertEqual(pu.resno_from_bw('numbering_fixture', '7.36'), 0)
            self.assertEqual(pu.resno_from_bw('numbering_fixture', '7.36+1'), 35)
            self.assertEqual(pu.resno_from_bw('numbering_fixture', '7.36+2'), 36)
            self.assertEqual(pu.resno_from_bw('numbering_fixture', '7.36+3'), 0)
            self.assertEqual(pu.resno_from_bw('numbering_fixture', '7.37'), 37)
            self.assertIn('REMARK 800 SITE BW_INSERTION 7.36 2', pu.bw_indel_metadata('numbering_fixture'))

    def test_metadata_refresh_is_idempotent_and_preserves_model_fingerprint(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'OR52E8.active.pdb'
            body = (ROOT / 'testdata/reec-small.pdb').read_bytes()
            path.write_bytes(body)
            status = {'valid': True, 'validated_ready': True,
                      'validated_sha256': hashlib.sha256(body).hexdigest()}
            Path(str(path)+'.json').write_text(json.dumps(status))
            self.assertTrue(pu.ensure_bw_metadata('OR52E8', path))
            self.assertFalse(pu.ensure_bw_metadata('OR52E8', path))
            after = path.read_bytes()
            self.assertEqual([l for l in after.splitlines() if l.startswith(b'ATOM  ')],
                             [l for l in body.splitlines() if l.startswith(b'ATOM  ')])
            self.assertEqual(json.loads(Path(str(path)+'.json').read_text())['validated_sha256'],
                             hashlib.sha256(after).hexdigest())

    def test_python_and_native_agree_across_the_receptor_database(self):
        executable = ROOT / 'test/bw_indel_test'
        self.assertTrue(executable.is_file(), 'Build test/bw_indel_test before running this suite')
        commands, expected = [], []
        for name, receptor in pu.prots.items():
            if not receptor.get('sequence') or not receptor.get('bw'):
                continue
            commands.append('R 0 0 0')
            for tag, anchor in receptor['bw'].items():
                commands.append(f'A {int(tag.split(".")[0])} {int(anchor)} 0')
            for field, code in [('deletion', 'D'), ('insertion', 'I')]:
                for tag, count in Counter(receptor.get(field, [])).items():
                    helix, position = map(int, tag.split('.'))
                    commands.append(f'{code} {helix} {position} {count}')
            for anchor_tag in receptor['bw']:
                helix = int(anchor_tag.split('.')[0])
                for position in range(20, 81):
                    tag = f'{helix}.{position}'
                    result = pu.resno_from_bw(name, tag)
                    deleted = tag in receptor.get('deletion', [])
                    if result or deleted:
                        commands.append(f'Q {helix} {position} 0')
                        expected.append(result)
                for tag, count in Counter(receptor.get('insertion', [])).items():
                    event_helix, position = map(int, tag.split('.'))
                    if event_helix != helix:
                        continue
                    for ordinal in range(1, count+1):
                        result = pu.resno_from_bw(name, f'{tag}+{ordinal}')
                        if result:
                            commands.append(f'Q {helix} {position} {ordinal}')
                            expected.append(result)
        result = subprocess.run([str(executable), '--probe'], input='\n'.join(commands)+'\n',
                                text=True, capture_output=True, timeout=60, cwd=ROOT)
        self.assertEqual(result.returncode, 0, result.stderr)
        actual = [int(value) for value in result.stdout.splitlines()]
        self.assertGreater(len(expected), 100000)
        self.assertEqual(actual, expected)

    def test_phew_uses_corrected_mapping_after_copy_and_reload(self):
        source = ROOT / 'pdbs/OR52/OR52E8.active.pdb'
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            annotated, legacy, saved = root / 'annotated.pdb', root / 'legacy.pdb', root / 'saved.pdb'
            annotated.write_bytes(source.read_bytes())
            pu.ensure_bw_metadata('OR52E8', annotated)
            legacy.write_text('\n'.join(line for line in source.read_text().splitlines()
                                        if not line.startswith(('REMARK 800 SITE BW_DELETION ',
                                                                'REMARK 800 SITE BW_INSERTION ')))+'\n')
            script = root / 'probe.phew'
            script.write_text(f'LOAD "{annotated}" A I\nLOAD "{legacy}" A A\nSTRAND A\n'
                              'BWCOPY I A\nECHO "ASN=" %5.42\nECHO "ILE=" %5.43\n'
                              f'ECHO "GAP=" %5.49\nSAVE "{saved}"\nLOAD "{saved}" A A\n'
                              'STRAND A\nECHO "RESTORED=" %5.42\n')
            result = subprocess.run([str(ROOT / 'bin/phew'), str(script)], cwd=ROOT,
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            for label, number in [('ASN', 208), ('ILE', 209), ('GAP', 0), ('RESTORED', 208)]:
                self.assertIn(f'{label}={number}', result.stdout)


if __name__ == '__main__':
    unittest.main()
