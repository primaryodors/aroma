"""Repair retries continue after generation and final-validation rejection."""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'hm'))
import fixfail
import modeling_backend
from reec_reports import report_path
from reec_backend import CandidateRejected


class FixfailRetryTest(unittest.TestCase):
    def test_native_rejection_is_offered_to_geometry_repair_with_contacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            candidate=Path(temporary)/'candidate.pdb'
            candidate.write_text('REMARK rejected candidate\n')
            contacts=[{'kind':'contact','hard':True}]
            Path(str(candidate)+'.json').write_text(json.dumps({'valid':False,'contacts':contacts}))
            options=SimpleNamespace(seed=4,iterations=12000)
            repaired=Path(temporary)/'repaired.pdb'
            with patch('hm_recovery.repair_existing',return_value=(repaired,Path(temporary))) as repair:
                result=modeling_backend.repair_rejected_candidate('TAAR5',candidate,options)
            self.assertEqual(result,repaired)
            repair.assert_called_once_with('TAAR5',candidate,options,modeling_backend.reec.ROOT,
                                           contact_report=contacts,orient=False)

    def test_geometry_repair_rejection_propagates_to_generation_retry(self):
        with tempfile.TemporaryDirectory() as temporary:
            candidate=Path(temporary)/'candidate.pdb'
            Path(str(candidate)+'.json').write_text('{"contacts":[]}')
            with patch('hm_recovery.repair_existing',side_effect=CandidateRejected('Repair failed')):
                with self.assertRaises(CandidateRejected):
                    modeling_backend.repair_rejected_candidate('TAAR5',candidate,SimpleNamespace())

    def exercise(self, phase, rejected, cap=None):
        previous=Path.cwd();argv=sys.argv[:]
        try:
            with tempfile.TemporaryDirectory() as temporary:
                root=Path(temporary);hm=root/'hm';hm.mkdir()
                destination=root/'pdbs/TAAR/TAAR5.active.pdb'
                destination.parent.mkdir(parents=True)
                prior=(ROOT/'testdata/reec-small.pdb').read_bytes()
                destination.write_bytes(prior)
                sidecar=Path(str(destination)+'.json');sidecar.write_text('{"valid":false}\n')
                calls=[]
                def construct(cls, env, **kwargs):
                    if len(calls)>rejected:raise AssertionError('Unexpected extra attempt')
                    if calls:
                        self.assertEqual(destination.read_bytes(),prior)
                        self.assertFalse(json.loads(sidecar.read_text())['valid'])
                    settings={'seed':env.options.seed,'iterations':env.options.iterations,
                              'clash_scale':env.options.clash_scale}
                    calls.append(settings)
                    model=SimpleNamespace(loop=SimpleNamespace(outputs=[]))
                    def make():
                        self.assertEqual(model.max_var_iterations,settings['iterations'])
                        if phase=='generation' and len(calls)<=rejected:
                            raise CandidateRejected('Simulated native rejection')
                        (hm/'candidate.pdb').write_bytes(prior)
                        model.loop.outputs=[{'failure':None,'molpdf':1,'name':'candidate.pdb'}]
                    model.make=make
                    return model
                def finalize(model, output, candidate, **kwargs):
                    output=Path(output)
                    if phase=='validation' and len(calls)<=rejected:
                        output.write_text('Rejected postprocessing candidate\n')
                        Path(str(output)+'.json').write_text('{"valid":false}\n')
                        raise CandidateRejected('Simulated final steric rejection')
                    output.write_bytes(b'REMARK accepted fresh model\n'+prior)
                    Path(str(output)+'.json').write_text('{"valid":true}\n')
                arguments=['TAAR5',str(destination),'--template',str(destination)]
                if cap is not None:arguments+=['--max-attempts',str(cap)]
                capture=io.StringIO()
                with patch.object(fixfail,'root_dir',str(root)), \
                     patch.object(fixfail,'script_dir',str(hm)), \
                     patch.object(fixfail,'__file__',str(hm/'fixfail.py')), \
                     patch.object(fixfail,'prepare'), \
                     patch.object(fixfail,'repair_reshape_script',return_value=None), \
                     patch.object(fixfail,'construct',side_effect=construct), \
                     patch.object(fixfail,'finalize',side_effect=finalize), \
                     patch.object(fixfail.subprocess,'run'), \
                     patch('modeling_backend.custom_alignment',return_value=(hm/'model.ali',['template='+str(destination)],'')), \
                     patch('data.globals.wait_cool_cpu'), \
                     patch('data.protutils.ensure_bw_metadata'), \
                     contextlib.redirect_stdout(capture),contextlib.redirect_stderr(capture):
                    code=fixfail.main(arguments)
                attempts=json.loads(next((root/'tmp').glob('TAAR5-fixfail-*/attempts.json')).read_text())
                accepted=json.loads(report_path(destination).read_text())['valid']
                return code,calls,attempts,accepted,capture.getvalue()
        finally:
            os.chdir(previous);sys.argv=argv

    def test_plain_repair_retries_native_rejections_until_a_new_model_passes(self):
        code,calls,attempts,accepted,_=self.exercise('generation',2)
        self.assertEqual(code,0);self.assertTrue(accepted)
        self.assertEqual([c['seed'] for c in calls],[0,1,2])
        self.assertEqual([c['iterations'] for c in calls],[3000,6000,12000])
        self.assertEqual([c['clash_scale'] for c in calls],[1,2,4])
        self.assertEqual([r['status'] for r in attempts],['rejected','rejected','accepted'])

    def test_final_validation_rejection_rolls_back_and_retries(self):
        code,calls,attempts,accepted,_=self.exercise('validation',1)
        self.assertEqual(code,0);self.assertTrue(accepted)
        self.assertEqual(len(calls),2)
        self.assertEqual([r['status'] for r in attempts],['rejected','accepted'])

    def test_explicit_cap_reports_failure_without_claiming_old_model_is_valid(self):
        code,calls,attempts,accepted,log=self.exercise('generation',1,cap=1)
        self.assertEqual(code,2);self.assertFalse(accepted)
        self.assertEqual(len(calls),1);self.assertIn('FAILED',log)
        self.assertEqual([r['status'] for r in attempts],['rejected'])


if __name__=='__main__':
    unittest.main()
