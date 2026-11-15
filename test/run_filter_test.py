import ast
import os
import io
import re
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
TREE = ast.parse((ROOT / "run.py").read_text())


class RunFilterTest(unittest.TestCase):
    def test_cli_option(self):
        parser = next(node for node in TREE.body
                      if isinstance(node, ast.If)
                      and ast.unparse(node.test) == "len(sys.argv) > 3")
        context = {"sys": SimpleNamespace(argv=["run.py", "all", "all", "onlyfailed"]),
                   "onlynew": False, "onlymissing": False, "onlyfailed": False}
        exec(compile(ast.Module(body=[parser], type_ignores=[]), "run.py", "exec"), context)
        self.assertTrue(context["onlyfailed"])
        self.assertFalse(context["onlynew"])
        self.assertFalse(context["onlymissing"])

    def test_onlyfailed_selection(self):
        selection = next(node for node in ast.walk(TREE)
                         if isinstance(node, ast.If)
                         and ast.unparse(node.test) == "onlynew or onlymissing or onlyfailed")
        wrapper = ast.parse("for suff in states:\n    pass\n    selected.append(suff)\n")
        wrapper.body[0].body[0] = selection
        code = compile(ast.fix_missing_locations(wrapper), "run.py", "exec")
        sizes = {"empty": 0, "failed": 500, "boundary": 100000, "success": 100001}
        reports = {
            "fallback": "Pose: 1\nDisqualified because: rejected\n0 pose(s) found.\n",
            "initial": "Pose: 1\nNOFAIL: output from initial placement\n0 pose(s) found.\n",
            "marker": "Pose: 1\nNOFAIL: output from initial placement\n",
            "accepted": "Pose: 1\n1 pose(s) found.\n",
            "mixed": "Disqualified because: rejected\nPose: 1\n2 pose(s) found.\n",
            "smallaccepted": "1 pose(s) found.\n",
        }
        sizes.update({name: 150000 for name in reports})
        sizes["smallaccepted"] = 50
        for refresh, resume in ((False, False), (True, False), (False, True), (True, True)):
            with self.subTest(refresh=refresh, resume=resume):
                context = dict(os=os, re=re, fam="OR1", rcpid="OR1A1", lignu="test",
                               states=["missing", *sizes], selected=[], onlyfailed=True,
                               onlynew=refresh, onlymissing=resume)
                def state(path):
                    return path.rsplit(".", 2)[-2]
                with patch("os.path.exists", side_effect=lambda path: state(path) in sizes), \
                     patch("os.path.getsize", side_effect=lambda path: sizes[state(path)]), \
                     patch("builtins.open", side_effect=lambda path, **kwargs:
                           io.StringIO(reports.get(state(path), ""))):
                    exec(code, context)
                self.assertEqual(context["selected"],
                                 ["empty", "failed", "boundary", "fallback", "initial", "marker"])


if __name__ == "__main__":
    unittest.main()
