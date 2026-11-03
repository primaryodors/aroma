import sys
import os
import os.path
import json
import subprocess
from pathlib import Path

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
sys.path.append(os.getcwd())

import data.globals
import data.protutils
import data.odorutils
from data.preparation import prepared_output_current

data.protutils.load_prots()
data.odorutils.load_odors()

with open('data/cavopts.json', 'r') as file:
    cavopts = json.load(file)

for protid in data.protutils.prots.keys():
    for suffix in ["active", "inactive"]:
        pdbfn = data.protutils.build_pdb_fname(protid, suffix)
        cavfn = pdbfn[0:-4]+".cvty"

        inputs = [pdbfn, "bin/cavity_search", "data/cavopts.json",
                  *sorted(Path("data").glob("*.dat"))]
        if prepared_output_current(cavfn, inputs):
            continue

        cmd = ["bin/cavity_search", "-p", pdbfn, "-o", cavfn]
        cmd.extend(cavopts)
        print(" ".join(cmd))
        subprocess.run(cmd)

