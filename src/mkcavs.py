import sys
import os
import os.path
import json
import subprocess

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
sys.path.append(os.getcwd())

import data.globals
import data.protutils
import data.odorutils

data.protutils.load_prots()
data.odorutils.load_odors()

with open('data/cavopts.json', 'r') as file:
    cavopts = json.load(file)

bincavfd = os.path.getmtime("bin/cavity_search")

for protid in data.protutils.prots.keys():
    for suffix in ["active", "inactive"]:
        pdbfn = data.protutils.build_pdb_fname(protid, suffix)
        cavfn = pdbfn[0:-4]+".cvty"

        if os.path.exists(cavfn):
            cavfd = os.path.getmtime(cavfn)
            pdbfd = os.path.getmtime(pdbfn)
            if cavfd > pdbfd and cavfd > bincavfd:
                continue

        cmd = ["bin/cavity_search", "-p", pdbfn, "-o", cavfn]
        cmd.extend(cavopts)
        print(" ".join(cmd))
        subprocess.run(cmd)

