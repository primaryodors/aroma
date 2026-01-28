
import sys
import re
import os
import os.path
import json
import subprocess

lsa = len(sys.argv)
if lsa < 3:
    print("Two residue Ballesteros-Weinstein numbers are required.")
    exit
elif lsa > 3:
    regex = sys.argv[3]
else:
    regex = False

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
sys.path.append(os.getcwd())

import data.globals
import data.protutils
import data.odorutils

data.protutils.load_prots()
data.odorutils.load_odors()

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")

for protid in data.protutils.prots:
    fam = data.protutils.family_from_protid(protid)
    if regex:
        if not re.match(regex, protid): continue
    elif fam[0:2] != "OR": continue
    resno1 = data.protutils.resno_from_bw(protid, sys.argv[1])
    resno2 = data.protutils.resno_from_bw(protid, sys.argv[2])
    distance = float(subprocess.run(["bin/phew", "data/measure.phew",
        f"pdbs/{fam}/{protid}.active.pdb",
        f"{resno1}", f"{resno2}"],
        capture_output=True, text=True).stdout.strip())
    if "best_agonist" in data.protutils.prots[protid]:
        oid = data.protutils.prots[protid]["best_agonist"]
        best = data.odorutils.odors[oid]["full_name"]
    else: best = ""
    print(f"{protid}: {resno1} {resno2} {distance} {best}")
