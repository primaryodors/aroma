import json
import re
import os
import sys
import subprocess

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
sys.path.append(os.getcwd())

import data.globals
import data.protutils
import data.odorutils

with open('data/binding_pocket.json', 'r') as file:
    bsrdat = json.load(file)

def get_pocket(protid, lig):
    o = data.odorutils.find_odorant(lig)
    if not "sdfname" in o:
        data.odorutils.ensure_sdf_exists(lig)
        o = data.odorutils.find_odorant(lig)
    if not o:
        raise Exception("Unknown ligand " + lig)
    matched = False
    for patt in bsrdat.keys():
        if patt == protid:
            matched = bsrdat[patt]
            break
        elif patt[-1:] == '*':
            pat = patt[0:-1]
            pl = len(pat)
            proti = protid[0:pl]
            if proti == pat:
                matched = bsrdat[patt]
                break
        elif re.search(patt, protid):
            matched = bsrdat[patt]
            break

    if matched:
        if "odorophores" in matched:
            for phore in matched["odorophores"].keys():
                cmd = ["test/moiety_test", o["sdfname"], phore]
                print(" ".join(cmd), "\n\n")
                proc = subprocess.run(cmd, stdout=subprocess.PIPE)
                # print(proc)
                for ln in proc.stdout.decode().split('\n'):
                    m = re.search("[0-9]+ times in molecule", ln)
                    if m:
                        inmol = int(re.sub("[^0-9]", "", m.string))
                        if inmol:
                            for key in matched["odorophores"][phore].keys():
                                matched[key] = matched["odorophores"][phore][key]
                            break

        return matched
    else: return ""
