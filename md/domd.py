
import sys
import re
import os
import os.path
import json
import shutil
import subprocess

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
sys.path.append(os.getcwd())
os.chdir("data")

import data.globals
import data.protutils
import data.odorutils
import data.dyncenter

data.protutils.load_prots()
data.odorutils.load_odors()

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")

argc = len(sys.argv)
if argc < 2:
    print("Both a protein ID and a ligand name are required.")
    exit()

mode="active"
protid = sys.argv[1]
if protid[-1] == 'a':
    mode="active"
    protid = protid[0:-1]
elif protid[-1] == 'i':
    mode="inactive"
    protid = protid[0:-1]
if not protid in data.protutils.prots.keys():
    print(f"Protein {protid} not found.")
    exit()
fam = data.protutils.family_from_protid(protid)

odor = data.odorutils.find_odorant(sys.argv[2])
if not odor:
    print(f"Odor {sys.argv[2]} not found.")
    exit()

ligname = odor["full_name"]
lignameu = ligname.replace(' ', '_')
origpdb = f"out/{fam}/{protid}/{protid}~{lignameu}.{mode}.model1.pdb"

if not os.path.exists(origpdb):
    print(f"File not found {origpdb}. Make sure you have already docked this protein + ligand pair", \
          "and that it generated a PDB file in out/{fam}/{protid} ending in .model1.pdb")
    exit()

cmd = ["atom2omd", "-ipdb", origpdb]
subprocess.run(cmd)

omdfname = f"out/{fam}/{protid}/{protid}~{ligname}.{mode}.model1.omd"
if not os.path.exists(omdfname):
    print(f"Failed to create omd file.")
    exit()