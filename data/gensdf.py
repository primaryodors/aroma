
import sys
import re
import os
import hashlib
import subprocess

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir('..')
sys.path.append(os.getcwd())

import data.globals
import data.protutils
import data.odorutils

data.odorutils.load_odors()

def echo_usage():
    print("Usage:\n\npython3 data/gensdf.py molecule_name SMILES\n\n")
    exit

name = sys.argv[1]
if not name: echo_usage()

smiles = sys.argv[2]
if not smiles: echo_usage()

canonical = subprocess.run(["obabel", "-:"+smiles, "-ocan"], capture_output=True, text=True).stdout.strip()
print("Canonical SMILES is "+canonical+"\n")

hash = hashlib.md5(canonical.encode()).hexdigest()
odors = data.odorutils.odors
if hash in odors.keys():
    print("Odor already exists.\n")
    exit()

odors[hash] = dict()
odors[hash]['full_name'] = name.replace("_", " ")
odors[hash]['smiles'] = canonical
odors[hash]['oid'] = hash

sortk = sorted(odors.keys())
odors1 = dict()
for k in sortk:
    odors1[k] = odors[k]

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir('..')
with open("data/odorant.json", "wb") as f:
    f.write(data.protutils.json_encode_pretty(odors1).encode())

nameu = name.replace(" ", "_")
output_file = f"sdf/{nameu}.sdf"
data.odorutils.smiles_to_sdf(canonical, output_file)

