
import sys
import re
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir('..')
sys.path.append(os.getcwd())

import data.globals
import data.protutils

data.protutils.load_prots()

f = open("data/alichk.dat", "r")
lines = f.read().split("\n")

agroups = \
{
    "aliphatic": "MAILVGP",
    "core_aromatic": "FWY",
    "disruptive": "GP",
    "hydrophilic": "DENQHKRSTY",
    "acids_amides": "DENQ",
    "acidic": "DE",
    "hydrochalcogenide": "STCY",
    "small": "ACGST",
}

def similarity(a, b):
    global agroups
    result = 0
    for aminos in agroups.values():
        if a in aminos and b in aminos: result += 1
    return result

alichk = dict()
label = False
for ln in lines:
    ln = ln.strip()
    if len(ln) < 50 and ln[-1:] == ':': label = ln.split(':')[0]
    elif label:
        alichk[label] = ln
        label = False

threshold = similarity('F', 'L')+1

for protid in data.protutils.prots.keys():
    label = False
    if protid[0:2] == "OR":
        fam = int(re.sub(r'[^0-9]', '', protid[2:4]))
        if fam < 50: label = "ORII"
        else: label = "ORI"
    elif protid[0:4] == "TAAR":
        label = "TAAR"
    else: continue

    prot = data.protutils.prots[protid]
    ali = prot["aligned"]
    chk = alichk[label]
    siti = min(len(ali), len(chk))

    for i in range(siti):
        cz = chk[i:i+1]
        if cz < 'A': continue
        c = ali[i:i+1]
        if c == cz: continue

        c0 = ali[i-1:i]
        c2 = ali[i+1:i+2]
        cz0 = chk[i-1:i]
        cz2 = chk[i+1:i+2]

        if (c0 == cz0): continue
        if (c2 == cz2): continue

        sim0 = similarity(c0, cz0)
        sim  = similarity(c,  cz )
        if sim >= threshold: continue

        siml = similarity(c0, cz )
        mismatch = False
        if siml > sim and siml > sim0: mismatch = True

        if mismatch: print(protid + " possible misalignment:\n" + chk[i-10:i+10] + "\n" + ali[i-10:i+10] + "\n")

