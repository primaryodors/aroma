
import sys
import re
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir('..')
sys.path.append(os.getcwd())

import data.globals
import data.protutils
import data.odorutils

data.protutils.load_prots()
data.odorutils.load_odors()

def sgn(num):
    if not num: return 0
    elif num < 0: return -1
    else: return 1

colors = \
    {
        "PAILV": "\033[38;5;247m",
        "G": "\033[38;5;243m",
        "M": "\033[38;5;185m",
        "C": "\033[38;5;220m",
        "STNQ": "\033[38;5;49m",
        "DE": "\033[38;5;196m",
        "KR": "\033[38;5;27m",
        "H": "\033[38;5;99m",
        "FWY": "\033[38;5;171m"
    }

orphan = "\033[2m\033[3m"
deorphan = ""; # "\033[1m"
reset = "\033[22m\033[23m\033[24m"

print("Legend: " + deorphan + "Deorphaned receptor" + reset + " " + orphan + "Orphan receptor" + reset + "\n")

bw = False
list = False
for a in sys.argv[1:]:
    a = a.split('=',2)
    if a[0] == 'bw': bw = a[1]
    elif a[0] == 'list': list = a[1]

if not bw:
    print("Error: no BW number.\n")
    exit

pieces = bw.split('.')
rgno = int(pieces[0])
ofst = int(pieces[1])

with open('data/sequences_aligned.txt') as f: c = f.read()
lines = c.split('\n')

var = {}
var["cls2"] = {}
var["cls1"] = {}
var["taar"] = {}
var["vn1r"] = {}
var_has = {}
var_has["cls2"] = {}
var_has["cls1"] = {}
var_has["taar"] = {}
var_has["vn1r"] = {}

col = 0
for ln in lines:
    if re.search("TMR[0-9]-*", ln):
        lookfor = "TMR" + str(rgno)[0:1]
        j = ln.find(lookfor)
        if j<0:
            lookfor = "HXR" + rgno[0:1]
            j = ln.find(lookfor)
        if j<15: continue

        if rgno < 10:
            col50 = ln.find('|', j)
        else:
            col50 = ln.find(' |', j) + 1
        col = False
    elif ln[0:2] == '% ':
        rel = ofst - 50
        j = 0
        i = 0
        while 1:
            i1 = col50 + i
            if i1 < 15 or i1 > len(ln): break
            c = ln[i1:i1+1]
            if c != ' ':
                if j == rel:
                    col = col50 + i
                    break
                j += sgn(rel)
            i += sgn(rel)
            if not i: break
    else:
        if col < 15: continue
        orid = ln[0:15].strip()
        if not orid: continue
        orid = orid.split(' ')[0]

        c = ln[col:col+1]
        if c == ' ': c = '-'

        vark = False
        if orid[0:2] == "OR":
            fam = int("".join(re.findall(r'\d', orid[2:4])))
            if fam < 50: vark = 'cls2'
            else: vark = 'cls1'
        else: vark = orid[0:4].lower()

        if not vark: continue
        if not vark in var: continue

        if not var[vark].get(c): var[vark][c] = 1
        else: var[vark][c] += 1

        ligands = len(data.odorutils.empirical_pairs(orid, True))
        if not c in var_has[vark]: var_has[vark][c] = ""
        if ligands: var_has[vark][c] += deorphan + orid + reset + " "
        else: var_has[vark][c] += orphan + orid + reset + " "

for vark in var.keys():
    vard = var[vark]
    if not len(vard): continue
    varh = var_has[vark]
    if vark == 'cls2': print("Class II")
    elif vark == 'cls1': print("Class I")
    elif vark == 'taar': print("TAAR")
    elif vark == 'vn1r': print("VN1R")

    vs = sorted(vard.items(), key=lambda x:x[1], reverse=True)

    ttl = 0
    for k in vard.keys():
        ttl += vard[k]
    thr = int(0.8 * ttl)
    vsum = 0
    fnd80 = 1
    for kv in vs:
        k = kv[0]
        vsum += vard[k]
        if vsum >= thr: break
        else: fnd80 += 1

    i = 0
    for cv in vs:
        i += 1
        c = cv[0]
        for líuon in colors.keys():
            if líuon.find(c) >= 0: print(colors[líuon], end="")
        print(c + " ", end="")
        pcnt = round(float(vard[c]) / ttl * 100, 3)
        print(pcnt, end="")
        print("%", end="")

        if i > fnd80 or c == list:
            if c in varh: print(" " + varh[c], end="")

        print(reset)

    print("\n")
