
#
# sequence_update.py
#
# Reads the sequences_aligned.txt file and updates all the PDBs with region info and binding site info.
#

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

os.chdir(os.path.dirname(os.path.abspath(__file__)))

c = ""
with open("sequences_aligned.txt", "r") as f:
    c = f.read()

if not c:
    print('NO INPUT FILE')
    exit()

lines = c.split("\n")
ells = dict()
sites = dict()
ballwein = dict()           # https://doi.org/10.1016/S1043-9471(05)80049-7

lln1 = ""
for i in range(len(lines)):
    ln = lines[i].upper().replace("\t", "    ")
    if not len(ln): continue
    if ln[0] == '#': continue
    if ln[0] == '%': continue

    if re.match("^[\\sLSM]+$", ln):
        ells.clear()
        j = -1
        jarr = dict()

        jarr = [(m.start(), ln[m.start()]) for m in re.finditer("[LSM]", ln)]
        for k, v in jarr: ells[k] = v
        # print(ells)

        continue

    if re.search("TMR[0-9][-]+", ln):
        ln += "     "
        rgnames = ["HEAD", "TMR1", "CYT1", "TMR2", "EXR1", "TMR3", "CYT2", "TMR4", "EXR2", "TMR5", "CYT3", "TMR6", "EXR3", "TMR7", "CYT4", "HXR8", "TAIL"]
        rgidx = 0
        k = 0
        rgns = dict()

        while True:
            lk = k
            j = ln.find("TMR", k)
            if j<0: j = ln.find("HXR", k)
            if j<0:
                j = len(ln)+100
                rgns[rgidx] = [lk, j-1, rgnames[rgidx]]
                break

            k = ln.find(' ', j)
            l = ln.find('|', lk)
            if l>=0 and l<j:
                ballwein[rgnames[rgidx]] = l

            l = ln.find('|', j)
            if l>=0:
                ballwein[rgnames[rgidx+1]] = l

            j1 = j-1
            rgns[rgidx] = [lk, j1, rgnames[rgidx]]
            rgidx += 1

            k1 = k-1
            rgns[rgidx] = [j, k1, rgnames[rgidx]]
            rgidx += 1

            if not k: break

        continue

    if re.search("^(OR[0-9]{1,2}[A-Z]{1,2}[0-9]{1,2}\\s+[. AC-IK-NP-TVWY-]+)|(TAAR[0-9]\\s+[. AC-IK-NP-TVWY-]+)|(VN1R[0-9]\\s+[. AC-IK-NP-TVWY-]+)", ln):
        regionse = dict()
        rcpid = ln.strip()[0:11].split(' ')[0]

        r = 0
        pos1 = 0
        posL = 0
        posM = 0
        posS = 0
        emmed = False
        spaced = False
        lastrgn = False

        rcpbw = dict()
        rcpbs = []

        for pos in range(len(ln)):
            char = ln[pos]
            if char == ' ': spaced = True
            if spaced and (char == 'M' or char == '.'): emmed = True
            if emmed:
                if char >= 'A' and char <= 'Y':
                    r += 1
                strgn = "????"
                for rgi in rgns.keys():
                    rg = rgns[rgi]
                    if pos >= rg[0] and pos <= rg[1]:
                        strgn = rg[2]
                        if not rgi in regionse:
                            regionse[rgi] = dict()
                        regionse[rgi][0] = strgn
                        if not 1 in regionse[rgi]:
                            regionse[rgi][1] = r
                        regionse[rgi][2] = r
                        break

                if strgn != lastrgn:
                    posL = 0
                    posM = 0
                    posS = 0

                if (strgn[0:3] == "TMR" or strgn[0:3] == "HXR") \
                    and pos == ballwein[strgn]:
                    tmrno = int(strgn[3:])
                    rcpbw[tmrno] = r
                    if rcpid in data.protutils.prots:
                        data.protutils.prots[rcpid]["bw"][f"{tmrno}.50"] = f"{r}"

                if strgn[0:3] == "EXR" and strgn in ballwein \
                    and pos == ballwein[strgn]:
                    exrno = int(strgn[3:])
                    a = exrno * 2
                    b = a+1
                    ab = int(f"{a}{b}")
                    rcpbw[ab] = r
                    if rcpid in data.protutils.prots:
                        data.protutils.prots[rcpid]["bw"][f"{ab}.50"] = f"{r}"

                if strgn[0:3] == "CYT" and strgn in ballwein \
                    and pos == ballwein[strgn]:
                    cytno = int(strgn[3:])
                    a = cytno * 2 - 1
                    b = a+1
                    ab = int(f"{a}{b}")
                    rcpbw[ab] = r
                    if rcpid in data.protutils.prots:
                        data.protutils.prots[rcpid]["bw"][f"{ab}.50"] = f"{r}"

                if pos in ells and ells[pos] != "":
                    rcpbs.append(r)

        for se in regionse.values():
            region = se[0]
            if region[0:3] == "TMR" or region[0:3] == "HXR":
                if rcpid in data.protutils.prots:
                    data.protutils.prots[rcpid]["region"][se[0]]["start"] = se[1]
                    data.protutils.prots[rcpid]["region"][se[0]]["end"] = se[2]

        # Read the PDB and then rewrite it with updated contents.
        subdir = rcpid[0:4]
        if rcpid[0:2] == "OR":
            lc = subdir[-1:]
            if lc < 'A' or lc > 'Z':
                subdir = subdir[0:3]

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        for suffix in ["upright", "active", "inactive", "apo", "bound"]:
            pdbname = f"../pdbs/{subdir}/{rcpid}.{suffix}.pdb"
            if not os.path.exists(pdbname):
                continue

            with open(pdbname, "r") as f:
                pdbdat = f.read().split("\n")

            remarks = []
            notrems = []

            for pdbln in pdbdat:
                if pdbln[0:6] == "REMARK":
                    remnostr = re.sub("[^0-9]", "", pdbln[7:10]).strip()
                    if not len(remnostr):
                        continue
                    remno = int(remnostr)
                    if remno < 600:
                        remarks.append(pdbln)
                else:
                    notrems.append(pdbln)

            notrems = "\n".join(notrems)
            notrems = notrems.replace("\n\nTER\nEND", "\nTER\nEND")

            remarks.append("REMARK 650")
            for se in regionse.values():
                region = se[0]
                if region[0:3] == "TMR" or region[0:3] == "HXR":
                    remarks.append(f"REMARK 650 HELIX {region} {se[1]} {se[2]}")
            remarks.append("REMARK 650")

            remarks.append("REMARK 800")
            for region in rcpbw:
                bw50 = rcpbw[region]
                remarks.append(f"REMARK 800 SITE BW {region}.50 {bw50}")

            for bsr in rcpbs:
                remarks.append(f"REMARK 800 SITE LIGAND_BINDING {bsr}")
            remarks.append("REMARK 800")

            pdbdat = "\n".join(remarks).strip() + f"\n{notrems}\n\n"

            with open(pdbname, "w") as f:
                f.write(pdbdat)

with open("../data/receptor.json", "w") as f:
    f.write(data.protutils.json_encode_pretty(data.protutils.prots))
