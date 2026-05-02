
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

# Get the atom numbers and residue numbers of all SG atoms.
# Then scan CONECT records for any SG-SG bonds.
# Save a list of disulfide-linked residue numbers.
sgatoms = dict()
disulfres = []
with open(origpdb, "r") as f:
    c = f.read().__str__()
    c = c.replace("\x0a", "\\n")
    lines = c.split("\\n")

    for ln in lines:
        if ln[0:4] == "ATOM":
            # 01234567890123456789012345678901234567890123456789012345678901234567890123456789
            # ATOM   4141  HD1 HIS A 261      7.697  19.264   2.085   1.00001.00           H 
            aname = ln[12:16].strip()
            if aname == "SG":
                atno = int(ln[5:11].strip())
                resno = int(ln[22:26].strip())
                sgatoms[atno] = resno
        if ln[0:6] == "CONECT":
            # 01234567890123456789
            # CONECT 1518 2823
            conect, atno1, atno2 = ln.split(' ', 2)
            atno1 = int(atno1)
            atno2 = int(atno2)
            if atno1 in sgatoms.keys() and atno2 in sgatoms.keys():
                disulfres.append(sgatoms[atno1])
                disulfres.append(sgatoms[atno2])

cmd = ["atom2omd", "-ipdb", origpdb]
print(" ".join(cmd))
subprocess.run(cmd)

omdfname = f"out/{fam}/{protid}/{protid}~{lignameu}.{mode}.model1.omd"
if not os.path.exists(omdfname):
    print(f"Failed to create omd file.")
    exit()

with open(omdfname, "r") as f:
    c = f.read().__str__()
    c = c.replace("\x0a", "\\n")
    lines = c.split("\\n")

first_atom = last_atom = 0
grpallow = [ # "ASN-OD1", "ASN-ND2", "ASN-HD1", "ASN-HD2",
             # "GLN-OE1", "GLN-NE2", "GLN-HE1", "GLN-HE2",
             # "THR-OG1", "THR-CG2", 
             "THR-HG1", "THR-HG2",
             # "HIS-ND1", "HIS-CD2", "HIS-CE1", "HIS-NE2", 
             "HIS-HD1", "HIS-HD2", "HIS-HE1", "HIS-HE2",
             # "ILE-CG1", "ILE-HG1", "ILE-CG2", "ILE-HG2",
             # "TRP-CD1", "TRP-CD2", "TRP-HD1", "TRP-HD2", 
             "TRP-NE1", "TRP-HE1", "TRP-CE2", "TRP-HE2", "TRP-CE3", "TRP-HE3",
             # "TRP-CZ2", "TRP-HZ2", "TRP-CZ3", "TRP-HZ3", "TRP-CH2", "TRP-HH2", 
            ]
translations = dict()
with open("md/translations", "r") as f:
    c = f.read().__str__()
    c = c.replace("\x0a", "\\n")
    xl8lines = c.split("\\n")
    for ln in xl8lines:
        ln = re.sub("\\s+", " ", ln).strip()
        if ln:
            key, value = ln.split(" ")
            translations[key] = value

xnew = dict()
for xl8 in translations.keys():
    if len(xl8) >= 4 and xl8[3] == '-':
        if not f"N{xl8}" in translations:
            xnew[f"N{xl8}"] = translations[xl8]
        if not f"C{xl8}" in translations:
            xnew[f"C{xl8}"] = translations[xl8]
translations.update(xnew)
# print(translations)

badatno = []
for i in range(len(lines)):
    ln = lines[i]
    if re.search("atom\\[[0-9]+\\]\\s+\\{\\s+type\\s+=\\s+\"[A-Z]{3}-[A-Z0-9]+\";\\s+position\\(\\s*[0-9.-]+,\\s+[0-9.-]+,\\s+[0-9.-]+\\);\\}", ln):
        if not first_atom:
            first_atom = i
        if "HXT" in ln:
            pieces = ln.split("]")
            atno = int(re.sub("[^0-9]", "", pieces[0]))
            badatno.append(atno)
            lines[i] = ""
            continue

    elif first_atom:
        last_atom = i-1


    if "members(" in ln:
        found = False
        for atno in badatno:
            if f"({atno}," in ln or f" {atno})" in ln:
                found = True
                break
        if found:
            lines[i] = ""
            continue

    if "</MetaData>" in ln:
        lines[i] = "\nforceField = \"gaff2\";\n" \
            + "ensemble = NVT;\n" \
            + "cutoffMethod = \"shifted_force\";\n" \
            + "electrostaticScreeningMethod = \"damped\";\n" \
            + "cutoffRadius = 10;\n" \
            + "dampingAlpha = 0.18;\n" \
            + "targetTemp = 310.2;\n" \
            + "tauThermostat = 1000;\n" \
            + "dt = 1.0;\n" \
            + "runTime = 1e4;\n" \
            + "tempSet = \"false\";\n" \
            + "sampleTime = 100;\n" \
            + "statusTime = 10;\n" \
            + ln

    if "Hmat:" in ln:
        ln = "        Hmat: {{ 300, 0, 0 }, { 0, 300, 0 }, { 0, 0, 300 }}"

apply_terminus_prefixes = ["N"] # , "CA", "C", "O", "OXT", "HN", "HA"]
for terminus_resaname in apply_terminus_prefixes:
    if not "XT" in terminus_resaname:
        for i in range(first_atom, last_atom-1):
            ln = lines[i]
            if re.search("[A-Z]{3}-"+terminus_resaname, ln):
                lines[i] = re.sub("([A-Z]{3}-"+terminus_resaname+"\")", "N\\1", ln)
                break
    if False:
        for i in range(last_atom, first_atom-1, -1):
            ln = lines[i]
            if re.search("[A-Z]{3}-"+terminus_resaname, ln):
                lines[i] = re.sub("([A-Z]{3}-"+terminus_resaname+"\")", "C\\1", ln)
                break

resno = 0
for i in range(len(lines)):
    ln = lines[i]

    if re.search("atom\\[[0-9]+\\]\\s+\\{\\s+type\\s+=\\s+\"[A-Z]{3,4}-[A-Z0-9]+\";\\s+position\\(\\s*[0-9.-]+,\\s+[0-9.-]+,\\s+[0-9.-]+\\);\\}", ln):
        m = re.search("[A-Z]{3,4}-[A-Z]+", ln)
        if m:
            grp = m.group()
            aa3let, aname = grp.split('-', 1)
            if aname == "N": resno += 1

            # Set the type of sulfur atoms of cystine cross links to ss.
            if resno in disulfres and aname == "SG":
                ln = lines[i] = ln[0:se[0]] + "ss" + ln[se[1]-1:]

    if re.search("atom\\[[0-9]+\\]\\s+\\{\\s+type\\s+=\\s+\"[A-Z]{3}-[A-Z0-9]+\";\\s+position\\(\\s*[0-9.-]+,\\s+[0-9.-]+,\\s+[0-9.-]+\\);\\}", ln):
        m = re.search("[A-Z]{3}-[A-Z]+[0-9]", ln)
        if m:
            grp = m.group()
            if not grp in grpallow:
                se = m.span()
                ln = lines[i] = ln[0:se[0]] + re.sub("[0-9]", "", grp) + ln[se[1]:]

    if "atom[" in ln:
        for xl8 in translations.keys():
            lines[i] = lines[i].replace(f"\"{xl8}\"", f"\"{translations[xl8]}\"")

with open(omdfname, "w") as f:
    for ln in lines:
        f.write(f"{ln}\n")

omdwarmname = omdfname.replace(".model1.omd", ".model1.warm.omd")
cmd = ["thermalizer", "-i", omdfname, "-o", omdwarmname, "-t", "310.2"]
print(" ".join(cmd))
subprocess.run(cmd)
if not os.path.exists(omdwarmname):
    print("Failed to thermalize.")
    exit()

cmd = ["openmd", omdwarmname]
print(" ".join(cmd))
subprocess.run(cmd)
