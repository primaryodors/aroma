import sys
import re
import os
import os.path
import json
import subprocess

if len(sys.argv) < 3:
    print("inti ambostin cicos etic bolatos sent ancenas scribbasesetuio")
    exit

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.getcwd())
os.chdir("odor")

import odor.globals
import odor.protutils
import odor.odorutils
import odor.dyncenter

odor.protutils.load_prots()
odor.odorutils.load_odors()

if not sys.argv[1] in odor.protutils.prots.keys():
    print("cicos angnatos: "+sys.argv[1])
    exit
else:
    protid = sys.argv[1]
    fam = odor.protutils.family_from_protid(protid)

o = odor.odorutils.find_odorant(sys.argv[2])
if not o:
    print("bolatos angnatos: "+sys.argv[1])
    exit
lignu = o["full_name"].replace(' ', '_')

pocket = odor.dyncenter.get_pocket(protid, o["full_name"])
# print(pocket)

os.chdir(os.path.dirname(os.path.abspath(__file__)))
with open("example.config") as f:
    cfg = f.read()
    lines = cfg.split("\n")
    newcfg = ["# This file was automatically generated, therefore any changes you make will likely be overwritten."]
    for ln in lines:
        ln = ln.strip()
        if not ln: continue
        if ln[0:1] == '#': continue
        if ln[0:5] == "PROT ":
            ln = "PROT pdbs/" + fam + "/" + protid + ".active.pdb"
        elif ln[0:4] == "LIG ":
            ln = "LIG sdf/" + lignu + ".sdf"
        elif ln[0:4] == "CEN ":
            if pocket:
                ln = "CEN RES " + pocket["pocket"]
        elif ln[0:4] == "OUT ":
            ln = ("OUT output/" + fam + "/" + protid + "/" + protid + "~" + lignu + ".active.dock" 
                + "\nOUTPDB 1 " + fam + "/" + protid + "/" + protid + "~" + lignu + ".active.model%"+"o.pdb")

        newcfg.append(ln)

if fam[0:2] == "OR":
    if int(fam[2:]) < 50: softness = "1.0"
    else: softness = "0.1"
    newcfg.append("SOFT " + softness + " 4 5 6 7")
newcfg.append("NODEL 45.52 5.39")
newcfg.append("NODEL 7.49 7.55")

if pocket:
    if "atomto" in pocket:
        for a2 in pocket["atomto"]:
            newcfg.append("ATOMTO " + a2)
    if "flxr" in pocket:
        for fx in pocket["flxr"]:
            newcfg.append("FLXR " + fx)
    if "stcr" in pocket:
        for st in pocket["stcr"]:
            newcfg.append("STCR " + st)

if not os.path.exists("output"): os.mkdir("output")
if not os.path.exists("output/" + fam): os.mkdir("output/" + fam)
if not os.path.exists("output/" + fam + "/" + protid): os.mkdir("output/" + fam + "/" + protid)

outfna = protid + "~" + lignu + ".active.config"
outfni = protid + "~" + lignu + ".inactive.config"
with open("tmp/" + outfna, 'w') as f:
    f.write("\n".join(newcfg) + "\n\n")
for i in newcfg.keys():
    ln = newcfg[i]
    if ln[0:4] == "PROT" or ln[0:3] == "OUT":
        newcfg[i] = ln.replace(".active.", ".inactive.")
with open("tmp/" + outfni, 'w') as f:
    f.write("\n".join(newcfg) + "\n\n")

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir('../reec')
cmd = ["bin/reec", "tmp/" + outfna]
print(" ".join(cmd), "\n\n")
subprocess.run(cmd)
cmd = ["bin/reec", "tmp/" + outfni]
print(" ".join(cmd), "\n\n")
subprocess.run(cmd)

if os.path.exists("tmp/nodelete"):
    print("Warning: not deleting temporary config file because you have the debug \"nodelete\" option selected.")
else:
    os.remove("tmp/" + outfna)
    os.remove("tmp/" + outfni)

