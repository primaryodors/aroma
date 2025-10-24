import sys
import re
import os
import os.path
import json
import subprocess

if len(sys.argv) < 3:
    print("Both a protein and a ligand are required.")
    exit

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.getcwd())
os.chdir("data")

import data.globals
import data.protutils
import data.odorutils
import data.dyncenter

data.protutils.load_prots()
data.odorutils.load_odors()

if not sys.argv[1] in data.protutils.prots.keys():
    protid = False
    popt = sys.argv[1]
else:
    protid = sys.argv[1]
    popt = "*"              # necessary for emp functionality

o = data.odorutils.find_odorant(sys.argv[2])
if not o:
    oid = False
    lopt = sys.argv[2]
else:
    oid = o["oid"]
    lopt = "*"

for rcpid in data.protutils.prots.keys():
    fam = data.protutils.family_from_protid(rcpid)
    if protid:
        if rcpid != protid: continue
    else:
        if popt != "*" and popt != "emp":
            if re.match("e[0-9]+", popt):
                if not "expression" in data.protutils.prots[rcpid]: continue
                minexpr = int(popt[1:])
                if int(data.protutils.prots[rcpid]["expression"]) < minexpr: continue
            elif not re.match(popt, rcpid): continue

    for ligid in data.odorutils.odors.keys():
        if oid:
            if ligid != oid: continue
        else:
            o = data.odorutils.odors[ligid]

        if popt == "emp" or lopt == "emp":
            if not "activity" in o: continue
            isemp = False
            for url in o["activity"].keys():
                acv = o["activity"][url]
                if rcpid in acv: isemp = True
            if not isemp: continue
        elif not oid:
            isnote = False
            if "aroma" in o:
                for url in o["aroma"]:
                    aromata = o["aroma"][url]
                    if lopt in aromata: isnote = True
            if not isnote:
                if not re.search(lopt, o["smiles"]) and not re.search(lopt, o["full_name"]):
                    # TODO: Moieties
                    continue

        lignu = o["full_name"].replace(' ', '_')
        isomers = data.odorutils.check_isomers(o["full_name"])
        data.odorutils.ensure_sdf_exists(o["full_name"])

        pocket = data.dyncenter.get_pocket(rcpid, o["full_name"])
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
                    ln = "PROT pdbs/" + fam + "/" + rcpid + ".active.pdb"
                elif ln[0:4] == "LIG ":
                    ln = "LIG sdf/" + lignu + ".sdf"
                    if isomers and len(isomers):
                        for iso in isomers:
                            ln += "\nISO sdf/" + iso.replace(' ', '_') + ".sdf"
                elif ln[0:4] == "CEN ":
                    if pocket:
                        if isinstance(pocket["pocket"], str):
                            ln = "CEN RES " + pocket["pocket"]
                        else:
                            ln = ""
                            for pkt in pocket["pocket"]:
                                ln = ln + "CEN RES " + pkt + "\n"
                elif ln[0:4] == "OUT ":
                    ln = ("OUT output/" + fam + "/" + rcpid + "/" + rcpid + "~" + lignu + ".active.dock" 
                        + "\nOUTPDB 1 output/" + fam + "/" + rcpid + "/" + rcpid + "~" + lignu + ".active.model%"+"o.pdb")

                newcfg.append(ln)

        if fam[0:2] == "OR":
            if int(fam[2:]) < 50: softness = "1.0"
            else: softness = "0.1"
            newcfg.append("SOFT " + softness + " 4 5 6 7")
        newcfg.append("NODEL 45.52 5.39")
        newcfg.append("NODEL 7.49 7.55")

        if pocket:
            if "atomto" in pocket:
                if isinstance(pocket["atomto"], str):
                    pocket["atomto"] = [pocket["atomto"]]
                for a2 in pocket["atomto"]:
                    newcfg.append("ATOMTO " + a2)
            if "flxr" in pocket:
                if isinstance(pocket["flxr"], str):
                    pocket["flxr"] = [pocket["flxr"]]
                for fx in pocket["flxr"]:
                    newcfg.append("FLXR " + fx)
            if "stcr" in pocket:
                if isinstance(pocket["stcr"], str):
                    pocket["stcr"] = [pocket["stcr"]]
                for st in pocket["stcr"]:
                    newcfg.append("STCR " + st)

        if not os.path.exists("output"): os.mkdir("output")
        if not os.path.exists("output/" + fam): os.mkdir("output/" + fam)
        if not os.path.exists("output/" + fam + "/" + rcpid): os.mkdir("output/" + fam + "/" + rcpid)

        outfna = rcpid + "~" + lignu + ".active.config"
        outfni = rcpid + "~" + lignu + ".inactive.config"
        with open("tmp/" + outfna, 'w') as f:
            f.write("\n".join(newcfg) + "\n\n")
        for i in range(len(newcfg)):
            ln = newcfg[i]
            if ln[0:4] == "PROT" or ln[0:3] == "OUT":
                newcfg[i] = ln.replace(".active.", ".inactive.")
        with open("tmp/" + outfni, 'w') as f:
            f.write("\n".join(newcfg) + "\n\n")

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        cmd = ["bin/aromadock", "tmp/" + outfna]
        print(" ".join(cmd), "\n\n")
        subprocess.run(cmd)
        cmd = ["bin/aromadock", "tmp/" + outfni]
        print(" ".join(cmd), "\n\n")
        subprocess.run(cmd)

        if os.path.exists("tmp/nodelete"):
            print("Warning: not deleting temporary config file because you have the debug \"nodelete\" option selected.")
        else:
            os.remove("tmp/" + outfna)
            os.remove("tmp/" + outfni)

