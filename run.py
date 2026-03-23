import sys
import re
import os
import os.path
import json
import subprocess
import pprint

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

onlynew = False
onlymissing = False
skipdock = False

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

if popt == "all": popt = "*"
if lopt == "all": lopt = "*"

g = ["L", "S2", "I2", "O1"]

resseek = False
if re.match("([A-Y]+[0-9]{1,2}[.][0-9]{2},?)+", popt):
    resseek = popt.split(",")
    popt = "res"

if (popt != "*" and popt != "emp" and popt != "ago" and popt != "top") and oid:
    if "ago" in sys.argv: lopt = "ago"
    elif "emp" in sys.argv: lopt = "emp"
    elif "top" in sys.argv: lopt = "top"

if len(sys.argv) > 3:
    for i in range(3, len(sys.argv)):
        if sys.argv[i] == "refresh": onlynew = True
        if sys.argv[i] == "resume": onlymissing = True

for rcpid in data.protutils.prots.keys():
    fam = data.protutils.family_from_protid(rcpid)
    sbf = data.protutils.subfamily_from_protid(rcpid)

    if protid:
        if rcpid != protid: continue
    else:
        if popt == "res":
            if not "bw" in data.protutils.prots[rcpid]: continue
            match = True
            for aabw in resseek:
                rn = data.protutils.resno_from_bw(rcpid, aabw)
                if not rn: match = False
            if not match: continue
        elif popt != "*" and popt != "emp" and popt != "ago" and popt != "top":
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

        if lopt != "*" or popt != "*":
            if popt == "emp" or lopt == "emp":
                if not "activity" in o: continue
                isemp = False
                for url in o["activity"].keys():
                    acv = o["activity"][url]
                    if rcpid in acv: isemp = True
                if not isemp: continue
            elif popt == "ago" or lopt == "ago":
                if not "activity" in o: continue
                isago = False
                for url in o["activity"].keys():
                    acv = o["activity"][url]
                    if rcpid in acv:
                        if "adjusted_curve_top" in acv[rcpid]:
                            if float(acv[rcpid]["adjusted_curve_top"]) > 0:
                                isago = True
                        elif "type" in acv[rcpid]:
                            if acv[rcpid]["type"] in ["vsa", "sa", "ma", "wa", "pa", "a"]:
                                isago = True
                        elif "ec50" in acv[rcpid]:
                            isago = True
                if not isago: continue
            elif lopt == "top":
                p = data.protutils.prots[rcpid]
                if not "best_agonist" in p: continue
                if ligid != p["best_agonist"]: continue
            elif lopt == "jk":
                skipdock = True
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
        runsuff = []
        if os.path.exists("cpllocal"):
            for suff in g:
                cpllocal = f"cpllocal/{fam}/{sbf}/{rcpid}.{suff}.pdb"
                if not os.path.exists(cpllocal):
                    if not os.path.exists(f"cpllocal/{fam}"): os.mkdir(f"cpllocal/{fam}")
                    if not os.path.exists(f"cpllocal/{fam}/{sbf}"): os.mkdir(f"cpllocal/{fam}/{sbf}")
                    cpl = f"coupled/{fam}/{sbf}/{rcpid}~hGNA{suff}.pdb"
                    if os.path.exists(cpl): data.protutils.prepare_coupled(cpl, cpllocal, rcpid)
                if os.path.exists(cpllocal):
                    runsuff.append(suff)
        if len(runsuff):
            pdbdir = f"cpllocal/{fam}/{sbf}"
        else:
            runsuff = [ "active", "inactive" ]
            pdbdir = f"pdbs/{fam}"

        lignu = o["full_name"].replace(' ', '_')
        isomers = data.odorutils.check_isomers(o["full_name"])
        forms = data.odorutils.check_forms(o["full_name"])
        data.odorutils.ensure_sdf_exists(o["full_name"])

        pocket = data.dyncenter.get_pocket(rcpid, o["full_name"])
        # pprint.pprint(pocket)

        for suff in runsuff:
            conffn = f"{rcpid}~{lignu}.{suff}.config"
            if not os.path.exists("out"): os.mkdir("out")
            if not os.path.exists("out/" + fam): os.mkdir("out/" + fam)
            if not os.path.exists("out/" + fam + "/" + rcpid): os.mkdir("out/" + fam + "/" + rcpid)

            if onlynew or onlymissing:
                outfn = f"out/{fam}/{rcpid}/{rcpid}~{lignu}.{suff}.dock"
                if os.path.exists(outfn) and os.path.getsize(outfn) > 1e5:
                    if onlymissing: continue
                    fmt = os.path.getmtime(outfn)
                    if fmt > os.path.getmtime("data/binding_pocket.json") \
                    and fmt > os.path.getmtime(f"sdf/{lignu}.sdf") \
                    and fmt > os.path.getmtime(f"{pdbdir}/{rcpid}.{suff}.pdb") \
                    and fmt > os.path.getmtime("bin/aromadock"):
                        continue

            print(f"Beginning {rcpid} ~ "+o["full_name"]+"...")
            os.chdir(os.path.dirname(os.path.abspath(__file__)))
            with open("example.config", "r") as f:
                cfg = f.read()
                lines = cfg.split("\n")
                newcfg = ["# This file was automatically generated, therefore any changes you make will likely be overwritten."]
                for ln in lines:
                    ln = ln.strip()
                    if not ln: continue
                    if ln[0:1] == '#': continue
                    if ln[0:5] == "PROT ":
                        ln = f"PROT {pdbdir}/{rcpid}.{suff}.pdb"
                    elif ln[0:4] == "LIG ":
                        ln = "LIG sdf/" + lignu + ".sdf"
                        if isomers and len(isomers):
                            for iso in isomers:
                                ln += "\nISO sdf/" + iso.replace(' ', '_') + ".sdf"
                        if forms and len(forms):
                            for form in forms:
                                ln += "\nFORM sdf/" + form.replace(' ', '_') + ".sdf"
                    elif ln[0:4] == "CEN ":
                        if pocket:
                            if isinstance(pocket["pocket"], str):
                                ln = "CEN RES " + pocket["pocket"]
                            else:
                                ln = ""
                                for pkt in pocket["pocket"]:
                                    ln = ln + "CEN RES " + pkt + "\n"
                    elif ln[0:4] == "OUT ":
                        ln = (f"OUT out/{fam}/{rcpid}/{rcpid}~{lignu}.{suff}.dock" 
                            + f"\nOUTPDB 1 out/{fam}/{rcpid}/{rcpid}~{lignu}.{suff}.model%"+"o.pdb")

                    newcfg.append(ln)

            if fam[0:2] == "OR":
                if int(fam[2:]) < 50: softness = "1.0"
                else: softness = "0.1"
                newcfg.append("SOFT " + softness + " 2 3 4 5 6 7")
            newcfg.append("NODEL 45.52 5.39")
            newcfg.append("NODEL 7.49 7.55")
            # newcfg.append("OUTBBP")
            # newcfg.append("OUTDISQ")
            # newcfg.append("OUTMC 1")
            newcfg.append("NORESWARN")
            newcfg.append("NOFAIL")
            # newcfg.append("MOVIE")

            if not skipdock:
                cmd = ["bin/ic", f"{pdbdir}/{rcpid}.{suff}.pdb", "-3.0", "nooil"]
                print(" ".join(cmd))
                proc = subprocess.run(cmd, stdout=subprocess.PIPE)
                for ln in proc.stdout.decode().split('\n'):
                    # Tyr35(1.43).OH-Ser75(2.55).OG: 3.56535 Å; -4.99529 kJ/mol.
                    # CNTCT 77 97
                    ln = ln.split(':')[0]
                    ln = re.sub("\\([0-9.]*\\)[.][A-Z0-9]+", "", ln)
                    ln = re.sub("[^0-9-]", "", ln)
                    if not "-" in ln: continue
                    ln = re.sub("-", " ", ln).strip()
                    if ln: newcfg.append("CNTCT "+ln)

            cavfn = f"{pdbdir}/{rcpid}.{suff}.cvty"
            if not os.path.exists(cavfn):
                with open('data/cavopts.json', 'r') as file:
                    cavopts = json.load(file)
                cmd = ["bin/cavity_search", "-p", f"{pdbdir}/{rcpid}.{suff}.pdb", "-o", cavfn]
                cmd.append(cavopts)
                print(" ".join(cmd))
                subprocess.run(cmd)
            newcfg.append(f"VCVTY {cavfn}")

            if skipdock: break

            if fam[0:2] == "OR":
                famno = int(re.sub("[^0-9]", "", fam))
                if famno < 50:
                    newcfg.append("CNTCT data/OR_ClassII_a.ic")
                else:
                    newcfg.append("CNTCT data/OR_ClassI_a.ic")

            if pocket:
                if "atomto" in pocket:
                    if isinstance(pocket["atomto"], str):
                        pocket["atomto"] = [pocket["atomto"]]
                    for a2 in pocket["atomto"]:
                        newcfg.append("ATOMTO " + a2)
                if "bridge" in pocket:
                    if isinstance(pocket["bridge"], str):
                        pocket["bridge"] = [pocket["bridge"]]
                    for st in pocket["bridge"]:
                        newcfg.append("BRIDGE " + st)
                if "search" in pocket:
                    if isinstance(pocket["search"], str):
                        newcfg.append("SEARCH " + pocket["search"])
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
                if "vest" in pocket:
                    print(pocket["vest"])
                    if isinstance(pocket["vest"], str):
                        pocket["vest"] = [pocket["vest"]]
                    for vb in pocket["vest"]:
                        newcfg.append("VESTIBULE " + vb)
                        print("VESTIBULE " + vb)
                if "mcoord" in pocket:
                    print(pocket["mcoord"])
                    newcfg.append("MCOORD " + pocket["mcoord"])

            newcfga = "\n".join(newcfg)

            if pocket:
                # TODO: bridge, flxr, and stcr
                if suff != "inactive" and "atomtoa" in pocket:
                    if isinstance(pocket["atomtoa"], str):
                        pocket["atomtoa"] = [pocket["atomtoa"]]
                    for a2 in pocket["atomtoa"]:
                        newcfga += "\n" + "ATOMTO " + a2
                if suff == "inactive" and "atomtoi" in pocket:
                    if isinstance(pocket["atomtoa"], str):
                        pocket["atomtoa"] = [pocket["atomtoa"]]
                    for a2 in pocket["atomtoa"]:
                        newcfga += "\n" + "ATOMTO " + a2

            with open("tmp/" + conffn, 'w') as f:
                f.write(newcfga + "\n\n")

            os.chdir(os.path.dirname(os.path.abspath(__file__)))
            cmd = ["bin/aromadock", "tmp/" + conffn]
            print(" ".join(cmd))
            subprocess.run(cmd)

        if os.path.exists("tmp/nodelete"):
            print("Warning: not deleting temporary config files because you have the debug \"nodelete\" option selected.")
        else:
            for suff in runsuff:
                rmvnm = "tmp/" + f"{rcpid}~{lignu}.{suff}.config"
                if os.path.exists(rmvnm): os.remove(rmvnm)

        print("Completed", rcpid, o["full_name"])

