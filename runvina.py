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

onlynew = False
onlymissing = False

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

if len(sys.argv) > 3:
    for i in range(3, len(sys.argv)):
        if sys.argv[i] == "refresh": onlynew = True
        if sys.argv[i] == "resume": onlymissing = True

for rcpid in data.protutils.prots.keys():
    fam = data.protutils.family_from_protid(rcpid)
    if protid:
        if rcpid != protid: continue
    else:
        if popt != "*" and popt != "emp" and popt != "ago":
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
        forms = data.odorutils.check_forms(o["full_name"])
        data.odorutils.ensure_sdf_exists(o["full_name"])

        pocket = data.dyncenter.get_pocket(rcpid, o["full_name"])

        if "odorophores" in pocket:
            for odorophore in pocket["odorophores"].keys():
                cmd = ["test/moiety_test", f"sdf/{lignu}.sdf", "odorophore"]
                print(" ".join(cmd))
                proc = subprocess.run(cmd, stdout=subprocess.PIPE)
                for ln in proc.stdout.decode().split('\n'):
                    if re.match(" occurs [1-9][0-9]* times", ln):
                        for key in pocket["odorophores"][odorophore]:
                            pocket[key] = pocket["odorophores"][odorophore][key]

        flxr = []
        stcr = []
        atomto = []
        if "flxr" in pocket:
            flxr = pocket["flxr"]
            if isinstance(flxr, str):
                flxr = [flxr]
        if "stcr" in pocket:
            stcr = pocket["stcr"]
            if isinstance(stcr, str):
                stcr = [stcr]
        if "atomto" in pocket:
            atomto = pocket["atomto"]
            if isinstance(atomto, str):
                atomto = [atomto]
        if "pocket" in pocket:
            pocket = pocket["pocket"]
        # print(atomto)
        # print(pocket)
        # exit()

        conffna = rcpid + "~" + lignu + ".active.config"
        conffni = rcpid + "~" + lignu + ".inactive.config"
        if not os.path.exists("output"): os.mkdir("output")
        if not os.path.exists("output/" + fam): os.mkdir("output/" + fam)
        if not os.path.exists("output/" + fam + "/" + rcpid): os.mkdir("output/" + fam + "/" + rcpid)

        if onlynew or onlymissing:
            outfna = f"output/{fam}/{rcpid}/{rcpid}~{lignu}.active.dock"
            if os.path.exists(outfna) and os.path.getsize(outfna) > 1e5:
                if onlymissing: continue
                fmt = os.path.getmtime(outfna)
                if fmt > os.path.getmtime("data/binding_pocket.json") \
                and fmt > os.path.getmtime(f"sdf/{lignu}.sdf") \
                and fmt > os.path.getmtime(f"pdbs/{fam}/{rcpid}.active.pdb") \
                and fmt > os.path.getmtime(f"pdbs/{fam}/{rcpid}.inactive.pdb"):
                    continue

        print(f"Beginning {rcpid} ~ "+o["full_name"]+"...")
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # TODO: Apply atomto directives

        # Determine which residues will be flexible and which will be rigid
        resflexy = [False] * len(data.protutils.prots[rcpid]["sequence"])
        for pkt in pocket.split(' '):
            aminos = re.search("^[A-Z]+", pkt)
            if aminos: aminos = aminos.group()
            bw = re.search("[0-9]{1,2}[.][0-9]{2}", pkt)
            if not bw: bw = re.search("[0-9]+", pkt)
            if bw: bw = bw.group()
            else:
                print(f"Bad pocket residue: {pkt}")
                exit()
            bang = re.match("[!]", pkt)
            resno = data.protutils.resno_from_bw(rcpid, bw)
            if not resno: continue
            if aminos:
                lett = data.protutils.aalet_at_resno(rcpid, resno)
                if not lett in aminos:
                    if bang: bang = False
                    else: continue
            resflexy[resno-1] = True
            # print(f"{resno} is bsr")

        for flx in flxr:
            for flxres in flx.split(' '):
                aminos = re.search("^[A-Z]+", flxres)
                if aminos: aminos = aminos.group()
                bw = re.search("[0-9]{1,2}[.][0-9]{2}", flxres)
                if not bw: bw = re.search("[0-9]+", flxres)
                if bw: bw = bw.group()
                else:
                    print(f"Bad flex residue: {flxres}")
                    exit()
                resno = data.protutils.resno_from_bw(rcpid, bw)
                if not resno: continue
                if aminos:
                    lett = data.protutils.aalet_at_resno(rcpid, resno)
                    if not lett in aminos: continue
                resflexy[resno-1] = True
                print(f"{resno} is flexy")

        for stc in stcr:
            for stcres in stc.split(' '):
                aminos = re.search("^[A-Z]+", stcres)
                if aminos: aminos = aminos.group()
                bw = re.search("[0-9]{1,2}[.][0-9]{2}", stcres)
                if not bw: bw = re.search("[0-9]+", stcres)
                if bw: bw = bw.group()
                else:
                    print(f"Bad flex residue: {stcres}")
                    exit()
                resno = data.protutils.resno_from_bw(rcpid, bw)
                if not resno: continue
                if aminos:
                    lett = data.protutils.aalet_at_resno(rcpid, resno)
                    if not lett in aminos: continue
                resflexy[resno-1] = False
                print(f"{resno} is sticky")

        rcpfn = f"pdbs/{fam}/{rcpid}.active.pdb"
        with open(rcpfn, 'r') as f:
            pdb = f.read()
            lines = pdb.split("\n")
            rigidpdb = ""
            flexpdb = ""
            for ln in lines:
                if ln[0:6] == "ATOM  ":
                    resno = re.sub("[^0-9]", "", ln[22:27])
                    if not resno: continue
                    resno = int(resno)
                    if not resno: continue
                    if resflexy[resno-1]:
                        flexpdb += ln + "\n"
                    else:
                        rigidpdb += ln + "\n"
        rigidfn = f"tmp/{rcpid}.rigid.pdb"
        flexfn = f"tmp/{rcpid}.flex.pdb"
        with open(rigidfn, 'w') as f:
            f.write(rigidpdb)
        with open(flexfn, 'w') as f:
            f.write(flexpdb)

        # Convert receptor to PDBQT
        cmdr = ["obabel", "-i", "pdb", rigidfn, "-xr", "-o", "pdbqt", "-O", f"tmp/{rcpid}.rigid.pdbqt"]
        cmdf = ["obabel", "-i", "pdb", flexfn, "-xs", "-o", "pdbqt", "-O", f"tmp/{rcpid}.flex.pdbqt"]
        print(" ".join(cmdr))
        subprocess.run(cmdr)
        print(" ".join(cmdf))
        subprocess.run(cmdf)

        # Convert ligand to PDBQT
        cmdl = ["obabel", "-i", "sdf", f"sdf/{lignu}.sdf", "-o", "pdbqt", "-O", f"tmp/{lignu}.pdbqt"]
        print(" ".join(cmdl))
        subprocess.run(cmdl)

        # Call Vina
        cmd = ["../AutoDock-Vina/build/linux/release/vina", "--receptor", f"tmp/{rcpid}.rigid.pdbqt",
            "--flex", f"tmp/{rcpid}.flex.pdbqt",
            "--ligand", f"tmp/{lignu}.pdbqt",
            "--center_x", "0", "--center_y", "15", "--center_z", "0",
            "--size_x", "20", "--size_y", "20", "--size_z", "20",
            "--exhaustiveness", "20", "--cpu", "1"]
        print(" ".join(cmd))
        subprocess.run(cmd)

        # TODO: interpret output

        if os.path.exists("tmp/nodelete"):
            print("Warning: not deleting temporary config file because you have the debug \"nodelete\" option selected.")
        else:
            os.remove("tmp/" + conffna)
            os.remove("tmp/" + conffni)

        print("Completed", rcpid, o["full_name"])

