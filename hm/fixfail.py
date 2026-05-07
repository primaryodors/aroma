import sys
import re
import os
import os.path
import json
import shutil
import subprocess
from modeller import *
from modeller.automodel import *

if "help" in sys.argv:
    print("Usage:\npython3 hm/fixfail.py [receptor] [docked PDB or ligand_name] ([options]...)\n")
    print("Available options:")
    print("  help\tShow this message.")
    print("  reset\tCreate a new homology model, using dohm.py, to use for processing, before performing any other function.")    
    print("  fit\tConduct a best-fit search for placing the ligand into one of the receptor's cavities. Requires a ligand name instead of a docked PDB.")
    print("  predock\tCreate a new docked PDB by attempting to dock a named ligand in the receptor before processing.")
    print("  loop\tAfter creating a fixed-fail model, attempt to dock the ligand and, if less that three output poses, go back and try the fixfail procedure again for a maximum of 100 attempts.")
    print("  dock\tAfter creating a fixed-fail model, attempt to dock the ligand and then exit.")

argc = len(sys.argv)
if argc < 2:
    # print("Both a protein ID and a PDB file are required.")
    print("A protein ID is required. Optionally, a PDB file may also be specified.")
    print("Example usage:")
    print("python3 hm/fixfail.py OR7D4 out/OR7/OR7D4/OR7D4~androstenone.active.model1.pdb")
    exit()

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

os.chdir("..")

aaletts = "ARNDCEQGHILKMFPSTWYV"
aacode3 = ["ALA", "ARG", "ASN", "ASP", "CYS",
           "GLU", "GLN", "GLY", "HIS", "ILE",
           "LEU", "LYS", "MET", "PHE", "PRO",
           "SER", "THR", "TRP", "TYR", "VAL"]

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
mode="active"
protid = sys.argv[1]
if protid[-1] == 'a':
    mode="active"
    protid = protid[0:-1]
elif protid[-1] == 'i':
    mode="inactive"
    protid = protid[0:-1]
if not protid in data.protutils.prots.keys():
    print("Protein", protid, "not found.")
    exit()
fam = data.protutils.family_from_protid(protid)
origpdb = f"pdbs/{fam}/{protid}.{mode}.pdb"

if argc > 2:
    inppdb = sys.argv[2]
    # Bypass file checks if the user provided a command keyword instead of a PDB file
    if inppdb in ["help"]:
        argc = 2  # Downgrade argument count so downstream logic ignores the missing PDB
    elif inppdb in ["reset", "fit", "predock", "loop", "dock"]:
        print("Second argument must be a path to a PDB file or a protein ID.", file=sys.stderr)
        exit()
    else:
        if not os.path.exists(inppdb):
            inppdb = f"out/{fam}/{protid}/{protid}~{sys.argv[2]}.{mode}.model1.pdb"
        if not os.path.exists(inppdb):
            print(f"Input file not found: {sys.argv[2]}")
            exit()

# If PDB files get too wonky, you can kill the fixfail.py process and rerun it with the reset argument.
doreset = False
if "reset" in sys.argv:
    # Rebuild the unmodified HM structure.
    cmd = ["python3", "hm/dohm.py", protid]
    data.globals.wait_cool_cpu()
    print(" ".join(cmd))
    subprocess.run(cmd)
    doreset = True

if argc > 2:
    odor = data.odorutils.find_odorant(inppdb.split('~')[1].split('.')[0])
else:
    odor = False

if odor:
    delcav = inppdb[0:-3]+"cvty"
    if os.path.exists(delcav):
        os.unlink(delcav)
    with open('data/cavopts.json', 'r') as file:
        cavopts = json.load(file)
    cmd = ["bin/cavity_search", "-p", origpdb if doreset else inppdb, "-o", delcav]
    cmd.extend(cavopts)
    print(" ".join(cmd))
    data.globals.wait_cool_cpu()
    subprocess.run(cmd)
else:
    inppdb = origpdb
    odor = False

delete_inppdb = False

# Ex.: python3 hm/fixfail.py OR5K1 hazelnut_pyrazine N7 279:OG1 C10 255:CG C1 104:CG
# Note: Thr must specify OG1, not OG2, or MODELLER will complain!
ligcontacts = []
if argc > 4:
    for i in range(3, argc-1):
        if re.match("[A-Za-z]+[0-9]+$", sys.argv[i]):
            j = i+1
            print(f"Matched {sys.argv[i]}")
            if re.match("[A-Za-z]{0,3}[0-9]+[:][0-9]?[A-Z]+[0-9]?$", sys.argv[j]):
                ligcontacts.append([sys.argv[i], sys.argv[j]])
                i=j

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")

if "fit" in sys.argv:
    # Perform a cavity fitting calculation on the ligand and receptor model.
    if not odor:
        print("An odorant is required for cavity prefit.")
        exit()
    origcvty = origpdb.replace(".pdb", ".cvty")
    sdfname = data.odorutils.ensure_sdf_exists(odor["full_name"])
    inppdb = f"tmp/{protid}.fit.pdb"
    delete_inppdb = True
    cmd = ["bin/cavity_fit", origpdb, origcvty, sdfname, "-a", "-c", "65536", "-n", "1", "-o", inppdb]
    data.globals.wait_cool_cpu()
    print(" ".join(cmd))
    subprocess.run(cmd)
    # TODO: If no result from cavity_fit, error out and exit.

if "predock" in sys.argv or ("reset" in sys.argv and not "fit" in sys.argv):
    # Perform an active-state dock on the new model...
    if odor:
        cmd = ["/bin/bash", "./dock.sh", protid, odor["full_name"], "noi"]
        data.globals.wait_cool_cpu()
        print(" ".join(cmd))
        subprocess.run(cmd)

tries = 0
while True:
    legal = ""
    tpls = ""
    with open(origpdb, "r") as f:
        c = f.read()
        for ln in c.split("\n"):
            if ln[0:11] == "REMARK   1 ":
                legal += ln + "\n"
            elif ln[0:23] == "REMARK 265 HM_TEMPLATES":
                tpls = ln[24:].strip()

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")

    cmd = ["make", "apps"]
    data.globals.wait_cool_cpu()
    subprocess.run(cmd)

    print(f"Protein: {protid}\nInput PDB: {inppdb}")

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")

    prot = data.protutils.prots[protid]
    seqlen = len(prot["sequence"])
    is_helix = '-' * seqlen
    for rgid in prot["region"]:
        if rgid[0:3] != "TMR" and rgid[0:3] != "HXR": continue
        rgnse = prot["region"][rgid]
        for resno in range(rgnse["start"], rgnse["end"]+1):
            is_helix = is_helix[0:resno-1] + '*' + is_helix[resno:]
    bw45_50 = data.protutils.resno_from_bw(protid, "45.50")
    if bw45_50:
        for resno in range(bw45_50, bw45_50+9):                         # include D/E45.51 in "helix".
            is_helix = is_helix[0:resno-1] + '*' + is_helix[resno:]

    startres = 0
    if (mode=="inactive"):  tplpdb = "hm/"+protid+"i_tpl.pdb"
    else:                   tplpdb = "hm/"+protid+"_tpl.pdb"

    rshpmft = "GPCR.ic"          # fake filename that doesn't exist so we can error out later

    if protid[0:2] == "OR":
        if int(re.sub(r"[^0-9]", "", protid[2:4])) >= 50:
            if (mode=="inactive"):      rshpmft = "OR_ClassI_i.ic"    # mota areberetus iagiion
            else:                       rshpmft = "OR_ClassI_a.ic"
        else:
            if (mode=="inactive"):      rshpmft = "OR_ClassII_i.ic"
            else:                       rshpmft = "OR_ClassII_a.ic"
    elif protid[0:2] == "TAAR":
        if (mode=="inactive"):          rshpmft = "TAAR_i.ic"
        else:                           rshpmft = "TAAR_a.ic"
    elif protid[0:2] == "VN1R":
        if (mode=="inactive"):          rshpmft = "VN1R_i.ic"
        else:                           rshpmft = "VN1R_a.ic"          # future expansion

    rshpmfn = "data/" + rshpmft
    if not os.path.exists(rshpmfn):
        print("No reshape file exists for this receptor.")
        exit()

    if len(ligcontacts):
        cmd = ["bin/poslig", inppdb, tplpdb]
        for lc in ligcontacts:
            cmd.extend(lc)
        print(" ".join(cmd))
        data.globals.wait_cool_cpu()
        subprocess.run(cmd)
    else:
        shutil.copyfile(inppdb, tplpdb)

    cmd = ["bin/ic", tplpdb, rshpmfn, "save"]
    print(" ".join(cmd))
    data.globals.wait_cool_cpu()
    subprocess.run(cmd)

    haslig = False
    with open(tplpdb, 'r') as ftpl:
        cout = ""
        c = ftpl.read()
        for ln in c.split("\n"):
            if ln[0:6] == "ATOM  ":
                resno = int(ln[22:28].strip())
                if resno < 1 or resno > len(is_helix): continue
                if is_helix[resno-1] == '-': continue
                if not startres or (resno and resno<startres): startres = resno
                cout += ln + "\n"
            if ln[0:6] == "HETATM":
                haslig = True

    with open(inppdb, 'r') as fin:
        cin = fin.read()
        for ln in cin.split("\n"):
            if ln[0:6] == "REMARK":
                cout = ln + "\n" + cout
            # if ln[0:6] == "ATOM  ":
              #  resno = int(ln[22:28].strip())
               # if resno < 1 or resno > len(is_helix): continue
                #if is_helix[resno-1] == '-': continue
                # if not startres or (resno and resno<startres): startres = resno
            if ln[0:6] == "HETATM":
                ln = ln[0:23] + "999" + ln[26:]
                cout += ln + "\n"
                haslig = True
        with open(tplpdb, "w") as fout:
            fout.write(cout)

        if delete_inppdb:
            os.remove(inppdb)

    env = Environ()
    if haslig:
        env.io.hetatm = True

    class AromaReceptorModel(DOPEHRLoopModel):
        def select_loop_atoms(self):
            return Selection(self.residue_range("1:A", f"{seqlen}:A"))

        def special_restraints(self, aln):
            rsr = self.restraints
            at = self.atoms

            prot = data.protutils.prots[protid]

            # TM helices and helix 8
            for rgid in prot["region"]:
                if rgid[0:3] != "TMR": continue # and rgid[0:3] != "HXR": continue
                rgnse = prot["region"][rgid]
                rsr.add(secondary_structure.Alpha(self.residue_range(str(rgnse["start"])+":A", str(rgnse["end"])+":A")))
                print("Helix "+str(rgnse["start"])+" - "+str(rgnse["end"]))

            # 45.52 - 45.58 helix
            bw45_50 = data.protutils.resno_from_bw(protid, "45.50")
            if bw45_50:
                rsr.add(secondary_structure.Alpha(self.residue_range(str(bw45_50+2)+":A", str(bw45_50+8)+":A")))

            # disulfide bond
            bw3_25 = data.protutils.resno_from_bw(protid, "3.25")
            if bw3_25 and bw45_50:
                if data.protutils.aalet_at_resno(protid, bw3_25) == 'C' and data.protutils.aalet_at_resno(protid, bw45_50) == 'C':
                    rsr.add(forms.Gaussian(group=physical.xy_distance,
                        feature=features.Distance(at["SG:"+str(bw3_25)+":A"],
                                                at["SG:"+str(bw45_50)+":A"]),
                                                mean=2.05, stdev=0.2))

            rigresno = self.residues[-1].intnum
            rigb = RigidBody(self.residue_range(f"{rigresno}:A", f"{rigresno}:A"))
            rsr.rigid_bodies.append(rigb)

            # Ligand-receptor contacts
            if len(ligcontacts):
                for lc in ligcontacts:
                    liga = lc[0]
                    resno, resa = lc[1].split(':')
                    rsr.add(forms.Gaussian(group=physical.xy_distance,
                        feature=features.Distance(at[f"{resa}:{resno}:A"],
                                                  at[f"{liga}:{rigresno}:A"]),
                                                  mean=3.5, stdev=0.5))

            # 5-7 tyrosine link
            bw5_58 = data.protutils.resno_from_bw(protid, "5.58")
            bw7_53 = data.protutils.resno_from_bw(protid, "7.53")
            if bw5_58 and bw7_53:
                if data.protutils.aalet_at_resno(protid, bw5_58) == 'Y' and data.protutils.aalet_at_resno(protid, bw7_53) == 'Y':
                    rsr.add(forms.Gaussian(group=physical.xy_distance,
                        feature=features.Distance(at["OH:"+str(bw5_58)+":A"],
                                                at["OH:"+str(bw7_53)+":A"]),
                                                mean=4.6, stdev=1.2))

            # 6-45 hydrogen bond
            bw6_55 = data.protutils.resno_from_bw(protid, "6.55")
            if bw6_55:
                atom6_55 = False
                if data.protutils.aalet_at_resno(protid, bw6_55) == 'Y':
                    atom6_55 = "OH"
                elif data.protutils.aalet_at_resno(protid, bw6_55) == 'H':
                    atom6_55 = "NE2"
                if atom6_55:
                    bw45_51 = bw45_50+1
                    atom45_51 = False
                    if data.protutils.aalet_at_resno(protid, bw45_51) == 'D':
                        atom45_51 = "OD1"
                    elif data.protutils.aalet_at_resno(protid, bw45_51) == 'E':
                        atom45_51 = "OE1"
                    elif data.protutils.aalet_at_resno(protid, bw45_51) == 'H':
                        atom45_51 = "NE2"
                    elif data.protutils.aalet_at_resno(protid, bw45_51) == 'N':
                        atom45_51 = "OD1"
                    elif data.protutils.aalet_at_resno(protid, bw45_51) == 'Q':
                        atom45_51 = "OE1"
                    if atom45_51:
                        rsr.add(forms.Gaussian(group=physical.xy_distance,
                            feature=features.Distance(at[atom6_55 +":"+str( bw6_55)+":A"],
                                                    at[atom45_51+":"+str(bw45_51)+":A"]),
                                                    mean=2.5, stdev=0.5))

        def special_patches(self, aln):
            # disulfide bond:
            bw3_25 = data.protutils.resno_from_bw(protid, "3.25")
            bw45_50 = data.protutils.resno_from_bw(protid, "45.50")
            if bw3_25 and bw45_50:
                if data.protutils.aalet_at_resno(protid, bw3_25) == 'C' and data.protutils.aalet_at_resno(protid, bw45_50) == 'C':
                    self.patch(residue_type='DISU', residues=(self.residues[str(bw3_25)+':A'],
                                                            self.residues[str(bw45_50)+':A']))

    # scan input pdb for sequence
    seq = ""
    lrno = 0
    with open(tplpdb, "r") as f:
        c = f.read()
        lines = c.split("\n")
        for ln in lines:
            if ln[0:5].strip() == "ATOM":
                resno = int(ln[22:28].strip())
                if resno:
                    if resno == lrno: continue
                    while lrno < resno-1:
                        seq += "-"
                        lrno += 1
                    aacode = ln[17:20]
                    if not aacode in aacode3:
                        seq += "."
                    else:
                        i = aacode3.index(aacode)
                        seq += aaletts[i]
                    lrno = resno

    # print(seq)

    os.chdir("hm")
    cmd = ["python3", "hm/build_alignment_file.py"]
    data.globals.wait_cool_cpu()
    subprocess.run(cmd)

    # scan ali file for protid
    alidat = ""
    with open("allgpcr.ali", "r") as f:
        c = f.read()
        lines = c.split("\n")
        reading = False
        lookfor = ">P1;"+protid
        for ln in lines:
            if ln.strip() == lookfor:
                reading = True
                continue
            elif reading:
                alidat += ln+"\n"
                if ln.find('*') >= 0:
                    break

    if haslig:
        alidat = alidat.replace("*", ".*")

    # duplicate the alidat var applying any gaps in seq
    alitpl = ""
    n = len(alidat)
    j = 0
    cryet = False
    for i in range(n):
        c = alidat[i]
        if not len(c): continue
        if ord(c) < ord(" "):
            cryet = True
            alitpl += c
        elif not cryet:
            alitpl += c
        elif ord(c) < ord('A') or ord(c) > ord('Z'):
            alitpl += c
        else:
            if j >= len(seq):
                d = '-'
            else:
                d = seq[j]
            j += 1
            if d == c or d == '-':
                alitpl += d
            else:
                alitpl += d         # some of the sequences do not match; use the PDB sequence.
                # print("Something went wrong:\n\n" + alitpl + "\n" + alidat + "\n" + c + "~" + d)
                # exit()

    if mode == "inactive":  tplfttl = protid+"i_tpl"
    else:                   tplfttl = protid+"_tpl"

    alitpl = alitpl.replace(protid, tplfttl)
    alitpl = alitpl.replace("sequence", "structure")
    startres_pad5 = f"{startres:<{5}}"
    alitpl = re.sub(":([0-9]+\\s+):([A-Z]):([0-9]+\\s+):([A-Z])", f":{startres_pad5}:\\2:999  :\\4", alitpl)

    with open("experimental.ali", "r") as f:
        expali = f.read().__str__()

    tmpalif = protid + "_tmp.ali"
    with open(tmpalif, "w") as f:
        f.write(f">P1;{protid}\n")
        f.write(alidat + "\n")
        f.write(f">P1;{tplfttl}\n")
        f.write(alitpl + "\n")
        f.write(f"\n\n{expali}\n\n")

    # directories for input atom files
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    env.io.atom_files_directory = ['.', './tpl']

    fix8uy0 = ""
    # t4hm = data.protutils.templates_for_hm(protid)
    # tpls += " " + " ".join(t4hm)
    # tpls = " ".join(list(set(tpls.split(" "))))
    # print(legal)
    # print(tpls)
    tplsfttls = (tplfttl,) # + tuple(t4hm)
    if "8uy0" in tplsfttls: fix8uy0 = "BRIDGE DE45.51 Y6.55"
    print(f"Templates: {tplsfttls}")
    a = AromaReceptorModel( env,
                            alnfile           = tmpalif,
                            knowns            = tplsfttls,
                            sequence          = protid
                        )
    a.starting_model = 0
    a.ending_model   = 9
    a.library_schedule = autosched.slow
    a.max_var_iterations = 300

    data.globals.wait_cool_cpu()
    a.make()

    # Find the best output file from MODELLER. If no output to use, exit the script.
    results = [x for x in a.loop.outputs if x['failure'] is None]
    if not len(results):
        print("FAIL.")
        exit()
    key = 'molpdf'
    i = 0
    j = 0
    best = 0.0
    for i in range(len(results)):
        ophis = results[i]
        if key in ophis:
            if not i or float(ophis[key]) < best:
                j = i
                best = float(ophis[key])
        i += 1
    # exit()
    model = results[j]
    print("Chose " + model['name'])

    phewcode = f"""
LET $rcpid = "{protid}"
LET $inpf = "pdbs/{fam}/{protid}.inactive.pdb"
LET $mdld = "hm/{model['name']}"

LOAD $inpf A I
LET %rcpseqln = %SEQLENI
LOAD $mdld A A

BWCOPY I A
STRAND I
UPRIGHT I
BWCENTER
STRAND A

{fix8uy0}

{legal}
REMARK 265 HM_TEMPLATES: {tpls}

UNCHAIN I
UNCHAIN O
STRAND A
UPRIGHT
BWCENTER
# IF $3.37 != "G" THEN ATOMTO %3.37 EXTENT @6.48

IF $3.25 != "C" OR $45.50 != "C" GOTO _not_disulfide
# DELATOM %3.25 HG
# DELATOM %45.50 HG
CONECT %3.25 SG %45.50 SG
_not_disulfide:
HYDRO

LET %y6 = 0
LET $atom6 = "OH"
IF $6.55 = "Y" THEN LET %y6 = %6.55
IF $6.55 = "H" THEN LET %y6 = %6.55
IF $6.55 = "H" THEN LET $atom6 = "NE2"
LET %de45 = 0
IF %45.51 = "D" OR %45.51 = "E" THEN LET %de45 = %45.51
IF NOT %y6 THEN GOTO _no_456_contact
IF NOT %de45 THEN GOTO _no_456_contact
BRIDGE %y6 %de45
_no_456_contact:

IF $5.42 != "C" OR $5.43 != "C" GOTO _not_Cu_binding_site
IF $5.39 != "M" GOTO _not_Cu_539
MEASURE %5.39 "SD" %5.42 "SG" &sdist
ECHO "5.39:SD - 5.42:SG distance: " &sdist
MEASURE %5.39 "SD" %5.43 "SG" &sdist
ECHO "5.39:SD - 5.43:SG distance: " &sdist
_not_Cu_539:
IF $5.46 != "M" GOTO _not_Cu_539
MEASURE %5.46 "SD" %5.42 "SG" &sdist
ECHO "5.46:SD - 5.42:SG distance: " &sdist
MEASURE %5.46 "SD" %5.43 "SG" &sdist
ECHO "5.46:SD - 5.43:SG distance: " &sdist
_not_Cu_546:
MEASURE %5.42 "SG" %5.43 "SG" &sdist
ECHO "5.42:SD - 5.43:SG distance: " &sdist
_not_Cu_binding_site:

LET $outf = "pdbs/{fam}/{protid}.{mode}.pdb"
SAVE $outf
"""

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")
    phewfn = f"hm/{protid}.{mode}.hm.phew"
    with open(phewfn, "w") as f:
        f.write(phewcode)
    cmd = ["bin/phew", phewfn]
    data.globals.wait_cool_cpu()
    subprocess.run(cmd)

    cmd = ["bin/ic", f"pdbs/{fam}/{protid}.{mode}.pdb", "minc", "save"]
    print(" ".join(cmd))
    data.globals.wait_cool_cpu()
    subprocess.run(cmd)

    if not "loop" in sys.argv and not "dock" in sys.argv:
        break
    else: tries += 1

    if odor:
        cmd = ["/bin/bash", "./dock.sh", protid, odor["full_name"], "noi"]
        data.globals.wait_cool_cpu()
        subprocess.run(cmd)

    dock_success = False
    dockfile = f"out/{fam}/{protid}/{protid}~{odor['full_name']}.{mode}.dock"
    if os.path.exists(dockfile):
        with open(dockfile, "r") as f:
            c = f.read()
            lines = c.split("\n")
            total = 0.0
            disqo = False
            poses = 0
            for ln in lines:
                if ln[0:7] == "Total: ":
                    if not total: total = float(ln[7:].strip())
                elif ln[0:22] == "Disqualified because: ":
                    disqo = True
                elif "pose(s) found." in ln:
                    poses = int(ln.split(' ')[0].strip())

                if poses >= 4 and total < 0 and not disqo:
                    print("SUCCESS")
                    dock_success = True
                    break

    if dock_success or not "loop" in sys.argv:
        break

    if tries >= 100:
        print("FAILED after 100 tries.")
        break
