import sys
import re
import os
import os.path
import json
import shutil
import subprocess
from modeller import *
from modeller.automodel import *

argc = len(sys.argv)
if argc < -3:
    print("Both a protein ID and a PDB file are required.")
    print("Example usages:")
    print("python3 hm/activate.py OR7A17 pdbs/OR7/accommodated/OR7A17.accommodated.pdb")
    print("python3 hm/activate.py OR7A17 pdbs/OR7/OR7A17.inactive.pdb sdf/vetynal.pdb")
    print("python3 hm/activate.py OR7A17 pdbs/OR7/OR7A17.inactive.pdb olfactophores/OR7/olfactophore_7A17.pdb")
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
protid = sys.argv[1]
if not protid in data.protutils.prots.keys():
    print("Protein", protid, "not found.")
    exit()
fam = data.protutils.family_from_protid(protid)

ligand = False
if argc > 2:
    inppdb = sys.argv[2]
else:
    inppdb = f"pdbs/{fam}/{protid}.inactive.pdb"                  # ALWAYS use the inactive model aas the input.
    member = protid[4:]
    if protid[0:2] == "OR": member = protid[2:]
    lookfor = f"olfactophores/{fam}/olfactophore_{member}.pdb"
    if os.path.exists(lookfor):
        ligand = lookfor

if not ligand and len(sys.argv) > 3:
    ligand = sys.argv[3]

print(f"Protein: {protid}\nInput PDB: {inppdb}")
if ligand:
    print(f"Ligand: {ligand}")
else:
    print("No ligand")

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
if ligand:
    cmd = ["make", "bin/olfactophore"]
    subprocess.run(cmd)
    outdir = "pdbs/"+fam+"/accommodated"
    outfile = outdir+"/"+protid+".accommodated.pdb"
    cmd = ["bin/olfactophore", inppdb, ligand, "-o", outfile]
    if not os.path.exists(outdir): os.mkdir(outdir)
    print(" ".join(cmd))
    subprocess.run(cmd)
    if not os.path.exists(outfile):
        print("Failed.")
        exit()
    inppdb = outfile

shutil.copyfile(inppdb, "hm/"+protid+"_tpl.pdb")

env = Environ()


class AromaReceptorModel(DOPEHRLoopModel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms

        prot = data.protutils.prots[protid]

        # TM helices and helix 8
        for rgid in prot["region"]:
            if rgid[0:3] != "TMR" and rgid[0:3] != "HXR": continue
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
                else:
                    bw45_52 = bw45_51+1
                    atom45_52 = False
                    if data.protutils.aalet_at_resno(protid, bw45_52) == 'D':
                        atom45_52 = "OD1"
                    elif data.protutils.aalet_at_resno(protid, bw45_52) == 'E':
                        atom45_52 = "OE1"
                    elif data.protutils.aalet_at_resno(protid, bw45_52) == 'H':
                        atom45_52 = "NE2"
                    elif data.protutils.aalet_at_resno(protid, bw45_52) == 'N':
                        atom45_52 = "OD1"
                    elif data.protutils.aalet_at_resno(protid, bw45_52) == 'Q':
                        atom45_52 = "OE1"
                    if atom45_52:
                        rsr.add(forms.Gaussian(group=physical.xy_distance,
                            feature=features.Distance(at[atom6_55 +":"+str( bw6_55)+":A"],
                                                      at[atom45_52+":"+str(bw45_52)+":A"]),
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
with open(inppdb, "r") as f:
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
cmd = ["php", "-f", "build_alignment_file.php"]
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
        d = seq[j]
        j += 1
        if d == c or d == '-':
            alitpl += d
        else:
            print("Something went wrong:\n\n" + alitpl + "\n" + alidat + "\n" + c + "~" + d)
            exit()

alitpl = alitpl.replace(protid, protid+"_tpl")
alitpl = alitpl.replace("sequence", "structure")

tmpalif = protid + "_tmp.ali"
with open(tmpalif, "w") as f:
    f.write(">P1;"+protid+"\n")
    f.write(alidat + "\n")
    f.write(">P1;"+protid+"_tpl\n")
    f.write(alitpl + "\n")

# directories for input atom files
os.chdir(os.path.dirname(os.path.abspath(__file__)))
env.io.atom_files_directory = ['.', '../atom_files']

a = AromaReceptorModel( env,
                        alnfile           = tmpalif,
                        knowns            = protid+"_tpl",
                        sequence          = protid,
                        assess_methods    = (assess.DOPE)
                        )
a.starting_model = 0
a.ending_model   = 9
a.library_schedule = autosched.slow
a.max_var_iterations = 300

a.make()

# TODO: Find the best output file from MODELLER, don't just assume #6. If no output to use, exit the script.
results = [x for x in a.outputs if x['failure'] is None]
if not len(results):
    print("FAIL.")
    exit()
key = 'DOPE score'
results.sort(key=lambda a: a[key], reverse=True)
model = results[0]

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

REMARK   1
REMARK   1 REFERENCE 1
REMARK   1  AUTH 1 B. Webb, A. Sali.
REMARK   1  TITL 1 Comparative Protein Structure Modeling Using Modeller.
REMARK   1  REF  1 Current Protocols in Bioinformatics 54, John Wiley & Sons, Inc., 5.6.1-5.6.37, 2016.
REMARK   1  DOI  1 10.1002/0471250953.bi0506s15
REMARK   1  
REMARK   1  AUTH 2 M.A. Marti-Renom, A. Stuart, A. Fiser, R. Sánchez, F. Melo, A. Sali.
REMARK   1  TITL 2 Comparative protein structure modeling of genes and genomes.
REMARK   1  REF  2 Annu. Rev. Biophys. Biomol. Struct. 29, 291-325, 2000.
REMARK   1  DOI  2 10.1146/annurev.biophys.29.1.291
REMARK   1
REMARK   1  AUTH 3 A. Sali & T.L. Blundell.
REMARK   1  TITL 3 Comparative protein modelling by satisfaction of spatial restraints.
REMARK   1  REF  3 J. Mol. Biol. 234, 779-815, 1993.
REMARK   1  DOI  3 10.1006/jmbi.1993.1626
REMARK   1  
REMARK   1  AUTH 4 A. Fiser, R.K. Do, & A. Sali.
REMARK   1  TITL 4 Modeling of loops in protein structures.
REMARK   1  REF  4 Protein Science 9. 1753-1773, 2000.
REMARK   1  DOI  4 10.1110/ps.9.9.1753
REMARK   1  
REMARK 265 HM_TEMPLATES: custom

# HYDRO

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

LET $outf = "pdbs/{fam}/{protid}.active.pdb"
SAVE $outf
"""

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir("..")
phewfn = f"hm/{protid}.hm.phew"
with open(phewfn, "w") as f:
    f.write(phewcode)
cmd = ["bin/phew", phewfn]
subprocess.run(cmd)


