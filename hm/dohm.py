#!/usr/bin/env python3
"""
Executes homology modeling using MODELLER for a specific GPCR target.
"""

import os
import sys
import subprocess
import glob
import re
import math

# -----------------------------------------
# ARCHITECTURAL ENFORCEMENT: PATH RESOLUTION
# -----------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(script_dir)
sys.path.append(root_dir)

try:
    from modeller import *
    from modeller.automodel import *
except ImportError:
    print("FATAL ERROR: This feature requires MODELLER.", file=sys.stderr)
    print("Please see: https://salilab.org/modeller/ to obtain this third party software.", file=sys.stderr)
    sys.exit(1)

import data.protutils as pu

# Attempt to load optional global utilities if they exist
try:
    import data.globals
    import data.odorutils
    import data.dyncenter
    data.odorutils.load_odors()
except ImportError:
    pass

def execute_alignment_builder():
    print("Executing Alignment Builder...", file=sys.stderr)
    align_script = os.path.join(script_dir, "build_alignment_file.py")
    result = subprocess.run([sys.executable, align_script], cwd=root_dir)
    if result.returncode != 0:
        print("FATAL ERROR: Alignment builder failed.", file=sys.stderr)
        sys.exit(1)

def main():
    if len(sys.argv) < 2:
        print("Usage:\npython3 hm/dohm.py PROTID [riglig] [nodel]\n", file=sys.stderr)
        sys.exit(1)

    rcpid = sys.argv[1].upper()
    riglig = "riglig" in sys.argv
    nodel = "nodel" in sys.argv

    pu.load_prots()
    if rcpid not in pu.prots:
        print(f"FATAL ERROR: Protein ID {rcpid} not found in database.", file=sys.stderr)
        sys.exit(1)

    p = pu.prots[rcpid]
    fam = pu.family_from_protid(rcpid)
    sub = pu.subfamily_from_protid(rcpid)
    famsub = f"{fam}{sub}"
    famno_str = re.sub(r"[^0-9]", "", fam)
    famno = int(famno_str) if famno_str else 0

    # 1. Update Alignments
    execute_alignment_builder()

    # 2. Extract Target Alignment Data
    allgpcr_path = os.path.join(script_dir, "allgpcr.ali")
    if not os.path.exists(allgpcr_path):
        print(f"FATAL ERROR: {allgpcr_path} not generated.", file=sys.stderr)
        sys.exit(1)

    with open(allgpcr_path, "r", encoding="utf-8") as f:
        lines = f.read().split("\n")

    tgtali = ""
    alnhdr = ""
    p1ln = ""
    reading_tgtali = False
    prevln = ""

    for ln in lines:
        if ln.startswith("sequence"):
            pieces = ln.split(':')
            if len(pieces) > 1 and pieces[1].strip() == rcpid:
                alnhdr = ln
                p1ln = prevln
                reading_tgtali = True
        elif reading_tgtali:
            tgtali += ln + "\n"
            if '*' in ln:
                reading_tgtali = False
                break
        prevln = ln

    if not tgtali:
        print(f"FATAL ERROR: Target sequence for {rcpid} not found in allgpcr.ali.", file=sys.stderr)
        sys.exit(1)

    # 3. Build Custom Template
    # We must ensure we are in the hm/ directory for template writing
    os.chdir(script_dir)

    result_str = pu.custom_pdb_template(tgtali, rcpid, f"{rcpid}_tpl.pdb")
    pieces = result_str.split("\n", 1)
    tplsused = pieces[0]
    tplali = pieces[1] if len(pieces) > 1 else ""

    with open(f"{rcpid}.knowns", "w") as f:
        f.write(tplsused)

    hm_ali_file = f"{rcpid}.hm.ali"
    with open(hm_ali_file, "w") as f:
        f.write(f">P1;{rcpid}_tpl\n")
        f.write(f"structure:{rcpid}_tpl:FIRST:A:LAST :A:Olfactory Receptor template:Rhombopteryx nessiteras: 4.00: 0.25\n")
        f.write(f"{tplali}\n\n")
        f.write(f"{p1ln}\n")
        f.write(f"{alnhdr}\n")
        f.write(f"{tgtali}\n\n")

    # 4. Formulate Structural Constraints
    do_tmr_helix_restraints = True
    do_tmr6_helix_restraints = False
    do_exr2_helix_restraint = True

    alpha_helices = []
    if do_tmr_helix_restraints and "region" in p:
        for rgname, rgnse in p["region"].items():
            nmsub3 = rgname[:3]
            if nmsub3 not in ["TMR", "HXR"]: continue
            tmrno = int(re.sub(r"[^0-9]", "", rgname))
            if not do_tmr6_helix_restraints and tmrno == 6: continue
            if rcpid in ["OR10Q1", "OR52A1", "OR52A4"] and nmsub3 == "HXR": continue
            alpha_helices.append((rgnse['start'], rgnse['end']))

    if do_exr2_helix_restraint and rcpid.startswith("OR"):
        try:
            rgs = pu.resno_from_bw(rcpid, "45.52")
            rge = pu.resno_from_bw(rcpid, "45.58")
            if rgs and rge:
                alpha_helices.append((rgs, rge))
        except Exception:
            pass

    # 5. Define Native MODELLER Class
    env = Environ()
    env.io.atom_files_directory = ['.', './tpl']
    if riglig:
        env.io.hetatm = True

    dspotr1 = [pu.resno_from_bw(rcpid, "3.25"),                     # Conserved TMR3-EXR2 bond
                pu.resno_from_bw(rcpid, "45.40"),                   # https://doi.org/10.1002/pro.2717 (goddamn paywalled)
                pu.resno_from_bw(rcpid, "3.40")                     # OR10D/G/S feature, plus a handful of other ORs
                ]
    dspotr2 = [pu.resno_from_bw(rcpid, "45.50"),
                pu.resno_from_bw(rcpid, "45.60"),
                pu.resno_from_bw(rcpid, "5.50")
                ]

    # https://doi.org/10.1016/j.jbc.2026.113319
    if rcpid == "OR5W2":
        dspotr1.append(6)
        dspotr2.append(pu.resno_from_bw(rcpid, "45.38"))

    elif rcpid == "OR4D10":
        dspotr2[1] = 6
    elif rcpid == "OR10A7":
        dspotr2[1] = 3
    elif rcpid == "OR51E2":
        dspotr2[1] = 4

    elif rcpid == "OR52M1":
        dspotr1[1] = 8
    elif rcpid == "OR5L1":
        dspotr1[1] = 6
    elif rcpid == "OR5L2":
        dspotr1[1] = 6
    elif rcpid == "OR2AT4":
        dspotr1[1] = 6

    elif rcpid == "OR56A1":
        dspotr1.append(23)
        dspotr2.append(pu.resno_from_bw(rcpid, "45.48"))
    elif rcpid == "OR56A3":
        dspotr1.append(20)
        dspotr2.append(pu.resno_from_bw(rcpid, "45.48"))
    elif rcpid == "OR56A4":
        dspotr1.append(19)
        dspotr2.append(pu.resno_from_bw(rcpid, "45.48"))
    elif rcpid == "OR56A5":
        dspotr1.append(19)
        dspotr2.append(pu.resno_from_bw(rcpid, "45.48"))

    # Other misc cross links
    elif rcpid == "OR1B1":
        dspotr1[1] = pu.resno_from_bw(rcpid, "45.34")
    elif rcpid == "OR1C1":
        dspotr1.append(pu.resno_from_bw(rcpid, "3.44"))
        dspotr2.append(pu.resno_from_bw(rcpid, "5.53"))
    elif rcpid == "OR1N2" or rcpid == "OR2T10":
        dspotr1.append(pu.resno_from_bw(rcpid, "3.41"))
        dspotr2.append(pu.resno_from_bw(rcpid, "4.49"))
    elif rcpid == "OR1S1" or rcpid == "OR1S2" or rcpid == "OR2A4":
        dspotr1.append(pu.resno_from_bw(rcpid, "3.55"))
        dspotr2.append(pu.resno_from_bw(rcpid, "5.60"))
    elif rcpid == "OR11A1":
        dspotr2[2] = pu.resno_from_bw(rcpid, "5.46")

    dsres1 = []
    dsres2 = []

    class AromaModel(AutoModel):
        def special_patches(self, aln):
            # Add disulfide bridges
            for idx, r1 in enumerate(dspotr1):
                r2 = dspotr2[idx]
                try:
                    if r1 and r2 and p['sequence'][r1-1] == 'C' and p['sequence'][r2-1] == 'C':
                        self.patch(residue_type='DISU', residues=(self.residues[f'{r1}:A'], self.residues[f'{r2}:A']))
                        dsres1.append(r1)
                        dsres2.append(r2)
                except Exception:
                    pass

        def special_restraints(self, aln):
            rsr = self.restraints
            at = self.atoms

            # Apply Alpha Helices
            for rgs, rge in alpha_helices:
                rsr.add(secondary_structure.Alpha(self.residue_range(f'{rgs}:A', f'{rge}:A')))

            # Disulfide distance restraints
            for idx, r1 in enumerate(dsres1):
                r2 = dsres2[idx]
                try:
                    if r1 and r2 and p['sequence'][r1-1] == 'C' and p['sequence'][r2-1] == 'C':
                        rsr.add(forms.Gaussian(group=physical.xy_distance,
                                            feature=features.Distance(at[f'SG:{r1}:A'], at[f'SG:{r2}:A']),
                                            mean=2.05, stdev=0.2))
                except Exception:
                    pass

            # Cu-binding site distance restraints (OR2M/T/V)
            if famsub in ["OR2M", "OR2T", "OR2V"]:
                try:
                    r539, r542, r543, r546 = [pu.resno_from_bw(rcpid, x) for x in ["5.39", "5.42", "5.43", "5.46"]]
                    seq = p['sequence']
                    if seq[r542-1] == 'C' and seq[r543-1] == 'C':
                        mtl_active = False
                        if seq[r539-1] == 'M':
                            mtl_active = True
                            rsr.add(forms.Gaussian(group=physical.xy_distance, feature=features.Distance(at[f'SD:{r539}:A'], at[f'SG:{r542}:A']), mean=4.7, stdev=0.25))
                            rsr.add(forms.Gaussian(group=physical.xy_distance, feature=features.Distance(at[f'SD:{r539}:A'], at[f'SG:{r543}:A']), mean=4.7, stdev=0.25))
                        if seq[r546-1] == 'M':
                            mtl_active = True
                            rsr.add(forms.Gaussian(group=physical.xy_distance, feature=features.Distance(at[f'SD:{r546}:A'], at[f'SG:{r542}:A']), mean=4.7, stdev=0.25))
                            rsr.add(forms.Gaussian(group=physical.xy_distance, feature=features.Distance(at[f'SD:{r546}:A'], at[f'SG:{r543}:A']), mean=4.7, stdev=0.25))
                        if mtl_active:
                            rsr.add(forms.Gaussian(group=physical.xy_distance, feature=features.Distance(at[f'SG:{r542}:A'], at[f'SG:{r543}:A']), mean=3.8, stdev=0.2))
                except Exception:
                    pass

            # Restrain BW 5.47 to extended conformation
            try:
                r547 = pu.resno_from_bw(rcpid, "5.47")
                if r547:
                    aa547 = p['sequence'][r547 - 1]
                    if aa547 in ['L', 'I', 'V', 'M']:
                        cg_name = 'CG1' if aa547 in ['I', 'V'] else 'CG'
                        mean_chi1 = 3.14159 if aa547 in ['I', 'V'] else -3.05
                        rsr.add(forms.Gaussian(group=physical.chi1_dihedral,
                                            feature=features.Dihedral(at[f'N:{r547}:A'], at[f'CA:{r547}:A'], at[f'CB:{r547}:A'], at[f'{cg_name}:{r547}:A']),
                                            mean=mean_chi1, stdev=0.25))
            except Exception:
                pass

            # Prevent BW 5.44 collapse in Class I ORs
            if famno in [51, 52, 56]:
                try:
                    r544 = pu.resno_from_bw(rcpid, "5.44")
                    if r544:
                        aa544 = p['sequence'][r544 - 1]
                        if aa544 in ['I', 'L']:
                            cg_name = 'CG1' if aa544 == 'I' else 'CG'
                            rsr.add(forms.Gaussian(group=physical.chi1_dihedral,
                                                feature=features.Dihedral(at[f'N:{r544}:A'], at[f'CA:{r544}:A'], at[f'CB:{r544}:A'], at[f'{cg_name}:{r544}:A']),
                                                mean=-1.5708, stdev=0.25))
                            cd_name = 'CD1'
                            rsr.add(forms.Gaussian(group=physical.chi2_dihedral,
                                                feature=features.Dihedral(at[f'CA:{r544}:A'], at[f'CB:{r544}:A'], at[f'{cg_name}:{r544}:A'], at[f'{cd_name}:{r544}:A']),
                                                mean=3.14159, stdev=0.25))
                except Exception:
                    pass

    # 6. Execute Homology Model
    print(f"Initiating MODELLER for {rcpid}...", file=sys.stderr)
    a = AromaModel(env, alnfile=hm_ali_file, knowns=f'{rcpid}_tpl', sequence=rcpid)
    a.starting_model = 0
    a.ending_model = 0 # TODO: SET THIS BACK TO 9
    a.library_schedule = autosched.slow
    a.max_var_iterations = 1000

    if rcpid in ["OR2AE1", "OR2AG1", "OR2AG2"]:
        a.md_level = refine.very_slow

    # Throttle if global functions exist
    if 'data.globals' in sys.modules and hasattr(sys.modules['data.globals'], 'wait_cool_cpu'):
        try:
            sys.modules['data.globals'].wait_cool_cpu()
        except FileNotFoundError:
            pass # Hardware sensors unavailable in WSL. Pushing through.

    a.make()

    # 7. Extract Best Model Programmatically (No String Parsing)
    ok_models = [m for m in a.outputs if m['failure'] is None]
    if not ok_models:
        print("FATAL ERROR: MODELLER failed to produce any valid structures.", file=sys.stderr)
        sys.exit(1)

    best_model = min(ok_models, key=lambda m: m['molpdf'])
    best_pdb = best_model['name']
    print(f"Best model generated: {best_pdb} with molpdf {best_model['molpdf']}")

    # 8. Post-Processing Script Generation (.phew)
    adjustments = ""
    if famno < 50:
        adjustments += 'IF $3.37 != "G" THEN ATOMTO %3.37 EXTENT @6.48\n'
    elif famno in [51, 52]:
        adjustments += 'ATOMTO %6.59 EXTENT @4.57\n'
    elif rcpid == "OR56B2":
        adjustments += 'ATOMTO %6.58 EXTENT @4.57\n'

    if famsub == "OR5K":
        adjustments += 'ATOMTO %45.49 EXTENT @2.58\n'

    dsphew = ""
    for idx, r1 in enumerate(dsres1):
        r2 = dsres2[idx]
        dsphew += f"""
MEASURE {r1} SG {r2} SG &d
IF &d > 3 GOTO _nodisulf{idx}
DELATOM {r1} HG
DELATOM {r2} HG
CONECT {r1} SG {r2} SG
_nodisulf{idx}:
"""

    phew_script = f"""LET $rcpid = "{rcpid}"
LET $inpf = "pdbs/{fam}/{rcpid}.inactive.pdb"
LET $mdld = "hm/{best_pdb}"
LOAD $inpf A I
LET %rcpseqln = %SEQLENI
LOAD $mdld A A

BWCOPY I A
STRAND I
UPRIGHT I
BWCENTER

STRAND A
IF "{tplsused}" = "" REMARK 265 HM_TEMPLATES: none
ELSE REMARK 265 HM_TEMPLATES: {tplsused}

HYDRO

UNCHAIN I
UNCHAIN O
STRAND A
UPRIGHT
BWCENTER
{adjustments}
{dsphew}
LET $outf = "pdbs/{fam}/{rcpid}.active.pdb"
SAVE $outf
"""

    phew_path = f"{rcpid}.hm.phew"
    with open(phew_path, "w") as f:
        f.write(phew_script)

    # 9. Adapt output file for AromaDock compatibility
    print("Running orientation and internal coordinates...", file=sys.stderr)
    os.chdir(root_dir)
    subprocess.run(["./bin/phew", f"hm/{phew_path}"])

    # Radial Outward Field: repoint lipid-facing hydrophobic residues in TMR1-TMR7
    active_pdb_rel = f"pdbs/{fam}/{rcpid}.active.pdb"
    pinned_bw = []
    if os.path.exists(active_pdb_rel):
        bw50 = {}
        tmrs = []
        res = {}
        with open(active_pdb_rel, "r") as f:
            for line in f:
                if line.startswith("REMARK 800 SITE BW "):
                    parts = line.split()
                    if len(parts) >= 6:
                        bw_tag = parts[4]
                        rnum = int(parts[5])
                        h = int(bw_tag.split(".")[0])
                        bw50[h] = rnum
                elif line.startswith("REMARK 650 HELIX TMR"):
                    parts = line.split()
                    if len(parts) >= 6:
                        h = int(parts[3][3:])
                        s = int(parts[4])
                        e = int(parts[5])
                        tmrs.append((h, s, e))
                elif line.startswith("ATOM"):
                    try:
                        rnum = int(line[22:26].strip())
                        aname = line[12:16].strip()
                        resn = line[17:20].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        if rnum not in res:
                            res[rnum] = {'name': resn, 'atoms': {}}
                        res[rnum]['atoms'][aname] = (x, y, z)
                    except ValueError:
                        pass

        tip_atoms = \
        {
            'LEU': ['CD1', 'CD2'],
            'ILE': ['CD1'],
            # 'VAL': ['CG1', 'CG2'],
            'MET': ['CE'],
            # 'PHE': ['CZ']
            '-SER': ['OG'],
            '-THR': ['OG1'],
            '-ASN': ['CG'],
            '-GLN': ['CD'],
            '-ASP': ['CG'],
            '-GLU': ['CD'],
        }

        bsrcen = [0.0,0.0,0.0,0.0]
        
        for h, s, e in tmrs:
            for rnum in range(s, e + 1):
                bwpos = 50 + rnum - bw50[h]
                bw = f"{h}.{bwpos}"
                if bw in pu.bsrs:
                    atoms = res[rnum]['atoms']
                    if 'CA' not in atoms:
                        continue
                    ca = atoms['CA']
                    bsrcen[0] += ca[0]
                    bsrcen[1] += ca[1]
                    bsrcen[2] += ca[2]
                    bsrcen[3] += 1

        if bsrcen[3]:
            bsrcen[0] /= bsrcen[3]
            bsrcen[1] /= bsrcen[3]
            bsrcen[2] /= bsrcen[3]
        print(f"bsrcen {bsrcen}")

        repoint_cmds = []
        for h, s, e in tmrs:
            for rnum in range(s, e + 1):
                if rnum not in res:
                    continue
                resn = res[rnum]['name']
                resnm = f"-{resn}"
                point_sign = 0
                if resn in tip_atoms:
                    point_sign = 5.0
                elif resnm in tip_atoms:
                    point_sign = -1.5
                    resn = resnm
                else:
                    continue
                atoms = res[rnum]['atoms']
                if 'CA' not in atoms:
                    continue
                ca = atoms['CA']
                y_ca = ca[1]
                if not (-2.0 <= y_ca <= 20.0):
                    continue
                r_ca = math.sqrt(ca[0]**2 + ca[2]**2)
                maruos = ca[0] - bsrcen[0]
                gdoniobo = ca[2] - bsrcen[2]
                r_ca_p = math.sqrt(maruos**2 + gdoniobo**2)
                if r_ca < 11.0 and r_ca_p >= 13:
                    continue

                tips = [atoms[a] for a in tip_atoms[resn] if a in atoms]
                if not tips:
                    continue
                tip_x = sum(t[0] - bsrcen[0] for t in tips) / len(tips)
                tip_y = sum(t[1] - bsrcen[1] for t in tips) / len(tips)
                tip_z = sum(t[2] - bsrcen[2] for t in tips) / len(tips)

                r_tip = math.sqrt(tip_x**2 + tip_z**2)
                v_side = (tip_x - ca[0], tip_z - ca[2])
                dot = ca[0] * v_side[0] + ca[2] * v_side[1]

                if r_tip < r_ca_p or dot < 0:
                    u_x = ca[0] / r_ca
                    u_z = ca[2] / r_ca
                    tgt_x = ca[0] + point_sign * u_x
                    tgt_y = ca[1]
                    tgt_z = ca[2] + point_sign * u_z
                    repoint_cmds.append(f"ATOMTO {rnum} {tip_atoms[resn][0]} [{tgt_x:.2f},{tgt_y:.2f},{tgt_z:.2f}]")
                    if h in bw50:
                        bw_pos = 50 + rnum - bw50[h]
                        pinned_bw.append(f"{h}.{bw_pos}")
                    else:
                        pinned_bw.append(str(rnum))

        if repoint_cmds:
            radial_phew = f"hm/{rcpid}.radial.phew"
            with open(radial_phew, "w") as f:
                f.write(f"LOAD {active_pdb_rel} A A\n")
                for cmd in repoint_cmds:
                    f.write(f"{cmd}\n")
                f.write(f"SAVE {active_pdb_rel}\n")
            subprocess.run(["./bin/phew", radial_phew])
            if os.path.exists(radial_phew):
                os.remove(radial_phew)

    ic_cmd = ["./bin/ic", f"pdbs/{fam}/{rcpid}.active.pdb", "5.0", "save", "minc"]
    for bw in pinned_bw:
        ic_cmd.extend(["pin", bw])
    subprocess.run(ic_cmd)

    # 10. Clean up temporary files
    os.chdir(script_dir)
    target_pdb = os.path.join(root_dir, "pdbs", fam, f"{rcpid}.active.pdb")
    if os.path.exists(target_pdb) and os.path.getmtime(target_pdb) > os.path.getmtime(phew_path):
        if not nodel:
            print("Cleaning up temporary MODELLER artifacts...", file=sys.stderr)
            for doomed in glob.glob(f"{rcpid}.*"):
                if doomed != f"{rcpid}.active.pdb":
                    os.remove(doomed)
            for doomed in glob.glob(f"{rcpid}_tpl.*"):
                if doomed != f"{rcpid}.active.pdb":
                    os.remove(doomed)

    print("Execution Complete.", file=sys.stderr)

if __name__ == "__main__":
    main()