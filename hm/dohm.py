#!/usr/bin/env python3
"""
Executes homology modeling using MODELLER for a specific GPCR target.
"""

import os
import sys
import subprocess
import glob
import re

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

    # 6. Execute Homology Model
    print(f"Initiating MODELLER for {rcpid}...", file=sys.stderr)
    a = AromaModel(env, alnfile=hm_ali_file, knowns=f'{rcpid}_tpl', sequence=rcpid)
    a.starting_model = 0
    a.ending_model = 9
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
    subprocess.run(["./bin/ic", f"pdbs/{fam}/{rcpid}.active.pdb", "5.0", "save", "minc"])

    # 10. Clean up temporary files
    os.chdir(script_dir)
    target_pdb = os.path.join(root_dir, "pdbs", fam, f"{rcpid}.active.pdb")
    if os.path.exists(target_pdb) and os.path.getmtime(target_pdb) > os.path.getmtime(phew_path):
        if not nodel:
            print("Cleaning up temporary MODELLER artifacts...", file=sys.stderr)
            for doomed in glob.glob(f"{rcpid}.*"):
                if doomed != f"{rcpid}.active.pdb":
                    os.remove(doomed)

    print("Execution Complete.", file=sys.stderr)

if __name__ == "__main__":
    main()