import sys
import re
import os
import os.path
import json
import shutil
import subprocess
from pathlib import Path

script_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(script_dir)
sys.path.insert(0, root_dir)
from modeling_backend import options as backend_options, api as backend_api, construct, forwarded_options, publication_guard, finalize
from reec_backend import receptor, directory

def main(argv=None):
    options = backend_options(sys.argv[1:] if argv is None else argv, repair=True)
    engine = backend_api(options)
    Environ, DOPEHRLoopModel = engine.Environ, engine.DOPEHRLoopModel
    Selection, RigidBody = engine.Selection, engine.RigidBody
    secondary_structure, forms, features, physical = engine.secondary_structure, engine.forms, engine.features, engine.physical
    autosched = engine.autosched
    sys.argv = [sys.argv[0], options.receptor] + options.legacy
    os.chdir(root_dir)
    import data.globals
    import data.protutils
    import data.odorutils
    import data.dyncenter
    from data.disulfides import pdb_disulfides

    data.protutils.load_prots()
    data.odorutils.load_odors()
    protid, mode = receptor(options.receptor)
    fam = data.protutils.family_from_protid(protid)
    origpdb = f"pdbs/{fam}/{protid}.{mode}.pdb"
    output_pdb = options.output or Path(root_dir) / origpdb
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    actions = {'reset', 'fit', 'predock', 'loop', 'dock'}
    if options.output and actions.intersection(options.legacy) & {'dock', 'loop'}:
        raise ValueError('Dock actions use the receptor model path; omit --output')
    if mode == 'inactive' and 'reset' in options.legacy:
        raise ValueError('Reset builds active models; omit reset for inactive repair')
    input_arg = options.legacy[0] if options.legacy and options.legacy[0] not in actions else None
    inppdb = origpdb
    odor = False
    if input_arg:
        candidate = Path(input_arg)
        if not candidate.exists():
            candidate = Path(f"out/{fam}/{protid}/{protid}~{input_arg}.{mode}.model1.pdb")
        odor_name = input_arg.split('~', 1)[-1].split(f'.{mode}.')[0] if '~' in input_arg else input_arg
        odor = data.odorutils.find_odorant(odor_name)
        if candidate.exists():
            inppdb = str(candidate.resolve())
        elif not odor or not actions.intersection(options.legacy) & {'fit', 'predock', 'reset'}:
            raise ValueError(f'Input file not found: {input_arg}')

    aaletts = "ARNDCEQGHILKMFPSTWYV"
    aacode3 = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    # If PDB files get too wonky, you can kill the fixfail.py process and rerun it with the reset argument.
    doreset = False
    if "reset" in sys.argv:
        # Rebuild the unmodified HM structure.
        cmd = [sys.executable, "hm/dohm.py", protid] + forwarded_options(options)
        data.globals.wait_cool_cpu()
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)
        doreset = True

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
        subprocess.run(cmd, check=True)
    else:
        inppdb = origpdb
        odor = False

    delete_inppdb = False

    tokens = [value for value in options.legacy[(1 if input_arg else 0):] if value not in actions]
    if len(tokens) % 2:
        raise ValueError('Ligand contacts use pairs: LIGAND_ATOM RES:ATOM')
    ligcontacts = []
    for ligand_atom, receptor_atom in zip(tokens[::2], tokens[1::2]):
        receptor_atom = re.sub(r'^[A-Za-z]{0,3}(?=[0-9])', '', receptor_atom)
        if not re.fullmatch(r'[0-9]+:[0-9]?[A-Z]+[0-9]?', receptor_atom):
            raise ValueError(f'Invalid receptor contact: {receptor_atom}')
        ligcontacts.append([ligand_atom, receptor_atom])

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
        subprocess.run(cmd, check=True)
        # TODO: If no result from cavity_fit, error out and exit.

    if "predock" in sys.argv or ("reset" in sys.argv and not "fit" in sys.argv):
        # Perform an active-state dock on the new model...
        if odor:
            cmd = ["/bin/bash", "./dock.sh", protid, odor["full_name"], "noi"]
            data.globals.wait_cool_cpu()
            print(" ".join(cmd))
            subprocess.run(cmd, check=True)

    source_disulfides = pdb_disulfides(origpdb)
    canonical = tuple(sorted((data.protutils.resno_from_bw(protid, "3.25") or 0,
                              data.protutils.resno_from_bw(protid, "45.50") or 0)))
    if all(canonical) and all(data.protutils.aalet_at_resno(protid, r) == 'C' for r in canonical):
        if canonical not in source_disulfides:
            source_disulfides.append(canonical)
    disulfide_code = "\n".join(f"CONECT {first} SG {second} SG" for first, second in source_disulfides)

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

        cmd = ["make", "bin/ic", "bin/phew"]
        if ligcontacts:
            cmd.append("bin/poslig")
        data.globals.wait_cool_cpu()
        subprocess.run(cmd, check=True)

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
        if not options.skip_reshape and not os.path.exists(rshpmfn):
            print("No reshape file exists for this receptor.")
            exit()

        if len(ligcontacts):
            cmd = ["bin/poslig", inppdb, tplpdb]
            for lc in ligcontacts:
                cmd.extend(lc)
            print(" ".join(cmd))
            data.globals.wait_cool_cpu()
            subprocess.run(cmd, check=True)
        else:
            shutil.copyfile(inppdb, tplpdb)

        if not options.skip_reshape:
            cmd = ["bin/ic", tplpdb, rshpmfn, "save"]
            print(" ".join(cmd))
            data.globals.wait_cool_cpu()
            subprocess.run(cmd, check=True)

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


        if ligcontacts and not haslig:
            raise ValueError('Ligand contacts require HETATM records')
        env = Environ()
        if haslig:
            env.io.hetatm = True

        class AromaReceptorModel(DOPEHRLoopModel):
            def select_loop_atoms(self):
                # Exclude most of each alpha helix from loop atoms.
                # Allow 4 residues at the beginning and end of each TM helix.
                # Exclude the entire EXR2 helix.
                prot = data.protutils.prots[protid]
                bleed = 4
                EXR2hxstart = data.protutils.resno_from_bw(protid, "45.50")
                EXR2hxend = data.protutils.resno_from_bw(protid, "45.57")
                return Selection(self.residue_range("1:A", f"{prot['region']['TMR1']['start']+bleed}:A"),
                                 self.residue_range(f"{prot['region']['TMR1']['end']-bleed}:A", f"{prot['region']['TMR2']['start']+bleed}:A"),
                                 self.residue_range(f"{prot['region']['TMR2']['end']-bleed}:A", f"{prot['region']['TMR3']['start']+bleed}:A"),
                                 self.residue_range(f"{prot['region']['TMR3']['end']-bleed}:A", f"{prot['region']['TMR4']['start']+bleed}:A"),
                                 self.residue_range(f"{prot['region']['TMR4']['end']-bleed}:A", f"{EXR2hxstart-1}:A"),
                                 self.residue_range(f"{EXR2hxend+1}:A", f"{prot['region']['TMR5']['start']+bleed}:A"),
                                 self.residue_range(f"{prot['region']['TMR5']['end']-bleed}:A", f"{prot['region']['TMR6']['start']+bleed}:A"),
                                 self.residue_range(f"{prot['region']['TMR6']['end']-bleed}:A", f"{prot['region']['TMR7']['start']+bleed}:A"),
                                 )

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

                for first, second in source_disulfides:
                    rsr.add(forms.Gaussian(group=physical.xy_distance,
                        feature=features.Distance(at[f"SG:{first}:A"], at[f"SG:{second}:A"]),
                        mean=2.05, stdev=0.2))

                if haslig:
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
                for first, second in source_disulfides:
                    self.patch(residue_type='DISU', residues=(self.residues[f'{first}:A'],
                                                            self.residues[f'{second}:A']))

        if options.template or options.alignment:
            from copy import copy
            from modeling_backend import custom_alignment
            supplied = copy(options)
            tplfttl = Path(tplpdb).stem
            supplied.template = [f'{tplfttl}={Path(tplpdb).resolve()}'] + options.template
            aligned, templates, _ = custom_alignment(supplied, protid, Path(script_dir))
            tmpalif = str(aligned)
            tplsfttls = tuple(value.split('=', 1)[0] for value in templates)
            os.chdir(script_dir)
        else:
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
            cmd = [sys.executable, "build_alignment_file.py"]
            data.globals.wait_cool_cpu()
            subprocess.run(cmd, check=True)

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

            if options.alignment:
                tmpalif = str(options.alignment)

            os.chdir(script_dir)
            tplsfttls = (tplfttl,)
        env.io.atom_files_directory = ['.', './tpl']

        fix8uy0 = ""
        # t4hm = data.protutils.templates_for_hm(protid)
        # tpls += " " + " ".join(t4hm)
        # tpls = " ".join(list(set(tpls.split(" "))))
        # print(legal)
        # print(tpls)
        if "8uy0" in tplsfttls: fix8uy0 = "BRIDGE DE45.51 Y6.55"
        print(f"Templates: {tplsfttls}")
        a = construct(AromaReceptorModel, env,
                                references        = (os.path.join(root_dir, origpdb),),
                                alnfile           = tmpalif,
                                knowns            = tplsfttls,
                                sequence          = protid
                            )
        a.starting_model = 0
        a.ending_model   = a.starting_model + options.models - 1 if options.models is not None else 9
        a.library_schedule = autosched.slow
        # REEC performs three minimization stages; its CLI default is 3000 total steps.
        a.max_var_iterations = options.iterations if options.iterations is not None else (3000 if options.backend == 'reec' else 300)

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

        with publication_guard(output_pdb):
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

        {disulfide_code}
        HYDRO

        LET %y6 = 0
        LET $atom6 = "OH"
        IF $6.55 = "Y" THEN LET %y6 = %6.55
        IF $6.55 = "H" THEN LET %y6 = %6.55
        IF $6.55 = "H" THEN LET $atom6 = "NE2"
        LET %de45 = 0
        IF $45.51 = "D" OR $45.51 = "E" THEN LET %de45 = %45.51
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
        IF $5.46 != "M" GOTO _not_Cu_546
        MEASURE %5.46 "SD" %5.42 "SG" &sdist
        ECHO "5.46:SD - 5.42:SG distance: " &sdist
        MEASURE %5.46 "SD" %5.43 "SG" &sdist
        ECHO "5.46:SD - 5.43:SG distance: " &sdist
        _not_Cu_546:
        MEASURE %5.42 "SG" %5.43 "SG" &sdist
        ECHO "5.42:SD - 5.43:SG distance: " &sdist
        _not_Cu_binding_site:

        LET $outf = "{output_pdb}"
        SAVE $outf
        """

            os.chdir(os.path.dirname(os.path.abspath(__file__)))
            os.chdir("..")
            phewfn = f"hm/{protid}.{mode}.hm.phew"
            with open(phewfn, "w") as f:
                # PHEW resolves labels against whole lines, including whitespace.
                f.write('\n'.join(line.strip() for line in phewcode.splitlines()) + '\n')
            cmd = ["bin/phew", phewfn]
            data.globals.wait_cool_cpu()
            subprocess.run(cmd, check=True)

            cmd = ["bin/ic", str(output_pdb), "minc", "save"]
            for resno in sorted({resno for pair in source_disulfides for resno in pair}):
                cmd.extend(["pin", str(resno)])
            print(" ".join(cmd))
            data.globals.wait_cool_cpu()
            subprocess.run(cmd, check=True)

            with directory(script_dir):
                finalize(a, output_pdb, model['name'])

        if not "loop" in sys.argv and not "dock" in sys.argv:
            break
        else: tries += 1

        if odor:
            cmd = ["/bin/bash", "./dock.sh", protid, odor["full_name"], "noi"]
            data.globals.wait_cool_cpu()
            subprocess.run(cmd, check=True)

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

    if delete_inppdb and os.path.exists(inppdb):
        os.remove(inppdb)


if __name__ == '__main__':
    main()
