#!/usr/bin/env python3
"""
Builds GPCR homology models with REEC by default; MODELLER is optional.
Rejected generations retry with fresh seeds and increasing relaxation budgets.
After three attempts, batch refreshes emit an explicitly labeled fallback model.
"""

import os
import sys
import subprocess
import glob
import re
import math
import shutil
import json
import hashlib
import tempfile

# -----------------------------------------
# ARCHITECTURAL ENFORCEMENT: PATH RESOLUTION
# -----------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(script_dir)
sys.path.append(root_dir)

from pathlib import Path
from reec_reports import report_path
from reec_backend import receptor, directory, CandidateRejected
from modeling_backend import options as backend_options, api as backend_api, prepare, construct, custom_alignment, publication_guard, finalize, pir_rows

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

def experimental_template(name, source, description, work, chain='A'):
    """Extract the PIR chain into a template with the engine's chain label."""
    fields = description.split(':')
    source_chain = fields[3].strip() or ' '
    if len(source_chain) != 1 or fields[5].strip() not in ('', source_chain.strip()):
        raise ValueError(f'Experimental template must describe one chain: {name}')
    records = []
    for line in source.read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith(('ATOM  ', 'HETATM', 'TER   ')) and len(line) > 21 and line[21] == source_chain:
            records.append(line[:21] + chain + line[22:])
    if not any(line.startswith('ATOM  ') and line[12:16].strip() == 'CA' for line in records):
        raise ValueError(f'Experimental template has no protein CA atoms in its PIR chain: {name}')
    template = work / f'{name}.pdb'
    template.write_text('\n'.join(records) + '\nEND\n')
    fields[3] = fields[5] = chain
    return template, ':'.join(fields)


def validated_subfamily_alignment(rcpid, work):
    """Combine trusted subfamily models with active experimental OR templates."""
    family = pu.family_from_protid(rcpid)
    subfamily = pu.subfamily_from_protid(rcpid)
    sources = []
    for source in sorted((Path(root_dir) / 'pdbs' / '_validated').glob('*.active.pdb')):
        name = source.name.removesuffix('.active.pdb')
        if (name in pu.prots and name != rcpid
                and pu.family_from_protid(name) == family
                and pu.subfamily_from_protid(name) == subfamily):
            sources.append((name, source))
    if not sources:
        return None

    rows = {}
    templates = []
    provenance = []
    codes = dict(zip(
        'ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split(),
        'ACDEFGHIKLMNPQRSTVWY'))
    for name, source in sources:
        sequence = pu.prots[name]['sequence']
        aligned = ''.join(pu.prots[name]['aligned'].split())
        if aligned.replace('-', '') != sequence:
            raise ValueError(f'GPCR alignment does not match the sequence for {name}')
        present = set()
        for line in source.read_text().splitlines():
            if not line.startswith('ATOM  ') or line[12:16].strip() != 'CA' or line[21] not in (' ', 'A'):
                continue
            number = int(line[22:26])
            if (line[26].strip() or number in present or not 1 <= number <= len(sequence)
                    or codes.get(line[17:20]) != sequence[number - 1]):
                raise ValueError(f'Validated {name} residue numbering does not match its sequence')
            present.add(number)
        if not present:
            raise ValueError(f'Validated {name} has no CA atoms in chain A')
        position = 0
        row = ''
        for amino in aligned:
            if amino != '-':
                position += 1
            row += amino if amino == '-' or position in present else '-'
        template_id = f'{rcpid}_validated_{name}'
        template = work / f'{template_id}.pdb'
        shutil.copyfile(source, template)
        rows[template_id] = (f'structureX:{template_id}:FIRST:A:LAST:A::::', row)
        templates.append(f'{template_id}={template}')
        provenance.append(str(source.relative_to(root_dir)))

    experimental = pir_rows(work / 'experimental.ali')
    experimental_ids = set(pu.templates_for_hm(rcpid)) | {
        '8uxv', '8hti', '9wg4', '8uxy', '8uy0', '8uyq', '9wpm',
        '9ldv', '9ldw', '9ldx', '9ldz'}
    for name, (description, row) in experimental.items():
        if name not in experimental_ids or any(tag in description.lower() for tag in ('inactive', 'apo ')):
            continue
        template = work / 'tpl' / f'{name}.pdb'
        if not template.is_file():
            raise FileNotFoundError(f'Experimental template is missing: {template}')
        template, description = experimental_template(name, template, description, work)
        rows[name] = (description, row)
        templates.append(f'{name}={template}')
        provenance.append(name)
    target_row = ''.join(pu.prots[rcpid]['aligned'].split())
    if target_row.replace('-', '') != pu.prots[rcpid]['sequence']:
        raise ValueError(f'GPCR alignment does not match the sequence for {rcpid}')
    rows[rcpid] = (f'sequence:{rcpid}:1:A::::::', target_row)
    width = max(len(row) for _, row in rows.values())
    alignment = work / f'{rcpid}.validated.ali'
    alignment.write_text(''.join(
        f'>P1;{name}\n{description}\n{row.ljust(width, "-")}*\n'
        for name, (description, row) in rows.items()))
    return alignment, templates, ' '.join(provenance)

def primary_experimental_alignment(rcpid, work):
    """Deliberate OR5V1 exception: retain 9lkb and use consOR5 for gaps."""
    if rcpid != 'OR5V1':
        return None
    return experimental_ensemble_alignment(rcpid, work)


def experimental_ensemble_alignment(rcpid, work):
    """Keep experimental donors separate for native similarity weighting."""
    rows = pir_rows(work / 'experimental.ali')
    if rcpid == 'OR5V1':
        identifiers = ['9lkb', '9wpm']
    else:
        candidates = set(pu.templates_for_hm(rcpid))
        if rcpid.startswith('OR'):
            family = int(re.sub(r'[^0-9]', '', pu.family_from_protid(rcpid)))
            candidates.update(pu.templates_for_hm('OR56B2' if family >= 51 else 'OR14A1'))
        identifiers = sorted(candidates)
    templates = []
    selected = {}
    for name in identifiers:
        if name not in rows:
            if rcpid == 'OR5V1':
                raise ValueError(f'Experimental alignment is missing: {name}')
            continue
        description, row = rows[name]
        if any(tag in description.lower() for tag in ('inactive', 'apo ')):
            continue
        source = work / 'tpl' / f'{name}.pdb'
        if not source.is_file():
            if rcpid == 'OR5V1':
                raise FileNotFoundError(f'Experimental template is missing: {source}')
            print(f'Skipping unavailable experimental template: {name}', file=sys.stderr)
            continue
        template, description = experimental_template(name, source, description, work)
        selected[name] = (description, row)
        templates.append(f'{name}={template}')
    if not templates:
        raise ValueError(f'No active experimental templates are available for {rcpid}')
    row = ''.join(pu.prots[rcpid]['aligned'].split())
    sequence = pu.prots[rcpid]['sequence']
    aligned_sequence = row.replace('-', '')
    if not sequence.startswith(aligned_sequence):
        raise ValueError(f'GPCR alignment does not match the sequence for {rcpid}')
    # Preserve canonical columns when the database alignment omits a terminal tail.
    row += sequence[len(aligned_sequence):]
    selected[rcpid] = (f'sequence:{rcpid}:1:A::::::', row)
    width = max(len(row) for _, row in selected.values())
    alignment = work / f'{rcpid}.experimental.ali'
    alignment.write_text(''.join(f'>P1;{name}\n{description}\n{row.ljust(width, "-")}*\n'
                                 for name, (description, row) in selected.items()))
    return alignment, templates, ' '.join(value.split('=', 1)[0] for value in templates)


def process_model(rcpid, fam, famno, famsub, output_pdb, best_pdb, tplsused,
                  dsres1, dsres2, dspotr1, dspotr2, model=None, coordinate_metals=False):
    """Publish modeled and validated structures through the same PHEW workflow."""
    model_path = Path(best_pdb)
    if not model_path.is_absolute():
        model_path = Path(script_dir) / model_path
    output_pdb = Path(output_pdb).resolve()
    pu.ensure_bw_metadata(rcpid, Path(root_dir) / f'pdbs/{fam}/{rcpid}.inactive.pdb')
    os.chdir(script_dir)
    with publication_guard(output_pdb):
        # 8. Post-Processing Script Generation (.phew)
        adjustments = ""
        if famno < 50:
            adjustments += 'IF $3.37 != "G" THEN ATOMTO %3.37 EXTENT @6.48\n'
            adjustments += 'IF $45.51 = "D" AND $6.55 = "Y" THEN BRIDGE %45.51 %6.55\n'
            adjustments += 'IF $45.51 = "E" AND $6.55 = "Y" THEN BRIDGE %45.51 %6.55\n'
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
    LET $mdld = "{model_path}"
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
    LET $outf = "{output_pdb}"
    SAVE $outf
    """

        phew_path = str(Path(script_dir) / f"{rcpid}.hm.phew")
        with open(phew_path, "w") as f:
            f.write(phew_script)

        # 9. Adapt output file for AromaDock compatibility
        print("Running orientation and internal coordinates...", file=sys.stderr)
        os.chdir(root_dir)
        subprocess.run(["./bin/phew", phew_path], check=True)

        # Radial Outward/Inward Field: flex pocket aliphatics out and hydrophilics in
        active_pdb_rel = str(output_pdb)
        pinned_bw = []
        if os.path.exists(active_pdb_rel):
            bw50 = {}
            tmrs = []
            res = {}
            atom_list = []
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
                                res[rnum] = \
                                {
                                    'name': resn,
                                    'atoms': {},
                                    'h': None,
                                    'bw': None,
                                }
                            res[rnum]['atoms'][aname] = (x, y, z)
                            atom_list.append((x, y, z, rnum, aname))
                        except ValueError:
                            pass

            for h, s, e in tmrs:
                for rnum in range(s, e + 1):
                    if rnum in res:
                        res[rnum]['h'] = h
                        if h in bw50:
                            bwpos = 50 + rnum - bw50[h]
                            res[rnum]['bw'] = f"{h}.{bwpos}"

            tip_atoms = \
            {
                'LEU': ['CD1', 'CD2'],
                'ILE': ['CD1'],
                'VAL': ['CG1', 'CG2'],
                'MET': ['CE'],
                '-SER': ['OG'],
                '-THR': ['OG1'],
                '-ASN': ['CG'],
                '-GLN': ['CD'],
                '-ASP': ['CG'],
                '-GLU': ['CD'],
            }

            # Identify BSRs within Y coordinate range
            bsr_rnums = []
            for h, s, e in tmrs:
                for rnum in range(s, e + 1):
                    bw = res.get(rnum, {}).get('bw')
                    if bw in pu.bsrs:
                        if rnum in res and 'CA' in res[rnum]['atoms']:
                            ca = res[rnum]['atoms']['CA']
                            if -2.0 <= ca[1] <= 20.0:
                                bsr_rnums.append(rnum)

            bsrcen = [0.0, 0.0, 0.0]
            for rnum in bsr_rnums:
                ca = res[rnum]['atoms']['CA']
                bsrcen[0] += ca[0]
                bsrcen[1] += ca[1]
                bsrcen[2] += ca[2]

            if bsr_rnums:
                bsrcen[0] /= len(bsr_rnums)
                bsrcen[1] /= len(bsr_rnums)
                bsrcen[2] /= len(bsr_rnums)
            print(f"bsrcen {bsrcen}")

            def dist_pt_to_seg(p, a, b):
                dx = b[0] - a[0]
                dy = b[1] - a[1]
                dz = b[2] - a[2]
                l2 = dx * dx + dy * dy + dz * dz
                if l2 == 0:
                    return math.sqrt((p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2 + (p[2] - a[2]) ** 2)
                t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy + (p[2] - a[2]) * dz) / l2))
                proj_x = a[0] + t * dx
                proj_y = a[1] + t * dy
                proj_z = a[2] + t * dz
                return math.sqrt((p[0] - proj_x) ** 2 + (p[1] - proj_y) ** 2 + (p[2] - proj_z) ** 2)

            sub_atoms = [a for a in atom_list if -5.0 <= a[1] <= 25.0]

            repoint_cmds = []
            for h, s, e in tmrs:
                for rnum in range(s, e + 1):
                    if rnum not in res:
                        continue
                    if rnum in dspotr1 or rnum in dspotr2:
                        continue
                    resn = res[rnum]['name']
                    resnm = f"-{resn}"
                    point_sign = 0
                    tip_atom_name = None
                    if resn in tip_atoms:
                        point_sign = 5.0
                        tip_atom_name = tip_atoms[resn][0]
                    elif resnm in tip_atoms:
                        point_sign = -5.0
                        tip_atom_name = tip_atoms[resnm][0]
                        resn = resnm
                    else:
                        continue

                    atoms = res[rnum]['atoms']
                    if 'CA' not in atoms or tip_atom_name not in atoms:
                        continue
                    ca = atoms['CA']
                    y_ca = ca[1]
                    if not (-2.0 <= y_ca <= 20.0):
                        continue

                    is_bsr = (rnum in bsr_rnums)
                    has_los = False
                    if not is_bsr:
                        for br in bsr_rnums:
                            if res[rnum]['h'] == res[br]['h']:
                                continue
                            ca_b = res[br]['atoms']['CA']
                            obstructed = False
                            for ax, ay, az, arnum, aname in sub_atoms:
                                if arnum == rnum or arnum == br:
                                    continue
                                if dist_pt_to_seg((ax, ay, az), ca, ca_b) < 1.5:
                                    obstructed = True
                                    break
                            if not obstructed:
                                if dist_pt_to_seg(bsrcen, ca, ca_b) < 10.0:
                                    has_los = True
                                    break

                    if not (is_bsr or has_los):
                        continue

                    dx = ca[0] - bsrcen[0]
                    dz = ca[2] - bsrcen[2]
                    r_ca_p = math.sqrt(dx * dx + dz * dz)
                    if r_ca_p < 0.1:
                        continue
                    ux = dx / r_ca_p
                    uz = dz / r_ca_p

                    if point_sign < 0:
                        tgt_x = bsrcen[0]
                        tgt_y = bsrcen[1]
                        tgt_z = bsrcen[2]
                    else:
                        tgt_x = ca[0] + point_sign * ux
                        tgt_y = ca[1]
                        tgt_z = ca[2] + point_sign * uz

                    repoint_cmds.append(f"ATOMTO {rnum} {tip_atom_name} [{tgt_x:.2f},{tgt_y:.2f},{tgt_z:.2f}]")
                    if res[rnum]['bw']:
                        pinned_bw.append(res[rnum]['bw'])
                    else:
                        pinned_bw.append(str(rnum))

            if repoint_cmds:
                radial_phew = f"hm/{rcpid}.radial.phew"
                with open(radial_phew, "w") as f:
                    f.write(f"LOAD {active_pdb_rel} A A\n")
                    for cmd in repoint_cmds:
                        f.write(f"{cmd}\n")
                    f.write(f"SAVE {active_pdb_rel}\n")
                subprocess.run(["./bin/phew", radial_phew], check=True)

        ic_cmd = ["./bin/ic", str(output_pdb), "5.0", "save", "minc"]
        for bw in pinned_bw:
            ic_cmd.extend(["pin", bw])
        for resno in sorted(set(dsres1 + dsres2)):
            ic_cmd.extend(["pin", str(resno)])
        subprocess.run(ic_cmd, check=True)

        metal_review = []
        if coordinate_metals:
            from data.metal_sites import configured_metal_sites, inspect_metal_sites
            sites = configured_metal_sites(rcpid, Path(root_dir) / 'data/binding_pocket.json')
            if sites:
                metal_script = Path(script_dir) / f'{rcpid}.metal.phew'
                commands = [f'LOAD "{output_pdb}" A A', 'STRAND A']
                commands.extend(f'MCOORD {site["element"]} {site["charge"]} '
                                + ' '.join(map(str, site['residues'])) for site in sites)
                commands.append(f'SAVE "{output_pdb}"')
                metal_script.write_text('\n'.join(commands)+'\n')
                try:
                    subprocess.run([str(Path(root_dir) / 'bin/phew'), str(metal_script)], check=True)
                    metal_review = inspect_metal_sites(output_pdb, sites)
                except (subprocess.CalledProcessError, ValueError) as exc:
                    raise CandidateRejected(f'Metal-site placement failed: {exc}') from exc

        with directory(script_dir):
            finalize(model, output_pdb, best_pdb, target=rcpid)
        if metal_review:
            report_path = Path(str(output_pdb)+'.json')
            report = json.loads(report_path.read_text())
            report['metal_sites'] = metal_review
            report['metal_placement'] = 'explicit opt-in after final IC processing'
            report_path.write_text(json.dumps(report, indent=2)+'\n')

    return phew_path

def _main(options):
    unknown = set(options.legacy) - {'nodel', 'riglig'}
    if unknown:
        raise ValueError(f'Unknown options: {sorted(unknown)}')
    prepare(options)
    sys.argv = [sys.argv[0], options.receptor] + options.legacy
    if len(sys.argv) < 2:
        print("Usage:\npython3 hm/dohm.py PROTID [riglig] [nodel]\n", file=sys.stderr)
        sys.exit(1)

    rcpid, state = receptor(options.receptor)
    if state != 'active':
        raise ValueError('dohm produces active models; use fixfail for inactive repair')
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

    output_pdb = (options.output or Path(root_dir) / 'pdbs' / fam / f'{rcpid}.active.pdb').resolve()
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    trusted_model = Path(root_dir) / 'pdbs' / '_validated' / f'{rcpid}.active.pdb'
    if trusted_model.is_file():
        pu.ensure_bw_metadata(rcpid, trusted_model)
        from data.disulfides import pdb_disulfides
        from dohm_recovery import review_source, repair_existing, copy_pair
        if trusted_model.resolve() == output_pdb.resolve():
            raise ValueError('Processed output must differ from the validated source')
        pairs = pdb_disulfides(trusted_model)
        dsres1 = [first for first, _ in pairs]
        dsres2 = [second for _, second in pairs]
        (Path(root_dir) / 'tmp').mkdir(exist_ok=True)
        review_work = Path(tempfile.mkdtemp(prefix=f'{rcpid}-validated-review-', dir=Path(root_dir) / 'tmp'))
        review = review_source(trusted_model, rcpid, len(p['sequence']), review_work)
        trusted_status_path = report_path(trusted_model)
        try:
            trusted_status = json.loads(trusted_status_path.read_text()) if trusted_status_path.is_file() else {}
        except ValueError:
            trusted_status = {}
        if (review['valid'] and trusted_status.get('validated_ready')
                and trusted_status.get('validated_sha256') == hashlib.sha256(trusted_model.read_bytes()).hexdigest()
                and not options.coordinate_metals):
            with publication_guard(output_pdb):
                copy_pair(trusted_model, output_pdb)
            print(f'Geometry-validated model copied to {output_pdb}; source conformation preserved.', file=sys.stderr)
            print('Execution Complete.', file=sys.stderr)
            return
        try:
            if not review['valid']:
                raise CandidateRejected(f'Validated source failed whole-protein geometry review: {review["failures"]}')
            print(f'Processing validated model for {rcpid} through PHEW...', file=sys.stderr)
            with directory(script_dir):
                phew_path = process_model(
                    rcpid, fam, famno, famsub, output_pdb, trusted_model,
                    str(trusted_model.relative_to(root_dir)), dsres1, dsres2, dsres1, dsres2,
                    coordinate_metals=options.coordinate_metals)
                if not nodel:
                    Path(phew_path).unlink(missing_ok=True)
                    (Path(script_dir) / f'{rcpid}.radial.phew').unlink(missing_ok=True)
        except CandidateRejected as exc:
            print(f'{exc}; repairing with REEC before publication.', file=sys.stderr, flush=True)
            repaired, repair_work = repair_existing(rcpid, trusted_model, options, root_dir)
            with publication_guard(output_pdb), publication_guard(trusted_model):
                copy_pair(repaired, trusted_model)
                copy_pair(repaired, output_pdb)
            print(f'Geometry-repaired model published and validated source updated. Work: {repair_work}', file=sys.stderr)
        print('Execution Complete.', file=sys.stderr)
        return

    engine = backend_api(options)
    Environ, AutoModel = engine.Environ, engine.AutoModel
    secondary_structure, forms, features, physical = engine.secondary_structure, engine.forms, engine.features, engine.physical
    autosched, refine = engine.autosched, engine.refine
    model_knowns = f'{rcpid}_tpl'
    validated = (custom_alignment(options, rcpid, Path(script_dir)) if options.template
                 else (primary_experimental_alignment(rcpid, Path(script_dir))
                       or validated_subfamily_alignment(rcpid, Path(script_dir))
                       or (experimental_ensemble_alignment(rcpid, Path(script_dir))
                           if options.backend == 'reec' else None)))
    if validated:
        alignment, templates, tplsused = validated
        model_knowns = tuple(value.split('=', 1)[0] for value in templates)
        hm_ali_file = str(alignment)
        os.chdir(script_dir)
    else:
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

    # 5. Define model callbacks shared by both engines
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

            # 5-7 tyrosine link
            bw5_58 = data.protutils.resno_from_bw(rcpid, "5.58")
            bw7_53 = data.protutils.resno_from_bw(rcpid, "7.53")
            if bw5_58 and bw7_53:
                if data.protutils.aalet_at_resno(rcpid, bw5_58) == 'Y' and data.protutils.aalet_at_resno(rcpid, bw7_53) == 'Y':
                    rsr.add(forms.Gaussian(group=physical.xy_distance,
                        feature=features.Distance(at["OH:"+str(bw5_58)+":A"],
                                                at["OH:"+str(bw7_53)+":A"]),
                                                mean=4.6, stdev=1.2))

            # 6-45 hydrogen bond
            bw45_50 = data.protutils.resno_from_bw(rcpid, "45.50")
            bw45_51 = bw45_50+1
            bw6_55 = data.protutils.resno_from_bw(rcpid, "6.55")
            print(f"Attempting Gaussian {bw45_51} ~ {bw6_55}")
            if bw6_55:
                atom6_55 = False
                if data.protutils.aalet_at_resno(rcpid, bw6_55) == 'Y':
                    atom6_55 = "OH"
                elif data.protutils.aalet_at_resno(rcpid, bw6_55) == 'H':
                    atom6_55 = "NE2"
                if atom6_55:
                    atom45_51 = False
                    atom45_51c = False
                    if data.protutils.aalet_at_resno(rcpid, bw45_51) == 'D':
                        atom45_51 = "OD1"
                        atom45_51c = "CG"
                    elif data.protutils.aalet_at_resno(rcpid, bw45_51) == 'E':
                        atom45_51 = "OE1"
                        atom45_51c = "CD"
                    elif data.protutils.aalet_at_resno(rcpid, bw45_51) == 'H':
                        atom45_51 = "NE2"
                        atom45_51c = "CD"
                    elif data.protutils.aalet_at_resno(rcpid, bw45_51) == 'N':
                        atom45_51 = "OD1"
                        atom45_51c = "CG"
                    elif data.protutils.aalet_at_resno(rcpid, bw45_51) == 'Q':
                        atom45_51 = "OE1"
                        atom45_51c = "CD"
                    if atom45_51:
                        print(f"Gaussian {bw45_51}:{atom45_51} ~ {bw6_55}:{atom6_55}")
                        rsr.add(forms.Gaussian(group=physical.xy_distance,
                            feature=features.Distance(at[atom6_55 +":"+str(bw6_55)+":A"],
                                                    at[atom45_51+":"+str(bw45_51)+":A"]),
                                                    mean=2.5, stdev=0.25))
                        rsr.add(forms.Gaussian(group=physical.xy_distance,
                            feature=features.Distance(at[atom6_55 +":"+str(bw6_55)+":A"],
                                                    at[atom45_51c+":"+str(bw45_51)+":A"]),
                                                    mean=3.0, stdev=0.5))

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

    # Retry complete generation and publication while preserving validation limits.
    base_iterations = options.iterations if options.iterations is not None else (3000 if options.backend == 'reec' else 1000)
    base_seed = options.seed
    base_annealing = options.anneal_steps
    attempt = 0
    last_rejection = None
    log_root = Path(root_dir) / 'tmp'
    log_root.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix=f'{rcpid}-dohm-', dir=log_root))
    attempts = []

    def record(entry, rejection=None):
        if rejection is not None:
            entry['reason'] = str(rejection)
            if rejection.report is not None and rejection.report.is_file():
                report = run_dir / f'attempt-{attempt}-event-{len(attempts)+1}.json'
                shutil.copyfile(rejection.report, report)
                entry['report'] = str(report)
        attempts.append(entry)
        (run_dir / 'attempts.json').write_text(json.dumps(attempts, indent=2) + '\n')

    while options.max_attempts == 0 or attempt < options.max_attempts:
        attempt += 1
        options.seed = (base_seed + attempt - 1) % 2147483648
        iterations = min(100000, base_iterations * 2**min(attempt-1, 16))
        options.clash_scale = min(8., 2.**min(attempt-1, 3))
        options.anneal_steps = (0 if attempt % 4 == 0 else
                               min(max(base_annealing, 512), base_annealing * (1+(attempt-1)//2)))
        settings = {'attempt': attempt, 'seed': options.seed, 'iterations': iterations,
                    'anneal_steps': options.anneal_steps, 'clash_scale': options.clash_scale}
        print(f'Generation attempt {attempt}: seed {options.seed}, {iterations} iterations, '
              f'clash scale {options.clash_scale}, {options.anneal_steps} annealing rounds',
              file=sys.stderr, flush=True)
        os.chdir(script_dir)
        dsres1.clear()
        dsres2.clear()
        a = construct(AromaModel, env, alnfile=hm_ali_file, knowns=model_knowns, sequence=rcpid)
        a.starting_model = 0
        a.ending_model = options.models - 1 if options.models is not None else 9
        a.library_schedule = autosched.slow
        a.max_var_iterations = iterations
        if rcpid in ['OR2AE1', 'OR2AG1', 'OR2AG2']:
            a.md_level = refine.very_slow
        if 'data.globals' in sys.modules and hasattr(sys.modules['data.globals'], 'wait_cool_cpu'):
            try:
                sys.modules['data.globals'].wait_cool_cpu()
            except FileNotFoundError:
                pass
        try:
            a.make()
        except CandidateRejected as exc:
            last_rejection = exc
            record(dict(settings, status='generation_rejected'), exc)
            print(f'Attempt {attempt} rejected: {exc}; increasing relaxation and changing seed.',
                  file=sys.stderr, flush=True)
            continue
        ok_models = sorted((m for m in a.outputs if m['failure'] is None), key=lambda m: m['molpdf'])
        if not ok_models:
            last_rejection = CandidateRejected('Modeling produced no valid candidates')
            record(dict(settings, status='generation_rejected'), last_rejection)
            continue
        published = False
        for best_model in ok_models:
            best_pdb = best_model['name']
            print(f'Processing candidate {best_pdb} with molpdf {best_model["molpdf"]}')
            os.chdir(script_dir)
            try:
                phew_path = process_model(rcpid, fam, famno, famsub, output_pdb, best_pdb,
                                         tplsused, dsres1, dsres2, dspotr1, dspotr2, a,
                                         coordinate_metals=options.coordinate_metals)
            except CandidateRejected as exc:
                last_rejection = exc
                options.last_rejected_model = Path(script_dir) / best_pdb
                record(dict(settings, status='processing_rejected', candidate=best_pdb), exc)
                print(f'Candidate {best_pdb} rejected after processing: {exc}', file=sys.stderr, flush=True)
                continue
            record(dict(settings, status='accepted', candidate=best_pdb))
            published = True
            break
        if published:
            print(f'Validated model published. Attempt log: {run_dir / "attempts.json"}', file=sys.stderr)
            break
    else:
        raise CandidateRejected(f'No accepted model after {attempt} generation attempts; '
                                f'attempt log: {run_dir / "attempts.json"}',
                                last_rejection.report if last_rejection is not None else None)

    # 10. Clean up temporary files
    os.chdir(script_dir)
    target_pdb = str(output_pdb)
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

def main(argv=None):
    options = backend_options(sys.argv[1:] if argv is None else argv)
    options.retain_rejected = True
    try:
        return _main(options)
    except (CandidateRejected, subprocess.CalledProcessError, OSError) as exc:
        from dohm_recovery import fallback_model
        rcpid, state = receptor(options.receptor)
        pu.load_prots()
        if rcpid not in pu.prots or state != 'active':
            raise
        family = pu.family_from_protid(rcpid)
        destination = (options.output or Path(root_dir) / f'pdbs/{family}/{rcpid}.active.pdb').resolve()
        if destination.exists():
            print(f'Model refresh rejected; preserving previous model: {destination}', file=sys.stderr)
            raise
        print(f'Model refresh could not pass: {exc}; emitting a labeled fallback.', file=sys.stderr)
        fallback_model(rcpid, destination, options, exc, root_dir)


if __name__ == "__main__":
    try:
        main()
    except CandidateRejected as exc:
        print(f'Model rejected: {exc}', file=sys.stderr)
        if exc.report is not None:
            print(f'Validation report: {exc.report}', file=sys.stderr)
        sys.exit(2)
