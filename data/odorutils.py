import sys
import json
import re
import os
import subprocess
import traceback
import tempfile

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("ERROR: RDKit not found in environment.", file=sys.stderr)
    exit()

odors = {}

def load_odors():
    """Extracts the odorant database and binds it to global memory."""
    global odors
    
    # 1. Determine the exact, absolute directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(current_dir, 'odorant.json')
    
    # 2. Validate existence
    if not os.path.exists(db_path):
        print(f"FATAL ERROR: Odorant database not found at {db_path}", file=sys.stderr)
        return
        
    # 3. Extract and bind globally
    try:
        with open(db_path, 'r', encoding='utf-8') as file:
            odor_data = json.load(file)
            
            # Lock the data into the global dictionary. No flattening.
            odors.clear()
            odors.update(odor_data)
            
    except Exception as e:
        print(f"FATAL ERROR: Failed to extract data from {db_path} - {e}", file=sys.stderr)

def ordered_ligand_ids(rcpid):
    """Use receptor-page order first, then retain unlisted database entries."""
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'receptor_ligand_order.php')
    result = subprocess.run(['php', helper, rcpid], check=True, capture_output=True, text=True)
    ordered = [oid for oid in json.loads(result.stdout) if oid in odors]
    return list(dict.fromkeys([*ordered, *odors]))


def empirical_pairs(rcpid, onedim=False, agonists_only=False):
    global odors
    result = {}
    for oid in odors.keys():
        o = odors[oid]
        if not o.get("activity"): continue
        for url in o["activity"]:
            acv = o["activity"][url]
            for arcp in acv.keys():
                if arcp != rcpid: continue
                data = acv[arcp]
                if agonists_only:
                    if data.get("type") in ("na", "ia"):
                        continue
                    if data.get("adjusted_curve_top") is not None and float(data["adjusted_curve_top"]) <= 0:
                        continue
                if onedim:
                    if data.get("type") == "na":
                        result[oid] = 0
                    elif data.get("type") == "ia":
                        result[oid] = 0
                    elif data.get("adjusted_curve_top") is not None and float(data["adjusted_curve_top"]) <= 0:
                        result[oid] = 0
                    elif data.get("adjusted_curve_top") is not None and float(data.get("adjusted_curve_top")) > 0:
                        result[oid] = 1
                    elif data.get("ec50") is not None and float(data.get("ec50")) < 0:
                        result[oid] = 1
                    else:
                        result[oid] = 0
                else:
                    result[oid] = data
    return result

def find_odorant(aroma):
    global odors
    if not aroma: return False

    if aroma in odors.keys():
        retval = odors[aroma]
        retval['oid'] = aroma
        return retval

    def normalize(name):
        return re.sub(r"[^a-z0-9+_-]", "", re.sub(r"[\s_]+", "_", name.strip().lower()))

    aroma1 = normalize(aroma)
    matches = []
    for oid, o in odors.items():
        names = [o['full_name']]
        if o.get('iupac'):
            names.append(o['iupac'])
        i = 1
        while "name"+str(i) in o:
            names.append(o["name"+str(i)])
            i += 1
        if o['smiles'] == aroma or any(normalize(name) == aroma1 for name in names):
            matches = [(oid, o)]
            break
    if not matches:
        compact = re.sub("[^a-z0-9+-]", "", aroma.lower())
        matches = [(oid, o) for oid, o in odors.items()
                   if re.sub("[^a-z0-9+-]", "", o['full_name'].lower()) == compact]
    if len(matches) != 1:
        return False
    oid, odor = matches[0]
    retval = dict(odor, oid=oid)
    fname = "sdf/" + odor['full_name'].replace(' ', '_') + ".sdf"
    if os.path.exists(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), fname)):
        retval["sdfname"] = fname
    return retval

def check_isomers(ligname, randomize=True):
    odor = find_odorant(ligname)
    if not odor: raise Exception("Odorant not found "+ligname)
    if not "isomers" in odor: return False

    result = []
    for iso in odor["isomers"].keys():
        if "preiso" in odor.keys():
            l = len(odor['preiso'])
            result.append((odor["full_name"][0:l] + iso + "-" + odor["full_name"][l:]).replace(' ', '_'))
        else:
            result.append(iso+"-"+odor["full_name"])
    return result

def check_forms(ligname, randomize=True):
    odor = find_odorant(ligname)
    if not odor: raise Exception("Odorant not found "+ligname)
    if not "forms" in odor: return False

    result = []
    for form in odor["forms"].keys():
        if "preform" in odor.keys():
            l = len(odor['preform'])
            result.append((odor["full_name"][0:l] + form + "-" + odor["full_name"][l:]).replace(' ', '_'))
        else:
            result.append(form+"-"+odor["full_name"])
    return result

def structure_molecule(smiles):
    """Read catalog atom labels and legacy branch-bond spelling without changing chemistry."""
    core = smiles.split('|',1)[0]
    labels = {}
    atom_index = -1
    tokens = r'(\[[^\]]+\]|Br|Cl|[A-Z][a-z]?|[bcnops*])|\{([^{}]+)\}'
    for match in re.finditer(tokens,core):
        if match.group(1):
            atom_index += 1
        else:
            if atom_index < 0:
                raise ValueError('An atom label precedes its atom')
            labels[atom_index] = match.group(2)
    core = re.sub(r'\{[^{}]+\}','',core)
    core = re.sub(r'([/\\])\(',lambda match:'('+match.group(1),core)
    options = Chem.SmilesParserParams()
    options.removeHs = False
    mol = Chem.MolFromSmiles(core,options)
    if mol is not None:
        for index,label in labels.items():
            Chem.SetAtomAlias(mol.GetAtomWithIdx(index),label)
    return mol


def _sdf_matches_smiles(path, smiles):
    if not os.path.isfile(path):
        return False
    expected = structure_molecule(smiles)
    if expected is None:
        return False
    expected = Chem.RemoveHs(expected)
    try:
        supplier = Chem.SDMolSupplier(path, removeHs=True)
        cached = supplier[0] if len(supplier) else None
        return (cached is not None
                and cached.GetNumAtoms() == expected.GetNumAtoms()
                and cached.GetNumBonds() == expected.GetNumBonds()
                and cached.HasSubstructMatch(expected, useChirality=True))
    except (OSError, ValueError, RuntimeError):
        return False


def _ensure_sdf_current(smiles, path):
    # A molecular graph cannot verify geometry directives such as ring flips.
    if "|" not in smiles and _sdf_matches_smiles(path, smiles):
        return
    smiles_to_sdf(smiles, path)
    if not _sdf_matches_smiles(path, smiles):
        raise RuntimeError(f"Generated SDF does not match the odorant structure: {path}")


def ensure_sdf_exists(odorant):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")
    o = find_odorant(odorant)
    if not o:
        raise ValueError(f'Unknown odorant: {odorant}')
    result = o.get("sdfname", "sdf/" + o['full_name'].replace(' ', '_') + ".sdf")
    _ensure_sdf_current(o['smiles'], result)
    o["sdfname"] = result

    isomers = check_isomers(o['full_name'])
    forms = check_forms(o['full_name'])
    if isomers:
        if len(isomers):
            for iso in o["isomers"].keys():
                if "preiso" in o.keys():
                    l = len(o['preiso'])
                    fname = ("sdf/" + (o["full_name"][0:l] + iso + "-" + o["full_name"][l:])).replace(' ', '_') + ".sdf"
                else:
                    fname = ("sdf/" + iso + "-"+o["full_name"]).replace(' ', '_') + ".sdf"
                _ensure_sdf_current(o["isomers"][iso], fname)
    if forms:
        if len(forms):
            for form in o["forms"].keys():
                if "preform" in o.keys():
                    l = len(o['preform'])
                    fname = ("sdf/" + (o["full_name"][0:l] + form + "-" + o["full_name"][l:])).replace(' ', '_') + ".sdf"
                else:
                    fname = ("sdf/" + form + "-"+o["full_name"]).replace(' ', '_') + ".sdf"
                _ensure_sdf_current(o["forms"][form], fname)
    return result

def smiles_to_sdf(smiles, output_file):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")
    pieces = smiles.split("|")
    ring_atoms = []
    for directive in pieces[1:]:
        match = re.fullmatch(r"rflp(?:\s*:\s*|\s+)(.+)", directive.strip())
        if match:
            ring_atoms.extend(match.group(1).split())
    if ring_atoms:
        # Direct Python callers also refresh the native writer after source changes.
        current = subprocess.run(["make", "-q", "bin/ringflip"])
        if current.returncode:
            subprocess.run(["make", "bin/ringflip"], check=True)
    destination = os.path.abspath(output_file)
    with tempfile.TemporaryDirectory(prefix=".aroma-sdf-", dir=os.path.dirname(destination)) as folder:
        staged = os.path.join(folder, "ligand.sdf")
        mol = structure_molecule(pieces[0])
        if mol is None:
            raise ValueError(f'Invalid ligand structure: {output_file}')
        else:
            mol = Chem.AddHs(mol)
            embedded = AllChem.EmbedMolecule(mol)
            for seed in (11973,11974,11975):
                if embedded >= 0:
                    break
                embedded = AllChem.EmbedMolecule(mol,randomSeed=seed,useRandomCoords=True)
            if embedded < 0:
                raise RuntimeError(f"Could not generate a 3D conformer: {output_file}")
            AllChem.UFFOptimizeMolecule(mol)
            with Chem.SDWriter(staged) as writer:
                writer.write(mol)
        if ring_atoms:
            cmd = ["bin/ringflip", staged, *ring_atoms]
            print(" ".join(cmd), flush=True)
            subprocess.run(cmd, check=True)
        if not _sdf_matches_smiles(staged, smiles):
            raise RuntimeError(f"Generated SDF does not match the odorant structure: {output_file}")
        os.replace(staged, destination)
    print(f"Saved {output_file}")

