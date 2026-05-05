import sys
import json
import re
import os
import subprocess
import traceback

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("WARNING: RDKit not found in environment. System defaulting to OpenBabel fallback.", file=sys.stderr)


def load_odors():
    """Extracts the odorant database with absolute path enforcement."""
    # 1. Determine the exact, absolute directory of THIS script (odorutils.py)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 2. Forge the unbreakable path to the database
    db_path = os.path.join(current_dir, 'odorant.json')
    
    # 3. Validate existence before execution
    if not os.path.exists(db_path):
        print(f"FATAL ERROR: Odorant database not found at {db_path}", file=sys.stderr)
        return None
        
    # 4. Extract and return the payload
    # 4. Extract and normalize the payload
    try:
        with open(db_path, 'r', encoding='utf-8') as file:
            odor_data = json.load(file)
            
            # --- NORMALIZATION ---
            # If the JSON is a dictionary (mapping IDs to molecules), flatten it.
            if isinstance(odor_data, dict):
                odor_data = list(odor_data.values())
                
            return odor_data
    except Exception as e:
        print(f"FATAL ERROR: Failed to extract data from {db_path} - {e}", file=sys.stderr)
        return None

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
                    if data.get("type") == "na":
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

    aroma1 = re.sub("[^a-z0-9+-]", "", aroma.lower().replace(' ', '_'))
    for oid, o in odors.items():
        namematch = False
        i = 1
        while "name"+str(i) in o.keys():
            if o["name"+str(i)] == aroma:
                namematch = True
            i += 1

        if o['smiles'] == aroma \
            or ("iupac" in o.keys() and o['iupac'] == aroma) \
            or namematch \
            or re.sub( "[^a-z0-9+-]", "", o['full_name'].lower().replace(' ', '_') ) == aroma1:
            retval = o
            retval['oid'] = oid

            pwd = os.getcwd()
            os.chdir(os.path.dirname(os.path.abspath(__file__)))
            os.chdir("..")
            fname = "sdf/" + o['full_name'].replace(' ', '_') + ".sdf"
            if os.path.exists(fname):
                retval["sdfname"] = fname

            os.chdir(pwd)
            return retval

    return False

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

def ensure_sdf_exists(odorant):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")
    o = find_odorant(odorant)
    if "sdfname" in o:
        result = o["sdfname"]
    else:
        output_file = "sdf/" + o['full_name'].replace(' ', '_') + ".sdf"
        smiles_to_sdf(o['smiles'], output_file)
        o["sdfname"] = result = output_file

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
                if not os.path.exists(fname):
                    smiles_to_sdf(o["isomers"][iso], fname)
    if forms:
        if len(forms):
            for form in o["forms"].keys():
                if "preform" in o.keys():
                    l = len(o['preform'])
                    fname = ("sdf/" + (o["full_name"][0:l] + form + "-" + o["full_name"][l:])).replace(' ', '_') + ".sdf"
                else:
                    fname = ("sdf/" + form + "-"+o["full_name"]).replace(' ', '_') + ".sdf"
                if not os.path.exists(fname):
                    smiles_to_sdf(o["forms"][form], fname)
    return result

def smiles_to_sdf(smiles, output_file):
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")
    try:
        pieces = smiles.split("|")
        smiles = pieces[0]

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            cmd = ["obabel", "--gen3d", "-osdf", f"-O{output_file}", f"-:{smiles}" ]
            print(" ".join(cmd))
            subprocess.run(cmd)
            return
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        if mol:
            with Chem.SDWriter(output_file) as writer:
                writer.removeHs = False
                writer.write(mol)
            print(f"Saved {output_file}")

            if len(pieces) > 1:
                pieces[1] = pieces[1].split(":")
                sub4 = pieces[1][0]
                rest = pieces[1][1]
                if sub4 == "rflp":
                    cmd = ["bin/ringflip", output_file]
                    for larg in rest.strip().split(" "):
                        cmd.append(larg)
                    print(" ".join(cmd), "\n")
                    subprocess.run(cmd)
        else:
            print(f"Invalid SMILES string: {smiles}")
    except Exception as e:
        print(f"Error occurred generating 3D structure: {e}")

