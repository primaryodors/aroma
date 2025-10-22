
import json
import re
import os
import subprocess
import data.globals

def load_odors():
    global odors
    with open('data/odorant.json', 'r') as file:
        odors = json.load(file)

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

    aroma1 = re.sub("/[^a-z0-9_+-]/", "", aroma.lower().replace(' ', '_'))
    for oid, o in odors.items():
        namematch = False
        i = 1
        while "name"+str(i) in o.keys():
            if o["name"+str(i)] == aroma:
                namematch = True
            i += 1

        if (o['smiles'] == aroma
            or re.sub( "/[^a-z0-9_+-]/", "", o['full_name'].lower().replace(' ', '_') ) == aroma1
            or ("iupac" in o.keys() and o['iupac'] == aroma)
            or namematch):
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
        result.append(iso+"-"+odor["full_name"])
    return result

def ensure_sdf_exists(odorant):
    o = find_odorant(odorant)
    if not "sdfname" in o:
        cmd = ["obabel", "-:"+o["smiles"], "--gen3D", "-osdf", "-Osdf/" + o['full_name'].replace(' ', '_') + ".sdf"]
        print(" ".join(cmd), "\n")
        subprocess.run(cmd)

    isomers = check_isomers(o['full_name'])
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    os.chdir("..")
    if len(isomers):
        for iso in o["isomers"].keys():
            fname = "sdf/" + iso + "-"+o["full_name"].replace(' ', '_') + ".sdf"
            if not os.path.exists(fname):
                cmd = ["obabel", "-:"+o["isomers"][iso], "--gen3D", "-osdf", "-O" + fname]
                print(" ".join(cmd), "\n")
                subprocess.run(cmd)

