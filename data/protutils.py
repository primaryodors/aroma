

import json
import re
import subprocess
import os
import os.path
import math
from natsort import natsorted
import data.globals
import data.geometry

def load_prots():
    global prots
    with open('data/receptor.json', 'r') as file:
        prots = json.load(file)

def bw_insdel(protid, tmrno, offset):
    global prots
    result = 0
    prot = prots[protid]

    if 'deletion' in prot:
        for dltn in prot['deletion']:
            pettia = dltn.split('.')
            dtmr = pettia[0]
            doff = pettia[1]

            if dtmr == tmrno:
                if doff <  50 and doff >= offset: result += 1
                if doff >= 50 and doff <  offset: result -= 1

    if 'insertion' in prot:
        for insn in prot['insertion']:
            pettia = insn.split('.')
            dtmr = pettia[0]
            doff = pettia[1]

            if dtmr == tmrno:
                if doff <  50 and doff >= offset: result -= 1
                if doff >= 50 and doff <  offset: result += 1

    return result

def resno_from_bw(protid, bw):
    global prots
    if not prots.get(protid): raise Exception("Protein not found: "+protid)
    bwraw = re.sub('[A-Z]', '', bw)
    aminos = bw[:len(bw)-len(bwraw)]
    bw = bwraw

    pieces = bwraw.split('.',2)
    if len(pieces) < 2: return int(bwraw)
    tmrno = int(pieces[0])
    offset = int(pieces[1])
    insdel = bw_insdel(protid, tmrno, offset)

    res50 = int(prots[protid]["bw"][str(tmrno)+".50"])
    result = res50 + offset - 50 + insdel

    if aminos:
        c = prots[protid]['sequence'][result-1]
        if not c in aminos: return 0

    return result

def aalet_at_resno(protid, resno):
    global prots
    if not prots.get(protid): raise Exception("Protein not found: "+protid)
    return prots[protid]["sequence"][resno-1]

def family_from_protid(protid):
    if protid[0:2] == "OR":
        return "OR" + str(int(re.sub("[^0-9]", "", protid[2:4])))
    else:
        return protid[0:4]

def subfamily_from_protid(protid):
    if protid[0:2] == "OR":
        fam = family_from_protid(protid)
        l = len(fam)
        return re.sub("[^A-Z]", "", protid[l:l+2])
    else:
        return ""

def build_pdb_fname(protid, suffix="active"):
    fam = family_from_protid(protid)
    return f"pdbs/{fam}/{protid}.{suffix}.pdb"

def templates_for_hm(protid):
    consOR51 = ["8uxv"]
    consOR52 = ["8hti"]
    OR51E2 = ["8f76"]
    consOR1 = ["8uxy"]
    OR1D2 = ["9w45"]
    consOR2 = ["8uy0"]
    consOR4 = ["8uyq"]
    consOR5 = ["9wpm"]
    OR5V1 = ["9lkb"]
    consOR6 = ["9ldv", "9ldw", "9ldx", "9ldz"]
    bmOR6A2 = ["9le0", "9le1", "9le2"]
    CLASSI = OR51E2 + consOR51 + consOR52
    CLASSII = consOR1 + consOR2 + consOR4 + consOR5 + consOR6
    TAAR1 = ["8jln", "8jlo", "8jlp", "8jlq", "8jlr", "8jso"]
    mTAAR = ["8iwe", "8iwm", "8itf", "8iw4", "8iw9", "8pm2"]
    ADORA2A = ["6gdg"]
    ADRB2 = ["7dhr", "8gej"]
    LPAR1 = ["7td0", "7yu3"]
    TAS2R = ["7xp6"]
    CB = ["5xr8", "5xra"]
    CHRM1 = ["6oij"]

    fam = family_from_protid(protid)
    if protid[0:2] == "OR":
        sub = subfamily_from_protid(protid)

        if fam == "OR1" or fam == "OR3" or fam == "OR7":
            return consOR1 + OR1D2

        if fam == "OR2" or fam == "OR13":
            return consOR2

        if fam == "OR4" or fam == "OR12":
            return consOR4

        if fam == "OR5" or fam == "OR8" or fam == "OR9":
            if protid == "OR5V1": return OR5V1
            else: return consOR5

        if fam == "OR6" or fam == "OR10" or fam == "OR11":
            if fam == "OR6" and (sub == 'A' or sub == 'B' or sub == 'P' or sub == 'Y'):
                return bmOR6A2
            else: return consOR6

        if fam == "OR14":
            return CLASSII

        if fam == "OR51":
            if sub == 'D' or sub == 'E': return OR51E2
            else: return consOR51

        if fam == "OR52":
            return consOR52

        if fam == "OR56":
            return CLASSI

    elif fam == "TAAR":
        return TAAR1 + mTAAR

    return ADORA2A + ADRB2 + LPAR1 + TAS2R + CB + CHRM1

def aa_similarity(a, b):
    if a == b: return 1.0
    apol = a in "DSTYHNQKER"
    bpol = b in "DSTYHNQKER"
    asml = a in "ACGST"
    bsml = b in "ACGST"
    aacd = a in "DE"
    bacd = b in "DE"
    aamd = a in "NQ"
    bamd = b in "NQ"
    abas = a in "HKR"
    bbas = b in "HKR"
    aaro = a in "FHWY"
    baro = b in "FHWY"
    aslf = a in "MC"
    bslf = b in "MC"

    # if both acidic
    if aacd and bacd: return 0.8

    # if both basic
    if abas and bbas:
        if aaro == baro: return 0.8
        return 0.5

    # if acid and amide
    if (aacd and bamd) or (bacd and aamd):
        return 0.75

    # if both polar or both nonpolar
    if apol == bpol:
        if asml and bsml: return 0.8
        if aslf and bslf: return 0.75
        if aaro == baro: return 0.75
        return 0.5

    # if both aromatic
    if aaro and baro: return 0.5

    # if both small
    if asml and bsml: return 0.333

    # if no match
    return 0

def aln_similarity(aln1, aln2):
    aln1 = aln1.strip().upper()
    aln2 = aln2.strip().upper()
    len1 = len(aln1)
    len2 = len(aln2)
    mlen = min(len1, len2)
    total = 0
    ttlsim = 0.0
    for i in range(0, mlen):
        a = aln1[i]
        b = aln2[i]
        if a < "A" or a > "Z" or b < "A" or b > "Z": continue
        sim = aa_similarity(a, b)
        ttlsim += sim
        total += 1

    if not total: return 0
    return ttlsim / total

def prepare_coupled(inpfn, outfn, rcpid, remarks):
    fam = family_from_protid(rcpid)
    phew  = f"LOAD \"{inpfn}\" A A\n"
    phew += f"LOAD \"pdbs/{fam}/{rcpid}.active.pdb\" A R\n"
    phew += "STRAND A\n"
    phew += "BWCOPY R A\n"
    phew += "UNCHAIN R\n"
    phew += "UPRIGHT\n"
    phew += "HYDRO\n"
    phew += "BWCENTER\n"
    if remarks: phew += remarks
    with open(f"pdbs/{fam}/{rcpid}.active.pdb", "r") as f:
        indat = f.read()
        lines = indat.split("\n")
        for ln in lines:
            if ln[0:16] == "REMARK 800 SITE ":
                phew += f"{ln}\n"
    phew += f"SAVE {outfn}\n"
    with open("tmp/cpl.phew", "w") as f:
        f.write(phew)
    cmd = [ "bin/phew", "tmp/cpl.phew" ]
    subprocess.run(cmd)

def ali_rel_resno(ali: str, helix: int, member: int):
    lines = ali.split("\n")

    if helix < 1 or helix > 8:
        raise Exception("Helix out of range")

    result = 0
    for i in range(helix-1):
        for c in lines[i]:
            if c >= 'A' and c <= 'Z': result += 1

    for i in range(100):
        c = lines[helix-1][i]
        if c >= 'A' and c <= 'Z': result += 1

    return result + member - 50

def custom_pdb_template(aln, rcpid, output_fname):                             # Writes a file and returns the templates and the .ali block.
    with open("../hm/experimental.ali", "r") as f:
        c = f.read().__str__()

    # Choose experimental structures by grouping.
    fam = family_from_protid(rcpid)
    exclude = []
    if fam == "OR1" or fam == "OR3" or fam == "OR7":
        exclude.extend(["8f76", "8hti", "8uxv", "8uy0", "8uyq", "9wpm", "9lkb", "9ldv", "9ldw", "9ldx", "9ldz", "9le0", "9le1", "9le2"])
    elif fam == "OR2" or fam == "OR13":
        exclude.extend(["8f76", "8hti", "8uxv", "9wpm", "9lkb", "9ldv", "9ldw", "9ldx", "9ldz", "9le0", "9le1", "9le2"])
    elif fam == "OR4" or fam == "OR12":
        exclude.extend(["8f76", "8hti", "8uxv", "9wpm", "9lkb", "9ldv", "9ldw", "9ldx", "9ldz", "9le0", "9le1", "9le2"])
    elif fam == "OR5" or fam == "OR8" or fam == "OR9":
        exclude.extend(["8f76", "8hti", "8uxv", "9ldv", "9ldw", "9ldx", "9ldz", "9le0", "9le1", "9le2"])
    elif fam == "OR6" or fam == "OR10" or fam == "OR11":
        exclude.extend(["8f76", "8hti", "8uxv", "8uy0", "9wpm", "9lkb"])
    elif fam == "OR14":
        exclude.extend(["8f76", "8hti", "8uxv", "9ldv", "9ldw", "9ldx", "9ldz", "9le0", "9le1", "9le2"])
    elif fam == "OR51" or fam == "OR52" or fam == "OR56":
        exclude.extend(["8uxy", "9w45", "8uy0", "8uyq", "9wpm", "9lkb", "9ldv", "9ldw", "9ldx", "9ldz", "9le0", "9le1", "9le2"])

    # Read in the alignments of the experimental structures.
    reading_aln = False
    alns = dict()
    strandids = dict()
    for ln in c.split("\n"):
        if ln[0:9] == "structure":
            if "inactive" in ln or "apo " in ln:
                continue
            reading_aln = True
            current_aln = ""
            pieces = ln.split(':')
            pdbid = pieces[1].strip()
            if pdbid in exclude: continue
            strandids[pdbid] = pieces[3].strip()
        elif reading_aln:
            current_aln += ln + "\n"
            if '*' in ln:
                reading_aln = False
                alns[pdbid] = current_aln
                current_aln = ""

    # Find the top few closest matches to the input sequence and weight them by similarity.
    closest_ids = []
    closest_sim = []
    for pdbid in alns.keys():
        if pdbid in exclude: continue
        simaln = aln_similarity(aln, alns[pdbid])
        # print(f"{pdbid}: {simaln}")
        for i in range(5):
            if i >= len(closest_ids):
                closest_ids.append(pdbid)
                closest_sim.append(simaln)
                break
            else:
                if simaln > closest_sim[i]:
                    for j in range(len(closest_ids)-1, i, -1):
                        closest_ids[j] = closest_ids[j-1]
                        closest_sim[j] = closest_sim[j-1]
                    closest_ids[i] = pdbid
                    closest_sim[i] = simaln
                    break

    print(closest_ids)
    weights = []

    # Make the highest weight 5 times the lowest weight, and make the weights add up to 1.
    # Since we want the highest to be 5 times as much as the lowest, that means we want it to be 4 times more.
    span = closest_sim[0] - closest_sim[4]
    tosub = closest_sim[4]; # - span/4
    for i in range(5):
        weights.append(closest_sim[i] - tosub)

    divisor = sum(weights)
    for i in range(5):
        weights[i] = weights[i] / divisor

    print(weights)

    # For each closest structure, grab the locations of 3.50:CA, 7.50:CA, and 6.50:CA.
    # Then determine the rotation and rescan the PDB, saving rotated atom locations for
    # backbone and CB atoms of all structures, and all atoms of the zeroth structure.
    frist = True
    atomxyz = dict()
    seq0 = dict()
    for pdbid in closest_ids:
        if pdbid in exclude: continue
        if not pdbid in strandids: continue
        relres3 = ali_rel_resno(alns[pdbid], 3, 50)
        relres6 = ali_rel_resno(alns[pdbid], 6, 50)
        relres7 = ali_rel_resno(alns[pdbid], 7, 50)
        pt3 = False
        pt6 = False
        pt7 = False

        atomxyz[pdbid] = dict()

        # Find the PDB and read through it to obtain the CA atom coordinates.
        pdbfname = f"../hm/tpl/{pdbid}.pdb"
        if not os.path.exists(pdbfname): raise Exception(f"File not found {pdbfname}")
        j = 0
        with open(pdbfname, "r") as f:
            c = f.read().__str__()
            lines = c.split("\n")
            for ln in lines:
                if ln[0:5] == "ATOM ":
                    if strandids[pdbid] == ln[21]:
                        aname = ln[12:16].strip()
                        if aname == "CA":
                            j += 1
                            if j == relres3 or j == relres6 or j == relres7:
                                x = float(ln[30:38])
                                y = float(ln[38:46])
                                z = float(ln[46:54])

                                if j == relres3: pt3 = [x,y,z]
                                if j == relres6: pt6 = [x,y,z]
                                if j == relres7: pt7 = [x,y,z]

            if not pt3: raise Exception(f"Residue 3.50 not found in {pdbid}")
            if not pt6: raise Exception(f"Residue 6.50 not found in {pdbid}")
            if not pt7: raise Exception(f"Residue 7.50 not found in {pdbid}")

            # Obtain the translation
            xlation = [-pt3[0], -pt3[1], -pt3[2]]
            pt3 = data.geometry.add_points(pt3, xlation)
            pt6 = data.geometry.add_points(pt6, xlation)
            pt7 = data.geometry.add_points(pt7, xlation)

            if frist:
                gpt3 = pt3
                gpt6 = pt6
                gpt7 = pt7
            else:
                # Obtain the rotations
                rot1 = data.geometry.align_points_3d(pt7, gpt7, gpt3)
                pt3 = data.geometry.rotate3D(pt3, gpt3, rot1, rot1[3])              # a trick since rotate3D() only uses the first three members of the axis argument
                pt6 = data.geometry.rotate3D(pt6, gpt3, rot1, rot1[3])
                pt7 = data.geometry.rotate3D(pt7, gpt3, rot1, rot1[3])
                # rot2 = data.geometry.align_points_3d(pt6, gpt6, gpt3)
                rot2 = data.geometry.subtract_points(pt7, pt3)
                rot2.append(data.geometry.find_angle_along_vector(pt6, gpt6, gpt3, rot2))

            # Rescan and fill in the atom positions, minding the alignments
            lresno = 0
            j = 0
            bwhelix = 1
            bwmember = -49
            for ln in lines:
                if ln[0:5] == "ATOM ":
                    if strandids[pdbid] == ln[21]:
                        resno = int(ln[22:26].strip())
                        aname = ln[12:16].strip()
                        a3let = ln[17:20]

                        if resno != lresno:
                            j += 1
                            if alns[pdbid][j] == "\n":
                                bwhelix += 1
                                bwmember = -49
                            else:
                                bwmember += 1

                            try:
                                while alns[pdbid][j] < 'A' or alns[pdbid][j] > 'Z':
                                    if alns[pdbid][j] == "\n":
                                        bwhelix += 1
                                        bwmember = -49
                                    else:
                                        bwmember += 1
                                    j += 1
                                    if j >= len(alns[pdbid]):
                                        alns[pdbid] += "*"
                            except:
                                if not pdbid in alns:
                                    print(f"No {pdbid} in alns.")
                                else:
                                    print("Index out of range")
                                    print(pdbid)
                                    print(alns[pdbid])
                                    print(len(alns[pdbid]))
                                print(j)
                                exit()

                            if frist:
                                seq0[f"{bwhelix}.{bwmember}"] = a3let

                            if bwmember == 50:
                                print(f"{pdbid} {bwhelix}.50 = {a3let}{resno}")

                        if frist or aname in ["N", "CA", "CB", "C", "O"]:
                            x = float(ln[30:38])
                            y = float(ln[38:46])
                            z = float(ln[46:54])

                            result = [x,y,z]
                            result = data.geometry.add_points(result, xlation)
                            if not frist:
                                result = data.geometry.rotate3D(result, gpt3, rot1, rot1[3])
                                result = data.geometry.rotate3D(result, gpt3, rot2, rot2[3])

                            atomxyz[pdbid][f"{bwhelix}.{bwmember}:{aname}"] = result
                            # if frist: print(f"atomxyz[{pdbid}][{a3let}{bwhelix}.{bwmember}:{aname}] = {result}")

                        lresno = resno

        frist = False

    # The OR5V1 cryo-EM is missing the EXR2 helix, but accurate predictions require that helix, so we'll fill it in from the consOR5 structure.
    pdbid0 = closest_ids[0]
    pdbid1 = closest_ids[1]
    outali = alns[closest_ids[0]]
    h = 5
    for m in range(19, 58):
        aname = f"{h}.{m}:CA"
        seqkey = f"{h}.{m}"
        if not aname in atomxyz[pdbid0]:
            if aname in atomxyz[pdbid1]:
                atomxyz[pdbid0][aname] = atomxyz[pdbid1][aname].copy()
                seq0[seqkey] = "GLY"

                mstr = m+49
                outalilns = outali.split("\n")
                outalilns[h-1] = outalilns[h-1][0:mstr] + "G" + outalilns[h-1][mstr+1:]
                outali = "\n".join(outalilns)
            else:
                print(f"WARNING: {aname} not found in {pdbid1}")

    atomxyz[pdbid0] = dict(natsorted(atomxyz[pdbid0].items()))
    seq0 = dict(natsorted(seq0.items()))

    # Generate a PDB of transposed and rotated 3D coordinates of the weighted average of the closest sequences.
    atno = 1
    resno = 0
    lresbw = ""
    with open(output_fname, "w") as f:
        for aname in atomxyz[pdbid0].keys():
            xyz = atomxyz[pdbid0][aname].copy()
            weight = weights[0]
            xyz[0] *= weight
            xyz[1] *= weight
            xyz[2] *= weight

            j = -1
            for lpdbid in closest_ids:
                j += 1
                if lpdbid == pdbid0: continue
                if aname in atomxyz[lpdbid]:
                    newxyz = atomxyz[lpdbid][aname].copy()
                    newxyz[0] *= weights[j]
                    newxyz[1] *= weights[j]
                    newxyz[2] *= weights[j]
                    xyz = data.geometry.add_points(xyz, newxyz)
                    weight += weights[j]

            if weight:
                xyz[0] /= weight
                xyz[1] /= weight
                xyz[2] /= weight

            ardata = aname.split(":")
            resbw = ardata[0]
            aname = ardata[1]
            if resbw != lresbw: resno += 1

            if resbw in seq0.keys():
                ln = "ATOM  "
                ln += str(atno).rjust(5)
                ln += "  " + aname.ljust(4)
                ln += seq0[resbw] + " A"
                ln += str(resno).rjust(4)
                ln += "    "
                x = xyz[0]
                ln += "-" if x < 0 else " "
                x = math.fabs(x)
                ln += f"{x:.3f}".zfill(7)
                y = xyz[1]
                ln += "-" if y < 0 else " "
                y = math.fabs(y)
                ln += f"{y:.3f}".zfill(7)
                z = xyz[2]
                ln += "-" if z < 0 else " "
                z = math.fabs(z)
                ln += f"{z:.3f}".zfill(7)
                ln += "  1.00  0.00           "
                ln += aname[0] + "  "
                f.write(ln+"\n")
                atno += 1

            lresbw = resbw
        f.write("TER\n")

        if 0:
            strord = ord('A')
            for pdbid in closest_ids:
                # if pdbid == pdbid0: continue
                resno = 0
                lresbw = ""
                strord += 1
                strand = chr(strord)
                for aname in atomxyz[pdbid].keys():
                    xyz = atomxyz[pdbid][aname]

                    ardata = aname.split(":")
                    resbw = ardata[0]
                    aname = ardata[1]
                    if resbw != lresbw: resno += 1

                    ln = "ATOM  "
                    ln += str(atno).rjust(5)
                    ln += "  " + aname.ljust(4)
                    ln += "ALA " + strand
                    ln += str(resno).rjust(4)
                    ln += "    "
                    x = xyz[0]
                    ln += "-" if x < 0 else " "
                    x = math.fabs(x)
                    ln += f"{x:.3f}".zfill(7)
                    y = xyz[1]
                    ln += "-" if y < 0 else " "
                    y = math.fabs(y)
                    ln += f"{y:.3f}".zfill(7)
                    z = xyz[2]
                    ln += "-" if z < 0 else " "
                    z = math.fabs(z)
                    ln += f"{z:.3f}".zfill(7)
                    ln += "  1.00  0.00           "
                    ln += aname[0] + "  "

                    f.write(ln+"\n")

                    atno += 1
                    lresbw = resbw
                f.write("TER\n")
        f.write("END\n")

    # The template's .ali data will be the same as those of the closest match.
    return " ".join(closest_ids) + "\n" + outali


def json_encode_pretty(array):
    return re.sub(r"(\s*)([^\s]*) ([{[(])\n", r"\1\2\n\1\3\n", json.dumps(array, indent = 4, default=lambda o: o.__dict__)).replace("\n\n", "\n")

bsrs = \
    [
        "2.53", "3.29", "3.32", "3.33", "3.36", "3.37", "3.40", "3.41",
        "4.53", "4.57", "4.60", "45.49", "45.51", "45.52", "45.53",
        "5.39", "5.42", "5.43", "5.46", "5.47",
        "6.48", "6.51", "6.55", "6.59", "7.38", "7.39", "7.42"
    ]
