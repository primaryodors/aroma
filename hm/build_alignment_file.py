#!/usr/bin/env python3
"""
Builds the Modeller alignment file (allgpcr.ali) by combining experimental structures,
coupled structures, and target receptor sequences.
"""

import os
import sys
import re
import urllib.request
import subprocess
from natsort import natsorted

# -----------------------------------------
# ARCHITECTURAL ENFORCEMENT: PATH RESOLUTION
# -----------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(script_dir)
sys.path.append(root_dir)

try:
    import data.protutils as pu
except ImportError as e:
    print(f"FATAL ERROR: Could not import data.protutils. Ensure you are running this from the correct environment. - {e}", file=sys.stderr)
    sys.exit(1)

# Initialize the global protein database
pu.load_prots()

def chunk_string(s: str, length: int) -> list:
    return [s[i:i+length] for i in range(0, len(s), length)]

def main():
    exp_ali_path = os.path.join(script_dir, "experimental.ali")
    out_ali_path = os.path.join(script_dir, "allgpcr.ali")
    tpl_dir = os.path.join(script_dir, "tpl")
    coupled_dir = os.path.join(root_dir, "coupled")
    tmp_dir = os.path.join(root_dir, "tmp")

    os.makedirs(tpl_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    if not os.path.exists(exp_ali_path):
        print(f"FATAL ERROR: {exp_ali_path} not found. You must place experimental.ali in the hm/ directory.", file=sys.stderr)
        sys.exit(1)

    with open(exp_ali_path, 'r', encoding='utf-8') as f:
        exp_content = f.read()

    deletions = {
        "5cxv": ("A", 1000, 1200)
    }

    already = set()
    prevln = ""
    
    print("Executing Phase 1: Parsing experimental alignments and syncing templates...", file=sys.stderr)
    for ln in exp_content.split("\n"):
        if prevln.startswith(">P1;"):
            rcpid = prevln[4:].strip()
            already.add(rcpid)
            
            pieces = ln.split(':')
            if len(pieces) > 1:
                pseqid = pieces[1].strip()
                fn = os.path.join(tpl_dir, f"{pseqid}.pdb")
                
                if not os.path.exists(fn) and re.match(r"^[0-9A-Za-z]{4}$", pseqid):
                    rcsbid = pseqid.upper()
                    url = f"https://files.rcsb.org/download/{rcsbid}.pdb"
                    print(f"Template missing. Downloading {rcsbid} from RCSB...", file=sys.stderr)
                    try:
                        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                        with urllib.request.urlopen(req, timeout=10) as response:
                            pdb_data = response.read().decode('utf-8')
                            
                        # Apply necessary hardcoded deletions to the PDB
                        if pseqid in deletions:
                            delstrand, delsr, deler = deletions[pseqid]
                            cleaned_lines = []
                            for line in pdb_data.split('\n'):
                                if line.startswith("ATOM  "):
                                    strand = line[21:22]
                                    resno_str = line[22:26].strip()
                                    if resno_str.isdigit():
                                        resno = int(resno_str)
                                        if strand == delstrand and delsr <= resno <= deler:
                                            continue
                                cleaned_lines.append(line)
                            pdb_data = "\n".join(cleaned_lines)
                            
                        with open(fn, 'w') as f_out:
                            f_out.write(pdb_data)
                        print(f"Secured {rcsbid}.", file=sys.stderr)
                    except Exception as e:
                        print(f"WARNING: Network execution failed for {rcsbid} - {e}", file=sys.stderr)
        prevln = ln

    print("Executing Phase 2: Processing coupled structures...", file=sys.stderr)
    cpl_content = ""
    coupled_ali_path = os.path.join(coupled_dir, "coupled.ali")
    
    if os.path.exists(coupled_dir):
        if os.path.exists(coupled_ali_path):
            with open(coupled_ali_path, 'r', encoding='utf-8') as f:
                cpl_content = f.read()
        else:
            lfam, lsub = "zero", "zero"
            
            for rcpid, p in pu.prots.items():
                fam = pu.family_from_protid(rcpid)
                sub = pu.subfamily_from_protid(rcpid)
                
                if fam == lfam and sub == lsub:
                    continue
                    
                path = os.path.join(coupled_dir, fam, sub)
                if not os.path.exists(path):
                    continue
                    
                files = natsorted([f for f in os.listdir(path) if f.endswith('.pdb') and '~' in f])
                
                for entry in files:
                    pieces = entry[:-4].split('~')
                    if len(pieces) < 2:
                        continue
                    rcpid_cpl, gprot = pieces[0], pieces[1]
                    
                    phew_script = f'LOAD "coupled/{fam}/{sub}/{entry}" A A\nECHO $SEQUENCEA\n'
                    phew_path = os.path.join(tmp_dir, "sequence.phew")
                    with open(phew_path, 'w') as f:
                        f.write(phew_script)
                        
                    # Execute bin/phew from the root directory safely
                    try:
                        result = subprocess.run(
                            ["bin/phew", "tmp/sequence.phew"],
                            cwd=root_dir,
                            capture_output=True,
                            text=True,
                            timeout=10
                        )
                        cplseq = re.sub(r"\s+", "", result.stdout)
                    except Exception as e:
                        print(f"ERROR executing phew binary on {entry}: {e}", file=sys.stderr)
                        continue
                            
                    if "aligned" in pu.prots.get(rcpid_cpl, {}):
                        cpllen = len(cplseq)
                        rcpseq = pu.prots[rcpid_cpl]["sequence"]
                        seqlen = len(rcpseq)
                        
                        if cpllen != seqlen:
                            print(f"WARNING: SEQUENCE LENGTH MISMATCH {entry}. Skipping.", file=sys.stderr)
                            continue
                            
                        cplaln = ""
                        rcpaln = pu.prots[rcpid_cpl]["aligned"]
                        j = 0
                        for char in rcpaln:
                            if 'A' <= char <= 'Z':
                                if j < len(cplseq):
                                    cplaln += cplseq[j]
                                    j += 1
                            else:
                                cplaln += char
                                
                        temp = cplaln[:-4] if cplaln.endswith('----') else cplaln
                        chunks = chunk_string(temp, 130)
                        caligned = "\n".join(chunks) + "\n" if chunks else ""
                        
                        p1row = f">P1;{rcpid_cpl}~{gprot}"
                        famno = int(re.sub(r"[^0-9]", "", fam)) if re.sub(r"[^0-9]", "", fam) else 0
                        memno = rcpid_cpl[len(fam)+len(sub):]
                        structrow = f"structure:{rcpid_cpl}~{gprot}:FIRST:A:LAST :A:Olfactory receptor family {famno} subfamily {sub} member {memno}:Homo sapiens: 2.00: 0.20"
                        
                        cpl_content += f"{p1row}\n{structrow}\n{caligned}---------------------------------------------------------------------------------------------------*\n\n"
                
                lfam, lsub = fam, sub

    print("Executing Phase 3: Constructing final alignment output...", file=sys.stderr)
    try:
        with open(out_ali_path, 'w', encoding='utf-8') as fp:
            fp.write(exp_content)
            fp.write("\n\n")
            
            for rcpid, p in pu.prots.items():
                if "aligned" not in p or rcpid in already:
                    continue
                
                temp = p["aligned"][:-4] if p["aligned"].endswith('----') else p["aligned"]
                chunks = chunk_string(temp, 130)
                paligned = "\n".join(chunks) + "\n" if chunks else ""
                
                fam = pu.family_from_protid(rcpid)
                mem = pu.member_from_protid(rcpid)
                pname = "(unspecified)"
                
                if fam == "TAAR":
                    pname = f"Trace amine-associated receptor {mem}"
                elif fam == "VN1R":
                    pname = f"Vomeronasal type 1 receptor number {mem}"
                elif fam == "MS4A":
                    pname = f"Membrane-spanning 4A receptor {mem}"
                else:
                    famn = fam[2:] if fam.startswith("OR") else fam
                    sub = pu.subfamily_from_protid(rcpid)
                    pname = f"Olfactory receptor family {famn} subfamily {sub} number {mem}"
                    
                seqlen = len(p["sequence"])
                deets = f"sequence:{rcpid}:1     :A:{seqlen}  :A:{pname}:Homo sapiens: 1.90: 0.19"
                
                fp.write(f">P1;{rcpid}\n{deets}\n{paligned}---------------------------------------------------------------------------------------------------*\n\n")
                
            fp.write("\n\n" + cpl_content)
            
        print("Execution complete. Alignment file forged.", file=sys.stderr)
    except Exception as e:
        print(f"FATAL ERROR: Failed to write {out_ali_path} - {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()