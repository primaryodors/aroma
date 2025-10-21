

import json
import re
import data.globals

def load_prots():
    global prots
    with open('data/receptor.json', 'r') as file:
        prots = json.load(file)

def bw_insdel(protid, tmrno, offset):
    #TODO:
    return 0

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

def family_from_protid(protid):
    if protid[0:2] == "OR":
        return "OR" + str(int(re.sub("[^0-9]", "", protid[2:4])))
    else:
        return protid[0:4]

def json_encode_pretty(array):
    return re.sub(r"(\s*)([^\s]*) ([{[(])\n", r"\1\2\n\1\3\n", json.dumps(array, indent = 4, default=lambda o: o.__dict__)).replace("\n\n", "\n")