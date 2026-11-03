from pathlib import Path


def pdb_disulfides(filename, chain="A"):
    lines = Path(filename).read_text().splitlines()
    sulfurs = {}
    for line in lines:
        if (line.startswith("ATOM  ") and line[12:16].strip() == "SG"
                and line[17:20] == "CYS" and line[21:22] == chain):
            sulfurs[int(line[6:11])] = int(line[22:26])
    residues = set(sulfurs.values())
    pairs = set()
    for line in lines:
        if line.startswith("CONECT"):
            serials = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()]
            if not serials or serials[0] not in sulfurs:
                continue
            for serial in serials[1:]:
                if serial in sulfurs and serial != serials[0]:
                    pairs.add(tuple(sorted((sulfurs[serials[0]], sulfurs[serial]))))
        elif line.startswith("SSBOND") and line[15:16] == chain and line[29:30] == chain:
            first, second = int(line[17:21]), int(line[31:35])
            if first in residues and second in residues and first != second:
                pairs.add(tuple(sorted((first, second))))
    return sorted(pairs)
