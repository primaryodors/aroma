
import sys
import re
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.chdir('..')
sys.path.append(os.getcwd())

import data.globals
import data.protutils
import data.odorutils

data.protutils.load_prots()
data.odorutils.load_odors()

colors = \
    {
        "PAILV": "\033[38;5;247m",
        "G": "\033[38;5;243m",
        "M": "\033[38;5;185m",
        "C": "\033[38;5;220m",
        "STNQ": "\033[38;5;49m",
        "DE": "\033[38;5;196m",
        "KR": "\033[38;5;27m",
        "H": "\033[38;5;99m",
        "FWY": "\033[38;5;171m"
    }

bsrs = \
    [
        "2.53", "3.29", "3.32", "3.33", "3.36", "3.37", "3.40", "3.41",
        "4.53", "4.57", "4.60", "45.49", "45.51", "45.52", "45.53",
        "5.39", "5.42", "5.43", "5.46", "5.47",
        "6.48", "6.51", "6.55", "6.59", "7.38", "7.39", "7.42"
    ]

colorless = "\033[0m"
orphan = "\033[2m\033[3m"
deorphan = ""; # "\033[1m"
reset = "\033[22m\033[23m\033[24m"

print("Legend: " + deorphan + "Deorphaned receptor" + reset + " " + orphan + "Orphan receptor" + reset + "\n")

if len(sys.argv) < 2:
    raise Exception("Error no receptors specified.")

rcppatt = sys.argv[1]
rcppl = len(rcppatt)
pattlc = rcppatt[-1:]

bwnos = []
for a in sys.argv[2:]:
    if a == 'bsr':
        bwnos = bwnos + bsrs
    else:
        bwnos.append(a)

print("         ", end="")
for bw in bwnos:
    print(" "+bw, end="")
print("")

for rcpid in data.protutils.prots.keys():
    p = data.protutils.prots[rcpid]
    if rcpid[:rcppl] != rcppatt and not re.match(rcppatt, rcpid):
        continue
    if data.odorutils.empirical_pairs(rcpid, True):
        print(rcpid.ljust(10, " "), end="" )
    else:
        print(orphan + rcpid.ljust(10, " ") + reset, end="")
    
    for bw in bwnos:
        bwl = len(bw)+1
        resno = data.protutils.resno_from_bw(rcpid, bw)
        if resno > 0:
            letter = p["sequence"][resno-1]
            for ckey in colors:
                if letter in ckey:
                    print(colors[ckey], end="")
            print(letter.ljust(bwl, " "), end="")
        else: print("-", end="")
        print(colorless, end="")
    print("")