import difflib
import math
import pathlib
import re
import sys


NUMBER = re.compile(r"(?<![\w.])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def equivalent(expected, actual):
    left = list(NUMBER.finditer(expected))
    right = list(NUMBER.finditer(actual))
    if len(left) != len(right):
        return False
    left_end = right_end = 0
    for a, b in zip(left, right):
        if expected[left_end:a.start()] != actual[right_end:b.start()]:
            return False
        av, bv = a.group(), b.group()
        if av != bv:
            if not all(any(c in value for c in ".eE") for value in (av, bv)):
                return False
            if not math.isclose(float(av), float(bv), rel_tol=5e-6, abs_tol=1e-9):
                return False
        left_end, right_end = a.end(), b.end()
    return expected[left_end:] == actual[right_end:]


if __name__ == "__main__":
    approved, received = map(pathlib.Path, sys.argv[1:])
    expected, actual = approved.read_text(), received.read_text()
    if not equivalent(expected, actual):
        sys.stdout.writelines(difflib.unified_diff(expected.splitlines(True), actual.splitlines(True),
                                                  fromfile=str(approved), tofile=str(received)))
        sys.exit(1)
