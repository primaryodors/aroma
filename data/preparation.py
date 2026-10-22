from collections import OrderedDict
from pathlib import Path
import subprocess


class PreparationCache:
    def __init__(self, capacity=8):
        self.capacity = capacity
        self.results = OrderedDict()

    def run(self, command, inputs, before_run=None):
        paths = [Path(command[0]), *(Path(path) for path in inputs)]
        fingerprint = tuple((str(path.resolve()), path.stat().st_mtime_ns, path.stat().st_size)
                            for path in paths)
        key = (tuple(command), fingerprint)
        if key in self.results:
            self.results.move_to_end(key)
            return self.results[key]
        if before_run:
            before_run()
        result = subprocess.run(command, stdout=subprocess.PIPE, check=True)
        self.results[key] = result
        while len(self.results) > self.capacity:
            self.results.popitem(last=False)
        return result
