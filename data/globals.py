
import time
import os
import os.path
import re
import subprocess

global prots
prots = []

global odors
odors = []

def wait_cool_cpu():
    cmd = ["sensors"]
    while 1:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE)
        waited = False
        for ln in proc.stdout.decode().split('\n'):
            temp = False
            high = False
            offset = -10

            # One style of `sensors` output.
            if re.search("Package id [0-9]+:\\s+[+-]?[0-9.]+", ln):
                condensed = re.sub("\\s+", ' ', ln).replace("Package id ", "Package_id_")
                pieces = condensed.split(' ')
                temp = float(re.sub("[^0-9.-]", "", pieces[1]))
                rhigh = re.search("(high|crit)\\s+=\\s+[+-]?[0-9.]+", ln)
                strhigh = rhigh.group()
                high = float(re.sub("[^0-9.-]", "", strhigh))
                if strhigh[0:4] == "high":
                    offset = -10
                elif strhigh[0:4] == "crit":
                    offset = -20

            # Another style.
            elif re.search("CPU Temperature:\\s+[+-]?[0-9.]+"):
                condensed = re.sub("\\s+", ' ', ln).replace("CPU Temperature", "CPU_Temperature")
                pieces = condensed.split(' ')
                temp = float(re.sub("[^0-9.-]", "", pieces[1]))

            # User setting takes precedence over whatever `sensors` says is too hot.
            if os.path.exists("cpu_temp_limit"):
                with open("cpu_temp_limit", "r") as f:
                    c = f.read()
                    high = float(re.sub("[^0-9.-]", "", c))

            # If processor is too hot, wait for external processes to end and CPU to cool.
            if temp:
                if not high:
                    high = 75
                diff = temp - high
                # print(f"Diff {diff} offset {offset}")
                if diff >= offset:
                    print("Processor too hot; waiting...")
                    waited = True
                    time.sleep(300)

        # If all clear, go for it.
        if not waited: break


