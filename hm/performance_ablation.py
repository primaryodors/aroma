#!/usr/bin/env python3
"""Build four-switch docking variants and time identical, isolated full docks."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
CONSTANTS = ROOT / "src/classes/constants.h"
SWITCHES = ("AROMA_OPT_SPATIAL", "AROMA_OPT_SHIELDING", "AROMA_OPT_HYDROGEN", "AROMA_OPT_CALLBACK")
CASES = (("all_enabled", None), ("spatial_off", SWITCHES[0]), ("shielding_off", SWITCHES[1]),
         ("hydrogen_off", SWITCHES[2]), ("callback_off", SWITCHES[3]))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def switches(disabled=None):
    return {name: int(name != disabled) for name in SWITCHES}


def set_switches(values):
    # Change just the value bytes, preserving every other part of the header.
    with CONSTANTS.open("r+b") as handle:
        content = handle.read()
        for name, value in values.items():
            match = re.search(rb"^#define " + name.encode() + rb" ([01])$", content, re.M)
            if not match:
                raise RuntimeError(f"Switch declaration absent: {name}")
            replacement = str(value).encode()
            if match.group(1) != replacement:
                handle.seek(match.start(1))
                handle.write(replacement)


def command(arguments, log, env=None):
    with Path(log).open("w") as output:
        return subprocess.run(arguments, cwd=ROOT, env=env, stdout=output, stderr=subprocess.STDOUT).returncode


def make_config(folder, source):
    result = []
    for line in source.splitlines():
        if line.lstrip().startswith("#"):
            continue
        fields = line.split()
        key = fields[0] if fields else ""
        if key == "OUT":
            line = f"OUT {folder / 'result.dock'}"
        elif key == "OUTPDB":
            line = f"OUTPDB 1 {folder / 'model%o.pdb'}"
        elif key == "PROGRESS":
            continue
        result.append(line)
    path = folder / "case.config"
    path.write_text("\n".join(result) + "\n")
    return path


def save_rows(folder, rows):
    (folder / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = ("case", "disabled", "wall_seconds", "user_seconds", "system_seconds", "max_rss_kib",
              "exit_code", "poses_found", "best_energy_kcal", "delta_seconds", "percent_change", "status")
    with (folder / "results.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_variants(folder, env):
    fingerprints = {}
    try:
        for name, disabled in CASES:
            case = folder / name
            case.mkdir(exist_ok=True)
            values = switches(disabled)
            set_switches(values)
            print(f"BUILD {name}: {values}", flush=True)
            code = command(["make", "-j2", "bin/aromadock", "test/residue_grid_test"], case / "build.log")
            if code:
                raise RuntimeError(f"Build failed for {name}; see {case / 'build.log'}")
            code = command(["test/residue_grid_test"], case / "spatial_test.log", env)
            if code:
                raise RuntimeError(f"Residue comparisons failed for {name}")
            match = re.search(r"fingerprint (\d+)", (case / "spatial_test.log").read_text())
            if not match:
                raise RuntimeError(f"Residue fingerprint absent for {name}")
            fingerprints[name] = match.group(1)
            if name != "all_enabled" and fingerprints[name] != fingerprints["all_enabled"]:
                raise RuntimeError(f"Residue selection differs for {name}: {fingerprints}")
            binary = case / "aromadock"
            shutil.copy2(ROOT / "bin/aromadock", binary)
            metadata = {"switches": values, "binary_sha256": digest(binary), "spatial_fingerprint": match.group(1)}
            (case / "build.json").write_text(json.dumps(metadata, indent=2) + "\n")
            print(f"BUILT {name}; residue selections match", flush=True)
    finally:
        set_switches(switches())
        print("RESTORE all-enabled applications", flush=True)
        code = command(["make", "-j2", "bin/aromadock", "bin/ic"], folder / "restore_build.log")
        if code:
            raise RuntimeError("All-enabled restoration build failed")
    (folder / "spatial_fingerprints.json").write_text(json.dumps(fingerprints, indent=2) + "\n")


def run_variants(folder, env, cpu, timeout):
    source = (ROOT / "tmp/OR1A1~d-limonene.active.config").read_text()
    rows = []
    for name, disabled in CASES:
        case = folder / name
        binary = case / "aromadock"
        metadata = json.loads((case / "build.json").read_text())
        if digest(binary) != metadata["binary_sha256"]:
            raise RuntimeError(f"Binary digest differs for {name}")
        config = make_config(case, source)
        timed = case / "time.txt"
        arguments = ["/usr/bin/time", "-o", str(timed), "-f",
                     "wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kib=%M\nexit_code=%x",
                     "taskset", "-c", str(cpu), str(binary), str(config), "--nosplash"]
        print(f"RUN {name}: full 10-pose, 50-iteration dock; CPU {cpu}", flush=True)
        began = time.monotonic()
        with (case / "run.log").open("w") as output:
            process = subprocess.Popen(arguments, cwd=ROOT, env=env, stdout=output,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            timed_out = False
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                import signal
                os.killpg(process.pid, signal.SIGTERM)
                code = process.wait()
                timed_out = True
        row = {"case": name, "disabled": disabled or "none", "exit_code": code,
               "wall_seconds": round(time.monotonic()-began, 2),
               "status": "timeout" if timed_out else ("complete" if code == 0 else "failed")}
        if timed.exists():
            for line in timed.read_text().splitlines():
                if "=" in line:
                    key, value = line.split("=", 1)
                    row[key] = float(value) if key.endswith("seconds") else int(value)
        log = (case / "run.log").read_text(errors="replace")
        poses = re.search(r"(\d+) pose\(s\) found", log)
        energy = re.search(r"Best candidate pose energy:\s*([-+\d.eE]+)", log)
        row["poses_found"] = int(poses.group(1)) if poses else 0
        row["best_energy_kcal"] = float(energy.group(1)) if energy else None
        if "Calculation time:" not in log and not timed_out:
            row["status"] = "incomplete"
        baseline = rows[0]["wall_seconds"] if rows else row["wall_seconds"]
        row["delta_seconds"] = round(row["wall_seconds"]-baseline, 2)
        row["percent_change"] = round((row["wall_seconds"]/baseline-1)*100, 2)
        row.update(metadata)
        rows.append(row)
        save_rows(folder, rows)
        print(f"DONE {name}: {row['wall_seconds']:.2f}s; {row['poses_found']} poses; {row['status']}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "out/performance_ablation_2026-11-11")
    parser.add_argument("--stage", choices=("all", "build", "run"), default="all")
    parser.add_argument("--seed", type=int, default=11973)
    parser.add_argument("--cpu", type=int, default=max(os.sched_getaffinity(0)))
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    folder = args.results.resolve()
    folder.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["AROMA_RANDOM_SEED"] = str(args.seed)
    metadata = {"seed": args.seed, "cpu": args.cpu, "platform": platform.platform(),
                "poses": 10, "iterations": 50, "runs_per_case": 1,
                "config_sha256": digest(ROOT / "tmp/OR1A1~d-limonene.active.config"),
                "protein_sha256": digest(ROOT / "pdbs/OR1/OR1A1.active.pdb"),
                "ligand_sha256": digest(ROOT / "sdf/d-limonene.sdf")}
    (folder / "protocol.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if args.stage in ("all", "build"):
        build_variants(folder, env)
    if args.stage in ("all", "run"):
        run_variants(folder, env, args.cpu, args.timeout)
    print(f"RESULTS {folder / 'results.csv'}", flush=True)


if __name__ == "__main__":
    main()
