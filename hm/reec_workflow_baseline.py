#!/usr/bin/env python3
"""Run real modeling workflows with isolated writable paths and input snapshots."""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


class RecipeCaptured(RuntimeError):
    pass


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def isolated_output(workspace, output):
    resolved = Path(output).resolve()
    if not resolved.is_relative_to(workspace):
        raise ValueError('Benchmark output must stay inside its isolated workspace')
    return resolved


def prepare_workspace(workspace):
    workspace.mkdir(parents=True, exist_ok=False)
    for name in ('tmp', 'out', 'output', 'sdf', 'results', 'coupled', 'hm', 'obj', 'bin', 'data'):
        (workspace / name).mkdir()
    for name in ('src', 'makefile'):
        (workspace / name).symlink_to(ROOT / name)
    for path in (ROOT / 'data').iterdir():
        if path.is_file():
            (workspace / 'data' / path.name).symlink_to(path)
    for path in (ROOT / 'hm').glob('*.py'):
        (workspace / 'hm' / path.name).symlink_to(path)
    for path in (ROOT / 'hm').glob('*.ali'):
        shutil.copy2(path, workspace / 'hm' / path.name)
    shutil.copytree(ROOT / 'hm/tpl', workspace / 'hm/tpl')
    shutil.copytree(ROOT / 'pdbs', workspace / 'pdbs')
    shutil.copytree(ROOT / 'sdf', workspace / 'sdf', dirs_exist_ok=True)
    for path in (ROOT / 'obj').glob('*.o'):
        shutil.copy2(path, workspace / 'obj' / path.name)
    for path in (ROOT / 'bin').iterdir():
        if path.is_file():
            shutil.copy2(path, workspace / 'bin' / path.name)
    canonical = {}
    for name in ('OR1A1', 'OR1B1', 'OR51E1', 'OR51E2', 'OR2AG1', 'OR5V1'):
        for path in (ROOT / 'pdbs').glob(f'*/{name}.*.pdb'):
            canonical[str(path.relative_to(ROOT))] = digest(path)
    manifest = {'source_root': str(ROOT), 'workspace': str(workspace),
                'canonical_pdb_hashes': canonical,
                'shared_read_only_inputs': ['data/*', 'src', 'makefile', 'hm/*.py'],
                'source_hashes': {str(path.relative_to(ROOT)): digest(path) for path in
                                  [ROOT / 'hm/dohm.py', ROOT / 'hm/fixfail.py',
                                   ROOT / 'hm/modeling_backend.py', ROOT / 'hm/reec_backend.py',
                                   ROOT / 'src/reec.cpp', ROOT / 'bin/reec', ROOT / 'bin/aromadock']}}
    (workspace / 'snapshot.json').write_text(json.dumps(manifest, indent=2) + '\n')


def bind_workspace(workspace, entry):
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / 'hm'))
    backend = importlib.import_module('reec_backend')
    backend.ROOT, backend.HM = workspace, workspace / 'hm'
    module = importlib.import_module(entry)
    module.__file__ = str(workspace / 'hm' / f'{entry}.py')
    module.script_dir, module.root_dir = str(workspace / 'hm'), str(workspace)
    for name in ('globals', 'protutils', 'odorutils', 'dyncenter'):
        shared = importlib.import_module('data.' + name)
        shared.__file__ = str(workspace / 'data' / f'{name}.py')
    os.chdir(workspace)
    return module, backend


def record_processing_stages(workspace):
    original = subprocess.run
    folder = workspace / 'results/stages'
    folder.mkdir(exist_ok=True)
    records = []

    def run(command, *positional, **keywords):
        capture = (isinstance(command, (list, tuple)) and len(command) > 1
                   and Path(str(command[0])).name == 'ic' and 'minc' in command and 'save' in command)
        if not capture:
            return original(command, *positional, **keywords)
        output = isolated_output(workspace, command[1])
        number = len(records) + 1
        record = {'command': [str(value) for value in command], 'output': str(output)}
        records.append(record)
        path = folder / f'{output.stem}.{number}.before-ic.pdb'
        shutil.copy2(output, path)
        record['before-ic'] = str(path)
        try:
            return original(command, *positional, **keywords)
        finally:
            if output.exists():
                path = folder / f'{output.stem}.{number}.after-ic.pdb'
                shutil.copy2(output, path)
                record['after-ic'] = str(path)
            (folder / f'{output.stem}.stages.json').write_text(json.dumps(records, indent=2) + '\n')

    subprocess.run = run
    return original


def capture_recipe(workspace, backend):
    def capture(receptor, args, work, templates, alignment, contacts, environment=None,
                fit_environment=False, environment_hydrogens=False):
        folder = workspace / 'results/recipes' / receptor
        folder.mkdir(parents=True, exist_ok=False)
        copied = []
        inputs = []
        for index, value in enumerate(templates):
            value = str(value)
            name, filename = value.split('=', 1) if '=' in value else ('', value)
            source = Path(filename)
            target = folder / f'template-{index}.pdb'
            shutil.copy2(source, target)
            copied.append(f'{name}={target}' if name else str(target))
            inputs.append({'source': str(source), 'sha256': digest(target)})
        target_alignment = folder / 'protein.ali'
        shutil.copy2(alignment, target_alignment)
        target_contacts = folder / 'receptor.contacts'
        text = '\n'.join(contacts) + '\n'
        for path in args.contacts:
            text += Path(path).read_text() + '\n'
        target_contacts.write_text(text)
        target_environment = None
        if environment:
            target_environment = folder / 'environment.pdb'
            shutil.copy2(environment, target_environment)
        recipe = {'receptor': receptor, 'alignment': str(target_alignment),
                  'contacts': str(target_contacts), 'templates': copied, 'template_inputs': inputs,
                  'environment': str(target_environment) if target_environment else None,
                  'fit_environment': fit_environment, 'environment_hydrogens': environment_hydrogens,
                  'chain': args.chain, 'anneal_residues': args.anneal_residues,
                  'alignment_sha256': digest(target_alignment), 'contacts_sha256': digest(target_contacts)}
        (folder / 'recipe.json').write_text(json.dumps(recipe, indent=2) + '\n')
        print(f'Captured workflow recipe: {folder / "recipe.json"}', flush=True)
        raise RecipeCaptured()

    original = backend.model
    backend.model = capture
    return original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run', 'process', 'verify', 'capture'])
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--entry', choices=['dohm', 'fixfail'], default='dohm')
    parser.add_argument('--target')
    parser.add_argument('--candidate', type=Path)
    parser.add_argument('--output', type=Path)
    args, workflow_args = parser.parse_known_args()
    workspace = args.workspace.resolve()
    if args.action == 'prepare':
        if workflow_args:
            parser.error('Unexpected preparation arguments')
        prepare_workspace(workspace)
        print(f'Isolated workspace: {workspace}')
        return 0
    if not (workspace / 'snapshot.json').is_file():
        parser.error('Workspace has no input snapshot manifest')
    if args.action in ('run', 'capture') and workflow_args[:1] != ['--']:
        parser.error('Place workflow arguments after --')
    if args.action in ('run', 'capture'):
        for index, value in enumerate(workflow_args):
            if value == '--output' and index + 1 < len(workflow_args):
                isolated_output(workspace, workflow_args[index+1])
            elif value.startswith('--output='):
                isolated_output(workspace, value.split('=', 1)[1])
    elif args.action == 'process' and args.output is not None:
        args.output = isolated_output(workspace, args.output)
    if args.action == 'verify':
        manifest = json.loads((workspace / 'snapshot.json').read_text())
        changed = [name for name, expected in manifest['canonical_pdb_hashes'].items()
                   if digest(ROOT / name) != expected]
        report = {'canonical_models_unchanged': not changed, 'changed': changed}
        (workspace / 'protection-check.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps(report))
        return 1 if changed else 0
    module, backend = bind_workspace(workspace, args.entry)
    original_run = record_processing_stages(workspace)
    original_model = capture_recipe(workspace, backend) if args.action == 'capture' else None
    started = time.perf_counter()
    initial_cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
    outcome = 0
    try:
        if args.action in ('run', 'capture'):
            if workflow_args[:1] == ['--']:
                workflow_args = workflow_args[1:]
            outcome = module.main(workflow_args) or 0
        else:
            if not args.target or not args.candidate or not args.output or args.entry != 'dohm':
                parser.error('Processing requires a target, candidate, output, and dohm entry')
            backend.pu.load_prots()
            from data.disulfides import pdb_disulfides
            family = backend.pu.family_from_protid(args.target)
            subfamily = backend.pu.subfamily_from_protid(args.target)
            family_no = int(''.join(c for c in family if c.isdigit()) or '0')
            pairs = pdb_disulfides(args.candidate)
            first, second = [a for a, b in pairs], [b for a, b in pairs]
            module.process_model(args.target, family, family_no, family + subfamily,
                                 args.output, args.candidate, str(args.candidate),
                                 first, second, first, second)
    except RecipeCaptured:
        return 0
    except backend.CandidateRejected as exc:
        print(f'Candidate rejected: {exc}', file=sys.stderr)
        return 2
    finally:
        subprocess.run = original_run
        if original_model is not None:
            backend.model = original_model
        final_cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_seconds = final_cpu.ru_utime + final_cpu.ru_stime - initial_cpu.ru_utime - initial_cpu.ru_stime
        print(f'Workflow elapsed seconds: {time.perf_counter()-started:.3f}', flush=True)
        print(f'Child CPU seconds: {cpu_seconds:.3f}', flush=True)
    return outcome


if __name__ == '__main__':
    raise SystemExit(main())
