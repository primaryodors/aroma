"""Keep completed receptor diagnostics outside the PDB collection."""
import hashlib
from pathlib import Path
import os


def archive_path(model):
    model = Path(model).resolve()
    storage = next((parent for parent in model.parents if parent.name == 'pdbs'), None)
    if storage is None or not model.is_file():
        return None
    digest = hashlib.sha256(str(model).encode()+b'\0'+model.read_bytes()).hexdigest()
    return storage.parent / 'tmp' / 'reec-reports' / f'{model.name}.{digest}.json'


def report_path(model):
    sidecar = Path(str(model)+'.json')
    if sidecar.is_file():
        return sidecar
    archived = archive_path(model)
    return archived if archived is not None and archived.is_file() else sidecar


def archive_report(model):
    """Remove a completed sidecar after retaining its coordinate-specific report."""
    sidecar = Path(str(model)+'.json')
    archived = archive_path(model)
    if archived is None or not sidecar.is_file():
        return
    archived.parent.mkdir(parents=True, exist_ok=True)
    temporary = archived.with_name(archived.name+f'.copy-{os.getpid()}')
    try:
        temporary.write_bytes(sidecar.read_bytes())
        os.replace(temporary, archived)
        sidecar.unlink()
    finally:
        temporary.unlink(missing_ok=True)
