"""File integrity and explicit validation helpers; no network access."""
import hashlib
import json
from pathlib import Path

MODULE = Path(__file__).resolve().parent
REPOSITORY = MODULE.parents[2]
MANIFEST = REPOSITORY / 'data/manuscript_2026/input_manifest.json'

def require(condition, message):
    if not condition:
        raise ValueError(message)

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def new_directory(path):
    path = Path(path).resolve()
    require(not path.is_relative_to(REPOSITORY), 'Keep patient-level inputs and generated outputs outside the public checkout.')
    require(not path.exists(), f'Output already exists: {path}; choose a new directory.')
    path.mkdir(parents=True)
    return path

def check_inputs(directory, group='analysis'):
    manifest = json.loads(MANIFEST.read_text())
    records = [r for r in manifest['files'] if r['group'] == group]
    for row in records:
        path = Path(directory) / row['file']
        require(path.is_file(), f'Missing input: {row["file"]}; see input_manifest.json.')
        require(sha256(path) == row['sha256'], f'Input checksum mismatch: {row["file"]}')
    return records

def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
