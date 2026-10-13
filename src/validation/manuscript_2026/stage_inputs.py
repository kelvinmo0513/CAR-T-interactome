"""Stage authorized external files from the October 7 reproducibility package."""
import argparse
import json
import shutil
from pathlib import Path
from common import MANIFEST, new_directory, require, sha256, write_json

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--include-reference', action='store_true')
    args = parser.parse_args()
    rows = json.loads(MANIFEST.read_text())['files']
    rows = [r for r in rows if r['group'] == 'analysis' or args.include_reference]
    for row in rows:
        source = args.bundle / row['source_path']
        require(source.is_file(), f'Missing bundle member: {row["source_path"]}')
        require(sha256(source) == row['sha256'], f'Source checksum mismatch: {row["source_path"]}')
    output = new_directory(args.output_dir)
    for row in rows:
        destination = output / row['file']
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.bundle / row['source_path'], destination)
        require(sha256(destination) == row['sha256'], 'Copy checksum mismatch')
    write_json(output/'input_integrity.json', {'status':'PASS','files':rows})
    print(f'Staged and verified {len(rows)} external files.')

if __name__ == '__main__':
    main()
