#!/usr/bin/env python3
"""Reproduce the selected manuscript analysis; never contacts CIBERSORTx."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from common import MODULE, check_inputs, new_directory, require, sha256, write_json

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--rscript',default='Rscript',help='Rscript executable on PATH or absolute path')
    parser.add_argument('--rebuild-reference',action='store_true',help='Also rebuild and validate the external annotated reference and recorded GUI outputs')
    args=parser.parse_args()
    inp=args.input_dir.resolve()
    records=check_inputs(inp)
    if args.rebuild_reference:
        records+=check_inputs(inp,'reference')
    require(shutil.which(args.rscript) is not None,'Rscript unavailable; install recorded R dependencies first')
    out=new_directory(args.output_dir)
    manifest=dict(status='RUNNING',python=platform.python_version(),dependencies={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','openpyxl']},inputs=records,
                  code_sha256={str(p.relative_to(MODULE)):sha256(p) for p in sorted(MODULE.rglob('*')) if p.suffix in ['.py','.R']},
                  selected_analysis='Four-cohort durability with oversampled two-cell AND Cologne reference; exploratory reference',service_submission=False,steps=[])
    log=out/'execution.log'
    def step(name,command):
        with log.open('a') as stream:
            stream.write('\n'+name+'\n');stream.flush()
            subprocess.run([str(x) for x in command],check=True,stdout=stream,stderr=subprocess.STDOUT)
        manifest['steps'].append(name)
        write_json(out/'run_manifest.json',manifest)
        print(name+': complete',flush=True)
    py=sys.executable
    try:
        step('Prepare external inputs',[py,MODULE/'prepare.py','--input-dir',inp,'--output-dir',out/'prepared'])
        step('Fit selected R models',[args.rscript,MODULE/'analyze.R',out/'prepared',out/'numerical'])
        step('Independent numerical verification',[py,MODULE/'verify.py','--prepared-dir',out/'prepared','--numerical-dir',out/'numerical','--output-dir',out/'verification'])
        step('Render from regenerated tables',[args.rscript,MODULE/'render.R',out/'numerical',out/'plots'])
        if args.rebuild_reference:
            ref=MODULE/'reference';dest=out/'reference'
            step('Rebuild annotated SCT reference',[args.rscript,ref/'prepare_reference.R','--rds',inp/'reference/annotated_ecotyper.rds','--selected-cells',inp/'reference/selected_cells_LOCAL_ONLY.csv','--output',dest/'pre_reference'])
            original=dest/'pre_reference/Kelvin_DP_SCT_reference.txt'
            require(sha256(original)==sha256(inp/'reference/frozen_AND_reference.txt'),'Rebuilt reference differs from approved reference')
            step('Exact two-cell oversampling',[py,ref/'oversample.py','--reference',original,'--selected-cells',inp/'reference/selected_cells_LOCAL_ONLY.csv','--output',dest/'oversampled'])
            derived=dest/'oversampled/AND_2cell_oversampled10x_exploratory.txt'
            require(sha256(derived)=='db37330fe9082caac34e4a4143d3d1139fa637ed0cb5d5a64c232bf1da7b88f3','Oversampled upload differs from approved upload')
            step('Independent reference readback',[py,ref/'validate.py','reference','--original',original,'--derived',derived,'--selected-cells',inp/'reference/selected_cells_LOCAL_ONLY.csv','--column-mapping',dest/'oversampled/column_mapping_LOCAL_ONLY.csv','--output',dest/'reference_readback.json'])
            cmd=[py,ref/'validate.py','service']
            for flag,filename in [('signature','job34_signature.tsv'),('classes','job34_classes.tsv'),('profiles','job34_profiles.tsv'),('mixture','Stahl_24_CPM.tsv'),('scores','job35_scores.tsv')]:
                cmd += ['--'+flag,inp/'cologne'/filename]
            step('Recorded service outputs',cmd+['--output',dest/'service_validation'])
            cmd=[py,ref/'validate.py','logs']
            for flag,filename in [('signature-log','job34_runtime_log.txt'),('deconvolution-log','job35_runtime_log.txt'),('signature-errors','job34_error_log.txt'),('deconvolution-errors','job35_error_log.txt')]:
                cmd += ['--'+flag,inp/'cologne'/filename]
            step('Recorded service logs',cmd+['--output',dest/'service_log_validation.json'])
        manifest['status']='PASS'
    except Exception as error:
        manifest['status']='FAIL'
        manifest['error']=str(error)
        raise
    finally:
        write_json(out/'run_manifest.json',manifest)
    print('Selected manuscript reproduction: PASS')

if __name__=='__main__':
    main()
