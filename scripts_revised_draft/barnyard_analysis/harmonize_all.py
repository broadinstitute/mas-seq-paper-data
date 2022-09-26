import os
import sys
import subprocess

long_h5_root_path = '/home/jupyter/mb-ml-data-disk/MAS-seq-analysis/data/t-cell-vdj/long/quant/revised/variants'
short_h5_path = '/home/jupyter/mb-ml-data-disk/MAS-seq-analysis/output/t-cell-vdj-cite-seq/M132TS_both.h5ad'
output_path = './barnyard_analysis'
file_list = './filelist'
input_prefix = 'M132TS_MAS_15x_'

with open(file_list) as of:
    for filename in of:
        
        # extract prefix
        output_prefix = filename[filename.find(input_prefix):].split('.')[0]
        print(f'Processing {output_prefix} ...')
        
        # bash cmd
        long_h5_path = os.path.join(long_h5_root_path, filename.strip())
        cmd = f"python harmonize_long_short_adata.py -s {short_h5_path} -l {long_h5_path} -o {output_path} -p {output_prefix}"
        
        proc = subprocess.Popen(cmd.split(), stderr=sys.stderr, stdout=sys.stdout)
        stdout, stderr = proc.communicate()
        proc_status = proc.wait()

