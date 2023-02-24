import os
import sys
import subprocess

long_h5_root_path = './barnyard_analysis'
short_h5_root_path = './barnyard_analysis'
output_path = './barnyard_analysis'
file_list = './filelist'
input_prefix = 'M132TS_MAS_15x_'

with open(file_list) as of:
    for filename in of:
        
        # extract prefix
        output_prefix = filename[filename.find(input_prefix):].split('.')[0]
        print(f'Processing {output_prefix} ...')
        
        # bash cmd
        short_h5_path = os.path.join(short_h5_root_path, f"{output_prefix}.harmonized.barnyard.short.h5ad")
        long_h5_path = os.path.join(long_h5_root_path, f"{output_prefix}.harmonized.barnyard.long.h5ad")
        cmd = f"python barnyard_analysis.py -s {short_h5_path} -l {long_h5_path} -o {output_path} -p {output_prefix}"
        
        proc = subprocess.Popen(cmd.split(), stderr=sys.stderr, stdout=sys.stdout)
        stdout, stderr = proc.communicate()
        proc_status = proc.wait()
