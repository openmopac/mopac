# Portable Python script to run basic regression tests of MOPAC
# Argument list: <path to testing directory> <path to MOPAC executable> <input file> <data file #1> ... <data file #N>

from shutil import copyfile
from sys import argv
import subprocess
import os

from compare_output import parse_mopac_out_file, compare_mopac_out_file

# make a local copy of the input & other necessary files
for file in argv[4:]:
   copyfile(os.path.join(argv[1],file),file)

# run MOPAC in the local directory
try:
    subprocess.run([argv[2],argv[4]], check=True)
except subprocess.SubprocessError as err:
    print("In attempting to run: ", err.cmd)
    print("stdout: ", err.stdout)
    print("stderr: ", err.stderr)
    raise

# only compare ".out" output files that have the same name as ".mop" or ".ent" input files
out_name = argv[4][:-3]+'out'
ref_path = os.path.join(argv[1],out_name)

# parse the 2 output files that we are comparing
ref_line, ref_list = parse_mopac_out_file(ref_path)
out_line, out_list = parse_mopac_out_file(out_name)

# Run the comparison
compare_mopac_out_file(out_line, out_list, ref_line, ref_list, float(argv[3]))
