# portable Python script for output file comparisons

from shutil import copyfile
from sys import argv
import subprocess
import difflib
import os
import re

# regular expression pattern for a time stamp or other signifier of timing output, "CLOCK" or "TIME" or "SECONDS", & system-dependent versioning
time_stamp = re.compile("([A-Z][a-z][a-z] [A-Z][a-z][a-z] [ 0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9] [0-9][0-9][0-9][0-9])"
                        "|(CLOCK)|(TIME)|(SECONDS)|(Version)")

# list of local files before the test
files_before = os.listdir(".")

# make a local copy of the input & other necessary files
for file in argv[3:]:
   copyfile(os.path.join(argv[1],file),file)

# run MOPAC in the local directory
subprocess.call([argv[2],argv[3]])
files_after = os.listdir(".")
files_change = set(files_after) - set(files_before)

for file in files_change:

    # load in both outputs for analysis
    out_new = open(file,"r")
    out_ref = open(os.path.join(argv[1],file),"r")

    file_test = True
    while file_test:

        # get a non-timing line from the reference output file
        line_test = True
        while line_test:
            bare_line = out_ref.readline()
            old_line = bare_line.rstrip('\r\n')
            if old_line != '' or bare_line == '':
                line_test = False
            if time_stamp.search(old_line):
                line_test = True

        # get a non-timing line from the test output file
        line_test = True
        while line_test:
            bare_line = out_new.readline()
            new_line = bare_line.rstrip('\r\n')
            if new_line != '' or bare_line == '':
                line_test = False
            if time_stamp.search(new_line):
                line_test = True

        # flag an error if the next two non-timing lines do not match
        if new_line != old_line:
            print("OLD: ",old_line)
            print("NEW: ",new_line)
            assert 0, "ERROR: mismatch between test & reference output files"

        # exit the loop if the most recent line
        if bare_line == '':
            file_test = False

    # close & erase file after comparisons
    out_new.close()
    os.remove(file)
