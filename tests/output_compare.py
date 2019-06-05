# portable Python script for output file comparisons

from shutil import copyfile
from sys import argv
import subprocess
import difflib
import os
import re
import math

# regular expression pattern for a time stamp or other signifier of timing output, "CLOCK" or "TIME" or "SECONDS", & system-dependent versioning
time_stamp = re.compile("([A-Z][a-z][a-z] [A-Z][a-z][a-z] [ 0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9] [0-9][0-9][0-9][0-9])"
                        "|(CLOCK)|(TIME)|(SECONDS)|(Version)")

# make a local copy of the input & other necessary files
for file in argv[3:]:
   copyfile(os.path.join(argv[1],file),file)

# run MOPAC in the local directory
subprocess.call([argv[2],argv[3]])

# only compare ".out" files
out_name = argv[3][:-3]+"out"

out_new = open(out_name,"r")
out_ref = open(os.path.join(argv[1],out_name),"r")

line_number = 0
max_error_location = 0
max_error_value = 0.0
file_test = True
while file_test:

    # get a non-timing line from the reference output file
    line_test = True
    while line_test:
        bare_line = out_ref.readline()
        old_line = bare_line.rstrip(' \r\n')
        line_number += 1
        if old_line != '' or bare_line == '':
            line_test = False
        if time_stamp.search(old_line):
            line_test = True

    # get a non-timing line from the test output file
    line_test = True
    while line_test:
        bare_line = out_new.readline()
        new_line = bare_line.rstrip(' \r\n')
        if new_line != '' or bare_line == '':
            line_test = False
        if time_stamp.search(new_line):
            line_test = True

    # flag an error if the next two non-timing lines do not match
    old_words = old_line.split()
    new_words = new_line.split()
    for old_word,new_word in zip(old_words,new_words):
        try:
            old_number = float(old_word)
            new_number = float(new_word)
            if math.isfinite(new_number):
                error = abs(old_number - new_number)
            if error > max_error_value:
                max_error_location = line_number
                max_error_value = error
        except:
            pass

    # exit the loop if the most recent line is blank (normal end-of-file behavior)
    if bare_line == '':
        file_test = False

#
print("largest error of",max_error_value,"on line",max_error_location)
if max_error_value > 1e-7:
    assert 0, "ERROR: numerical deviation between test & reference is larger than threshold"
