# portable Python script for numerical output file comparisons of MOPAC

from shutil import copyfile
from sys import argv
from itertools import zip_longest
import subprocess
import difflib
import os
import re
import math

# TODO:
# - detect eigenvector blocks by waiting for "Root No." or "ROOT NO."
# - parse eigenvector #, eigenvalues, & eigenvector components, ignore leading characters for electronic case
# - ignore possibly degenerate blocks at the edge of the spectrum (unless root # starts at 1)

# thresholds of acceptable errors
NUMERIC_THRESHOLD = 1e-7
EIGVEC_THRESHOLD = 1e-4

# regular expression pattern for a time stamp or other signifier of timing output, "CLOCK" or "TIME" or "SECONDS", & system-dependent versioning
skip_criteria = re.compile('([A-Z][a-z][a-z] [A-Z][a-z][a-z] [ 0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9] [0-9][0-9][0-9][0-9])'
                           '|(CLOCK)|(TIME)|(SECONDS)|(Version)')

# make a local copy of the input & other necessary files
for file in argv[3:]:
   copyfile(os.path.join(argv[1],file),file)

# run MOPAC in the local directory
subprocess.call([argv[2],argv[3]])

# only compare ".out" output files, which have the same name as ".mop" or ".ent" input files
out_name = argv[3][:-3]+'out'
ref_path = os.path.join(argv[1],out_name)

# open the pair of output files being compared
with open(ref_path,'r') as out_ref, open(out_name,"r") as out_new:

    # loop over pairs of output file lines to be compared
    for line_number, (ref_line, new_line) in enumerate(zip_longest(out_ref, out_new, fillvalue='')):

        # decide if the lines need to be analyzed further
        ref_skip = skip_criteria.search(ref_line)
        new_skip = skip_criteria.search(new_line)
        if (ref_skip is not None and new_skip is None) or (ref_skip is None and new_skip is not None):
            assert 0, f'''ERROR: skip mismatch between output files on line {line_number}
                          REF: {ref_line}
                          NEW: {new_line}'''
        elif ref_skip is None:

            # engage/parse eigenvector mode

            # compare words & approximately compare numbers
            for ref_word, new_word in zip(ref_line.split(), new_line.split()):

                # try to convert to floats
                try:
                    ref_float = float(ref_word)
                    new_float = float(new_word)

                    error = abs(ref_float - new_float)
#            if math.isfinite(new_number):
                    if error > NUMERIC_THRESHOLD:
                        assert 0, f'''ERROR: float mismatch between output files on line {line_number}
                                      REF: {ref_line}
                                      NEW: {new_line}'''

                # otherwise assume a word
                except ValueError:
                    if ref_word != new_word:
                        assert 0, f'''ERROR: word mismatch between output files on line {line_number}
                                      REF: {ref_line}
                                      NEW: {new_line}'''
