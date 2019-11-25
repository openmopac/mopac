# Portable Python script for numerical output file comparisons of MOPAC
# Argument list: <path to testing directory> <path to MOPAC executable> <input file> <data file #1> ... <data file #N>

from shutil import copyfile
from sys import argv
from itertools import zip_longest
import subprocess
import difflib
import os
import re
import math
import numpy as np

# MOPAC testing has historically been based on checking that the output files corresponding
# to a test set of input files remain consistent with the output of past MOPAC versions.
# This Python script automates such a test, assuming that we allow for small deviations
# between numerical outputs (with a uniform threshold for simplicity).
# All version/system-dependent output (timing & version info) is ignored.
# Eigenvectors are only stable in the sense of distance between degenerate subspaces,
# which requires more careful identification and analysis of eigenvector blocks.
# We cannot test the stability of all eigenvector information (if we cannot guarantee
# completeness of a degenerate subspace), and untestable data is ignored.
# In principle, we could interpret the symmetry labels to use some of the edge data when
# we can independently determine the size of the subspace, but this is way more trouble than it is worth.

# TODO:
# - pre-filter skipped lines to avoid inconsistency in their number
# - figure out why "IS LESS THAN CUTOFF" lines aren't being flagged as inconsistent
# - debug & wrap up initial testing phase

# thresholds of acceptable errors
NUMERIC_THRESHOLD = 1e-3
EIGVEC_THRESHOLD = 1e-3

# regular expression pattern for a time stamp or other signifier of timing output, "CLOCK" or "TIME" or "SECONDS", & system-dependent versioning
skip_criteria = re.compile('([A-Z][a-z][a-z] [A-Z][a-z][a-z] [ 0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9] [0-9][0-9][0-9][0-9])'
                           '|(CLOCK)|(TIME)|(SECONDS)|(Version)')

# regular expression pattern for an eigenvector block
eigen_criteria = re.compile('(Root No.)|(ROOT NO.)')

# make a local copy of the input & other necessary files
for file in argv[3:]:
   copyfile(os.path.join(argv[1],file),file)

# run MOPAC in the local directory
subprocess.call([argv[2],argv[3]])
print(argv[2],argv[3])

# only compare ".out" output files, which have the same name as ".mop" or ".ent" input files
out_name = argv[3][:-3]+'out'
ref_path = os.path.join(argv[1],out_name)

# initialize state variables of the parsing process
mode = 'standard'
empty_line_counter = 0
ref_modes = []
new_modes = []
ref_eigenvalues = []
new_eigenvalues = []
ref_eigenvectors = []
new_eigenvectors = []

# extract non-skipped lines of output files
with open(ref_path,'r') as ref_file:
    ref_list = [ line for line in ref_file if skip_criteria.search(line) is None ]

with open(out_name,'r') as out_file:
    out_list = [ line for line in out_file if skip_criteria.search(line) is None ]

# loop over pairs of output file lines to be compared
for (ref_line, new_line) in zip_longest(ref_list, out_list, fillvalue=''):

    # count empty lines (used to navigate eigenvalue & eigenvector blocks)
    assert (ref_line == '' and new_line == '') or (ref_line != '' and new_line != ''), \
        f'''ERROR: empty mismatch between output files
            REF: {ref_line}
            NEW: {new_line}'''
    if ref_line == '':
        empty_line_counter += 1
        continue

    # switch to eigenvalue mode if eigenvector blocks are detected
    ref_eigen = eigen_criteria.search(ref_line)
    new_eigen = eigen_criteria.search(new_line)
    assert (ref_eigen is None and new_eigen is None) or (ref_eigen is not None and new_eigen is not None), \
        f'''ERROR: eigen mismatch between output files
            REF: {ref_line}
            NEW: {new_line}'''
    if ref_eigen is not None:
        mode = 'eigenvalue'
        empty_line_counter = 0
        ref_modes += [ int(val) for val in ref_line.split()[2:] ]
        new_modes += [ int(val) for val in new_line.split()[2:] ]

    # parse line in eigenvalue mode
    if mode == 'eigenvalue':

        # NOTE: we are not looking at the symmetry labels for now, the numerical stability of their ordering is unclear

        # read the list of eigenvalues after 2 blank lines
        if empty_line_counter == 2:
            ref_eigenvalues += [ float(val) for val in ref_line.split() ]
            new_eigenvalues += [ float(val) for val in new_line.split() ]

        # check for the end of an eigenvalue block & switch to eigenvector mode after 3 blank lines
        if empty_line_counter >= 3:
            mode = 'eigenvector'
            empty_line_counter = 0

    # parse line in eigenvector mode
    if mode == 'eigenvector':

        # analysis at the end of an eigenvector block that isn't immediately followed by another one
        if empty_line_counter == 2:

            # regroup eigenvector information into a proper matrix layout
            nrow = len(ref_eigenvectors) / len(ref_eigenvalues)
            ncol = len(ref_eigenvalues)
            ref_eigenmatrix = np.empty((nrow,ncol))
            new_eigenmatrix = np.empty((nrow,ncol))
            for block_offset in range(0,ncol,8):
                block_end = min(block_offset+8,nrow)
                block_size = block_end - block_offset
                ref_eigenmatrix[:,block_offset:block_end] = \
                    np.reshape(ref_eigenvectors[block_offset*ncol:block_end*ncol],(nrow,block_size),order='C')

            # identify degenerate eigenvector blocks
            if ref_modes[0] == 1:
                edge_list = []
            else:
                edge_list = [0]
            edge_list += [ index for index in range(ncol-1) if np.abs(ref_eigenvalues[i] - ref_eigenvalues[i+1]) < NUMERIC_THRESHOLD ]
            if ref_modes[-1] == nrow:
                edge_list += [ncol-1]

            # calculate the distance between degenerate subspaces
            for i in range(len(edge_list)-1):
                begin = edge_list[i]
                end = edge_list[i+1]
                overlap = new_eigenvectors[:,begin:end].T @ ref_eigenvectors[:,begin:end]
                singular_values = np.linalg.svd(overlap)
                assert (max(singular_values) < 1.0 + EIGVEC_THRESHOLD) and (min(singular_values) > 1.0 - EIGVEC_THRESHOLD), \
                    f'''ERROR: degenerate subspace mismatch between output files
                        REF: {ref_eigenvectors[:,begin:end]}
                        NEW: {new_eigenvectors[:,begin:end]}'''

            # clean up the eigenvector data at the end of the analysis
            mode = 'standard'
            ref_modes = []
            new_modes = []
            ref_eigenvalues = []
            new_eigenvalues = []
            ref_eigenvectors = []
            new_eigenvectors = []
        else:
            # accumulate eigenvector information in the as-input order
            index = -len(ref_modes)-1
            ref_eigenvectors += [ float(val) for val in ref_line.split()[index:] ]
            new_eigenvectors += [ float(val) for val in new_line.split()[index:] ]

    # compare words & approximately compare numbers in standard mode
    if mode == 'standard':
        for ref_word, new_word in zip(ref_line.split(), new_line.split()):

            # try to convert to floats
            try:
                ref_float = float(ref_word)
                new_float = float(new_word)

                error = abs(ref_float - new_float)
#                    if math.isfinite(new_number):
                if error > NUMERIC_THRESHOLD:
                    assert 0, f'''ERROR: float mismatch between output files
                                  REF: {ref_line}
                                  NEW: {new_line}'''

            # otherwise assume a word
            except ValueError:
                if ref_word != new_word:
                    assert 0, f'''ERROR: word mismatch between output files
                                  REF: {ref_line}
                                  NEW: {new_line}'''
