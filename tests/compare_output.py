# Portable Python script for numerical output file comparisons of MOPAC
# Argument list: <output file #1> <output file #2>

from shutil import copyfile
from sys import argv
import subprocess
import os
import re
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
# This comparison is insensitive to differences in whitespace and number of empty lines.
# Some input files used for testing contain reference data in comments, which are ignored here.

# More fine-grained numerical tests using PyTest are planned after development of a Python interface for MOPAC,
# which will make it easier to assign different numerical tolerances to different quantities

NUMERIC_THRESHOLD = 1e-4 # large because of numerical errors in unoccupied orbital energies, energy gradients, & relaxed geometries
HEAT_THRESHOLD = 1e-4
DEGENERACY_THRESHOLD = 1e-3
EIGVEC_THRESHOLD = 1e-3

# regular expression pattern for a time stamp or other signifier of timing output, "CLOCK" or "TIME" or "SECONDS", & system-dependent versioning
skip_criteria = re.compile('([A-Z][a-z][a-z] [A-Z][a-z][a-z] [ 0-9][0-9] [0-9][0-9]:[0-9][0-9]:[0-9][0-9] [0-9][0-9][0-9][0-9])'
                           '|(CLOCK)|(TIME)|(SECONDS)|(Version)|(THE VIBRATIONAL FREQUENCY)|(ITERATION)|(SCF CALCULATIONS)|(Stewart)'
                           '|(remaining)|(\*  THREADS)|(\*  ISOTOPE)|(\*  DENOUT)|(\*  OLDENS)|(\*  SETUP)|(ITER.)|(\*\* )|(web-site)|(MOPAC)|(GRADIENT NORM)')

# regular expression pattern for an eigenvector block
eigen_criteria = re.compile('(Root No.)|(ROOT NO.)')

def is_float(string):
    '''check if a string contains a float'''
    try:
        float(string.replace('D','E'))
        return True
    except ValueError:
        return False

def to_float(string):
    '''check if a string contains a float'''
    try:
        return float(string.replace('D','E'))
    except ValueError:
        return False

def parse_mopac_out_file(path):
    '''parse a MOPAC output file at a given path into a list of basic elements (strings, numbers, matrices)'''
    parse_line = []
    parse_list = []
    mode = 'standard'
    with open(path,'r') as file:
        for line_num, line in enumerate(file):

            # this hack separates floats that don't have a space between them because of a minus sign & trailing comma's
            word_list = line.replace('-',' -').replace('E -','E-').replace('D -','D-').replace('=',' = ').replace(',',' , ').split()

            # ignore molecular dimensions block
            if 'MOLECULAR DIMENSIONS (Angstroms)' in line:
                mode = 'dim'
                continue
            elif mode == 'dim':
                if 'SCF CALCULATIONS' in line:
                    mode = 'standard'
                else:
                    continue

            # switch to or continue iter mode
            if 'RHF CALCULATION' in line or 'UHF CALCULATION' in line or 'Geometry optimization using BFGS' in line:
                mode = 'iter'
                continue
            elif mode == 'iter':
                if 'SCF FIELD WAS ACHIEVED' in line or 'THERE IS NOT ENOUGH TIME FOR ANOTHER CYCLE' in line:
                    mode = 'standard'
                else:
                    continue

            # skip lines as necessary
            if skip_criteria.search(line):
                continue

            # switch to or continue geo mode
            if 'ATOM    CHEMICAL      BOND LENGTH      BOND ANGLE     TWIST ANGLE' in line:
                mode = 'geo'
            elif mode == 'geo':
                if len(word_list) == 0:
                    mode = 'standard'
                else:
                    continue

            # switch to or continue lmo mode
            if 'NUMBER OF CENTERS  LMO ENERGY     COMPOSITION OF ORBITALS' in line:
                mode = 'lmo'
            elif mode == 'lmo':
                if 'LOCALIZED ORBITALS' in line:
                    mode = 'standard'
                else:
                    continue

            # switch to or continue grad mode
            if 'LARGEST ATOMIC GRADIENTS' in line:
                mode = 'grad'
                blank_count = 0
            # simple-minded skipping based on counting blank lines
            elif mode == 'grad':
                if len(word_list) == 0:
                    blank_count += 1
                if blank_count == 3:
                    mode = 'standard'
                else:
                    continue

            # switch to or continue vibe mode
            if 'DESCRIPTION OF VIBRATIONS' in line:
                mode = 'vibe'
            elif mode == 'vibe':
                if 'FORCE CONSTANT IN INTERNAL COORDINATES' in line or 'SYMMETRY NUMBER FOR POINT-GROUP' in line:
                    mode = 'standard'
                else:
                    continue

            # switch to or continue eigen mode
            if eigen_criteria.search(line):
                if mode != 'eigen':
                    eigen_line_num = line_num+1
                    mode = 'eigen'
                    label_list = []
                    value_list = []
                    vector_list = []
                    num_eigen = []
                label_list += [ int(word) for word in word_list[2:] ]
                num_eigen.append(len(word_list) - 2)

            # eigen parsing
            elif mode == 'eigen':

                # save eigenvalues in a list
                if len(word_list) == num_eigen[-1] and len(value_list) < len(label_list):

                    # check if the list of numbers is just another label
                    label_check = True
                    try:
                        for word,label in zip(word_list,label_list[-len(word_list):]):
                            if int(word) != label:
                                label_check = False
                    except ValueError:
                        label_check = False

                    if label_check == False:
                        value_list += [ float(word) for word in word_list ]

                # ignore symmetry labels
                elif len(word_list) == 2*num_eigen[-1] and is_float(word_list[-2]) and not is_float(word_list[-1]):
                    pass

                # save eigenvectors in a matrix
                elif len(word_list) > num_eigen[-1] and all([is_float(word) for word in word_list[-num_eigen[-1]:]]):
                    vector_list += [ float(word) for word in word_list[-num_eigen[-1]:] ]

                # ignore blank lines
                elif len(word_list) == 0:
                    pass

                # switch back to standard mode & reformat eigenvectors
                else:
                    mode = 'standard'

                    # reshape into a matrix
                    nrow = len(vector_list) // len(label_list)
                    ncol = len(label_list)
                    eigenmatrix = np.empty((nrow,ncol))

                    offset = 0
                    for num in num_eigen:
                        eigenmatrix[:,offset:offset+num] = np.reshape(vector_list[offset*nrow:(offset+num)*nrow],(nrow,num),order='C')
                        offset += num

                    # renormalize the eigenvectors (MOPAC uses a variety of normalizations)
                    for col in eigenmatrix.T:
                        col /= np.linalg.norm(col)

                    # output eigenvalue (if known) and eigenvectors
                    if len(value_list) == len(label_list):
                        parse_list.append((value_list,eigenmatrix,label_list[0] == 1,label_list[-1] == nrow))
                    else:
                        parse_list.append((label_list,eigenmatrix,label_list[0] == 1,label_list[-1] == nrow))
                    parse_line.append(eigen_line_num)

            # standard parsing
            if mode == 'standard':
                for word in word_list:
                    if is_float(word):
                        if 'FINAL HEAT OF FORMATION =' in line and word is word_list[5]:
                            parse_list.append(('HOF',to_float(word)))
                        else:
                            parse_list.append(to_float(word))
                    else:
                        parse_list.append(word)
                    parse_line.append(line_num+1)

    return parse_line, parse_list

def compare_mopac_out_file(out_line, out_list, ref_line, ref_list, heat_error_threshold):
    '''Compares the output to the given reference'''

    if len(ref_list) != len(out_list):
        assert len(ref_list) == len(out_list), f'ERROR: output file size mismatch, {len(ref_list)} vs. {len(out_list)}'
        #print(f'WARNING: output file size mismatch, {len(ref_list)} vs. {len(out_list)}')

    for (out_line0, out, ref_line0, ref) in zip(out_line, out_list, ref_line, ref_list):
    #    print(ref, "vs.", out)
        # check that types match
        assert type(ref) == type(out), f'ERROR: type mismatch between {ref} on reference line {ref_line0} and {out} on output line {out_line0}'

        # compare strings
        if type(ref) is str:
            assert ref == out, f'ERROR: string mismatch between {ref} on reference line {ref_line0} and {out} on output line {out_line0}'

        # compare floats
        elif type(ref) is float:
    #        assert abs(ref - out) < NUMERIC_THRESHOLD, f'ERROR: numerical mismatch between {ref} and {out} on output line {line}'
            if abs(ref - out) > NUMERIC_THRESHOLD:
                print(f'WARNING: numerical mismatch between {ref} on reference line {ref_line0} and {out} on output line {out_line0}')

        # compare heats of formation
        elif len(ref) == 2:
            assert abs(ref[1] - out[1]) < heat_error_threshold, f'ERROR: numerical heat mismatch between {ref[1]} on reference line {ref_line0} and {out[1]} on output line {out_line0}'
            if abs(ref[1] - out[1]) > HEAT_THRESHOLD:
                print(f'WARNING: numerical heat mismatch between {ref[1]} on reference line {ref_line0} and {out[1]} on output line {out_line0}')
 
        # compare eigenvalues & eigenvectors
        elif len(ref) == 4:
            ref_val, ref_vec, ref_begin, ref_end = ref
            out_val, out_vec, ref_begin, ref_end = out

            for refv, outv in zip(ref_val,out_val):
    #            assert abs(refv - outv) < NUMERIC_THRESHOLD, f'ERROR: numerical mismatch between {refv} and {outv} on output line {line}'
                if abs(refv - outv) > NUMERIC_THRESHOLD:
                    print(f'WARNING: eigenvalue mismatch between {refv} on reference line {ref_line0} and {outv} on output line {out_line0}')

                # build list of edges denoting degenerate subspaces
                if ref_begin:
                    edge_list = [0]
                else:
                    edge_list = []
                edge_list += [ i+1 for i in range(len(ref_val)-1) if np.abs(ref_val[i] - ref_val[i+1]) > DEGENERACY_THRESHOLD ]
                if ref_end:
                    edge_list += [len(ref_val)]

                # test the distance between each pair of degenerate subspaces
                for i in range(len(edge_list)-1):
                    overlap = ref_vec[:,edge_list[i]:edge_list[i+1]].T @ out_vec[:,edge_list[i]:edge_list[i+1]]
    #                print("overlap = ",overlap)
                    sval = np.linalg.svd(overlap, compute_uv=False)
                    assert (sval[0] < 1.0 + EIGVEC_THRESHOLD) and (sval[-1] > 1.0 - EIGVEC_THRESHOLD), \
                        f'ERROR: degenerate subspace mismatch between reference line {ref_line0} and output line {out_line0}, overlap range in [{min(sval)},{max(sval)}]'

# stub main for command-line comparisons
if __name__ == "__main__":

    # parse the 2 output files that we are comparing
    ref_line, ref_list = parse_mopac_out_file(argv[1])
    out_line, out_list = parse_mopac_out_file(argv[2])

    # Run the comparison
    compare_mopac_out_file(out_line, out_list, ref_line, ref_list)
