# MOPAC source code

This directory contains all Fortran source code associated with MOPAC and accompanying programs
(PARAM, MAKPOL, & BZ). Historically, MOPAC was developed with all source code contained in a single,
unorganized source directory. The organization of source code into subdirectories here is an attempt
to identify and reverse-engineer the substructure of distinct features and contributions to MOPAC,
to make it easier for future open-source developers to understand the codebase.

Most source files are named after the first subroutine in the file, and all data modules are contained
in separate source files with the suffix `*_C.F90`.
