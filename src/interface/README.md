# MOPAC API

This directory contains the application programming interface (API) to MOPAC, which presently has
very limited functionality. MOPAC's primary interface is disk-based rather than API-based, with
input passed in through input files and output read from output files (the AUX file being
machine-readable). The only functionality available right now is access to version information
and an ability to enable/disable expanded disk-based output for GUI support.

To allow MOPAC to function with Windows DLL's, the API subroutines must be marked by `dllexport` preprocessor statements that are used by the Intel Fortran compiler to build a DLL symbol table.
