# MOPAC include files

For OS-level ABI-compatibility, the MOPAC API is bound to C. The data structures and function prototypes of this C-like API
are provided in the `mopac.h` header file. It is possible but inconvenient to use the C-like API in other Fortran software,
so a Fortran convenience wrapper is also provided here. Because it needs to be compiled with the software that is calling the
MOPAC API, the Fortran wrapper source files are being provided with a permissive MIT license to allow for broad use.
