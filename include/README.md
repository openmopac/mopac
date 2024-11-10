# MOPAC include files

For OS-level ABI-compatibility, the MOPAC API is bound to C. The data structures and function prototypes of this C-like API
are provided in the `mopac_api.h` header file. It is possible but inconvenient to use the C-like API in other Fortran software,
so a Fortran convenience wrapper is also provided here. Because this wrapper needs to be compiled with the software using the
MOPAC API, it is being provided with a permissive MIT license to allow broad use.
