This is a minimal wrapper library for packaging the small fraction of the
Intel MKL library that is needed by MOPAC as a self-contained shared library.
The complete set of shared MKL libraries required by MOPAC are quite large,
so we are repackaging a subset of them to facilitate a smaller binary distribution.
This solution is made possible by the position-independent code (PIC) compatibility
of MKL's static libraries, which allows them to be repackaged in this manner.

Conceptually, it should be possible to strip unused functionality directly from a
shared library to reduce its size [https://doi.org/10.1145/3359789.3359823], but this
is not yet a readily available multi-platform capability.
