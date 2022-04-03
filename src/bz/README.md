# BZ

This directory contains the Windows BZ program for visualizing band structures and Fermi
surfaces from periodic calculations performed by MOPAC. It generates a Brillouin zone by
downfolding the supercells created by MAKPOL to recover the original unit cell. This program
is presently available for Windows only as it uses the QuickWin library for visualization.
There is some overlap between the source code of MOPAC and BZ, but they are being kept
separate for now (merge efforts are welcome).

`blas.F90` is from the BLAS library of basic linear algebra subroutines. `minv.F90` is from
the IBM-SSP library of mathematical and statistical subroutines. `rsp.F90` is derived and
adapted from the EISPACK library of eigenvalue problem solvers.
