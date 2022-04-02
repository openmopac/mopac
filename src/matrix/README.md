# Matrix operations

This directory contains all of the files in MOPAC that have been identified as mostly containing
dense linear algebra operations. They have been centrally collected in preparation for their
consolidation into a more unified subsystem for dense linear algebra within MOPAC. In particular,
the prior support for GPUs in MOPAC, which is not presently functioning, needs to be re-introduced
in a more uniform and maintainable way. Some of the files in this directory may be relocated when
their matrix operations have been appropriately encapsulated.

`ddiag.F90` is from the EISPACK library of eigensolvers. `linpack.F90` is from LINPACK library of
linear solvers. `minv.F90` is from the IBM-SSP library of mathematical and statistical subroutines.
`esp_utilities.F90` contains subroutines from LINPACK and BLAS. `interp.F90` was derived and adapted
from software written by R. Nicholas Camp and Harry F. King.

There are also cases where dense matrix storage and operations are directly embedded in other files
that have not been relocated to this directory. The present list of such subroutines and files is:

- `densf @ polar.F90`
- `bmakuf @ polar.F90`
- `bdenup @ polar.F90`
- `bdenin @ polar.F90`
- `aval @ polar.F90`
- `trugud @ polar.F90`
- `trudgu @ polar.F90`
- `trugdu @ polar.F90`
- `trsub @ polar.F90`
- `transf @ polar.F90`
- `tf @ polar.F90`
- `coscl1 @ cosmo.F90`
- `coscl2 @ cosmo.F90`
- `determinant @ charst.F90`
- `polar @ polar.F90`
