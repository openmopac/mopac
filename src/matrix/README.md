dense matrix operations

NOTE: remove diag_for_GPU & re-introduce diag, either as a historical method w/ a special keyword, or a new form of pseudo-diagonalization

NOTE: MKL-specific routines have been removed from the rest of the code, but some remain here.
The code had been forcing the maximum number of threads using:

    call mkl_set_num_threads(mkl_get_max_threads())

but this is not a good idea and has been removed. It prevents single-thread use of MOPAC, and MKL
has a sensible default (1 thread per core) if the number of threads isn't specified by MKL_NUM_THREADS.

NOTE: some files in this folder should be migrated back to other components once their matrix operations have been isolated
and implemented in a separate subroutine.

NOTE: there are subroutines inside larger files that need to be extracted & re-implemented w/ BLAS & LAPACK calls.
They will be listed here as they are identified and before they are fixed:

densf @ polar.F90
bmakuf @ polar.F90
bdenup @ polar.F90
bdenin @ polar.F90
aval @ polar.F90
trugud @ polar.F90
trudgu @ polar.F90
trugdu @ polar.F90
trsub @ polar.F90
transf @ polar.F90
tf @ polar.F90

NOTE: there are large temporary matrix workspaces throughout MOPAC that perhaps should be consolidated somehow.
They will be listed here for consideration:

polar @ polar.F90
