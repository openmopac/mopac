dense matrix operations

NOTE: remove diag_for_GPU & re-introduce diag, either as a historical method w/ a special keyword, or a new form of pseudo-diagonalization

NOTE: MKL-specific routines have been removed from the code, but it may be worth considering the
reintroduction of commands to set the number of threads:

    num_threads = mkl_get_max_threads()
    call mkl_set_num_threads(num_threads)

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
