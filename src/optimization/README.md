# Geometry optimization

This directory contains various optimizers that are primarily used for geometry optimization. These
are all established optimization algorithms, and many of the implementations in MOPAC were derived
from standard, publicly-available implementations.

`powsq.F90` was originally integrated into the pre-MOPAC MINDO/3 codebase within the Dewar group by
one of its original authors, Andrew Komornicki. `flepo.F90` and `nllsq.F90` seem to have been from
the original MNDO program, but this provenance is unclear. `lbfgs.F90` was derived and adapted from
Algorithm 778 of the ACM Transactions on Mathematical Software [https://doi.org/10.1145/279232.279236].
