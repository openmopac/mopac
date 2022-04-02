# MOZYME

This directory contains the MOZYME feature of MOPAC, which performs reduced-scaling electronic
structure calculations based directly on localized molecular orbitals (LMOs) instead of canonical
orbitals, which allows for the use of sparse linear algebra rather than dense linear algebra. Many
of the core computational steps are then reduced from the cubic-scaling cost of dense linear algebra
to a linear-scaling cost in the number of atoms, but there are still steps and prefactors that cause
overall super-linear scaling, notably electrostatics are still performed with brute-force pair
summation rather than any fast hierarchical method and the SCF cycle and geometric relaxation cycle
inevitably retains some scaling with system size, depending on the degree of electronic and atomic
polarizability in a system.

This feature of MOPAC was formerly covered by a US patent owned by Fujitsu [5,604,686],
which has now expired.
