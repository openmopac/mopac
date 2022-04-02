==================
MOPAC Contributors
==================

The development of MOPAC began in 1981 within the Dewar group at the University of Texas at Austin.
Prof. Michael Dewar and his group had been developing semiempirical models since the late 1960's,
and this activity had produced a substantial amount of internal software as well as two programs,
MINDO/3 and MNDO, released through the Quantum Chemistry Program Exchange (QCPE). The main developer
of MOPAC, James (Jimmy) J. P. Stewart, was visting the Dewar group on sabbatical leave from the
University of Strathclyde and took on the task of refactoring and consolidating this collection of
software into a more cohesive and user-friendly program with a unified system of input and output
files that was released as QCPE Program #455 in 1983.

The overwhelming majority of development and maintenance of MOPAC from 1981-2021 has been performed
by Jimmy Stewart. The longstanding policy of MOPAC has been to make its source code freely available
to interested developers and to accept the donation of new features, provided either in a modified
version of MOPAC or a self-contained reference program to be integrated into MOPAC. One of the
stipulations of these donations was the prior academic publication of the new feature, and the
historic contributions listed below include a DOI link to the relevant publications in these cases.

As MOPAC has an old codebase, attribution will inevitably be imperfect and have error bars. This
contributor list is an attempt at being maximally inclusive, and requests for corrections are welcome.
From a minimally inclusive perspective, Jimmy Stewart was the pre-open-source copyright owner of
MOPAC (donations to the code were made through a copyright transfer agreement), wrote the majority
of the code, and has borne all responsibility for maintenance, support, and promoting the interests
of MOPAC. Between this minimum and maximum, an earnest attempt at a formal author list is contained
in `CITATION.cff`, which signifies major contributions that persistent in the present MOPAC codebase
that are not part of some other original work of which MOPAC is a derivative work.

Pre-MOPAC Contributors
======================

MINDO/3 Authors
---------------

The 1975 methodology paper [`DOI:10.1021/ja00839a001 <https://doi.org/10.1021/ja00839a001>`_]
was authored by Richard C. Bingham, Michael J. S. Dewar, and Donald H. Lo.
The MINDO/3 program [QCPE Program #279 (1975)] authors according to the QCPE listing were
M. J. S. Dewar, H. Metiu, P. J. Student, A. Brown, R. C. Bingham, D. H. Lo, C. A. Ramsden,
H. Kollmar, P. Werner, and P. K. Bischof.

MNDO Authors
------------

The 1977 methodology paper [`DOI:10.1021/ja00457a004 <https://doi.org/10.1021/ja00457a004>`_]
was authored by Michael J. S. Dewar and Walter Thiel.
The MNDO program [QCPE Program #428 (1981)] authors according to the QCPE listing were
W. Thiel, P. Weiner, J. Stewart, and M. J. S. Dewar.

*NOTE: As with much of the software previously distributed through the QCPE, the original MINDO/3
and MNDO programs are no longer distributed and are likely to be permanently lost. If anyone has
a copy of either of these programs, then please consider sharing them on GitHub for the purpose
of historical preservation and better understanding of the origins of MOPAC.*

Historic MOPAC Contributors
===========================

Andrew Komornicki
   adapted the POWSQ geometry optimizer
   [`DOI:10.1021/ja00763a011 <https://doi.org/10.1021/ja00763a011>`_]
   to the MINDO/3 program after its QPCE release
   [`DOI:10.1021/ja00444a012 <https://doi.org/10.1021/ja00444a012>`_]
   but before the development of MOPAC itself

Santiago Olivella
   semiempirical energy partitioning (ENPART subroutine)
   [`DOI:10.1002/jhet.5570180625 <https://doi.org/10.1002/jhet.5570180625>`_]

Peter Pulay
   design & optimization of pseudodiagonalization
   [`DOI:10.1002/jcc.540030214 <https://doi.org/10.1002/jcc.540030214>`_]

Harry King & R. Nicholas Camp
   original implementation of the Camp-King SCF converger
   [`DOI:10.1063/1.441834 <https://doi.org/10.1063/1.441834>`_]

John McKelvey
   adaptation of the Camp-King converger for MOPAC & improved output formatting

Roger Sargent, Dimitris Agrafiotis, & Henry Rzepa
   implementation of the Broyden-Fletcher-Goldfarb-Shanno (BFGS) optimizer.

Larry Davis & Larry Burggraf
   design of the dynamic reaction coordinate (DRC) & intrinsic reaction coordinate (IRC) features
   [`DOI:10.1002/jcc.540080808 <https://doi.org/10.1002/jcc.540080808>`_]

Frank Jensen
   efficiency improvements to the Eigenvector-Following (EF) method
   [`DOI:10.1063/1.469144 <https://doi.org/10.1063/1.469144>`_]

Juan Carlos Paniagua
   improvements to the orbital localization procedure
   [`DOI:10.1002/qua.560260307 <https://doi.org/10.1002/qua.560260307>`_]

Jorge Medrano
   expanded bonding analysis
   [`DOI:10.1002/jcc.540060205 <https://doi.org/10.1002/jcc.540060205>`_]

Tsuneo Hirano
   revision of energy partitioning & thermodynamic corrections

James Friedheim
   testing & bug hunting

Eamonn Healy
   testing & feature validation

James Ritchie
   bug fixes (SCF restarting)

Masamoto Togashi, Jerzy Rudzinski, Zdenek Slanina, & Eiji Osawa
   bug fixes (vibrational analysis)

Michael Frisch
   bug fixes (density matrix)

Patrick Redington
   bug fixes (heavy atom matrix elements)

Ernest Davidson
   improvements to the 2-electron matrix elements

Daniel Liotard
   partial analytical derivatives of the density matrix & 2-electron matrix elements
   [`DOI:10.1016/0166-1280(90)85012-C <https://doi.org/10.1016/0166-1280(90)85012-C>`]

Yukio Yamaguchi
   partial support for analytical derivatives
   [`DOI:10.1016/0097-8485(78)80005-9 <https://doi.org/10.1016/0097-8485(78)80005-9>`_]

George Purvis III
   expanded STO-6G orbital implementation up to principal quantum number 6
   for use in analytical derivatives

Henry Kurtz
   implementation of polarizability and hyperpolarizability
   [`DOI:10.1002/jcc.540110110 <https://doi.org/10.1002/jcc.540110110>`_]

Prakashan Korambath
   frequency dependence of hyperpolarizability
   [`DOI:10.1021/bk-1996-0628.ch007 <https://doi.org/10.1021/bk-1996-0628.ch007>`_]

David Danovich
   implementation of point-group symmetry & Green's function corrections to ionization potentials
   [`DOI:10.1039/P29930000321 <https://doi.org/10.1039/P29930000321>`_]

Michael Coolidge
   use of symmetry to accelerate vibrational analysis
   [`DOI:10.1002/jcc.540120807 <https://doi.org/10.1002/jcc.540120807>`_]

Andreas Klamt
   implementation of the COSMO solvation model
   [`DOI:10.1039/P29930000799 <https://doi.org/10.1039/P29930000799>`_]

Anna Stewart
   copyediting of MOPAC documentation

Victor Danilov
   edited the MOPAC7 manual & identified bugs in the MECI feature

John Simmie
   conversion of the MOPAC7 manual to LaTeX

Walter Thiel & Alexander Voityuk
   reference implementation of semiempirical models with d orbitals
   [`DOI:10.1007/BF01134863 <https://doi.org/10.1007/BF01134863>`_]

Kenneth Merz, Jr.
   implementation of atomic charge model for electrostatic potentials (ESP)
   [`DOI:10.1002/jcc.540110404 <https://doi.org/10.1002/jcc.540110404>`_]

Bingze Wang
   implementation of parametric electrostatic potentials (PMEP)
   [`DOI:10.1002/jcc.540150210 <https://doi.org/10.1002/jcc.540150210>`_]

Stephan Grimme
   reference implementation of the D3 dispersion model
   [`DOI:10.1063/1.3382344 <https://doi.org/10.1063/1.3382344>`_]

Jan Rezac
   expanded implementation of classical energy corrections (hydrogen bonding, halogen bonding, dispersion)
   [`DOI:10.1021/ct200751e <https://doi.org/10.1021/ct200751e>`_]

Gerd Rocha
   expanded BLAS/LAPACK support, Intel MKL for multi-threading, & cuBLAS/MAGMA for GPU acceleration
   [`DOI:10.1021/ct3004645 <https://doi.org/10.1021/ct3004645>`_]

Rebecca Gieseking
   implementation of the INDO/S spectroscopy model
   [`DOI:10.1002/jcc.26455 <https://doi.org/10.1002/jcc.26455>`_]

Jonathan Moussa
   open-source transition: reorganization & clean-up of the codebase, portability testing & debugging,
   minor performance tuning, transition to CMake-based build system, automation of continuous integration & deployment

Open-Source MOPAC Contributors
==============================

One of the major benefits of modern open-source software development is that contributions are passively recorded by
git version control. As such, refer to the commit records for a complete list of contributions. Major new feature
contributions will continue to be added to this list over time.
