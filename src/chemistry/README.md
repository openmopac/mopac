# Chemical detection & modification

This directory contains subroutines for determining and manipulating Lewis dot structures of
organic molecules (with some support for organometallic molecules). Its functionality includes
an ability to add hydrogen atoms, which are often missing from experimental structures determined
by x-ray crystallography. It also identifies standard functional groups and can label them according
to the PDB file format.

The optimizer of water molecule orientation in `orient_water.F90` is derived and adapted from the Constrained Optimization BY Linear Approximation (COBYLA) algorithm and software originally written
by M. J. D. Powell [https://www.zhangzk.net/software.html].
