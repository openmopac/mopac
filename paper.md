---
title: 'MOPAC: An open-source semiempirical molecular orbital program'
tags:
  - Fortran
  - physical chemistry
  - quantum chemistry
  - electronic structure
  - semiempirical models
authors:
  - name: Jonathan E. Moussa
    orcid: 0000-0003-3701-1830
    corresponding: true
    affiliation: 1
  - name: James J. P. Stewart
    affiliation: 2
affiliations:
 - name: Molecular Sciences Software Institute, Virginia Tech, Blacksburg, VA 24060, United States
   index: 1
   ror: 02smfhw86
 - name: Stewart Computational Chemistry, Colorado Springs, CO 80921, United States
   index: 2
date: 21 February 2025
bibliography: paper.bib

---

# Summary

The Molecular Orbital PACkage (MOPAC) is a Fortran program that calculates chemical and physical
properties of molecules, crystals, and nanostructures. MOPAC primarily operates as a command-line
program that takes an input file defining a molecule by its approximate atomic coordinates and
produces an output file containing the molecule's heat of formation and optimized coordinates
alongside other useful properties. Keywords can then be added to the top of the input file to
adjust MOPAC's behavior and request the use of other features. MOPAC is similar in function to
*ab initio* quantum chemistry software, but it uses semiempirical models such as PM7 [@PM7] to
reduce computational costs substantially while preserving as much accuracy as possible.

# Statement of need

In the broader context of chemistry, physics, and materials science, the computer simulation of
atomistic systems is intended as a partner to experiment -- to help interpret experimental results,
to extend their reach beyond what can be observed directly, and even to replace experiments
altogether in some cases. Much of the rationale for this is motivated by the reduction of costs,
that some physical properties ought to be fundamentally cheaper to understand and predict through
computer simulations rather than by experiment. However, an unfortunate reality is that *ab initio*
quantum chemistry has extremely high computational costs that grow very rapidly with molecular
size, which has long been an impediment to this partnership with experiment. Historically,
semiempirical models of thermochemistry were developed to mitigate these costs while preserving
enough quantum mechanical structure to retain a high degree of model transferability.
Semiempirical models of thermochemistry were pioneered by Pople [@Pople1; @Pople2] and later
championed by Dewar [@Dewar], culminating in the development of the Modified Neglect of Diatomic
Overlap (MNDO) model form [@MNDO], the Austin Model 1 (AM1) parameterization [@AM1], and the
MOPAC program [@MOPAC1]. Historically, MOPAC has been the development platform for the MNDO-family
semiempirical models, but these models are also implemented in other software such as
Gaussian [@Gaussian], CP2K [@CP2K], Sparrow [@Sparrow], and ULYSSES [@ULYSSES].

In the four decades since MOPAC was first released, the computing power of personal computers
has increased a million fold. As a result, the practical cost considerations of quantum chemistry
have changed significantly. Practical *ab initio* quantum chemistry calculations were once limited
to Hartree-Fock calculations in a minimal atomic-orbital basis set and small molecules with a few
atoms, but now density functional theory (DFT) calculations in moderately-sized basis sets (split
valence and polarization) of large molecules with hundreds of atoms are routine. Relative to such
routine DFT calculations, semiempirical calculations with MOPAC are roughly a thousand times faster
but half as accurate. Because of its distinct balance between cost and accuracy,
semiempirical quantum chemistry software such as MOPAC is well-suited for simulation tasks that
are highly sensitive to cost such as interactive simulations for chemical exploration and education
[@interactive] and high-throughput virtual screening of molecules [@screening]. More generally,
it is useful to perform calculations at a semiempirical level to estimate results and check for
problems before committing to a much more expensive calculation at an *ab initio* level. For some
expensive applications such as protein modeling, only the semiempirical level may be affordable to
quantum chemistry users who do not have access to supercomputers. Semiempirical quantum mechanical
models continue to be a middle ground in atomistic simulation between *ab initio* quantum mechanics
and molecular mechanics based on classical force fields [@SQM_perspective].

# Features

The original 1990 software paper for MOPAC [@MOPAC2] contains a detailed description of its
features and architecture. At the time of its publication, MOPAC was public-domain software that
was distributed by the Quantum Chemistry Program Exchange (QCPE) [@QCPE]. In 1993, MOPAC was
acquired by Fujitsu Limited and became commercial, closed-source software until its open-source
release in 2022. During its nearly three decades as commercial software and its recent open-source
existence, many new features have been developed for MOPAC. While it was originally developed for
thermochemistry calculations of organic molecules in vacuum containing chemical elements without
valence *d* orbitals, it has since been extended to solids [@solids] and molecules in solution
[@solvent], most elements of the periodic table [@PM6], and more physical properties such as
electronic spectroscopy [@INDO].

A major focus of the last two decades of MOPAC development has been the modeling of proteins and
enzymes. Proteins usually contain hundreds or thousands of atoms, which are challenging to study
with quantum chemistry software because of their rapid growth of costs with molecular size. The
MOZYME solver [@MOZYME] can perform calculations in a localized molecular orbital (LMO) basis,
which has costs that grow only linearly with molecular size. Most protein work is based on
structures deposited in the Protein Data Bank (PDB) [@PDB], and MOPAC can handle PDB structure
files and add missing hydrogen atoms that are not reliably resolved by X-ray crystallography.
To study proton transfer reactions in enzymes, MOZYME is able to locate transition states near
an active site within a larger enzyme structure [@LOCATETS]. While most of the semiempirical
models in MOPAC are intended for general use, it also has a model that is specifically optimized
for biomolecular applications such as proteins and enzymes [@PM6ORG].

As an open-source project, the focus of MOPAC development has shifted to increasing its
accessibility for important use cases. For example, high-throughput applications of MOPAC have
been limited by its need for disk-based input and output on computers with a large number of
compute cores but a slow file system. Since the 23.0.0 release of MOPAC, the most commonly used
functionality of MOPAC is now available through an application programming interface (API) that
avoids all disk usage and captures all essential state information from MOPAC in the input
and output data of its API calls. This API has C bindings to avoid the application binary
interface (ABI) incompatibility problems of Fortran, which allows it to serve as a base API layer
for a Python API wrapper [@mopactools] and increases the accessibility of MOPAC to the Python
software ecosystem.

# Acknowledgements

The Molecular Sciences Software Institute is supported by grant CHE-2136142 from the National Science Foundation.

# References
