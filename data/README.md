This directory contains the reference data that has been used to train and test the semiempirical models contained within MOPAC.

Each piece of reference data is contained in a MOPAC input file with special syntax on the 3rd line of the input file to
denote what observables have reference data for the given molecule, what their values are, and what the source of the reference
data is. The bibliography for the reference data is contained in `references.txt`. PARAM parses this extra line of input
in addition to the standard first line of MOPAC keywords.

Because most of the reference data is heats of formation and most geometries are not fully constrained by experiment,
the reference geometries are relaxed to the local minimum of the semiempirical model being trained. The initial values
of the reference geometries serve to place the molecules in a particular basin of convergence in the potential energy
surface.
