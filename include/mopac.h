/* Molecular Orbital PACkage (MOPAC)
 * Copyright (C) 2021, Virginia Polytechnic Institute and State University
 *
 * MOPAC is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MOPAC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Diskless/stateless Application Programming Interface (API) to core MOPAC operations */
#ifndef MOPAC_H
#define MOPAC_H

/*
This C-like MOPAC API consists of 4 data structures - mopac_system, mopac_properties, mopac_state,
and mozyme_state - and 6 main functions - mopac_scf, mopac_relax, mopac_vibe, mozyme_scf, mozyme_relax,
and mozyme_vibe - that are organized by the solver they use - conventional MOPAC or reduced-scaling
MOZYME - and the calculation they perform - a single self-consistent field (SCF) cycle, a geometric
relaxation, or a vibrational calculation. For each function, mopac_system defines the atomistic
system and various solver and model options and mopac_properties contains physical properties
produced by the calculation. All data in mopac_system must be specified on input, while
mopac_properties does not need to be initialized. Memory in mopac_properties is allocated by MOPAC,
and it should be deallocated by using the destroy_mopac_properties function provided by the API.

The mopac_state and mozyme_state structures contain information that specifies the electronic ground
state. They must either be given a trivial initialization (mpack = 0 or numat = 0) to start a
calculation from an atomic guess (MOPAC) or a Lewis-structure guess (MOZYME), or the output state of a
previous calculation can be used to initialize the new calculation with the final state of the previous
calculation, which is useful for running molecular dynamics simulations. The output states should be
deallocated using functions provided by the API - destroy_mopac_state or destroy_mozyme_state. It is
also possible to construct custom initial states, and the memory for these states must be allocated by
using the create_mopac_state or create_mozyme_state API functions with the sizes of all arrays
specified on input (mpack for MOPAC; numat, noccupied, nvirtual, icocc_dim, icvir_dim, cocc_dim, and
cvir_dim for MOZYME).

The API can also be used to run MOPAC calculations from an input file using the run_mopac_from_input
function. This is equivalent to running MOPAC from the command line, except that MOPAC and all of its
dependent math libraries can stay resident in memory between calculations when using API calls, which
may improve performance in high-throughput settings for which disk access is not a bottleneck. In
contrast, the other API functions run MOPAC without any disk access, which may be beneficial for
high-throughput settings with file systems that perform poorly when large numbers of small files are
being read and written simultaneously. If MOPAC encounters an error, then run_mopac_from_input returns
1, otherwise it returns 0 for calculations that have completed successfully.
*/

/* data that defines the atomistic system and MOPAC job options */
struct mopac_system {
  /* number of atoms */
  int natom;
  /* number of atoms that are allowed to move (first natom_move atoms in array) */
  int natom_move;
  /* net charge */
  int charge;
  /* number of spin excitations, floor[(number of alpha electrons)/2 - (number of beta electrons)/2] */
  int spin;
  /* semiempirical model: PM7 = 0, PM6-D3H4 = 1, PM6-ORG = 2, PM6 = 3, AM1 = 4, RM1 = 5 */
  int model;
  /* dielectric constant for COSMO implicit solvent, must be 1 (no solvent) for nlattice > 0 */
  double epsilon;
  /* atomic number of each atom */
  int *atom; /* [natom] */
  /* (x,y,z) coordinates of each atom (Angstroms) */
  double *coord; /* [3*natom] */
  /* number of lattice vectors / translation vectors / periodic dimensions */
  int nlattice;
  /* number of lattice vectors that are allowed to move (first nlattice_move vectors in array) */
  int nlattice_move;
  /* external hydrostatic pressure (Gigapascals) */
  double pressure;
  /* (x,y,z) coordinates of each lattice vectors (Angstroms) */
  double *lattice; /* [3*nlattice] */
  /* relative numerical tolerances (equivalent to GNORM and RELSCF keyword values) */
  double tolerance;
  /* time limit for a MOPAC calculation (seconds) */
  int max_time;
};

/* calculated ground-state properties of an atomistic system and MOPAC job info */
struct mopac_properties {
  /* heat of formation (kcal/mol) */
  double heat;
  /* dipole moment vector (Debye) */
  double dipole[3];
  /* atomic partial charges */
  double *charge; /* [natom] */
  /* (x,y,z) coordinates of each moveable atom (Angstroms) */
  double *coord_update; /* [3*natom_move] */
  /* (x,y,z) heat gradients for each moveable atom (kcal/mol/Angstrom) */
  double *coord_deriv; /* [3*natom_move] */
  /* vibrational frequencies of normal modes (1/cm) */
  double *freq; /* [3*natom_move], NULL if unavailable */
  /* (x,y,z) displacement vectors of normal modes */
  double *disp; /* [3*natom_move,3*natom_move], NULL if unavailable */
  /* bond-order matrix in compressed sparse column (CSC) matrix format (0-based indexing)
   * with insignificant bond orders (<0.01) truncated
   * diagonal matrix entries are atomic valencies */
  /* > first index of each atom in CSC bond-order matrix */
  int *bond_index; /* [natom+1] */
  /* > list of atoms bonded to each atom in CSC format */
  int *bond_atom; /* [bond_index[natom]] */
  /* > bond order of atoms bonded to each atom in CSC format */
  double *bond_order; /* [bond_index[natom]] */
  /* (x,y,z) coordinates of each moveable lattice vectors (Angstroms) */
  double *lattice_update; /* [3*nlattice_move] */
  /* (x,y,z) heat gradients for each moveable lattice vector (kcal/mol/Angstrom) */
  double *lattice_deriv; /* [3*nlattice_move] */
  /* stress tensor (Gigapascals) in Voigt form (xx, yy, zz, yz, xz, xy) */
  double stress[6]; /* 0 if unavailable */
  /* number of MOPAC error messages (negative value indicates that allocation of error_msg failed) */
  int nerror;
  /* text of MOPAC error messages */
  char **error_msg; /* [nerror][*] */
};

/* data that describes the electronic state using standard molecular orbitals */
struct mopac_state {
  /* MOPAC data format is adapted from molkst_C and Common_arrays_C modules */
  /* > number of matrix elements in packed lower triangle matrix format */
  int mpack; /* 0 if state is unavailable */
  /* > flag for unrestricted Hartree-Fock ground state (0 == restricted, 1 == unrestricted) */
  int uhf;
  /* > alpha density matrix */
  double *pa; /* [mpack] */
  /* > beta density matrix */
  double *pb; /* [mpack], NULL if uhf == 0 */
};

/* data that describes the electronic state using localized molecular orbitals */
struct mozyme_state {
  /* MOZYME data format is adapted from molkst_C, Common_arrays_C, and MOZYME_C modules */
  /* > number of real atoms */
  int numat; /* 0 if state is unavailable */
  /* > number of Lewis bonds per real atom */
  int *nbonds; /* [numat] */
  /* > list of Lewis-bonded real atoms for each real atom */
  int *ibonds; /* [9,numat] */
  /* > number of orbitals per real atom */
  int *iorbs; /* [numat] */
  /* > number of occupied molecular orbitals */
  int noccupied;
  /* > number of atoms in each occupied LMO */
  int *ncf; /* [noccupied] */
  /* > number of virtual molecular orbitals */
  int nvirtual;
  /* > number of atoms in each virtual LMO */
  int *nce; /* [nvirtual] */
  /* > size of array icocc */
  int icocc_dim;
  /* > index of each real atom in the occupied LMOs */
  int *icocc; /* [iccoc_dim] */
  /* > size of array icvir */
  int icvir_dim;
  /* > index of each real atom in the virtual LMOs */
  int *icvir; /* [icvir_dim] */
  /* > size of array cocc */
  int cocc_dim;
  /* > atomic orbital coefficients of the occupied LMOs */
  double *cocc; /* [cocc_dim] */
  /* > size of array cvir */
  int cvir_dim;
  /* > atomic orbital coefficients of the virtual LMOs */
  double *cvir; /* [cvir_dim] */
};

/* MOPAC electronic ground state calculation */
void mopac_scf(struct mopac_system *system,
               struct mopac_state *state,
               struct mopac_properties *properties);

/* MOPAC geometry relaxation */
void mopac_relax(struct mopac_system *system,
                 struct mopac_state *state,
                 struct mopac_properties *properties);

/* MOPAC vibrational calculation */
void mopac_vibe(struct mopac_system *system,
                struct mopac_state *state,
                struct mopac_properties *properties);

/* MOZYME electronic ground state calculation */
void mozyme_scf(struct mopac_system *system,
                struct mozyme_state *state,
                struct mopac_properties *properties);

/* MOZYME geometry relaxation */
void mozyme_relax(struct mopac_system *system,
                  struct mozyme_state *state,
                  struct mopac_properties *properties);

/* MOZYME vibrational calculation */
void mozyme_vibe(struct mopac_system *system,
                 struct mozyme_state *state,
                 struct mopac_properties *properties);

/* allocate memory for mopac_state */
void create_mopac_state(struct mopac_state *state);

/* allocate memory for mozyme_state */
void create_mozyme_state(struct mozyme_state *state);

/* deallocate memory in mopac_properties */
void destroy_mopac_properties(struct mopac_properties *properties);

/* deallocate memory in mopac_state */
void destroy_mopac_state(struct mopac_state *state);

/* deallocate memory in mozyme_state */
void destroy_mozyme_state(struct mozyme_state *state);

/* run MOPAC conventionally from an input file */
int run_mopac_from_input(char *path_to_file);

#endif