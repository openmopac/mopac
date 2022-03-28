! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module Param_global_C
!
!  All the data that are global to PARAM, and are not used at all in MOPAC
!
    implicit none
                                            !       Maximum number of
!
    integer, parameter :: maxmol = 10000    !          molecules
    integer, parameter :: maxpms = 4000     !          parameter
    integer, parameter :: maxfns = 15000    !          reference data
    integer, parameter :: maxrfs = 4000     !          journal references
    integer, parameter :: maxatm = 190000   !          atoms in all systems
    integer, parameter :: maxsym = 20000    !          symmetry functions
!
!  Dependent parameters
!
    integer, parameter :: maxpat = 3*maxatm !          geometric parameters
!
!  Integer scalars
!
    logical ::  large
!
!
!
    integer ::    &
    nmols,        & ! Number of molecules in PARAM
    atom_pKa,     & ! Atom number of hydrogen atom used in pKa work
    atom_pKa_O,   & ! Atom number of oxygen atom that the hydrogen atom in pKa work is attached to
    nfns,         & ! Number of reference data in current run
    numvar,       & ! Number of parameters in PARAM
    ifiles_8,     & ! output channel for PARAM
    molnum,       & ! Molecule number of the current molecule
    nref,         & ! Number of reference directories to be used
    maxpab,       & ! Number of density matrix elements
    aaaint
!
!  Real scalars
!
    double precision :: &    ! Reference value of
      refdip,     & ! Reference value of Dipole
      refhof,     & ! Reference value of Heat of formation
      refips,     & ! Reference value of Ionization potential
      wtdip,      & ! Weighting factor for dipole
      wtgeo,      & ! Weighting factor for geometry
      wthof,      & ! Weighting factor for heat of formation
      wtips,      & ! Weighting factor for ionization potential
      wtz,        & ! Number of mers in a cluster in a solid
      refpKa,     & ! Reference value of pKa
      wtpKa,      & ! Weighting factor for pKa
      refcer,     & ! Error in calculated ionization potential
      refder,     & ! Error in calculated dipole moment
      refger,     & ! Error in calculated geometry
      refher,     & ! Error in calculated heat of formation
      reftot,     & ! Total error for current system
      parab,      & ! Restraining force constant (on values of parameters)
      power,      & ! Power for error function (normally 2)
      penalty,    & ! Penalty for exceeding parameter boundaries = penalty*error**2
    phase = 1.d0, & ! Alternate step directions when building derivatives to minimize bias
      xxxxxx
!
!  Character scalars
!
      character (len=2000) :: contrl = " "  ! Keywords for PARAM
!

      character (len=40) :: molnam         ! Names of file of current system
!
!   Arrays
!
    logical :: &
    is_a_ref(maxmol), &   ! "true" if the H.o.F. reference datum is used only as a reference to
                          ! another reference H.O.F.
    save_parameters,  &   ! "true" if parameters are to be saved
    locked(maxpms),   &   ! "false" if parameter is selected for optimization using "(all)"
    dorbs(107)
!
!
!
    integer :: &
    ihrefs(6, maxmol),      & ! Number of, and indices for, all I.P. levels for all molecules
    molele(10, maxmol),     & ! Number of elements in, and atomic numbers of elements, in each molecule
    locvar(2, maxpms),      & ! (1,*) : Element number for parameter, (2,*) : parameter number
    qn(maxmol),             & ! Quantum number of electronic state use in "root=" option
    spin_state(maxmol),     & ! Spin-state of electronic state use in "root=" option
    lions(maxmol),          & ! Eigenvalue number for ionization potential
    Atom_pKas(maxmol,2),    & ! Atom numbers for O and H in oKa calculation, for all molecules
    natmss(maxmol),         & ! Number of atoms (real plus others) in each molecule
    numats(maxmol),         & ! Number of real atoms in each molecule
    lablss(maxatm),         & ! Labels for every atom in all molecules
    nas(maxatm),            & ! Bond-length connectivity for all atoms in all molecules
    nbs(maxatm),            & ! Bond-angle connectivity for all atoms in all molecules
    ncs(maxatm),            & ! Dihedral connectivity for all atoms in all molecules
    nats(maxatm),           & ! Atomic numbers of all atoms in all molecules
    nfirss(maxatm),         & ! NFIRST for all atoms in all molecules (see MOPAC)
    nlasts(maxatm),         & ! NLAST for all atoms in all molecules (see MOPAC)
    ncloss(maxmol),         & ! NCLOSE for all molecules (see MOPAC)
    nlecss(maxmol),         & ! NELECS for all molecules (see MOPAC)
    nlm61(maxmol),          & ! LM61 for all molecules (see MOPAC)
    nnmos(maxmol),          & ! NMOS for all molecules (see MOPAC)
    nopens(maxmol),         & ! NOPEN for all molecules (see MOPAC)
    norbss(maxmol),         & ! NORBS for all molecules (see MOPAC)
    nnalpha_open(maxmol),   &
    nnbeta_open(maxmol),    &
    n2elecs(maxmol),        & ! N2ELEC for all molecules (see MOPAC)
    ndeps(maxmol),          & ! Number of symmetry-dependent functions for all molecules
    idepfs(maxsym),         & ! IDEP symmetry functions for all molecules
    locdes(maxsym),         & ! LOCDEP symmetry functions for all molecules
    locpas(maxsym),         & ! LOCDER symmetry functions for all molecules
    nvars(maxmol),          & ! Number of geometric variables in all molecules
    locs(2, maxpat),        & ! Indices of all geometric variables in all molecules
    msdels(maxmol),         & ! Magnetic component of spin for all molecules
    l123s(3,maxmol),        & ! Number of unit cells (for solids)
    asasasa
!
!
!
    double precision :: &
    geos(3, maxatm),    & ! All geometries for all systems
    dipls(maxmol),      & ! Reference dipoles for all systems
    heats(maxmol),      & ! Reference H.o.F. for all systems
    hions(maxmol),      & ! Reference ionization potentials for all systems
    pKas(maxmol),       & ! Reference pKa for all systems
    weight(6, maxmol),  & ! The weights for all reference data for all molecules
    xparamp(maxpms),    & ! Values of all parameters used in PARAM
    hofcal(maxmol),     & ! Calculated H.o.F. for all molecules
    fns(maxfns),        & ! Calculated values of reference data
    fnsnew(maxfns),     & ! Predicted values of reference data (using RAPID)
    error(maxfns),      & ! Errors in calculated and predicted values of reference data
    factor(maxfns),     & ! Weighting factors, used in reducing significance of inaccurate data
    fracts(maxmol),     & ! FRACT for all molecules (see MOPAC)
    depmuls(maxsym),    & ! DEPMUL symmetry functions for all molecules
    botlim(maxpms),     & ! Bottom limit for all parameters
    toplim(maxpms),     & ! Top limit for all parameters
    valold(maxpms),     & ! Values of all parameters in previous cycle
    valvar(maxpms),     & ! Values of all parameters
    asaaaa
!
!
!
    character (len=12) :: &
    geotxt(300, maxmol),  & ! Descriptions of all geometric reference data in all molecules
    refgeo(300)             ! Descriptions of each geometric reference datum in a molecule
    character (len=600) :: &
    refdir(20)              ! Names of reference directories to be used
    character (len=80) :: &
    names(maxmol)           ! File-names of reference data, one per system
    character (len=63) :: &
    tfns(maxfns)            ! Text describing each reference datum (mol. name + parameter type)
    character (len=8)  :: &
    refers(maxmol,7)        ! Reference mnemonics for all reference data in all molecules
    character (len=4)  :: &
    i_r(maxmol)             ! Irreducible representation of electronic state use in "root=" option
    character (len=241):: &
    keys(maxmol),         & ! Keywords for all molecules
    titls(maxmol)           ! Titles for all molecules
    character (len=80) :: &
    comnts(maxmol)          ! Comments for all molecules
    character (len=120) :: &
    source(3, 400)          ! Original references for all reference data


!
!
!
    double precision, dimension (:,:), allocatable :: &
    diffns     !  Differentials of all reference data by all parameters
    double precision, dimension (:), allocatable :: &
    pas,          & ! Alpha density matrix for all molecules
    pbs,          & ! Beta density matrix for all molecules
    asasa
!
    save
end module Param_global_C
