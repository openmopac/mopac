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

module atomradii_C
    implicit none
    !
    ! Tables of atomic radii
    !
    ! A value of 999.d0 always means that no appropriate datum is available
    !
    !.. Global Arrays ..
    integer, parameter :: periodic_table_length = 107
    logical, dimension (periodic_table_length) :: is_metal
    double precision, dimension (periodic_table_length) :: atom_radius_vdw
    double precision, dimension (periodic_table_length) :: atom_radius_covalent
    double precision, dimension (periodic_table_length) :: atom_radius_cosmo
    double precision, dimension (periodic_table_length) :: atom_radius_cosmo_oldcav
    double precision, dimension (periodic_table_length) :: atom_radius_estat
    double precision, allocatable :: radius(:)  !  This will hold the atomic radii for each atom
    !
    !.. Data Declarations ..
    !
    ! ********* van der Waals radii *********
    ! Values taken from
    !  A. Bondi, J.Phys.Chem., 68 (1964) 441.
    ! Where these are missing, taken from
    ! table of Andrey Bliznyuk
    ! ***************************************
    data atom_radius_vdw / &
    !      H       He
         & 1.20d0, 1.40d0, &
    !
    !      Li      Be      B       C       N       O       F       Ne
         & 1.82d0, 1.87d0, 1.69d0, 1.70d0, 1.55d0, 1.52d0, 1.47d0, 1.54d0, &
    !
    !      Na      Mg      Al      Si      P       S       Cl      Ar
         & 2.27d0, 1.73d0, 2.06d0, 2.10d0, 1.80d0, 1.80d0, 1.75d0, 1.88d0, &
    !
    !      K       Ca      Sc      Ti      V       Cr      Mn      Fe
         & 2.75d0, 2.17d0, 2.26d0, 2.26d0, 2.15d0, 2.05d0, 2.10d0, 2.06d0, &
    !      Co      Ni      Cu      Zn      Ga      Ge      As      Se
         & 2.05d0, 1.63d0, 1.40d0, 1.39d0, 1.87d0, 2.10d0, 1.85d0, 1.90d0, &
    !      Br      Kr
         & 1.85d0, 2.02d0, &
    !
    !      Rb      Sr      Y       Zr      Nb      Mo      Tc      Ru
         & 3.23d0, 2.94d0, 2.90d0, 2.85d0, 2.80d0, 2.20d0, 2.20d0, 2.20d0, &
    !      Rh      Pd      Ag      Cd      In      Sn      Sb      Te
         & 2.20d0, 1.63d0, 1.72d0, 1.58d0, 1.93d0, 2.17d0, 2.16d0, 2.06d0, &
    !      I       Xe
         & 1.98d0, 2.16d0, &
    !
    !      Cs      Ba      La      Ce      Pr      Nd      Pm      Sm
         & 3.42d0, 2.97d0, 2.40d0, 2.40d0, 2.40d0, 2.40d0, 2.40d0, 2.40d0, &
    !      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
         & 2.40d0, 2.40d0, 2.40d0, 2.40d0, 2.40d0, 2.40d0, 2.40d0, 2.40d0, &
    !      Lu      Hf      Ta      W       Re      Os      Ir      Pt
         & 2.40d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, 1.75d0, &
    !      Au      Hg      Tl      Pb      Bi      Po      At      Rn
         & 1.66d0, 1.55d0, 1.96d0, 2.02d0, 2.26d0, 2.26d0, 2.25d0, 2.30d0, &
    !
    !      Fr      Ra      Ac      Th      Pa      U       Np      Pu
         & 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, 1.86d0, 2.20d0, 2.20d0, &
    !      Am      Cm      Bk      Cf      Es      Fm      Md      No
         & 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0, &
    !      Lr      Rf      Db      Sg      Tv
         & 2.20d0, 2.20d0, 2.20d0, 2.20d0, 2.20d0 /
    !
    ! ********* Covalent radii *********
    ! Values taken from
    !  http://www.webelements.com
    ! **********************************
    data atom_radius_covalent / &
    !      H       He
         & 0.37d0, 0.32d0, &
    !
    !      Li      Be      B       C       N       O       F       Ne
         & 1.34d0, 0.90d0, 0.82d0, 0.77d0, 0.75d0, 0.73d0, 0.71d0, 0.69d0, &
    !
    !      Na      Mg      Al      Si      P       S       Cl      Ar
         & 1.54d0, 1.30d0, 1.18d0, 1.11d0, 1.06d0, 1.02d0, 0.99d0, 0.97d0, &
    !
    !      K       Ca      Sc      Ti      V       Cr      Mn      Fe
         & 1.96d0, 1.74d0, 1.44d0, 1.36d0, 1.25d0, 1.27d0, 1.39d0, 1.25d0, &
    !      Co      Ni      Cu      Zn      Ga      Ge      As      Se
         & 1.26d0, 1.21d0, 1.38d0, 1.31d0, 1.26d0, 1.22d0, 1.19d0, 1.16d0, &
    !      Br      Kr
         & 1.14d0, 1.10d0, &
    !
    !      Rb      Sr      Y       Zr      Nb      Mo      Tc      Ru
         & 2.11d0, 1.92d0, 1.62d0, 1.48d0, 1.37d0, 1.45d0, 1.56d0, 1.26d0, &
    !      Rh      Pd      Ag      Cd      In      Sn      Sb      Te
         & 1.35d0, 1.31d0, 1.53d0, 1.48d0, 1.44d0, 1.41d0, 1.38d0, 1.35d0, &
    !      I       Xe
         & 1.33d0, 1.30d0, &
    !
    !      Cs      Ba      La      Ce      Pr      Nd      Pm      Sm
         & 2.25d0, 1.98d0, 1.69d0, 1.d0  , 1.d0  , 1.d0  , 1.d0  , 1.d0  , &
    !      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
         & 1.d0  , 1.d0  , 1.d0  , 1.d0  , 1.d0  , 1.d0  , 1.d0  , 1.d0   , &
    !      Lu      Hf      Ta      W       Re      Os      Ir      Pt
         & 1.60d0, 1.50d0, 1.38d0, 1.46d0, 1.59d0, 1.28d0, 1.37d0, 1.28d0, &
    !      Au      Hg      Tl      Pb      Bi      Po      At      Rn
         & 1.44d0, 1.49d0, 1.48d0, 1.47d0, 1.46d0, 999.d0, 999.d0, 1.45d0, &
    !
    !      Fr      Ra      Ac      Th      Pa      U       Np      Pu
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 2.d0,   999.d0, &
    !      Am      Cm      Bk      Cf      Es      Fm      Md      No
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 2.d0,   2.d0,   2.d0,   &
    !      Lr      Rf      Db      Sg      Tv
         & 2.d0,   2.d0,   2.d0,   2.d0, 2.d0 /
    !
    ! ********* Radii for cosmo solvation model *********
    ! These are the original Mopac values.  They are used when the old
    ! cavity option, "OLDCAV", is used
    ! ***************************************************
    data atom_radius_cosmo_oldcav / &
    !      H       He
         & 1.08d0, 1.00d0, &
    !
    !      Li      Be      B       C       N       O       F       Ne
         & 1.80d0, 999.d0, 999.d0, 1.53d0, 1.48d0, 1.36d0, 1.30d0, 999.d0, &
    !
    !      Na      Mg      Al      Si      P       S       Cl      Ar
         & 2.30d0, 999.d0, 2.05d0, 2.10d0, 1.75d0, 1.70d0, 1.65d0, 999.d0, &
    !
    !      K       Ca      Sc      Ti      V       Cr      Mn      Fe
         & 2.80d0, 2.75d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Co      Ni      Cu      Zn      Ga      Ge      As      Se
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Br      Kr
         & 1.80d0, 999.d0, &
    !
    !      Rb      Sr      Y       Zr      Nb      Mo      Tc      Ru
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Rh      Pd      Ag      Cd      In      Sn      Sb      Te
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      I       Xe
         & 2.05d0, 999.d0, &
    !
    !      Cs      Ba      La      Ce      Pr      Nd      Pm      Sm
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lu      Hf      Ta      W       Re      Os      Ir      Pt
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Au      Hg      Tl      Pb      Bi      Po      At      Rn
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !
    !      Fr      Ra      Ac      Th      Pa      U       Np      Pu
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Am      Cm      Bk      Cf      Es      Fm      Md      No
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lr      Rf      Db      Sg      Tv
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0 /
    ! ********* Radii for cosmo solvation model *********
    ! These are the VDW radii for use with the new COSMO
    ! cavity.
    ! ***************************************************
           data atom_radius_cosmo / &
    !      H       He
         & 1.3d0, 999.d0, &
    !
    !      Li      Be      B       C       N       O       F       Ne
         & 999.d0, 999.d0, 2.0475d0,2.0d0, 1.83d0, 1.72d0, 1.72d0, 999.d0, &
    !
    !      Na      Mg      Al      Si      P       S       Cl      Ar
         & 999.d0, 999.d0, 999.d0, 2.457d0,2.106d0,2.16d0, 2.05d0, 999.d0, &
    !
    !      K       Ca      Sc      Ti      V       Cr      Mn      Fe
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Co      Ni      Cu      Zn      Ga      Ge      As      Se
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 2.223d0,2.223d0, &
    !      Br      Kr
         & 2.16d0, 999.d0, &
    !
    !      Rb      Sr      Y       Zr      Nb      Mo      Tc      Ru
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Rh      Pd      Ag      Cd      In      Sn      Sb      Te
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      I       Xe
         & 2.32d0, 999.d0, &
    !
    !      Cs      Ba      La      Ce      Pr      Nd      Pm      Sm
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lu      Hf      Ta      W       Re      Os      Ir      Pt
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Au      Hg      Tl      Pb      Bi      Po      At      Rn
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !
    !      Fr      Ra      Ac      Th      Pa      U       Np      Pu
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Am      Cm      Bk      Cf      Es      Fm      Md      No
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lr      Rf      Db      Sg      Tv
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0 /
    !
    !
    ! ********* Whether the element should be *********
    !           treated as a metal in determining
    !           bonding.
    ! ***************************************************
    data is_metal / &
    !      H       He
         & .false.,.false., &
    !
    !      Li      Be      B       C       N       O       F       Ne
         & .true. ,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
    !
    !      Na      Mg      Al      Si      P       S       Cl      Ar
         & .true. ,.true. ,.false.,.false.,.false.,.false.,.false.,.false., &
    !
    !      K       Ca      Sc      Ti      V       Cr      Mn      Fe
         & .true. ,.true. ,.false.,.false.,.false.,.false.,.false.,.false., &
    !      Co      Ni      Cu      Zn      Ga      Ge      As      Se
         & .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
    !      Br      Kr
         & .false.,.false., &
    !
    !      Rb      Sr      Y       Zr      Nb      Mo      Tc      Ru
         & .true.,.true.,.false.,.false.,.false.,.false.,.false.,.false., &
    !      Rh      Pd      Ag      Cd      In      Sn      Sb      Te
         & .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
    !      I       Xe
         & .false.,.false., &
    !
    !      Cs      Ba      La      Ce      Pr      Nd      Pm      Sm
         & .true.,.true.,.false.,.false.,.false.,.false.,.false.,.false., &
    !      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
         & .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
    !      Lu      Hf      Ta      W       Re      Os      Ir      Pt
         & .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
    !      Au      Hg      Tl      Pb      Bi      Po      At      Rn
         & .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
    !
    !      Fr      Ra      Ac      Th      Pa      U       Np      Pu
         & .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
    !      Am      Cm      Bk      Cf      Es      Fm      Md      No
         & .false.,.false.,.false.,.false.,.true. ,.true. ,.true. ,.true. , &
    !      Lr      -       -       -       -
         & .true. ,.true. ,.true. ,.true., .true.  /
    !
    ! ********* Radii for electrostatics ****************
    ! These are the original Mopac values.
    ! ***************************************************
    data atom_radius_estat / &
    !      H       He
         & 1.20d0, 1.20d0, &
    !
    !      Li      Be      B       C       N       O       F       Ne
         & 1.37d0, 1.45d0, 1.45d0, 1.50d0, 1.50d0, 1.40d0, 1.35d0, 1.30d0, &
    !
    !      Na      Mg      Al      Si      P       S       Cl      Ar
         & 1.57d0, 1.36d0, 1.24d0, 1.17d0, 1.80d0, 1.75d0, 1.70d0, 999.d0, &
    !
    !      K       Ca      Sc      Ti      V       Cr      Mn      Fe
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Co      Ni      Cu      Zn      Ga      Ge      As      Se
         & 999.d0, 999.d0, 999.d0, 1.00d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Br      Kr
         & 2.30d0, 999.d0, &
    !
    !      Rb      Sr      Y       Zr      Nb      Mo      Tc      Ru
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Rh      Pd      Ag      Cd      In      Sn      Sb      Te
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      I       Xe
         & 999.d0, 999.d0, &
    !
    !      Cs      Ba      La      Ce      Pr      Nd      Pm      Sm
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lu      Hf      Ta      W       Re      Os      Ir      Pt
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Au      Hg      Tl      Pb      Bi      Po      At      Rn
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !
    !      Fr      Ra      Ac      Th      Pa      U       Np      Pu
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Am      Cm      Bk      Cf      Es      Fm      Md      No
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, 999.d0, &
    !      Lr      Rf      Db      Sg       Tv
         & 999.d0, 999.d0, 999.d0, 999.d0, 999.d0 /
    save
end module atomradii_C
