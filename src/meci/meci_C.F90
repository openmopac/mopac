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

module meci_C
  integer, parameter :: mmci = 70
  integer ::         &
  &  nmos,           & !  Type        Number of M.O.s in active space in C.I.
                       !  Definition  Number of M.O.s involved in a MECI calculation
                       !  Min. value  0
                       !  Max. value  norbs
  &  lab,            & !  Number of microstates used in the C.I.
  &  labsiz,         & !  Number of states to be calculated (a subset of lab)
  &  root_requested, & !  Number of the root that was requested
  &  nstate,         & !  Type        Number of States used in defining the system
                       !  Definition  The complete manifold of staet vectors defining the state.
                       !  Min. value  1 (for a non-degenerate state)
                       !  Max. value  Lesser of norbs and 22
  &  nelec =-20,     & !  Minus the maximum number of electrons in the active space
  &  nmeci = 22,     & !  Maximum number of M.O.s in active space
  &  maxci,          & !  Largest number of configurations allowed
  &  is,             &
  &  iiloop,         &
  &  jloop,          &
  &  k,              &
  &  dummy

  double precision ::    &
  &  cif1,           & !
  &  cif2,           & !
  &  cdiagi


  integer :: nbo(3)

  integer, dimension(:), allocatable :: &
  &  nalmat, & !
  &  ispin     ! Type        Spins of each State (Singlet = 1)
               ! Definition  "S" quantum number for each State
               ! Min. value   0,  0,  0, etc.
               ! Max. value  10, 10, 10, etc  (9 = Nonet, 10 = "???????")

  integer, dimension(:,:), allocatable :: ispqr

  double precision, dimension(:), allocatable :: &
  &  cdiag,  & !
  &  occa,   & ! Type        Initial M.O. occupancy
               ! Definition  Number of electrons in each M.O. of the SCF in the active space
               ! Units       Electrons
               ! Min. value  0.0, 0.0, etc.
               ! Max. value  2.0, 2.0, etc
  &  conf,   & !  State vectors
  &  eig,    & !  State eigenvalues
  &  vectci    ! Type        State vector of interest
               ! Definition  nstate State vectors, each of length lab, if not degenerate, nstate = 1

  double precision, dimension(:,:), allocatable :: &
  &  rjkaa,  & ! Type
  &  rjkab,  & !
  &  deltap, & ! Type        Change in molecular orbital occupancy as a result of C.I.
               ! Definition  deltap(i,i) = Number of electrons gained or lost in forming the State
               ! Units       Electrons
  &  cxy       !

  double precision, dimension(:,:,:), allocatable :: &
  &  dijkl
  double precision, dimension(:,:,:,:), allocatable :: &
  &  xy

  integer, dimension(:,:), allocatable :: &
  &  microa, & ! Type       microstates
               ! Definition occupancy of alpha M.O.s in microstates
               ! Units      Electrons
  &  microb    ! Type       Beta microstates
               ! Definition occupancy of beta M.O.s in microstates
               ! Units      Electrons
end module meci_C
