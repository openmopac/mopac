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

      module maps_C

      double precision, dimension(:), allocatable  :: &
   &  surf,    &
   &  react

      integer :: &
  &  ijlp    , & !
  &  ilp,      & !
  &  jlp,      & !
  &  jlp1,     & !
  &  ione,     & !
  &  latom,    & !
  &  lparam,   & !
  &  kloop,    &
  &  latom1,   &
  &  lpara1,   &
  &  latom2,   &
  &  lpara2
      double precision :: &
  &  rxn_coord1 = 0.d0,   &
  &  rxn_coord2 = 0.d0,   &
  &  rxn_coord,           &
  &  rc_escf,             &  !  Reaction coordinate Heat of Formation
  &  rc_dipo,             &  !  Reaction coordinate dipole
  &  ekin,                &  !  Kinetic energy
  &  dummy
      end module maps_C
