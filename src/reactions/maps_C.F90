! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
