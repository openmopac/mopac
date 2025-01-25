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

module to_screen_C
!
!  This module contains all the quantities relating to the system being calculated
!  that are likely to be used in the subroutine "to_screen".
!
  double precision :: &
  rot(3), &    ! Type          Rotational constants
               ! Definition    Frequencies of the three rotational constants
               ! Units         cm**(-1)
               ! Min inclusive 0.0
               ! Max inclusive no limit
               !
  xyzmom(3), & ! Type          Principal moments of inertia
               ! Definition    Moments of inertia for the system.
               ! Units         10**(-40)*gram-cm**2
               ! Min inclusive 0.0
               ! Max inclusive no limit
               !
  dip(4,3)     ! Type          Dipole moment array
               ! Definition    Point-charge and hybrid components of dipole in X, Y, and Z, and totals
               ! Units         Debye
               !
  double precision, dimension (:,:), allocatable :: &
  fcint,     & ! Type          Internal force constants
               ! Definition    Force constants for all coordinates of all atoms, using the coordinate system supplied
               ! Units         Millidynes/Angstrom and millidynes per radian
               !
  to_a,      & ! Type          Generic two-dimensional matrix
               ! Definition    Used for transferring matrix "a" in matou1 to here
               ! Units         Unknown
  redmas       ! Type          Reduced masses
               ! Definition    Effective mass of vibrational frequencies
               ! Units         Atomic mass units (amu)
               !
  double precision, dimension (:), allocatable :: &
  dipt,      & ! Type          Transition dipoles for vibrational frequencies
               ! Definition    <ground state|eR(x)|vibrational state>
               ! Units         electrons
force_const, & ! Type          Force constant
               ! Definition    Force constant for of vibrational frequencies
               ! Units         millidynes per Angstrom
               !
  travel,    & ! Type          Distance traveled during a vibration
               ! Definition
               ! Units         Angstroms
  freq,      & ! Type          Vibrational frequencies
               ! Definition    Normal mode frequencies
               ! Units         cm(-1)

  cnorml,    & ! Type          Normal modes of vibration
               ! Definition    Orthonormal coordinates of vibration
               ! Units         (None) Normalized to unity
  to_b         ! Type          Generic one-dimensional matrix

end module to_screen_C
