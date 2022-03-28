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
