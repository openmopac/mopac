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

  module conref_C
    double precision, dimension(2,10) :: fpcref
!
!   Fundamental physical constants
!
!   fpcref(1,*) = 2018 CODATA  fpcref(2,8) = old constants
!
!   In the 2019 revision of SI units, 7 fundamental constants were switched from
!   uncertain measured values to exact defined values. Similarly, the non-SI calorie
!   has an exact defined value (4.184 J) that is used instead of old measured values.
!
    data fpcref(1, 1) / 1.602176634d0 /         ! Elementary charge, in Coulombs*10**(-19) (exact)
    data fpcref(2, 1) / 1.60217733d0 /
    data fpcref(1, 2) / 14.399645478456d0 /     ! a0*(one au in eV) = ref(3)*ref(4)
    data fpcref(2, 2) / 14.399d0 /
    data fpcref(1, 3) / 0.529177210903d0 /      ! a0 (+/- 1.5e-10 RSD)
    data fpcref(2, 3) / 0.529167d0 /
    data fpcref(1, 4) / 27.211386245988d0 /     ! one au in eV (+/- 1.9e-12 RSD)
    data fpcref(2, 4) / 27.21d0 /
    data fpcref(1, 5) / 1.987204258640832d0 /   ! "R", the gas constant, in calories (exact)
    data fpcref(2, 5) / 1.98726d0 /
    data fpcref(1, 6) / 6.62607015d-27 /        ! Planck's constant in erg-seconds (exact)
    data fpcref(2, 6) / 6.626d-27 /
    data fpcref(1, 7) / 1.380649d-16 /          ! Boltzmann's constant, in J/K (exact)
    data fpcref(2, 7) / 1.3807d-16 /
    data fpcref(1, 8) / 2.99792458d+10 /        ! Speed of light, in cm/sec (exact)
    data fpcref(2, 8) / 2.99776d+10 /
    data fpcref(1, 9) / 23.060547830619029d0 /  ! 1eV in kcal/mol = (ref(1)*ref(10)*10**(-3))/4.184
    data fpcref(2, 9) / 23.061d0 /
    data fpcref(1, 10) / 6.02214076d+23 /       ! Avogadro's number, in 1/mol (exact)
    data fpcref(2, 10) / 6.02205d+23 /
  end module conref_C
