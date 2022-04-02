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

  subroutine alpb_and_xfac_mndod
    use parameters_C, only : xfac, alpb
    xfac = 0.d0
    alpb = 0.d0
    alpb(11, 1) =  1.05225212D0  !     Sodium - Hydrogen
    xfac(11, 1) =  1.00000000d0  !     Sodium - Hydrogen
    alpb(11, 6) =  1.05225212D0  !     Sodium - Carbon
    xfac(11, 6) =  1.00000000d0  !     Sodium - Carbon
    alpb(12, 1) =  1.35052992D0  !     Magnesium - Hydrogen
    xfac(12, 1) =  1.00000000d0  !     Magnesium - Hydrogen
    alpb(12, 6) =  1.48172071D0  !     Magnesium - Hydrogen
    xfac(12, 6) =  1.00000000d0  !     Magnesium - Hydrogen
    alpb(16,12) =  1.48172071D0  !     Sulfur - Magnesium
    xfac(16,12) =  1.00000000d0  !     Sulfur - Magnesium
    alpb(13, 1) =  1.38788000D0  !     Aluminum - Hydrogen
    xfac(13, 1) =  1.00000000d0  !     Aluminum - Hydrogen
    alpb(13, 6) =  1.38788000D0  !     Aluminum - Carbon
    xfac(13, 6) =  1.00000000d0  !     Aluminum - Carbon
    alpb(13,13) =  1.38788000D0  !     Aluminum - Aluminum
    xfac(13,13) =  1.00000000d0  !     Aluminum - Aluminum
  end subroutine alpb_and_xfac_mndod
