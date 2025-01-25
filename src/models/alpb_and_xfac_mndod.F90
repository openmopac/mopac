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
