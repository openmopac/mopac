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

submodule (mopac_api) mopac_api_initialize
  implicit none

contains

  ! intialize geometry and data structures for a MOPAC/MOZYME calculation
  subroutine mopac_initialize(system, status)
    type(mopac_system), intent(in) :: system
    integer, intent(out) :: status
    ! TO DO: implement this
    ! Common_arrays_C
    ! molkst_C
    ! molmec_C
    ! overlaps_C ???
    ! MOZYME_C (conditional)
    ! chanel_C (redirect all to /dev/null)
    ! COSMO init ???
    ! parameters_C
  end subroutine mopac_initialize

end submodule mopac_api_initialize