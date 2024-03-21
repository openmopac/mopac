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

! Diskless Application Programming Interface to core MOPAC functionality
module mopac_api
!dec$ attributes dllexport :: mopac_api

  implicit none
  private

  public :: mopac_system, mopac_state, mopac_properties, mozyme_state
  public :: mopac_scf, mopac_relax, mopac_vibe, mozyme_scf, mozyme_relax

contains

  subroutine mopac_scf(system, state, properties)
  !dec$ attributes dllexport :: mopac_scf
  end subroutine mopac_scf

  subroutine mopac_relax(system, state, properties)
  !dec$ attributes dllexport :: mopac_relax
  end subroutine mopac_relax

  subroutine mopac_vibe(system, state, properties)
  !dec$ attributes dllexport :: mopac_vibe
  end subroutine mopac_vibe

  subroutine mozyme_scf(system, state, properties)
  !dec$ attributes dllexport :: mozyme_scf
  end subroutine mozyme_scf

  subroutine mozyme_relax(system, state, properties)
  !dec$ attributes dllexport :: mozyme_relax
  end subroutine mozyme_relax

end module diskless_api