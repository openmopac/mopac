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

module interface_C
!dec$ attributes dllexport :: gui, iw0
!
!  This module contains the data exposed by the MOPAC shared library API
!
    logical :: gui = .true.   !  By default, output information for a Graphical User Interface
    integer :: iw0 = -1       !  Abbreviated output channel, for GUI (By default, not used)
end module interface_C