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

program mopac
  use chanel_C, only : iw0
  use molkst_C, only : gui, line, verson
  implicit none
  call getdatestamp(line, verson)
#ifdef MOPAC_OS
  verson(7:7) = MOPAC_OS
#else
  verson(7:7) = "X"
#endif
  gui = .false.
  iw0 = -1
  call run_mopac 
end program mopac
  
