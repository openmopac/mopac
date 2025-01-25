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

subroutine upcase (keywrd)
  implicit none
  character(len=*), intent(inout) :: keywrd
!
  integer :: i, icapa, iline, ilowa, ilowz
  icapa = Ichar ("A")
  ilowa = Ichar ("a")
  ilowz = Ichar ("z")
  do i = 1, Len (keywrd)
    iline = Ichar (keywrd(i:i))
    if (iline>=ilowa .and. iline<=ilowz) then
      keywrd(i:i) = Char(iline+icapa-ilowa)
    end if
  end do
end subroutine upcase
