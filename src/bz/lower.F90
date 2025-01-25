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

subroutine lower (a, n)
  implicit none
  character(len=*), intent(inout) :: a
  integer, intent(in) :: n
  integer :: i, j, lowa, lupa, lupz
  intrinsic Char, Ichar
  lowa = Ichar ("a")
  lupa = Ichar ("A")
  lupz = Ichar ("Z")
  do i = 1, n
    if (Ichar (a(i:i))>=lupa .and. Ichar (a(i:i))<=lupz) then
      j = Ichar (a(i:i)) + lowa - lupa
      a(i:i) = Char(j)
    end if
  end do
end subroutine lower
