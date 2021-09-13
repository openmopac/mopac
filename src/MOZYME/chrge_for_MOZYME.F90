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

subroutine chrge_for_MOZYME (p, q)
!
! Evaluate the atomic partial charge.
!
    use molkst_C, only: mpack, numat
    use MOZYME_C, only : iorbs
    implicit none
    double precision, dimension (mpack), intent (in) :: p
    double precision, dimension (numat), intent (out) :: q
!
    integer :: i, ii, j
    double precision :: sum
    integer, external :: ijbo
!
    do i = 1, numat
      ii = ijbo (i, i)
      sum = 0.d0
      do j = 1, iorbs(i)
        ii = ii + j
        sum = sum + p(ii)
      end do
      q(i) = sum
    end do
end subroutine chrge_for_MOZYME
