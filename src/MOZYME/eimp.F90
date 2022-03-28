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

subroutine eimp ()
!
!   Store in the s-s location of array p the value of the Fock terms relating
!   atom A and atom B, for all pairs of atoms.  This quantity will be used by
!   the diagonalizer
!
  use molkst_C, only: numat
  use common_arrays_C, only : p, f
  use MOZYME_C, only : iorbs
  implicit none
  integer :: i, j, k, l, m
  double precision :: sum
  integer, external :: ijbo
  do i = 1, numat
    do j = 1, i - 1
      k = ijbo (i, j)
      if (k >= 0) then
        l = iorbs(i) * iorbs(j)
        if (l /= 0) then
          l = k + l
          sum = 0.d0
          do m = k + 1, l
            sum = sum + f(m) ** 2
          end do
          p(k+1) = sum
        end if
      end if
    end do
  end do
end subroutine eimp
