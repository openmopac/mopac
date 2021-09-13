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

    integer function ijbo (ii, jj)
    !**********************************************************************
    !
    !  Given atom numbers II and JJ,
    !  IJBO = -1 if the atoms are separated by more than SQRT(CUTOF1), else
    !  IJBO = -2 if the atoms are separated by more than SQRT(CUTOF2), else
    !  IJBO = N, where "N+1" in the address of the starting location in the
    !            H, P, and F arrays.
    !
    !**********************************************************************
      use overlaps_C, only: cutof1, cutof2
      use common_arrays_C, only: coord
      use MOZYME_C, only : lijbo, nijbo, iij, numij, ijall, iijj
      implicit none
      integer, intent (in) :: ii, jj
      integer :: ik, il, iu, j, ind_i, ind_j, jm1, jm2
      double precision :: r
      if (lijbo) then
        ijbo = nijbo(ii, jj)
        return
      end if
!
      r =      (coord(1, ii)-coord(1, jj)) ** 2 &
           & + (coord(2, ii)-coord(2, jj)) ** 2 &
           & + (coord(3, ii)-coord(3, jj)) ** 2
      if (r > cutof1) then
        ijbo = -1
        return
      end if
      if (r > cutof2) then
        ijbo = -2
        return
      end if
    !
    !   use lookup array to find ijbo
    !
      if (ii > jj) then
        ind_i = ii
        ind_j = jj
      else
        ind_i = jj
        ind_j = ii
      end if
      il = iij(ind_i)
      iu = numij(ind_i)
      j = (il+iu+1) / 2
      jm1 = 0
      jm2 = 0
      do
        ik = ijall(j)
        if (ik < ind_j) then
          il = j
          j = (j+iu+1) / 2
        else
          if (ik == ind_j) exit
          iu = j
          j = (j+il) / 2
          if (j == jm2) then
            ijbo = -2
            return
          end if
          jm1 = j
          jm2 = jm1
        end if
      end do
      ijbo = iijj(j)
      return
    end function ijbo
