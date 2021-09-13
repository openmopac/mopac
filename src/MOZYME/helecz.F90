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

double precision function helecz ()
!
!   Return the total electronic energy,
!   from e = 0.5P(H + F)
!
    use molkst_C, only: numat
    use MOZYME_C, only : iorbs
    use common_arrays_C, only : p, h, f
    implicit none
    integer :: i, j, k, l, l1, l2, m
    double precision :: ed, ee
    integer, external :: ijbo
    ed = 0.0d00
    ee = 0.0d00
    do i = 1, numat
      do j = 1, i - 1
          k = ijbo (i, j)
          if (k >= 0) then
            l = k + iorbs(i) * iorbs(j)
            do m = k + 1, l
              ee = ee + p(m) * (h(m)+f(m))
            end do
          end if
      end do
        k = ijbo (i, i)
        do l1 = 1, iorbs(i)
          do l2 = 1, l1 - 1
            k = k + 1
            ee = ee + p(k) * (h(k)+f(k))
          end do
          k = k + 1
          ed = ed + p(k) * (h(k)+f(k))
        end do
    end do
    ee = ee + .5d00 * ed
    helecz = ee
end function helecz
