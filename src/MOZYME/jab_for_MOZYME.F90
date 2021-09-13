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

subroutine jab_for_MOZYME (ia, ja, pja, pjb, w, f1, f2)
    implicit none
    integer, intent (in) :: ia, ja
    double precision, dimension (16), intent (in) :: pja, pjb
    double precision, dimension (100), intent (in) :: w
    double precision, dimension (*), intent (inout) :: f1, f2
    integer :: i, i5, i6, iia, ija, ioff, j, joff
    integer, dimension (10, 16) :: jindx, jjndx
    double precision, dimension (10) :: suma, sumb
    data jindx / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, &
   & 18, 19, 20, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 61, 62, 63, 64, 65, &
   & 66, 67, 68, 69, 70, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, &
   & 24, 25, 26, 27, 28, 29, 30, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 71, &
   & 72, 73, 74, 75, 76, 77, 78, 79, 80, 31, 32, 33, 34, 35, 36, 37, 38, 39, &
   & 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, &
   & 58, 59, 60, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 61, 62, 63, 64, 65, &
   & 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, &
   & 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100 /
    data jjndx / 1, 11, 21, 31, 41, 51, 61, 71, 81, 91, 2, 12, 22, 32, 42, 52, &
   & 62, 72, 82, 92, 4, 14, 24, 34, 44, 54, 64, 74, 84, 94, 7, 17, 27, 37, 47, &
   & 57, 67, 77, 87, 97, 2, 12, 22, 32, 42, 52, 62, 72, 82, 92, 3, 13, 23, 33, &
   & 43, 53, 63, 73, 83, 93, 5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 8, 18, 28, &
   & 38, 48, 58, 68, 78, 88, 98, 4, 14, 24, 34, 44, 54, 64, 74, 84, 94, 5, 15, &
   & 25, 35, 45, 55, 65, 75, 85, 95, 6, 16, 26, 36, 46, 56, 66, 76, 86, 96, 9, &
   & 19, 29, 39, 49, 59, 69, 79, 89, 99, 7, 17, 27, 37, 47, 57, 67, 77, 87, &
   & 97, 8, 18, 28, 38, 48, 58, 68, 78, 88, 98, 9, 19, 29, 39, 49, 59, 69, 79, &
   & 89, 99, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 /
    do i = 1, 10
      suma(i) = 0.0d0
      sumb(i) = 0.0d0
    end do
    do j = 1, 16
      do i = 1, 10
        suma(i) = suma(i) + pja(j) * w(jindx(i, j))
        sumb(i) = sumb(i) + pjb(j) * w(jjndx(i, j))
      end do
    end do
   !
    i = 0
    do i5 = 1, 4
      iia = ia + i5 - 1
      ija = ja + i5 - 1
      ioff = (iia*(iia-1)) / 2 + ia - 1
      joff = (ija*(ija-1)) / 2 + ja - 1
      do i6 = 1, i5
        ioff = ioff + 1
        joff = joff + 1
        i = i + 1
        f1(ioff) = f1(ioff) + sumb(i)
        f2(joff) = f2(joff) + suma(i)
      end do
    end do
end subroutine jab_for_MOZYME
