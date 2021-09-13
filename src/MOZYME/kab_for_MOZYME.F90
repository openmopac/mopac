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

subroutine kab_for_MOZYME (ia, ja, pk, w, f)
    implicit none
    integer, intent (in) :: ia, ja
    double precision, dimension (16), intent (in) :: pk
    double precision, dimension (100), intent (in) :: w
    double precision, dimension (*), intent (inout) :: f
    integer :: i, j, j1, j2, j3, m
    integer, dimension (16, 16) :: kkind
    double precision, dimension (16) :: sum
    data kkind / 1, 2, 4, 7, 11, 12, 14, 17, 31, 32, 34, 37, 61, 62, 64, 67, &
   & 2, 3, 5, 8, 12, 13, 15, 18, 32, 33, 35, 38, 62, 63, 65, 68, 4, 5, 6, 9, &
   & 14, 15, 16, 19, 34, 35, 36, 39, 64, 65, 66, 69, 7, 8, 9, 10, 17, 18, 19, &
   & 20, 37, 38, 39, 40, 67, 68, 69, 70, 11, 12, 14, 17, 21, 22, 24, 27, 41, &
   & 42, 44, 47, 71, 72, 74, 77, 12, 13, 15, 18, 22, 23, 25, 28, 42, 43, 45, &
   & 48, 72, 73, 75, 78, 14, 15, 16, 19, 24, 25, 26, 29, 44, 45, 46, 49, 74, &
   & 75, 76, 79, 17, 18, 19, 20, 27, 28, 29, 30, 47, 48, 49, 50, 77, 78, 79, &
   & 80, 31, 32, 34, 37, 41, 42, 44, 47, 51, 52, 54, 57, 81, 82, 84, 87, 32, &
   & 33, 35, 38, 42, 43, 45, 48, 52, 53, 55, 58, 82, 83, 85, 88, 34, 35, 36, &
   & 39, 44, 45, 46, 49, 54, 55, 56, 59, 84, 85, 86, 89, 37, 38, 39, 40, 47, &
   & 48, 49, 50, 57, 58, 59, 60, 87, 88, 89, 90, 61, 62, 64, 67, 71, 72, 74, &
   & 77, 81, 82, 84, 87, 91, 92, 94, 97, 62, 63, 65, 68, 72, 73, 75, 78, 82, &
   & 83, 85, 88, 92, 93, 95, 98, 64, 65, 66, 69, 74, 75, 76, 79, 84, 85, 86, &
   & 89, 94, 95, 96, 99, 67, 68, 69, 70, 77, 78, 79, 80, 87, 88, 89, 90, 97, &
   & 98, 99, 100 /
    do j = 1, 16
      sum(j) = 0.0d0
      do i = 1, 16
        sum(j) = sum(j) + pk(i) * w(kkind(i, j))
      end do
    end do
   !
    if (ia > ja) then
      !
      !   IA IS GREATER THAN JA, THEREFORE USE HALF OF TRIANGLE
      !
      m = 0
      do j1 = ia, ia + 3
        j = (j1*(j1-1)) / 2
        do j2 = ja, ja + 3
          m = m + 1
          j3 = j + j2
          f(j3) = f(j3) - sum(m)
        end do
      end do
    else if (ia < ja) then
      !
      !   IA IS LESS THAN JA, THEREFORE USE OTHER HALF OF TRIANGLE
      !
      m = 0
      do j1 = ia, ia + 3
        do j2 = ja, ja + 3
          m = m + 1
          j3 = (j2*(j2-1)) / 2 + j1
          f(j3) = f(j3) - sum(m)
        end do
      end do
    else
      !
      !   MOZYME - type - array elements are stored simply.
      !   (The factor of 0.5 is because MOZYME uses only the
      !    TOTAL density matrix.)
      !
      do i = 1, 16
        f(i) = f(i) - sum(i) * 0.5d0
      end do
    end if
end subroutine kab_for_MOZYME
