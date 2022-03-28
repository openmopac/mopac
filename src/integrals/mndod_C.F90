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

      module mndod_C
      integer, dimension(9,9) :: indexd, indx
      integer, dimension(243) :: intij, intkl, intrep
      integer, dimension(107) :: iii, iiid
      integer, dimension(3,3) :: indpp
      integer, dimension(5,3) :: inddp
      integer, dimension(5,5) :: inddd
      integer, dimension(500) :: iaf, ial
      integer, dimension(45,45) :: ind2
      integer, dimension(491) :: isym
      integer :: nalp
      double precision, dimension(45,0:2,-2:2) :: ch
      double precision, dimension(500) :: alpb, xfac
      double precision, dimension(6,107) :: aij
      double precision, dimension(52,107) :: repd
      double precision, dimension(10,2) :: cored
      double precision, dimension(3,3) :: sp
      double precision, dimension(5,5) :: sd
      double precision, dimension(6,3,3) :: pp
      double precision, dimension(15,5,3) :: dp
      double precision, dimension(15,5,5) :: d_d
      double precision, dimension(30) :: fx   !   Factorials:  fx(i) = i!
      double precision, dimension(30,30) :: b  !   Binomial expansion (Pascal's triangle): b(1,1) = 1

      data iii/ 2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 21*0/
      data iiid/ 30*3, 18*4, 32*5, 6*6, 21*0/
      data intij/ 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4&
        , 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10&
        , 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, &
        13, 13, 13, 13, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16&
        , 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 19, 19, 19, 19, 19, &
        20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22&
        , 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, &
        24, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28&
        , 28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 31, &
        31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 34, 34, 34, 34&
        , 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, &
        37, 37, 37, 37, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 40, 40, 40, 41&
        , 42, 42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 44, 45, 45, 45, &
        45, 45, 45, 45, 45, 45, 45/
      data intkl/ 15, 21, 28, 36, 45, 12, 19, 23, 39, 11, 15, 21, 22, 26, 28, &
        36, 45, 13, 24, 32, 38, 34, 37, 43, 11, 15, 21, 22, 26, 28, 36, 45, 17&
        , 25, 31, 16, 20, 27, 44, 29, 33, 35, 42, 15, 21, 22, 28, 36, 45, 3, 6&
        , 11, 21, 26, 36, 2, 12, 19, 23, 39, 4, 13, 24, 32, 38, 14, 17, 31, 1, &
        3, 6, 10, 15, 21, 22, 28, 36, 45, 8, 16, 20, 27, 44, 7, 14, 17, 25, 31&
        , 18, 30, 40, 2, 12, 19, 23, 39, 8, 16, 20, 27, 44, 1, 3, 6, 10, 11, 15&
        , 21, 22, 26, 28, 36, 45, 3, 6, 10, 15, 21, 22, 28, 36, 45, 2, 12, 19, &
        23, 39, 4, 13, 24, 32, 38, 7, 17, 25, 31, 3, 6, 11, 21, 26, 36, 8, 16, &
        20, 27, 44, 1, 3, 6, 10, 15, 21, 22, 28, 36, 45, 9, 29, 33, 35, 42, 18&
        , 30, 40, 7, 14, 17, 25, 31, 4, 13, 24, 32, 38, 9, 29, 33, 35, 42, 5, &
        34, 37, 43, 9, 29, 33, 35, 42, 1, 3, 6, 10, 11, 15, 21, 22, 26, 28, 36&
        , 45, 5, 34, 37, 43, 4, 13, 24, 32, 38, 2, 12, 19, 23, 39, 18, 30, 40, &
        41, 9, 29, 33, 35, 42, 5, 34, 37, 43, 8, 16, 20, 27, 44, 1, 3, 6, 10, &
        15, 21, 22, 28, 36, 45/
      data intrep/ 1, 1, 1, 1, 1, 3, 3, 8, 3, 9, 6, 6, 12, 14, 13, 7, 6, 15, 8&
        , 3, 3, 11, 9, 14, 17, 6, 7, 12, 18, 13, 6, 6, 3, 2, 3, 9, 11, 10, 11, &
        9, 16, 10, 11, 7, 6, 4, 5, 6, 7, 9, 17, 19, 32, 22, 40, 3, 33, 34, 27, &
        46, 15, 33, 28, 41, 47, 35, 35, 42, 1, 6, 6, 7, 29, 38, 22, 31, 38, 51&
        , 9, 19, 32, 21, 32, 3, 35, 33, 24, 34, 35, 35, 35, 3, 34, 33, 26, 34, &
        11, 32, 44, 37, 49, 1, 6, 7, 6, 32, 38, 29, 21, 39, 30, 38, 38, 12, 12&
        , 4, 22, 21, 19, 20, 21, 22, 8, 27, 26, 25, 27, 8, 28, 25, 26, 27, 2, &
        24, 23, 24, 14, 18, 22, 39, 48, 45, 10, 21, 37, 36, 37, 1, 13, 13, 5, &
        31, 30, 20, 29, 30, 31, 9, 19, 40, 21, 32, 35, 35, 35, 3, 42, 34, 24, &
        33, 3, 41, 26, 33, 34, 16, 40, 44, 43, 50, 11, 44, 32, 39, 10, 21, 43, &
        36, 37, 1, 7, 6, 6, 40, 38, 38, 21, 45, 30, 29, 38, 9, 32, 19, 22, 3, &
        47, 27, 34, 33, 3, 46, 34, 27, 33, 35, 35, 35, 52, 11, 32, 50, 37, 44, &
        14, 39, 22, 48, 11, 32, 49, 37, 44, 1, 6, 6, 7, 51, 38, 22, 31, 38, 29&
        /
      end module mndod_C
