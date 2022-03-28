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

      subroutine kab(ia, ja, pk, w, f)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!#      SUBROUTINE KAB(IA,JA, PK, W, KINDEX, F)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ia
      integer , intent(in) :: ja
      double precision , intent(in) :: pk(*)
      double precision , intent(in) :: w(*)
      double precision , intent(inout) :: f(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: m, j1, j, j2, j3
      double precision, dimension(16) :: sum
!-----------------------------------------------
      sum(1) = pk(1)*w(1) + pk(2)*w(2) + pk(3)*w(4) + pk(4)*w(7) + pk(5)*w(11)&
         + pk(6)*w(12) + pk(7)*w(14) + pk(8)*w(17) + pk(9)*w(31) + pk(10)*w(32)&
         + pk(11)*w(34) + pk(12)*w(37) + pk(13)*w(61) + pk(14)*w(62) + pk(15)*w&
        (64) + pk(16)*w(67)
      sum(2) = pk(1)*w(2) + pk(2)*w(3) + pk(3)*w(5) + pk(4)*w(8) + pk(5)*w(12)&
         + pk(6)*w(13) + pk(7)*w(15) + pk(8)*w(18) + pk(9)*w(32) + pk(10)*w(33)&
         + pk(11)*w(35) + pk(12)*w(38) + pk(13)*w(62) + pk(14)*w(63) + pk(15)*w&
        (65) + pk(16)*w(68)
      sum(3) = pk(1)*w(4) + pk(2)*w(5) + pk(3)*w(6) + pk(4)*w(9) + pk(5)*w(14)&
         + pk(6)*w(15) + pk(7)*w(16) + pk(8)*w(19) + pk(9)*w(34) + pk(10)*w(35)&
         + pk(11)*w(36) + pk(12)*w(39) + pk(13)*w(64) + pk(14)*w(65) + pk(15)*w&
        (66) + pk(16)*w(69)
      sum(4) = pk(1)*w(7) + pk(2)*w(8) + pk(3)*w(9) + pk(4)*w(10) + pk(5)*w(17)&
         + pk(6)*w(18) + pk(7)*w(19) + pk(8)*w(20) + pk(9)*w(37) + pk(10)*w(38)&
         + pk(11)*w(39) + pk(12)*w(40) + pk(13)*w(67) + pk(14)*w(68) + pk(15)*w&
        (69) + pk(16)*w(70)
      sum(5) = pk(1)*w(11) + pk(2)*w(12) + pk(3)*w(14) + pk(4)*w(17) + pk(5)*w(&
        21) + pk(6)*w(22) + pk(7)*w(24) + pk(8)*w(27) + pk(9)*w(41) + pk(10)*w(&
        42) + pk(11)*w(44) + pk(12)*w(47) + pk(13)*w(71) + pk(14)*w(72) + pk(15&
        )*w(74) + pk(16)*w(77)
      sum(6) = pk(1)*w(12) + pk(2)*w(13) + pk(3)*w(15) + pk(4)*w(18) + pk(5)*w(&
        22) + pk(6)*w(23) + pk(7)*w(25) + pk(8)*w(28) + pk(9)*w(42) + pk(10)*w(&
        43) + pk(11)*w(45) + pk(12)*w(48) + pk(13)*w(72) + pk(14)*w(73) + pk(15&
        )*w(75) + pk(16)*w(78)
      sum(7) = pk(1)*w(14) + pk(2)*w(15) + pk(3)*w(16) + pk(4)*w(19) + pk(5)*w(&
        24) + pk(6)*w(25) + pk(7)*w(26) + pk(8)*w(29) + pk(9)*w(44) + pk(10)*w(&
        45) + pk(11)*w(46) + pk(12)*w(49) + pk(13)*w(74) + pk(14)*w(75) + pk(15&
        )*w(76) + pk(16)*w(79)
      sum(8) = pk(1)*w(17) + pk(2)*w(18) + pk(3)*w(19) + pk(4)*w(20) + pk(5)*w(&
        27) + pk(6)*w(28) + pk(7)*w(29) + pk(8)*w(30) + pk(9)*w(47) + pk(10)*w(&
        48) + pk(11)*w(49) + pk(12)*w(50) + pk(13)*w(77) + pk(14)*w(78) + pk(15&
        )*w(79) + pk(16)*w(80)
      sum(9) = pk(1)*w(31) + pk(2)*w(32) + pk(3)*w(34) + pk(4)*w(37) + pk(5)*w(&
        41) + pk(6)*w(42) + pk(7)*w(44) + pk(8)*w(47) + pk(9)*w(51) + pk(10)*w(&
        52) + pk(11)*w(54) + pk(12)*w(57) + pk(13)*w(81) + pk(14)*w(82) + pk(15&
        )*w(84) + pk(16)*w(87)
      sum(10) = pk(1)*w(32) + pk(2)*w(33) + pk(3)*w(35) + pk(4)*w(38) + pk(5)*w&
        (42) + pk(6)*w(43) + pk(7)*w(45) + pk(8)*w(48) + pk(9)*w(52) + pk(10)*w&
        (53) + pk(11)*w(55) + pk(12)*w(58) + pk(13)*w(82) + pk(14)*w(83) + pk(&
        15)*w(85) + pk(16)*w(88)
      sum(11) = pk(1)*w(34) + pk(2)*w(35) + pk(3)*w(36) + pk(4)*w(39) + pk(5)*w&
        (44) + pk(6)*w(45) + pk(7)*w(46) + pk(8)*w(49) + pk(9)*w(54) + pk(10)*w&
        (55) + pk(11)*w(56) + pk(12)*w(59) + pk(13)*w(84) + pk(14)*w(85) + pk(&
        15)*w(86) + pk(16)*w(89)
      sum(12) = pk(1)*w(37) + pk(2)*w(38) + pk(3)*w(39) + pk(4)*w(40) + pk(5)*w&
        (47) + pk(6)*w(48) + pk(7)*w(49) + pk(8)*w(50) + pk(9)*w(57) + pk(10)*w&
        (58) + pk(11)*w(59) + pk(12)*w(60) + pk(13)*w(87) + pk(14)*w(88) + pk(&
        15)*w(89) + pk(16)*w(90)
      sum(13) = pk(1)*w(61) + pk(2)*w(62) + pk(3)*w(64) + pk(4)*w(67) + pk(5)*w&
        (71) + pk(6)*w(72) + pk(7)*w(74) + pk(8)*w(77) + pk(9)*w(81) + pk(10)*w&
        (82) + pk(11)*w(84) + pk(12)*w(87) + pk(13)*w(91) + pk(14)*w(92) + pk(&
        15)*w(94) + pk(16)*w(97)
      sum(14) = pk(1)*w(62) + pk(2)*w(63) + pk(3)*w(65) + pk(4)*w(68) + pk(5)*w&
        (72) + pk(6)*w(73) + pk(7)*w(75) + pk(8)*w(78) + pk(9)*w(82) + pk(10)*w&
        (83) + pk(11)*w(85) + pk(12)*w(88) + pk(13)*w(92) + pk(14)*w(93) + pk(&
        15)*w(95) + pk(16)*w(98)
      sum(15) = pk(1)*w(64) + pk(2)*w(65) + pk(3)*w(66) + pk(4)*w(69) + pk(5)*w&
        (74) + pk(6)*w(75) + pk(7)*w(76) + pk(8)*w(79) + pk(9)*w(84) + pk(10)*w&
        (85) + pk(11)*w(86) + pk(12)*w(89) + pk(13)*w(94) + pk(14)*w(95) + pk(&
        15)*w(96) + pk(16)*w(99)
      sum(16) = pk(1)*w(67) + pk(2)*w(68) + pk(3)*w(69) + pk(4)*w(70) + pk(5)*w&
        (77) + pk(6)*w(78) + pk(7)*w(79) + pk(8)*w(80) + pk(9)*w(87) + pk(10)*w&
        (88) + pk(11)*w(89) + pk(12)*w(90) + pk(13)*w(97) + pk(14)*w(98) + pk(&
        15)*w(99) + pk(16)*w(100)
      if (ia > ja) then
        m = 0
        do j1 = ia, ia + 3
          j = (j1*(j1 - 1))/2
          f(j+ja:3+j+ja) = f(j+ja:3+j+ja) - sum(m+1:4+m)
          m = 4 + m
        end do
      else
!
!   IA IS LESS THAN JA, THEREFORE USE OTHER HALF OF TRIANGLE
!
        m = 0
        do j1 = ia, ia + 3
          do j2 = ja, ja + 3
            m = m + 1
            j3 = (j2*(j2 - 1))/2 + j1
            f(j3) = f(j3) - sum(m)
          end do
        end do
      end if
      return
      end subroutine kab
