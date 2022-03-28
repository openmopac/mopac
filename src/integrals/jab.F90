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

      subroutine jab(ia, ja, pja, pjb, w, f)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!#      SUBROUTINE JAB(IA,JA,LLPERM,JINDEX, JJNDEX,PJA,PJB,W, F)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ia
      integer , intent(in) :: ja
      double precision , intent(in) :: pja(16)
      double precision , intent(in) :: pjb(16)
      double precision , intent(in) :: w(*)
      double precision , intent(inout) :: f(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, i5, iia, ija, ioff, joff, i6
      double precision, dimension(10) :: suma, sumb
!-----------------------------------------------

      suma(1) = pja(1)*w(1) + pja(2)*w(11) + pja(3)*w(31) + pja(4)*w(61) + pja(&
        5)*w(11) + pja(6)*w(21) + pja(7)*w(41) + pja(8)*w(71) + pja(9)*w(31) + &
        pja(10)*w(41) + pja(11)*w(51) + pja(12)*w(81) + pja(13)*w(61) + pja(14)&
        *w(71) + pja(15)*w(81) + pja(16)*w(91)
      suma(2) = pja(1)*w(2) + pja(2)*w(12) + pja(3)*w(32) + pja(4)*w(62) + pja(&
        5)*w(12) + pja(6)*w(22) + pja(7)*w(42) + pja(8)*w(72) + pja(9)*w(32) + &
        pja(10)*w(42) + pja(11)*w(52) + pja(12)*w(82) + pja(13)*w(62) + pja(14)&
        *w(72) + pja(15)*w(82) + pja(16)*w(92)
      suma(3) = pja(1)*w(3) + pja(2)*w(13) + pja(3)*w(33) + pja(4)*w(63) + pja(&
        5)*w(13) + pja(6)*w(23) + pja(7)*w(43) + pja(8)*w(73) + pja(9)*w(33) + &
        pja(10)*w(43) + pja(11)*w(53) + pja(12)*w(83) + pja(13)*w(63) + pja(14)&
        *w(73) + pja(15)*w(83) + pja(16)*w(93)
      suma(4) = pja(1)*w(4) + pja(2)*w(14) + pja(3)*w(34) + pja(4)*w(64) + pja(&
        5)*w(14) + pja(6)*w(24) + pja(7)*w(44) + pja(8)*w(74) + pja(9)*w(34) + &
        pja(10)*w(44) + pja(11)*w(54) + pja(12)*w(84) + pja(13)*w(64) + pja(14)&
        *w(74) + pja(15)*w(84) + pja(16)*w(94)
      suma(5) = pja(1)*w(5) + pja(2)*w(15) + pja(3)*w(35) + pja(4)*w(65) + pja(&
        5)*w(15) + pja(6)*w(25) + pja(7)*w(45) + pja(8)*w(75) + pja(9)*w(35) + &
        pja(10)*w(45) + pja(11)*w(55) + pja(12)*w(85) + pja(13)*w(65) + pja(14)&
        *w(75) + pja(15)*w(85) + pja(16)*w(95)
      suma(6) = pja(1)*w(6) + pja(2)*w(16) + pja(3)*w(36) + pja(4)*w(66) + pja(&
        5)*w(16) + pja(6)*w(26) + pja(7)*w(46) + pja(8)*w(76) + pja(9)*w(36) + &
        pja(10)*w(46) + pja(11)*w(56) + pja(12)*w(86) + pja(13)*w(66) + pja(14)&
        *w(76) + pja(15)*w(86) + pja(16)*w(96)
      suma(7) = pja(1)*w(7) + pja(2)*w(17) + pja(3)*w(37) + pja(4)*w(67) + pja(&
        5)*w(17) + pja(6)*w(27) + pja(7)*w(47) + pja(8)*w(77) + pja(9)*w(37) + &
        pja(10)*w(47) + pja(11)*w(57) + pja(12)*w(87) + pja(13)*w(67) + pja(14)&
        *w(77) + pja(15)*w(87) + pja(16)*w(97)
      suma(8) = pja(1)*w(8) + pja(2)*w(18) + pja(3)*w(38) + pja(4)*w(68) + pja(&
        5)*w(18) + pja(6)*w(28) + pja(7)*w(48) + pja(8)*w(78) + pja(9)*w(38) + &
        pja(10)*w(48) + pja(11)*w(58) + pja(12)*w(88) + pja(13)*w(68) + pja(14)&
        *w(78) + pja(15)*w(88) + pja(16)*w(98)
      suma(9) = pja(1)*w(9) + pja(2)*w(19) + pja(3)*w(39) + pja(4)*w(69) + pja(&
        5)*w(19) + pja(6)*w(29) + pja(7)*w(49) + pja(8)*w(79) + pja(9)*w(39) + &
        pja(10)*w(49) + pja(11)*w(59) + pja(12)*w(89) + pja(13)*w(69) + pja(14)&
        *w(79) + pja(15)*w(89) + pja(16)*w(99)
      suma(10) = pja(1)*w(10) + pja(2)*w(20) + pja(3)*w(40) + pja(4)*w(70) + &
        pja(5)*w(20) + pja(6)*w(30) + pja(7)*w(50) + pja(8)*w(80) + pja(9)*w(40&
        ) + pja(10)*w(50) + pja(11)*w(60) + pja(12)*w(90) + pja(13)*w(70) + pja&
        (14)*w(80) + pja(15)*w(90) + pja(16)*w(100)
      sumb(1) = pjb(1)*w(1) + pjb(2)*w(2) + pjb(3)*w(4) + pjb(4)*w(7) + pjb(5)*&
        w(2) + pjb(6)*w(3) + pjb(7)*w(5) + pjb(8)*w(8) + pjb(9)*w(4) + pjb(10)*&
        w(5) + pjb(11)*w(6) + pjb(12)*w(9) + pjb(13)*w(7) + pjb(14)*w(8) + pjb(&
        15)*w(9) + pjb(16)*w(10)
      sumb(2) = pjb(1)*w(11) + pjb(2)*w(12) + pjb(3)*w(14) + pjb(4)*w(17) + pjb&
        (5)*w(12) + pjb(6)*w(13) + pjb(7)*w(15) + pjb(8)*w(18) + pjb(9)*w(14)&
         + pjb(10)*w(15) + pjb(11)*w(16) + pjb(12)*w(19) + pjb(13)*w(17) + pjb(&
        14)*w(18) + pjb(15)*w(19) + pjb(16)*w(20)
      sumb(3) = pjb(1)*w(21) + pjb(2)*w(22) + pjb(3)*w(24) + pjb(4)*w(27) + pjb&
        (5)*w(22) + pjb(6)*w(23) + pjb(7)*w(25) + pjb(8)*w(28) + pjb(9)*w(24)&
         + pjb(10)*w(25) + pjb(11)*w(26) + pjb(12)*w(29) + pjb(13)*w(27) + pjb(&
        14)*w(28) + pjb(15)*w(29) + pjb(16)*w(30)
      sumb(4) = pjb(1)*w(31) + pjb(2)*w(32) + pjb(3)*w(34) + pjb(4)*w(37) + pjb&
        (5)*w(32) + pjb(6)*w(33) + pjb(7)*w(35) + pjb(8)*w(38) + pjb(9)*w(34)&
         + pjb(10)*w(35) + pjb(11)*w(36) + pjb(12)*w(39) + pjb(13)*w(37) + pjb(&
        14)*w(38) + pjb(15)*w(39) + pjb(16)*w(40)
      sumb(5) = pjb(1)*w(41) + pjb(2)*w(42) + pjb(3)*w(44) + pjb(4)*w(47) + pjb&
        (5)*w(42) + pjb(6)*w(43) + pjb(7)*w(45) + pjb(8)*w(48) + pjb(9)*w(44)&
         + pjb(10)*w(45) + pjb(11)*w(46) + pjb(12)*w(49) + pjb(13)*w(47) + pjb(&
        14)*w(48) + pjb(15)*w(49) + pjb(16)*w(50)
      sumb(6) = pjb(1)*w(51) + pjb(2)*w(52) + pjb(3)*w(54) + pjb(4)*w(57) + pjb&
        (5)*w(52) + pjb(6)*w(53) + pjb(7)*w(55) + pjb(8)*w(58) + pjb(9)*w(54)&
         + pjb(10)*w(55) + pjb(11)*w(56) + pjb(12)*w(59) + pjb(13)*w(57) + pjb(&
        14)*w(58) + pjb(15)*w(59) + pjb(16)*w(60)
      sumb(7) = pjb(1)*w(61) + pjb(2)*w(62) + pjb(3)*w(64) + pjb(4)*w(67) + pjb&
        (5)*w(62) + pjb(6)*w(63) + pjb(7)*w(65) + pjb(8)*w(68) + pjb(9)*w(64)&
         + pjb(10)*w(65) + pjb(11)*w(66) + pjb(12)*w(69) + pjb(13)*w(67) + pjb(&
        14)*w(68) + pjb(15)*w(69) + pjb(16)*w(70)
      sumb(8) = pjb(1)*w(71) + pjb(2)*w(72) + pjb(3)*w(74) + pjb(4)*w(77) + pjb&
        (5)*w(72) + pjb(6)*w(73) + pjb(7)*w(75) + pjb(8)*w(78) + pjb(9)*w(74)&
         + pjb(10)*w(75) + pjb(11)*w(76) + pjb(12)*w(79) + pjb(13)*w(77) + pjb(&
        14)*w(78) + pjb(15)*w(79) + pjb(16)*w(80)
      sumb(9) = pjb(1)*w(81) + pjb(2)*w(82) + pjb(3)*w(84) + pjb(4)*w(87) + pjb&
        (5)*w(82) + pjb(6)*w(83) + pjb(7)*w(85) + pjb(8)*w(88) + pjb(9)*w(84)&
         + pjb(10)*w(85) + pjb(11)*w(86) + pjb(12)*w(89) + pjb(13)*w(87) + pjb(&
        14)*w(88) + pjb(15)*w(89) + pjb(16)*w(90)
      sumb(10) = pjb(1)*w(91) + pjb(2)*w(92) + pjb(3)*w(94) + pjb(4)*w(97) + &
        pjb(5)*w(92) + pjb(6)*w(93) + pjb(7)*w(95) + pjb(8)*w(98) + pjb(9)*w(94&
        ) + pjb(10)*w(95) + pjb(11)*w(96) + pjb(12)*w(99) + pjb(13)*w(97) + pjb&
        (14)*w(98) + pjb(15)*w(99) + pjb(16)*w(100)
      i = 0
      do i5 = 1, 4
        iia = ia + i5 - 1
        ija = ja + i5 - 1
        ioff = (iia*(iia - 1))/2 + ia - 1
        joff = (ija*(ija - 1))/2 + ja - 1
        do i6 = 1, i5
          ioff = ioff + 1
          joff = joff + 1
          i = i + 1
          f(ioff) = f(ioff) + sumb(i)
          f(joff) = f(joff) + suma(i)
        end do
      end do
      return
      end subroutine jab
