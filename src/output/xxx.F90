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

      subroutine xxx(type, i, j, k, l, r)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i
      integer , intent(in) :: j
      integer , intent(in) :: k
      integer , intent(in) :: l
      character , intent(in) :: type
      character , intent(out) :: r*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(4) :: ijk
      integer :: m, loop, ii, i2
!-----------------------------------------------
!***********************************************************************
!
!    XXX WILL FORM A UNIQUE STRING LABEL 'R' FOR GAUSSIAN-TYPE INPUT
!    THE LABEL WILL BE LETTER (EITHER R, P, OR F, NORMALLY), FOLLOWED
!    BY THE CONNECTIVITY, IN THE ORDER I, J, K, L.
!    'R' IS 13 CHARACTERS LONG IN ORDER TO ACCOMMODATE 3 DIGITS PER
!    LABEL, WHEN NECESSARY
!
!***********************************************************************
      r = type
      ijk(1) = i
      ijk(2) = j
      ijk(3) = k
      ijk(4) = l
      m = 1
      do loop = 1, 4
        ii = ijk(loop)
        if (ii == 0) cycle
!
!   IF LABELS GREATER THAN 99 ARE USED, UNCOMMENT THE FOLLOWING CODE
!
!#         I2=II/100
!#         IF(I2.NE.0) THEN
!#            M=M+1
!#            R(M:M)=CHAR(ICHAR('0')+I2)
!#            II=II-I2*100
!#         ENDIF
        i2 = ii/10
        if (i2 /= 0) then
          m = m + 1
          r(m:m) = char(ichar('0') + i2)
          ii = ii - i2*10
        end if
        m = m + 1
        r(m:m) = char(ichar('0') + ii)
      end do
      return
      end subroutine xxx
