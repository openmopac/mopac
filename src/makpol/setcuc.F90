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

subroutine setcuc (coord, natoms, tvec, id, xtoc)
    use common_keywrd, only : keywrd
    implicit none
    integer, intent (in) :: natoms
    double precision, dimension (3, natoms), intent (inout) :: coord
    double precision, dimension (3, 3), intent (in) :: tvec
    double precision, dimension (3, 3) :: xtoc
    integer, intent (in) :: id
    integer :: i, j, j1, j2, j3, k, l, store_l(3)
    double precision :: sum, summax, summin, storen(3), xyz(3)
    logical :: l_bcc
    l_bcc = (index(keywrd, " BCC") /= 0)
!
!  Make normalized translation vectors
!
    do i = 1, id
      sum = 0.d0
      do j = 1, 3
        sum = sum + tvec(j, i)**2
      end do
      storen(i) = Sqrt (sum)
      sum = 1.d0 / storen(i)
      do j = 1, 3
        xtoc(j, i) = tvec(j, i) * sum
      end do
    end do
!
!   Construct vectors perpendicular to known vectors,
!   if system is of dimension 1 or 2
!
    if (id == 1) then
!
!    Make a vector orthogonal to the first vector
!
      summax = 0.d0
      summin = 1.d0
      do j = 1, 3
        if (Abs (summax) < Abs (xtoc(j, 1))) then
          j1 = j
          summax = xtoc(j, 1)
        end if
        if (Abs (summin) >= Abs (xtoc(j, 1))) then
          j2 = j
          summin = xtoc(j, 1)
        end if
      end do
      j3 = 6 - j1 - j2
      xtoc(j3, 2) = 0.d0
      xtoc(j2, 2) = summax
      xtoc(j1, 2) = -summin
    end if
    if (id /= 3) then
!
!  Make a vector orthogonal to the first and second vectors.
!
      summin = 1.d0
      do j = 1, 3
        sum = xtoc(j, 1)**2 + xtoc(j, 2)**2
        if (summin > sum) then
          j2 = j
          summin = sum
        end if
      end do
      xtoc(j2, 3) = 1.d0
      do j = 1, 3
        if (j /= j2) then
          sum = 0.d0
          do i = 1, 2
            sum = sum + xtoc(j, i) * xtoc(j2, i)
          end do
          xtoc(j, 3) = -sum
        end if
      end do
      sum = 1.d0 / Sqrt (xtoc(1, 3)**2+xtoc(2, 3)**2+xtoc(3, 3)**2)
      do j = 1, 3
        xtoc(j, 3) = xtoc(j, 3) * sum
      end do
    end if
!
!   Expand normalized vectors to their full size
!
    do i = 1, id
      do j = 1, 3
        xtoc(j, i) = xtoc(j, i) * storen(i)
      end do
    end do
!
!   Make unused vectors of length 10
!
    do i = id + 1, 3
      do j = 1, 3
        xtoc(j, i) = xtoc(j, i) * 10.d0
      end do
    end do
!
!  Invert translation vector matrix
!
    call minv (xtoc, 3, sum)
!
!  FORCE ALL ATOMS IN THE CENTRAL UNIT CELL INTO THE DOMAIN 0-1
!
    store_l = 0
    do i = 1, natoms
    xyz = 0.d0
      do k = 1, 3        
        do j = 1, id
          xyz(j) = xyz(j) + coord(k, i) * xtoc(j, k)
        end do
      end do
!
!   XYZ holds the crystal fractional coordinates
!
      k = 0
      do j = 1, id
        l = Int (xyz(j)+20.001d0) - 20
        k = k + l
        store_l(j) = l
      end do
      if (.not. l_bcc .or. mod(k,2) == 0) then
        continue
        do j = 1, id
          l = store_l(j)
          do k = 1, 3
            coord(k, i) = coord(k, i) - l * tvec(k, j)
          end do
        end do
      end if
      continue
    end do
    return
end subroutine setcuc
