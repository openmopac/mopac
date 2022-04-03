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

subroutine minv (a, n, d)
  implicit none
  double precision, dimension(*), intent(inout) :: a
  integer, intent(in) :: n
  double precision, intent(out) :: d
!
  integer :: i, ij, ik, iz, j, ji, jk, jp, jq, jr, k, ki, kj, kk, nk
  double precision :: biga, hold
  integer, dimension(n) :: l, m
  !
  !     SEARCH FOR LARGEST ELEMENT
  !
  d = 1.0d0
  nk = -n
  do k = 1, n
    nk = nk + n
    l(k) = k
    m(k) = k
    kk = nk + k
    biga = a(kk)
    do j = k, n
      iz = n * (j-1)
      do i = k, n
        ij = iz + i
        if (Abs (biga) < Abs (a(ij))) then
          biga = a(ij)
          l(k) = i
          m(k) = j
        end if
      end do
    end do
    !
    !     INTERCHANGE ROWS
    !
    j = l(k)
    if (j > k) then
      ki = k - n
      do i = 1, n
        ki = ki + n
        hold = -a(ki)
        ji = ki - k + j
        a(ki) = a(ji)
        a(ji) = hold
      end do
    end if
    !
    !     INTERCHANGE COLUMNS
    !
    i = m(k)
    if (i > k) then
      jp = n * (i-1)
      do j = 1, n
        jk = nk + j
        ji = jp + j
        hold = -a(jk)
        a(jk) = a(ji)
        a(ji) = hold
      end do
    end if
    !
    !     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
    !     CONTAINED IN BIGA)
    !
    if (biga == 0.0d0) then
      d = 0.d0
      return
    end if
    do i = 1, n
      if (i /= k) then
        ik = nk + i
        a(ik) = a(ik) / (-biga)
      end if
    end do
    !  REDUCE MATRIX
    do i = 1, n
      ik = nk + i
      hold = a(ik)
      ij = i - n
      do j = 1, n
        ij = ij + n
        if (i/=k .and. j/=k) then
          kj = ij - i + k
          a(ij) = hold*a(kj) + a(ij)
        end if
      end do
    end do
    !
    !     DIVIDE ROW BY PIVOT
    !
    kj = k - n
    do j = 1, n
      kj = kj + n
      if (j /= k) then
        a(kj) = a(kj) / biga
      end if
    end do
    !
    !     PRODUCT OF PIVOTS
    !
    d = Max (-1.d25, Min (1.d25, d))
    d = d * biga
    !
    !     REPLACE PIVOT BY RECIPROCAL
    !
    a(kk) = 1.d0 / biga
  end do
  !
  !     FINAL ROW AND COLUMN INTERCHANGE
  !
  k = n
  do
    k = (k-1)
    if (k <= 0) exit
    i = l(k)
    if (i > k) then
      jq = n * (k-1)
      jr = n * (i-1)
      do j = 1, n
        jk = jq + j
        hold = a(jk)
        ji = jr + j
        a(jk) = -a(ji)
        a(ji) = hold
      end do
    end if
    j = m(k)
    if (j > k) then
      ki = k - n
      do i = 1, n
        ki = ki + n
        hold = a(ki)
        ji = ki - k + j
        a(ki) = -a(ji)
        a(ji) = hold
      end do
    end if
  end do
end subroutine minv
