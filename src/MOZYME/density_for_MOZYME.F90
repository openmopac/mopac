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

subroutine density_for_MOZYME (p, mode, nclose_loc, partp)
   !***********************************************************************
   !
   !   DENSIT COMPUTES THE DENSITY MATRIX GIVEN THE EIGENVECTOR MATRIX, AND
   !          INFORMATION ABOUT THE M.O. OCCUPANCY.
   !
   !  INPUT:  COCC  = COMPRESSED OCCUPIED EIGENVECTOR MATRIX
   !
   !   ON EXIT: P   = DENSITY MATRIX
   !
   !***********************************************************************
    use molkst_C, only: numat, mpack, keywrd
    use MOZYME_C, only : lijbo, nijbo, ncf, ncocc, &
      nncf, iorbs, cocc, icocc
    use chanel_C, only: iw
    implicit none
    integer, intent (in) :: mode, nclose_loc
    double precision, dimension (mpack), intent (in) :: partp
    double precision, dimension (mpack), intent (inout) :: p
    logical :: first = .true.
    logical, save :: prnt
    integer :: i, j, j1, ja, jj, k, k1, k2, ka, kk, l, loop, nj
    double precision :: spinfa, sum
    integer, external :: ijbo
    if (first) then
      first = .false.
      prnt = Index (keywrd, " DIAG ") /= 0
    end if
    if (mode == 0) then
      !
      !  INITIALIZE P:  FULL DENSITY CALCULATION
      !
      p(:) = 0.d0
    else if (mode ==-1) then
      !
      !   FULL DENSITY MATRIX IS IN PARTP, REMOVE DENSITY DUE TO
      !   OLD LMO's USED IN SCF
      !
      p(:) = -0.5d0 * partp(:)
    else
      !
      !   PARTIAL DENSITY MATRIX IS IN PARTP, BUILD THE REST OF P
      !
      p(:) = 0.5d0 * partp(:)
    end if

    do i = 1, nclose_loc
      loop = ncocc(i)
      ja = 0
      if (lijbo) then
        do jj = nncf(i) + 1, nncf(i) + ncf(i)
          j = icocc(jj)
          nj = iorbs(j)
          ka = loop
          do kk = nncf(i) + 1, nncf(i) + ncf(i)
            k = icocc(kk)
            if (j == k) then
              l = nijbo (j, k)
              do j1 = 1, nj
                sum = cocc(ja+j1+loop)
                do k1 = 1, j1
                  k2 = ka + k1
                  l = l + 1
                  p(l) = p(l) + cocc(k2) * sum
                end do
              end do
            else if (j > k .and. nijbo (j, k) >= 0) then
              l = nijbo (j, k)
              do j1 = 1, nj
                sum = cocc(ja+j1+loop)
                do k1 = 1, iorbs(k)
                  k2 = ka + k1
                  l = l + 1
                  p(l) = p(l) + cocc(k2) * sum
                end do
              end do
            end if
            ka = ka + iorbs(k)
          end do
          ja = ja + nj
        end do
      else
        do jj = nncf(i) + 1, nncf(i) + ncf(i)
          j = icocc(jj)
          nj = iorbs(j)
          ka = loop
          do kk = nncf(i) + 1, nncf(i) + ncf(i)
            k = icocc(kk)
            l = ijbo (j, k)
            if (j == k) then
              do j1 = 1, nj
                sum = cocc(ja+j1+loop)
                do k1 = 1, j1
                  k2 = ka + k1
                  l = l + 1
                  p(l) = p(l) + cocc(k2) * sum
                end do
              end do
            else if (j > k .and. l >= 0) then
              do j1 = 1, nj
                sum = cocc(ja+j1+loop)
                do k1 = 1, iorbs(k)
                  k2 = ka + k1
                  l = l + 1
                  p(l) = p(l) + cocc(k2) * sum
                end do
              end do
            end if
            ka = ka + iorbs(k)
          end do
          ja = ja + nj
        end do
      end if
    end do
    if (mode == 0 .or. mode == 1) then
      !
      !    FULL DENSITY CALCULATION.  MULTIPLY BY 2 FOR SPIN
      !
      spinfa = 2.d0
    else if (mode ==-1) then
      !
      !   MAKING PARTIAL DENSITY MATRIX. REVERSE SIGN ONCE MORE
      !
      spinfa = -2.d0
    else
      spinfa = 1.d0
    end if
   !
    if (Abs(spinfa - 1.d0) > 1.d-10) then
      p(:) = spinfa * p(:)
    end if
    if (prnt) then
      sum = 0.d0
      if (lijbo) then
        do i = 1, numat
          j = nijbo(i, i) + 1
          if (iorbs(i) > 0) then
            sum = sum + p(j)
          end if
          if (iorbs(i) == 4) then
            sum = sum + p(j+2) + p(j+5) + p(j+9)
          end if
        end do
      else
        do i = 1, numat
          j = ijbo (i, i) + 1
          if (iorbs(i) > 0) then
            sum = sum + p(j)
          end if
          if (iorbs(i) == 4) then
            sum = sum + p(j+2) + p(j+5) + p(j+9)
          end if
        end do
      end if
      write (iw,*) " COMPUTED NUMBER OF ELECTRONS:", sum
    end if
end subroutine density_for_MOZYME
