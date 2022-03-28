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

      subroutine deri22(c, b, work, foc2, ab, minear, fci, w, diag, &
        scalar, ninear)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : nopen, numat, norbs, n2elec
      use common_arrays_C, only : nfirst, nlast
      use meci_C, only : nmos, nbo, nelec
      implicit none
      integer , intent(in) :: minear
      integer , intent(in) :: ninear
      double precision :: c(norbs,norbs), work(norbs,norbs), foc2(norbs**2)
      double precision, allocatable  :: dp(:), dpa(:)
      double precision  :: b(minear)
      double precision  :: ab(minear)
      double precision  :: fci(ninear)
      double precision  :: w(n2elec)
      double precision , intent(in) :: diag(*)
      double precision , intent(in) :: scalar(*)
!
      integer :: i, l, j, nend, loop, ninit, n1, n2, icount, ncol
      double precision :: wj(1) = 0.d0, wk(1) = 0.d0  ! Dummys - wj and wk are used by solids only.
      double precision, external :: ddot, dot
!***********************************************************************
!  1) BUILD THE 2-ELECTRON FOCK MATRIX DEPENDING ON B AS FOLLOWS :
!     DP = C * SCALE*B * C' ...  DP DENSITY MATRIX 'DERIVATIVE',
!     FOC2 = 0.5 * TRACE ( DP * (2<J>-<K>) ) DONE IN FOCK2 & FOCK1.
!  2) HALF-TRANSFORM ONTO M.O. BASIS : DPT =  FOC2 * C
!     AND COMPUTE DIAGONAL BLOCKS ELEMENTS OF C' * FOC2, EXTRACTING
!     IN FCI ELEMENTS OVER C.I-ACTIVE M.O ONLY.
!  3) COMPUTE SUPERVECTOR AB = (DIAG + A) * B DEFINED BY THE MATRIX :
!     AB(I,J)= ( DIAG(I,J)*B(I,J)+DPT(I,J) )*SCALAR(I,J)  WITH I.GT.J,
!     DIAG(I,J)=(EIGS(I)-EIGS(J))/(O(J)-O(I)) >0, O OCCUPANCY NUMBERS,
!     EIGS EIGENVALUES OF FOCK OPERATOR WITH EIGENVECTORS C IN A.O.
!
!   INPUT
! C(NORBS,NORBS)   : M.O. EIGENVECTORS (COLUMNWISE).
! B(*)             : B SUPERVECTOR PACKED BY OFF-DIAGONAL BLOCKS, SCALED
! WORK(*)          : WORK AREA OF SIZE N*N.
! NORBS            : NUMBER OF M.O.S
! NELEC,NMOS       : LAST FROZEN CORE M.O. , C.I-ACTIVE BAND LENGTH.
!           IN COMMON
! DIAG,SCALAR AS DEFINED IN 'DERI0'.
!   OUTPUT
! FOC2(*)       : 2-ELECTRON FOCK MATRIX, PACKED CANONICAL.
! AB(*)         : ANTISYMMETRIC MATRIX PACKED IN SUPERVECTOR FORM WITH
!                 THE CONSECUTIVE FOLLOWING BLOCKS:
!              1) OPEN-CLOSED  I.E. B(IJ)=B(I,J) WITH I OPEN & J CLOSED
!                 AND I RUNNING FASTER THAN J,
!              2) VIRTUAL-CLOSED SAME RULE OF ORDERING,
!              3) VIRTUAL-OPEN   SAME RULE OF ORDERING.
! FCI(*)        : FOCK DIAGONAL BLOCKS ELEMENTS OVER C.I-ACTIVE M.O.
!            FOC2 CAN BE EQUIVALENCED WITH WORK IN THE CALLING SEQUENCE.
!***********************************************************************
!
!  NOTE: NORBS AND NORD ARE THE SAME ADDRESS.  THE NAME NORBD IS NOT
!        USED HERE.
!
!
!     DERIVATIVE OF THE DENSITY MATRIX IN DP (PACKED,CANONICAL).
!     ----------------------------------------------------------
!     DP = C * B * C' .
!
!     STEP 0 : UNSCALE VECTOR B.
      n2 = 0
      allocate (dp(norbs**2), dpa(norbs**2), stat=i)
      if (i /= 0) then
        call memory_error ("deri22")
        return
      end if
      b = b*scalar(:minear)
!
!     STEP 1 : WORK = C * B    .  DP TEMPORARY ARRAY.
      l = 1
      if (nbo(2)/=0 .and. nbo(1)/=0) then
!        OPEN-CLOSED
        call mxm (c(1,nbo(1)+1), norbs, b(l), nbo(2), work, nbo(1))
!       CLOSED-OPEN
        call mxmt (c, norbs, b(l), nbo(1), work(1,nbo(1)+1), nbo(2))
        l = l + nbo(2)*nbo(1)
      end if
      if (nbo(3)/=0 .and. nbo(1)/=0) then
!       VIRTUAL-CLOSED
        if (l > 1) then
          call mxm (c(1,nopen+1), norbs, b(l), nbo(3), dp, nbo(1))
          icount = 0
        do j = 1, nbo(1)
          do i = 1, norbs
            icount = icount + 1
            work(i, j) = work(i, j) + dp(icount)
          end do
        end do
        else
          call mxm (c(1,nopen+1), norbs, b(l), nbo(3), work, nbo(1))
        end if
!       CLOSED-VIRTUAL
        call mxmt (c, norbs, b(l), nbo(1), work(1,nopen+1), nbo(3))
        l = l + nbo(3)*nbo(1)
      end if
      if (nbo(3)/=0 .and. nbo(2)/=0) then
!       VIRTUAL-OPEN
        call mxm (c(1,nopen+1), norbs, b(l), nbo(3), dp, nbo(2))
         icount = 0
      do j = nbo(1) + 1, nbo(1) + nbo(2)
        do i = 1, norbs
          icount = icount + 1
          work(i, j) = work(i, j) + dp(icount)
        end do
      end do
!       OPEN-VIRTUAL
        call mxmt (c(1,nbo(1)+1), norbs, b(l), nbo(2), dp, nbo(3))
     icount = 0
      do j = nopen + 1, nopen + nbo(3)
        do i = 1, norbs
          icount = icount + 1
          work(i, j) = work(i, j) + dp(icount)
        end do
      end do
      end if
!
!     STEP 2 : DP= WORK * C'   WITH DP PACKED,CANONICAL.
      l = 0
      do i = 1, norbs
        do j = 1, i
          l = l + 1
          dp(l) = ddot(norbs,work(i,1),norbs,c(j,1),norbs)
        end do
      end do
!
!     2-ELECTRON FOCK MATRIX BUILD WITH THE DENSITY MATRIX DERIVATIVE.
!     ----------------------------------------------------------------
!     RETURNED IN FOC2 (PACKED CANONICAL).
      foc2 = 0.D0
      dpa(:l) = 0.5D0*dp(:l)
        call fock2 (foc2, dp, dpa, w, wj, wk, numat, nfirst, nlast, 2)
!
!     BUILD DP AND EXTRACT FCI.
!     --------------------------
!
!     DP(NORBS,NEND) = FOC2(NORBS,NORBS) * C(NORBS,NEND).
      nend = max(nopen,nelec + nmos)
      l = 1
      do i = 1, nend
        call supdot (dp(l), foc2, c(1,i), norbs)
        l = l + norbs
      end do
!     EXTRACT FCI
      l = 1
      nend = 0
      do loop = 1, 3
        ninit = nend + 1
        nend = nend + nbo(loop)
        n1 = max(ninit,nelec + 1)
        n2 = min(nend,nelec + nmos)
        if (n2 < n1) cycle
        do i = n1, n2
          if (i <= ninit) cycle
          ! TODO: GBR future modifications
          call mxm (c(1,i), 1, dp(norbs*(ninit-1)+1), norbs, fci(l), i - ninit)
          l = l + i - ninit
        end do
      end do
      ncol = n2 - ninit + 1
      if (ncol > 0 .and. n2 < norbs) then
        ! TODO: GBR future modifications
        call mtxm (c(1, n2+1), norbs-n2, dp(norbs*(ninit-1)+1), norbs, fci(l), ncol)
        l = l + ncol * (norbs-n2)
      end if
      do i = nelec + 1, nelec + nmos
        fci(l) = -dot(c(1,i),dp(norbs*(i-1)+1),norbs)
        l = l + 1
      end do
!
!     NEW SUPERVECTOR AB = (DIAG + C'* FOC2 * C) * B , SCALED.
!     --------------------------------------------------------
!
!     PART 1 : AB(I,J) = (C' * DP)(I,J) DONE BY BLOCKS.
      l = 1
      if (nbo(2)/=0 .and. nbo(1)/=0) then
        ! TODO: GBR future modifications
        call mtxm (c(1,nbo(1)+1), nbo(2), dp, norbs, ab(l), nbo(1))
        l = l + nbo(2)*nbo(1)
      end if
      if (nbo(3)/=0 .and. nbo(1)/=0) then
        call mtxm (c(1,nopen+1), nbo(3), dp, norbs, ab(l), nbo(1))
        l = l + nbo(3)*nbo(1)
      end if
      if (nbo(3)/=0 .and. nbo(2)/=0) call mtxm (c(1,nopen+1), nbo(3), &
        dp(norbs*nbo(1)+1), norbs, ab(l), nbo(2))
!
!     PART 2 : AB = SCALE * (D * B + AB) AND RESCALE BASIS VECTOR B.
!
      ab = (diag(:minear)*b+ab)*scalar(:minear)
      b = b/scalar(:minear)
      return
      end subroutine deri22
