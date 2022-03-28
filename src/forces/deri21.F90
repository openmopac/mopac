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

      subroutine deri21(a, nvar_nvo, minear, ifirst, vnert, pnert, b, ncut)
      implicit none
      integer  :: nvar_nvo
      integer  :: minear
      integer , intent(in) :: ifirst
      integer  :: ncut
      double precision  :: a(minear,nvar_nvo)
      double precision  :: vnert(nvar_nvo)
      double precision  :: pnert(nvar_nvo)
      double precision  :: b(minear,nvar_nvo)
!
      integer :: i, l
      double precision, dimension(4) :: work
      double precision :: cutoff, sum2, sum
!***********************************************************************
!
!     LEAST-SQUARE ANALYSIS OF A SET OF ANALYTICAL POINTS {A} :
!
!     PRODUCE A SUBSET OF NCUT ORTHONORMALIZED VECTORS B, OPTIMUM IN A
!     LEAST-SQUARE SENSE WITH RESPECT TO THE INITIAL SPACE {A}.
!     EACH NEW HIERARCHIZED VECTOR B EXTRACTS A MAXIMUM PERCENTAGE FROM
!     THE REMAINING DISPERSION OF THE SET {A} OUT OF THE PREVIOUS
!     {B} SUBSPACE.
!   INPUT
!     A(MINEAR,nvar_nvo): ORIGINAL SET {A}.
!     ifirst        : MAXIMUM ALLOWED SIZE OF THE BASIS B.
!   OUTPUT
!     VNERT(nvar_nvo)   : LOWEST EIGENVECTOR OF A'* A.
!     PNERT(nvar_nvo)     : SQUARE ROOT OF THE ASSOCIATED EIGENVALUES
!                     IN DECREASING ORDER.
!     B(MINEAR,NCUT): OPTIMUM ORTHONORMALIZED SUBSET {B}.
!
!***********************************************************************
!
!     VNERT = A' * A
      cutoff = 0.85D0
      sum2 = 0.D0
      call mtxmc (a, nvar_nvo, a, minear, work)
      work(:(nvar_nvo*(nvar_nvo + 1))/2) = -work(:(nvar_nvo*(nvar_nvo + 1))/2)
!     DIAGONALIZE IN DECREASING ORDER OF EIGENVALUES
      if (abs(work(1))<1.D-28 .and. nvar_nvo==1) then
        pnert(1) = sqrt((-work(1)))
        work(1) = 1.D15
        vnert(1) = 1.D0
        ncut = 1
        go to 50
      else
        call rsp (work, nvar_nvo, pnert, vnert)
!     FIND NCUT ACCORDING TO CUTOFF, BUILD WORK = VNERT * (PNERT)**-0.5
        sum = 0.D0
        do i = 1, nvar_nvo
          sum = sum - pnert(i)
        end do
        l = 1
        do i = 1, ifirst
          sum2 = sum2 - pnert(i)/sum
          pnert(i) = sqrt(Abs(pnert(i)))
          if (nvar_nvo > 0) then
            work(l:nvar_nvo-1+l) = vnert(l:nvar_nvo-1+l)/pnert(i)
            l = nvar_nvo + l
          end if
          if (sum2 < cutoff) cycle
          ncut = i
          go to 50
        end do
        ncut = ifirst
!     ORTHONORMALIZED BASIS
!     B(MINEAR,NCUT) = A(MINEAR,nvar_nvo)*WORK(nvar_nvo,NCUT)
      end if
   50 continue
      call mxm (a, minear, work, nvar_nvo, b, ncut)
      return
      end subroutine deri21
