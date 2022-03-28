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

      subroutine jcarin(xparam, step, preci, b, ncol, il, iu)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numat, nvar, ndep, id, l1u, l2u, l3u, l123, &
      natoms
      use common_arrays_C, only : geo, coord, loc, tvec
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: ncol
      integer , intent(in) :: il
      integer , intent(in) :: iu
      double precision , intent(in) :: step
      logical , intent(in) :: preci
      double precision , intent(in) :: xparam(3*natoms)
      double precision , intent(inout) :: b(iu - il + 1,3*numat*l123)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ivar, jvar, j, ij, ii, im, jl, kl, j1, j2
      double precision, dimension(:), allocatable :: coold
!-----------------------------------------------
!     JACOBIAN dCARTESIAN/dINTERNAL, WORKED OUT BY FINITE DIFFERENCE.
!  INPUT
!     XPARAM(*) : INTERNAL COORDINATES
!     STEP      : STEP SIZE FOR FINITE DIFFERENCE METHOD
!     PRECI     : .TRUE. IF 2-POINTS FINITE DIFFERENCES TO BE USED,
!                 .FALSE. OTHERWISE.
!  OUTPUT
!     B(NVAR,NCOL) : JACOBIAN, STEP TIME TOO LARGE.
!
!
       ncol = 3*numat*l123
       allocate(coold(3*numat*l123))
!
!     INTERNAL OF CENTRAL POINT
      do ivar = 1, nvar
        geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar)
      end do
!
      if (id == 0) then
!
!        MOLECULAR SYSTEM
!        ----------------
        jvar = 0
        do ivar = il, iu
          jvar = jvar + 1
!        STEP FORWARD
          geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar) + step
          if (ndep /= 0) call symtry
          call gmetry (geo, coord)
          j = 0
          do j1 = 1, numat
            do j2 = 1, 3
              j = j + 1
              b(jvar, j) = coord(j2, j1)
            end do
          end do
          geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar)
        end do
        if (preci) then
          jvar = 0
          do ivar = il, iu
            jvar = jvar + 1
!           STEP BACKWARD
            geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar) - step
            if (ndep /= 0) call symtry
            call gmetry (geo, coord)
            j = 0
            do j1 = 1, numat
              do j2 = 1, 3
                j = j + 1
                b(jvar, j) = b(jvar, j) - coord(j2, j1)
              end do
            end do
            geo(loc(2, ivar), loc(1, ivar)) = xparam(ivar)
          end do
        else
!           CENTRAL POINT
          if (ndep /= 0) call symtry
          call gmetry (geo, coord)
          jvar = 0
          do ivar = il, iu
            jvar = jvar + 1
            j = 0
            do j1 = 1, numat
              do j2 = 1, 3
                j = j + 1
                b(jvar, j) = b(jvar, j) - coord(j2, j1)
              end do
            end do
          end do
        end if
      else
!
!        SOLID STATE
!        -----------
        jvar = 0
        do ivar = il, iu
          jvar = jvar + 1
!        STEP FORWARD
          geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar) + step
          if (ndep /= 0) call symtry
          call gmetry (geo, coord)
          ij = 0
          do ii = 1, numat
            do im = -l1u, l1u
              do jl = -l2u, l2u
                do kl = -l3u, l3u
                  b(jvar,ij+1:3+ij) = coord(:,ii) + tvec(:,1)*im + tvec(:,2)*jl + tvec(:,3)*kl
                  ij = 3 + ij
                end do
              end do
            end do
          end do
          geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar)
        end do
        if (preci) then
          jvar = 0
          do ivar = il, iu
            jvar = jvar + 1
!           STEP BACKWARD
            geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar) - step
            if (ndep /= 0) call symtry
            call gmetry (geo, coord)
            ij = 0
            do ii = 1, numat
              do im = -l1u, l1u
                do jl = -l2u, l2u
                  do kl = -l3u, l3u
                    b(jvar,ij+1:3+ij) = b(jvar,ij+1:3+ij) - &
                    coord(:,ii) - tvec(:,1)*im - tvec(:,2)*jl - tvec(:,3)*kl
                    ij = 3 + ij
                  end do
                end do
              end do
            end do
            geo(loc(2,ivar),loc(1,ivar)) = xparam(ivar)
          end do
        else
!           CENTRAL POINT
          if (ndep /= 0) call symtry
          call gmetry (geo, coord)
          ij = 0
          do ii = 1, numat
            do im = -l1u, l1u
              do jl = -l2u, l2u
                do kl = -l3u, l3u
                  ij = ij + 1
                  coold(3*ij - 2) = coord(1,ii) + tvec(1,1)*im + tvec(1,2)*jl + tvec(1,3)*kl
                  coold(3*ij - 1) = coord(2,ii) + tvec(2,1)*im + tvec(2,2)*jl + tvec(2,3)*kl
                  coold(3*ij - 0) = coord(3,ii) + tvec(3,1)*im + tvec(3,2)*jl + tvec(3,3)*kl
                end do
              end do
            end do
          end do
          jvar = 0
          do ivar = il, iu
            jvar = jvar + 1
            b(jvar,:ncol) = b(jvar,:ncol) - coold(:ncol)
          end do
        end if
      end if
      return
      end subroutine jcarin
