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

      subroutine matout(a, b, nc, nr, ndim)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numat, method_indo
      use common_arrays_C, only : nfirst, nlast, nat
      use chanel_C, only : iw
      use elemts_C, only : elemnt
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nc
      integer , intent(inout) :: nr
      integer , intent(in) :: ndim
      double precision , intent(in) :: a(ndim,ndim)
      double precision , intent(in) :: b(ndim)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(:), allocatable :: natom
      integer :: i, jlo, jhi, l, j, ka, kc, kb, la, lc, lb
      character :: ref_atorbs(9,2)*2, atorbs(9)*2
      character, dimension(:), allocatable :: itext*2, jtext*2

      save atorbs
!-----------------------------------------------
!**********************************************************************
!
!      MATOUT PRINTS A SQUARE MATRIX OF EIGENVECTORS AND EIGENVALUES
!
!    ON INPUT A CONTAINS THE MATRIX TO BE PRINTED.
!             B CONTAINS THE EIGENVALUES.
!             NC NUMBER OF MOLECULAR ORBITALS TO BE PRINTED.
!             NR IS THE SIZE OF THE SQUARE ARRAY TO BE PRINTED.
!             NDIM IS THE ACTUAL SIZE OF THE SQUARE ARRAY "A".
!             NFIRST AND NLAST CONTAIN ATOM ORBITAL COUNTERS.
!             NAT = ARRAY OF ATOMIC NUMBERS OF ATOMS.
!
!
!***********************************************************************
      data ref_atorbs/ ' S', 'PX', 'PY', 'PZ', 'X2', 'XZ', 'Z2', 'YZ', 'XY', &
                       ' S', 'PX', 'PY', 'PZ', 'Z2', 'X2', 'XY', 'XZ', 'YZ'/
      if (method_indo) then
        atorbs(:) = ref_atorbs(:,2)
      else
        atorbs(:) = ref_atorbs(:,1)
      end if
      i = max(ndim, 1)
      allocate (itext(i), jtext(i), natom(i))
      if (nlast(numat) /= nr) go to 30
      do i = 1, numat
        jlo = nfirst(i)
        jhi = nlast(i)
        l = nat(i)
        itext(jlo:jhi) = atorbs(:jhi-jlo+1)
        jtext(jlo:jhi) = elemnt(l)
        natom(jlo:jhi) = i
      end do
      go to 50
   30 continue
      nr = abs(nr)
      do i = 1, nr
        itext(i) = '  '
        jtext(i) = '  '
        natom(i) = i
      end do
   50 continue
      ka = 1
      kc = 6
   60 continue
      kb = min0(kc,nc)
      write (iw, 100) (i,i=ka,kb)
      if (b(1) /= 0.D0) write (iw, 110) (b(i),i=ka,kb)
      write (iw, 120)
      la = 1
      lc = 40
   70 continue
      lb = min0(lc,nr)
      do i = la, lb
        if (itext(i) == ' S') write (iw, 120)
        write (iw, 130) itext(i), jtext(i), natom(i), (a(i,j),j=ka,kb)
      end do
      if (lb == nr) go to 90
      la = lc + 1
      lc = lc + 40
      go to 70
   90 continue
      if (kb == nc) return
      ka = kc + 1
      kc = kc + 6
      go to 60
!
  100 format(/,/,/,/,3x,' ROOT NO.',i5,9i12)
  110 format(/,8x,10f12.5)
  120 format('  ')
  130 format(2(1x,a2),i4,f10.5,10f12.5)
!
      end subroutine matout
